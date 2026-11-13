#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import sys

import json
import logging
import shutil
from contextlib import contextmanager
from pathlib import Path

from inputs import INPUTS, canonical_chrom, open_text, resolve_input, sha256
from outputs import emit_outputs, verify_outputs, verify_receipt, write_output_hashes, write_table
from preparation import (LOG, match_coding, match_intergenic, pool, prepare, select_lncrnas)

VERSION = '1.0.0'


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='Prepare human GENCODE/GRCh38 paired null-model FASTAs. '
                                'Does not run FragIncRNA or download inputs.')
    p.add_argument('--input-dir', type=Path, default=Path('input'))
    p.add_argument('--output-dir', type=Path, required=True)
    for key in INPUTS:
        p.add_argument('--' + key.replace('_', '-'), type=Path,
                       help=f'Explicit {key} override (otherwise resolved in --input-dir)')
    p.add_argument('--n-transcripts', type=int, default=5000)
    p.add_argument('--shuffles-per-transcript', type=int, default=10)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--threads', type=int, default=1, help='Preparation workers; does not launch searches')
    p.add_argument('--k-values', default='15,20,32', help='Window QC only; never changes FASTAs')
    p.add_argument('--chromosomes', default=','.join(map(str, range(1, 23))),
                   help='Explicit primary nuclear chromosome list, default autosomes; X/Y opt-in')
    p.add_argument('--max-per-gene', type=int, default=1, help='lncRNA isoform cap; 0 disables')
    p.add_argument('--lncrna-transcript-types', default='lncRNA,retained_intron,processed_transcript',
                   help='Accepted transcript biotypes within strictly lncRNA genes')
    p.add_argument('--coding-isoforms', action='store_true', help='Allow multiple controls from a coding gene')
    p.add_argument('--missing-coding-annotation', choices=('exclude', 'error'), default='exclude',
                   help='ALL-regions coding FASTA can exceed CHR GTF scope; list and exclude unresolved IDs, '
                   'or require every coding record to have an annotation')
    p.add_argument('--flank', type=int, default=10000, help='All-gene exclusion flank in bases')
    p.add_argument('--length-tolerance', type=float, default=.15, help='Maximum relative coding length error')
    p.add_argument('--gc-tolerance', type=float, default=.05, help='Maximum absolute GC fraction difference')
    p.add_argument('--repeat-tolerance', type=float, help='Optional maximum absolute repeat fraction error')
    p.add_argument('--max-splice-mismatch', type=float, default=0,
                   help='Maximum spliced-reference substitution fraction; exon lengths must always agree')
    p.add_argument('--gencode-version', type=int, default=49, help='Expected annotation version; 0 disables check')
    p.add_argument('--intergenic-candidates', type=int, default=512, help='Bounded draws per parent')
    p.add_argument('--shuffle-retries', type=int, default=32, help='Maximum shuffle attempts per replicate')
    p.add_argument('--cache-dir', type=Path, help='Reusable, provenance-checked disk cache; ~4 GB for human inputs')
    p.add_argument('--allow-partial', action='store_true', help='Permit explicit unmatched pairs/counts')
    p.add_argument('--dry-run', action='store_true', help='Read input headers and report estimates; no writes')
    group = p.add_mutually_exclusive_group()
    group.add_argument('--force', action='store_true', help='Atomically replace a completed output owned by this tool')
    group.add_argument('--resume', action='store_true', help='Verify completed run or restart compatible interrupted run')
    return p


def configuration(args: argparse.Namespace) -> dict:
    config = vars(args).copy()
    config['k_values'] = sorted(set(int(v) for v in args.k_values.split(',')))
    config['chromosomes'] = sorted(set(canonical_chrom(c.strip()) for c in args.chromosomes.split(',')))
    config['lncrna_transcript_types'] = sorted(set(v.strip() for v in args.lncrna_transcript_types.split(',') if v.strip()))
    if not config['lncrna_transcript_types']:
        raise ValueError('--lncrna-transcript-types must not be empty')
    primary = {str(i) for i in range(1, 23)} | {'X', 'Y'}
    if not config['chromosomes'] or set(config['chromosomes']) - primary:
        raise ValueError('--chromosomes must explicitly name primary nuclear chromosomes 1..22,X,Y')
    if not config['k_values'] or min(config['k_values']) < 1 or max(config['k_values']) > 32:
        raise ValueError('--k-values must be in 1..32')
    for key in ('n_transcripts', 'shuffles_per_transcript', 'threads', 'intergenic_candidates', 'shuffle_retries'):
        if config[key] < 1:
            raise ValueError(f'--{key.replace("_", "-")} must be positive')
    if args.flank < 0 or args.max_per_gene < 0:
        raise ValueError('Flank and isoform cap must be nonnegative')
    for key in ('length_tolerance', 'gc_tolerance', 'repeat_tolerance', 'max_splice_mismatch'):
        value = config[key]
        if value is not None and not 0 <= value <= 1:
            raise ValueError(f'{key} must be a fraction in [0,1]')
    # Runtime paths/worker counts do not change the scientific output fingerprint.
    omitted = {'input_dir', 'output_dir', 'cache_dir', 'threads', 'dry_run', 'force', 'resume'} | set(INPUTS)
    return {key: value for key, value in config.items() if key not in omitted}


def code_hashes() -> dict:
    return {path.name: sha256(path) for path in sorted(Path(__file__).parent.glob('*.py'))}


def input_hashes(paths: dict) -> dict:
    result = {}
    for key, path in paths.items():
        LOG.info('SHA-256 %s', path.name)
        result[key] = sha256(path)
    return result


def dry_run(paths: dict, config: dict) -> None:
    probes = {}
    for key, path in paths.items():
        with open_text(path) as handle:
            lines = []
            for _, line in zip(range(30), handle):
                lines.append(line.rstrip()[:350])
        if not lines:
            raise ValueError(f'Empty input: {path}')
        if key.endswith('fasta') and not lines[0].startswith('>'):
            raise ValueError(f'Expected FASTA: {path}')
        if key.endswith('gtf'):
            if any('GRCh37' in line for line in lines):
                raise ValueError(f'Incompatible annotation assembly: {path}')
        probes[key] = {'path': str(path), 'bytes': path.stat().st_size, 'first_lines': lines[:2]}
    print(json.dumps({'status': 'dry_run', 'validation_scope': 'path resolution and header probes only; '
                      'complete CRC, ID, coordinate and sequence validation occurs during preparation',
                      'inputs': probes, 'config': config, 'estimated_cache_disk_bytes': 5 * 1024 ** 3,
                      'estimated_ram_bytes': 3 * 1024 ** 3,
                      'expected_records': {'real_lncrna': config['n_transcripts'],
                                           'dinuc_shuffle': config['n_transcripts'] * config['shuffles_per_transcript'],
                                           'intergenic_matched': config['n_transcripts'],
                                           'coding_positive': config['n_transcripts'],
                                           'repeat_masked_query': config['n_transcripts']}}, indent=2))


@contextmanager
def exclusive_lock(path: Path):
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise ValueError(f'Run/cache lock exists: {path}. Check the recorded PID before removing a stale lock.') from error
    try:
        with os.fdopen(fd, 'w') as handle:
            handle.write(f'{os.getpid()}\n')
        yield
    finally:
        path.unlink(missing_ok=True)


def save_toml(path: Path, config: dict) -> None:
    # None represents an unset optional constraint and is omitted (TOML has no null).
    lines = []
    for key, value in sorted(config.items()):
        if value is None:
            continue
        if isinstance(value, bool):
            text = str(value).lower()
        else:
            text = json.dumps(value)
        lines.append(f'{key} = {text}')
    path.write_text('\n'.join(lines) + '\n')



def run(args: argparse.Namespace) -> dict:
    config = configuration(args)
    paths = {key: resolve_input(args.input_dir, name, getattr(args, key)) for key, name in INPUTS.items()}
    output = args.output_dir.resolve()
    if output.is_symlink() or args.output_dir.is_symlink():
        raise ValueError('Output directory may not be a symlink')
    if any(path.is_relative_to(output) for path in paths.values()):
        raise ValueError('Output would contain immutable inputs; choose a separate output directory')
    if args.dry_run:
        dry_run(paths, config)
        return {'status': 'dry_run'}
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.with_name('.' + output.name + '.preparing')
    cache = (args.cache_dir.resolve() if args.cache_dir else stage / '_cache')
    if args.cache_dir and (cache == output or cache.is_relative_to(output)
                           or output.is_relative_to(cache) or any(p.is_relative_to(cache) for p in paths.values())):
        raise ValueError('--cache-dir must be separate from output and must not contain input files')
    with exclusive_lock(output.with_name('.' + output.name + '.lock')):
        hashes = input_hashes(paths)
        fingerprint = {'config': config, 'inputs': hashes, 'code': code_hashes(), 'version': VERSION}
        if output.exists():
            manifest_path = output / 'provenance/run_manifest.json'
            if not manifest_path.is_file():
                raise ValueError(f'Refusing to overwrite directory not owned by this tool: {output}')
            manifest = json.loads(manifest_path.read_text())
            if manifest.get('tool') != 'FragIncRNA-null-models' or manifest.get('status') != 'complete':
                raise ValueError('Existing output lacks a completed pipeline manifest')
            if args.resume:
                if manifest.get('fingerprint') != fingerprint:
                    raise ValueError('Resume provenance mismatch: inputs, code or scientific configuration changed')
                verify_receipt(output)
                LOG.info('Completed output hashes verified; no regeneration required')
                return json.loads((output / 'qc/summary.json').read_text())
            if not args.force:
                raise ValueError(f'Output exists: {output}; use --resume to verify, --force to replace, or a new directory')
        if stage.exists():
            stage_manifest = stage / 'preparing.json'
            if not args.resume or not stage_manifest.exists() or json.loads(stage_manifest.read_text()) != fingerprint:
                raise ValueError(f'Existing incomplete output {stage}; --resume requires matching provenance')
            # Restart compatible incomplete outputs, preserving only this run's cache.
            for child in stage.iterdir():
                if child.name not in {'_cache', 'preparing.json'}:
                    shutil.rmtree(child) if child.is_dir() else child.unlink()
        else:
            stage.mkdir()
            (stage / 'preparing.json').write_text(json.dumps(fingerprint, sort_keys=True))
        cache.mkdir(parents=True, exist_ok=True)
        genome = db = None
        try:
            with exclusive_lock(cache / 'preparation.lock'):
                working_config = dict(config, threads=args.threads)
                genome, repeats, db, info = prepare(working_config, paths, hashes, cache, fingerprint['code'])
                parents, strata = select_lncrnas(pool(db, 'lncrna'), config['n_transcripts'], config['seed'],
                                               config['max_per_gene'], config['allow_partial'])
                LOG.info('Selected %d stratified lncRNA transcripts; matching coding controls', len(parents))
                coding, coding_diag = match_coding(parents, pool(db, 'coding'), working_config)
                LOG.info('Coding matches %d/%d; sampling intergenic controls', len(coding), len(parents))
                intergenic, intergenic_diag, universe = match_intergenic(
                    parents, genome, repeats, info['genes'], info['eligible_chromosomes'], working_config)
                LOG.info('Intergenic matches %d/%d', len(intergenic), len(parents))
                missing = {'coding': sorted(p['transcript_id'] for p in parents if p['transcript_id'] not in coding),
                           'intergenic': sorted(p['transcript_id'] for p in parents if p['transcript_id'] not in intergenic)}
                if any(missing.values()) and not config['allow_partial']:
                    error_path = output.with_name(output.name + '.matching_failure.json')
                    error_path.write_text(json.dumps({'fingerprint': fingerprint, 'unmatched': missing,
                                                      'coding_diagnostics': coding_diag,
                                                      'intergenic_diagnostics': intergenic_diag}, indent=2))
                    raise ValueError(f'Infeasible matching: coding={len(coding)}/{len(parents)}, '
                                     f'intergenic={len(intergenic)}/{len(parents)}. Details: {error_path}. '
                                     'No tolerances were relaxed. Explicitly change tolerances/candidate budget '
                                     'or use --allow-partial in a new run.')
                fasta_paths, summary = emit_outputs(stage, db, genome, repeats, parents, coding, intergenic,
                                                    working_config)
                LOG.info('Independently reparsing output FASTAs, metadata and gzip members')
                verification = verify_outputs(stage, fasta_paths, config['k_values'])
                summary.update(requested_transcripts=config['n_transcripts'], selected_transcripts=len(parents),
                               seed=config['seed'], partial=any(missing.values()) or len(parents) != config['n_transcripts'],
                               unmatched=missing, fasta_paths=fasta_paths, verification=verification,
                               intergenic_universe=universe, input_validation={k: v for k, v in info.items() if k != 'genes'},
                               matching_tolerances={k: config[k] for k in ['length_tolerance', 'gc_tolerance', 'repeat_tolerance']},
                               #limitations=LIMITATIONS
                )
                write_table(stage / 'metadata/matching_diagnostics.tsv', coding_diag + intergenic_diag)
                write_table(stage / 'metadata/exclusion_counts.tsv',
                            [dict(reason=k, count=v) for k, v in sorted(info['exclusions'].items())])
                write_table(stage / 'qc/strata.tsv', strata)
                (stage / 'qc/summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
                (stage / 'qc/report.md').write_text(
                    '# Null-model preparation QC\n\n'
                    + '\n'.join(f'- {kind}: {count}' for kind, count in summary['counts'].items())
                    + f'\n\nSeed: {config["seed"]}. Complete output CRC/sequence/window verification passed.\n\n')
                    #+ '\n'.join('- ' + text for text in LIMITATIONS) + '\n')
                (stage / 'provenance').mkdir(exist_ok=True)
                write_table(stage / 'provenance/inputs.sha256.tsv', [dict(input=key, path=str(paths[key]),
                            sha256=hashes[key], size_bytes=paths[key].stat().st_size) for key in sorted(paths)])
                save_toml(stage / 'provenance/resolved_config.toml',
                          dict(config, threads=args.threads, input_dir=str(args.input_dir.resolve()),
                               **{key: str(path) for key, path in paths.items()}))
                manifest = {'tool': 'FragIncRNA-null-models', 'status': 'complete', 'fingerprint': fingerprint,
                            'python_version': sys.version, 'preparation_threads': args.threads,
                            'counts': summary['counts'], 'fasta_paths': fasta_paths}
                (stage / 'provenance/run_manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
                db.close()
                db = None
                genome.close()
                genome = None
            if not args.cache_dir:
                shutil.rmtree(cache)
            (stage / 'preparing.json').unlink()
            write_output_hashes(stage)
            backup = output.with_name('.' + output.name + '.previous')
            if backup.exists():
                raise ValueError(f'Prior publication backup exists: {backup}; inspect before proceeding')
            if output.exists():
                output.rename(backup)
            try:
                stage.rename(output)
            except BaseException:
                if backup.exists():
                    backup.rename(output)
                raise
            if backup.exists():
                shutil.rmtree(backup)
            LOG.info('Published validated dataset: %s', output)
            return summary
        except BaseException:
            if db is not None:
                db.close()
            if genome is not None:
                genome.close()
            # Keep explicit caches; remove only the unpublished directory owned by this run.
            if stage.exists():
                shutil.rmtree(stage)
            raise


def main() -> int:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
    args = parser().parse_args()
    try:
        summary = run(args)
    except (ValueError, OSError, KeyError, csv.Error) as error:
        LOG.error('%s', error)
        return 2
    except KeyboardInterrupt:
        LOG.error('Interrupted; no incomplete dataset published')
        return 130
    if not args.dry_run:
        print(json.dumps({'output_dir': str(args.output_dir.resolve()), 'counts': summary['counts'],
                          'partial': summary['partial'], 'verification': summary['verification']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
