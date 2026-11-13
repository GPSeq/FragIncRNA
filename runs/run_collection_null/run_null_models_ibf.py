#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import sys

import hashlib
import json
import logging
import math
import shutil
import signal
import tarfile
import threading
import time
import tomllib
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from null_ibf_archives import analyze_results, create_archive, verify_completion
from null_ibf_data import (LABELS, campaign_lock, compatible_json, discover_datasets, discover_genomes,
                           no_symlinks, regular_file, remove_private, safe_tree, write_tsv)
from run_fraglncrna_grid import iso_now, render_config, run_executable, sha256_file, slug_for, write_json

FORMAT = 'fraglncrna-null-ibf-campaign-v1'
OUTPUT_OPTIONS = ('positions', 'calibration_bins', 'fragments', 'repeat_classes', 'reference_hit_classes',
                  'discover_loci', 'expected_loci_output', 'window_output', 'kmer_fasta_output')
INPUT_OPTIONS = ('threshold_table', 'expected_loci_bed', 'exon_blocks', 'window_query_ids')
STATUS_FIELDS = ['genome', 'dataset', 'k', 'hash_functions', 'status', 'elapsed_seconds', 'archive_path',
                 'query_count_if_known', 'reference_sha256', 'queries_sha256', 'job_signature', 'message']
SUMMARY_FIELDS = ['genome', 'dataset', 'k', 'hash_functions', 'aggregation_status', 'metric', 'query_count',
                  'numeric_recovery_count', 'mean_recovery', 'median_recovery', 'pass_basic_count',
                  'pass_basic_fraction', 'pass_strict_count', 'pass_strict_fraction', 'at_least_basic_count',
                  'at_least_basic_fraction', 'valid_windows', 'zero_valid_queries', 'shards_built', 'shards_reused']


def signature(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def argument_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description='Run named null-model query sets through shared sharded IBFs. '
                                'Candidate-local exact confirmation is retained; no global-exact jobs.',
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--genomes-dir', type=Path, required=True)
    p.add_argument('--genomes', nargs='+', help='Each filename/quoted glob must match exactly one reference, or use absolute paths')
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--datasets-manifest', type=Path)
    source.add_argument('--datasets-dir', type=Path)
    p.add_argument('--datasets', default=','.join(LABELS), help='Comma-separated logical dataset subset')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--threads', type=int, required=True, help='Threads per C++ process', default=25)
    p.add_argument('--parallel-genomes', type=int, default=4)
    p.add_argument('--k-values', default='15,20,32')
    p.add_argument('--hash-functions', nargs='+', type=int, default=[3])
    p.add_argument('--binary', type=Path, default=Path('build/lncrna_mers_v2'))
    p.add_argument('--template', type=Path, default=Path(__file__).with_name('default_null_ibf.toml'))
    p.add_argument('--dry-run', action='store_true', help='Validate and print the grid; no writes or executable calls')
    p.add_argument('--max-runs', type=int, default=0, help='First N new dataset query executions in plan order; 0 means unlimited')
    p.add_argument('--continue-on-error', action='store_true')
    p.add_argument('--keep-work', action='store_true', help='Keep successful dataset configs/results; completed group caches are still removed')
    p.add_argument('--min-free-disk-gb', type=float, default=10,
                   help='Minimum free GiB at preflight and before each query/archive publication')
    return p


def check_template(source: str, template: Path) -> tuple[dict, dict]:
    parsed = tomllib.loads(source)
    if set(parsed) != {'v2'}:
        raise ValueError('Null IBF template must contain only [v2]')
    config = parsed['v2']
    required = dict(mode='ibf', scoring_method='per_bin_verified', per_bin_exact_verify=True,
                    fdr_control_per_bin=True, ibf_layout='sharded', ibf_shard_storage='disk',
                    ibf_shard_retention='keep', ibf_shard_reuse=True, ibf_resume=True,
                    ambiguity='skip', repeat='include')
    for key, value in required.items():
        if config.get(key) != value:
            raise ValueError(f'IBF-only verified template requires {key}={value!r}')
    for key in ('ibf_build_only', 'ibf_dry_run'):
        if config.get(key, False):
            raise ValueError(f'{key}=true is not a query campaign')
    for key in ('reference_dir', 'output_dir', 'save_index', 'load_index', 'reference_repeat_bed',
                'query_repeat_bed', 'ibf_local_cache_dir'):
        if config.get(key):
            raise ValueError(f'{key} is not supported in this single-reference, query-N-masked campaign; '
                             'use a separate reviewed reference-mask/index experiment')
    if config.get('calibration_bins'):
        raise ValueError('Calibration reports are not part of this IBF-only campaign')
    auxiliary = {}
    for key in INPUT_OPTIONS:
        if config.get(key):
            path = Path(config[key]).expanduser()
            path = regular_file(path if path.is_absolute() else template.parent / path)
            config[key] = str(path)
            auxiliary[key] = {'path': str(path), 'sha256': sha256_file(path)}
    # Catch duplicate TOML keys/incompatible renderer without creating anything.
    render_config(source, {})
    return config, auxiliary


def nearest_existing(path: Path) -> Path:
    while not path.exists():
        path = path.parent
    return path


def disk_check(root: Path, minimum: float) -> float:
    available = shutil.disk_usage(nearest_existing(root)).free / 1024 ** 3
    if available < minimum:
        raise ValueError(f'Insufficient free disk: {available:.2f} GiB available, '
                         f'--min-free-disk-gb={minimum}. Private work/caches are preserved.')
    return available


def code_identity() -> dict:
    directory = Path(__file__).parent
    return {name: sha256_file(directory / name) for name in (
        'run_null_models_ibf.py', 'null_ibf_data.py', 'null_ibf_archives.py', 'run_fraglncrna_grid.py')}


def preflight(args: argparse.Namespace) -> dict:
    if (args.threads < 1 or args.parallel_genomes < 1 or args.max_runs < 0
            or not math.isfinite(args.min_free_disk_gb) or args.min_free_disk_gb < 0):
        raise ValueError('Threads/concurrency must be positive; max-runs/disk reserve must be nonnegative')
    ks = [int(k.strip()) for k in args.k_values.split(',')]
    if len(set(ks)) != len(ks) or not ks or min(ks) < 1 or max(ks) > 32:
        raise ValueError('--k-values must be unique integers in [1,32]')
    ks.sort()
    hashes = args.hash_functions
    if len(set(hashes)) != len(hashes) or not hashes or min(hashes) < 1 or max(hashes) > 5:
        raise ValueError('--hash-functions must contain unique counts in [1,5]')
    hashes = sorted(hashes)
    requested = args.datasets.split(',')
    if len(set(requested)) != len(requested) or set(requested) - set(LABELS):
        raise ValueError(f'--datasets must contain distinct labels from {LABELS}')
    selected = [label for label in LABELS if label in requested]
    binary = regular_file(args.binary)
    if not os.access(binary, os.X_OK):
        raise ValueError(f'Binary is not executable: {binary}')
    template = regular_file(args.template)
    source = template.read_text()
    settings, auxiliary = check_template(source, template)
    genomes = discover_genomes(args.genomes_dir, args.genomes)
    datasets, metadata, identifiers = discover_datasets(args.datasets_manifest, args.datasets_dir, selected, ks)
    root = no_symlinks(args.output_dir)
    inputs = [binary, template, *genomes, *[Path(v['path']) for v in datasets.values()],
              *[Path(v['path']) for v in metadata.values()], *[Path(v['path']) for v in auxiliary.values()]]
    if args.datasets_manifest:
        inputs.append(regular_file(args.datasets_manifest))
    genome_dir = no_symlinks(args.genomes_dir)
    if root.is_relative_to(genome_dir) or any(path.is_relative_to(root) for path in inputs):
        raise ValueError('Output root must be separate from genome inputs and cannot contain any input file')
    if root.exists():
        no_symlinks(root / 'runner.log')
        safe_tree(root / '_work', root)
        safe_tree(root / 'archives', root)
        safe_tree(root / 'manifests', root)
        safe_tree(root / 'summaries', root)
        if not (root / 'manifests/campaign.json').is_file() and any(root.iterdir()):
            raise ValueError(f'Refusing to adopt an unrelated nonempty output root: {root}')
    free = disk_check(root, args.min_free_disk_gb)
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    workers = min(args.parallel_genomes, len(genomes))
    if args.threads * workers > cpus:
        print(f'[WARN] {args.threads * workers} requested CPU threads exceed {cpus} available CPUs', file=sys.stderr)
    if workers > 1:
        print(f'[WARN] {workers} independent genome caches/builds can run concurrently; '
              f'{free:.1f} GiB free is a reserve check, not an index-size or RAM estimate.', file=sys.stderr)
    if not metadata:
        print('[WARN] No source metadata supplied; FASTA IDs remain unchanged, but parent/replicate linkage '
              'must be supplied separately for paired analysis.', file=sys.stderr)
    if 'repeat_masked_lncrna' in selected:
        masked = datasets['repeat_masked_lncrna']
        if not masked['n_bases']:
            print('[WARN] Repeat-masked dataset contains no N bases; sensitivity may be identical to baseline.', file=sys.stderr)
        if masked['lowercase_acgt_bases']:
            print('[WARN] Lowercase ACGT is counted as ordinary sequence; only ambiguous/N windows are skipped.', file=sys.stderr)
        print('[MASK QC] repeat-masked valid-window denominators: ' + json.dumps(masked['valid_windows']), flush=True)
    registry = {'format': FORMAT, 'datasets': {label: {k: v for k, v in value.items()
                 if k not in ('valid_windows', 'original_windows', 'zero_valid_queries')}
                 for label, value in datasets.items()}, 'metadata': metadata}
    campaign = {'format': FORMAT, 'binary': {'path': str(binary), 'sha256': sha256_file(binary)},
                'template': {'path': str(template), 'sha256': sha256_file(template)},
                'settings': settings, 'auxiliary_inputs': auxiliary, 'threads': args.threads,
                'dataset_registry_signature': signature(registry), 'runner_code': code_identity()}
    # Dry runs validate identities of existing roots without writing them.
    for path, value in [(root / 'manifests/campaign.json', campaign),
                        (root / 'manifests/dataset_registry.json', registry)]:
        if path.exists() and json.loads(path.read_text()) != value:
            raise ValueError(f'Provenance mismatch: {path}; choose a new output root')
    return dict(root=root, binary=binary, template=template, source=source, settings=settings,
                auxiliary=auxiliary, genomes=genomes, datasets=datasets, metadata=metadata,
                identifiers=identifiers, ks=ks, hashes=hashes, selected=selected,
                registry=registry, campaign=campaign)


def generate_config(plan: dict, args: argparse.Namespace, job: dict) -> str:
    results = job['work'] / 'results'
    cache = job['group_work'] / 'cache'
    overrides = {
        'reference': str(job['genome']), 'queries': plan['datasets'][job['dataset']]['path'],
        'output': str(results / 'lncrna.tsv'), 'metrics': str(results / 'metrics.json'),
        'threads': args.threads, 'k': job['k'], 'hash_functions': job['h'],
        'ibf_shard_dir': str(cache / 'ibf'), 'reference_cache': str(cache / 'reference.bases'),
        'log_file': str(results / 'tool.log'), 'ibf_tuning_report': str(results / 'ibf_tuning.json'),
        'per_bin_output': str(results / 'lncrna.tsv.per_bin.tsv'),
        'ibf_build_only': False, 'ibf_dry_run': False, 'save_index': '', 'load_index': '',
        'ibf_shard_retention': 'keep', 'ibf_shard_reuse': True, 'ibf_resume': True,
    }
    for key in INPUT_OPTIONS:
        if key in plan['auxiliary']:
            overrides[key] = plan['auxiliary'][key]['path']
    for key in OUTPUT_OPTIONS:
        value = plan['settings'].get(key)
        if value:
            suffix = ''.join(Path(value).suffixes[-2:]) or '.tsv'
            overrides[key] = str(results / (key + suffix))
    rendered = render_config(plan['source'], overrides)
    return rendered


def build_jobs(plan: dict, args: argparse.Namespace, hash_references: bool) -> list[dict]:
    jobs = []
    for genome in plan['genomes']:
        if hash_references:
            print(f'[HASH REFERENCE] {genome.name}', flush=True)
        ref_hash = sha256_file(genome) if hash_references else ''
        for h in plan['hashes']:
            for k in plan['ks']:
                group_key = f'{slug_for(genome)}/hash{h}/k{k}'
                group_work = plan['root'] / '_work' / group_key
                # Conservative superset of the C++ reference/index identity. The C++
                # manifest independently checks semantic SHA, layout, compiler and ABI.
                group_data = {'format': FORMAT, 'genome': str(genome), 'reference_sha256': ref_hash,
                              'k': k, 'hash_functions': h, 'settings': plan['settings'],
                              'binary': plan['campaign']['binary'], 'auxiliary': plan['auxiliary']}
                group_signature = signature(group_data)
                for label in plan['selected']:
                    job_data = dict(group_signature=group_signature, dataset=label,
                                    queries=plan['registry']['datasets'][label],
                                    metadata=plan['metadata'], campaign_signature=signature(plan['campaign']))
                    job_signature = signature(job_data)
                    archive = plan['root'] / 'archives' / group_key / (label + '.tar.gz')
                    jobs.append(dict(key=group_key + '/' + label, group_key=group_key, group_work=group_work,
                                     group_data=dict(group_data, group_signature=group_signature),
                                     dataset=label, genome=genome, k=k, h=h,
                                     work=group_work / 'datasets' / label, archive=archive,
                                     receipt=archive.with_name(archive.name + '.complete.json'),
                                     signature=job_signature, job_data=job_data,
                                     reference_sha256=ref_hash))
    return jobs


class Campaign:
    def __init__(self, plan: dict, args: argparse.Namespace, jobs: list[dict], logger: logging.Logger):
        self.plan, self.args, self.jobs, self.log = plan, args, jobs, logger
        self.root = plan['root']
        self.stop = threading.Event()
        self.lock = threading.Lock()
        self.completed: dict[str, dict] = {}
        self.rows = {}
        self.failure_count = 0
        self.attempted = 0
        self.interrupted = False
        self.allowed: set[str] = set()
        self.archive_stamps = {}

    def stamp(self, path: Path) -> tuple:
        no_symlinks(path)
        info = path.stat()
        return info.st_size, info.st_mtime_ns, info.st_ino

    def save_summaries(self) -> None:
        # Caller holds self.lock, or no workers have started yet.
        ordered = [self.rows[job['key']] for job in self.jobs]
        write_tsv(self.root / 'summaries/job_status.tsv', ordered, STATUS_FIELDS)
        recovery = []
        for job in self.jobs:
            if job['key'] in self.completed:
                receipt = self.completed[job['key']]
                recovery.append(dict(genome=job['genome'].name, dataset=job['dataset'], k=job['k'],
                                     hash_functions=job['h'], **receipt['result_summary'], **{
                                         key: receipt['index_counts'][key] for key in ('shards_built', 'shards_reused')}))
        write_tsv(self.root / 'summaries/recovery_by_dataset.tsv', recovery, SUMMARY_FIELDS)

    def initialize(self) -> None:
        # Validate every existing artifact before any executable is launched.
        for job in self.jobs:
            safe_tree(job['group_work'], self.root)
            group_file = job['group_work'] / 'group.json'
            if job['group_work'].exists() and not group_file.exists() and any(job['group_work'].iterdir()):
                raise ValueError(f'Unmanaged nonempty group directory: {job["group_work"]}; preserve and inspect')
            if group_file.exists() and json.loads(group_file.read_text()) != job['group_data']:
                raise ValueError(f'Incompatible retained group cache: {group_file}; choose a new output root')
            receipt = verify_completion(job['archive'], job['receipt'], job['signature'], repair=True)
            if receipt:
                self.completed[job['key']] = receipt
                self.archive_stamps[job['key']] = (self.stamp(job['archive']), self.stamp(job['receipt']))
            previous_state = job['work'] / 'state.json'
            state = json.loads(previous_state.read_text()) if previous_state.exists() else {}
            self.rows[job['key']] = dict(
                genome=job['genome'].name, dataset=job['dataset'], k=job['k'], hash_functions=job['h'],
                status='complete' if receipt else state.get('status', 'pending'),
                elapsed_seconds=receipt['elapsed_seconds'] if receipt else state.get('elapsed_seconds', ''),
                archive_path=str(job['archive'].relative_to(self.root)),
                query_count_if_known=self.plan['datasets'][job['dataset']]['query_count'],
                reference_sha256=job['reference_sha256'], queries_sha256=self.plan['datasets'][job['dataset']]['sha256'],
                job_signature=job['signature'], message='verified completion' if receipt else state.get('message', ''))
        pending = [j for j in self.jobs if j['key'] not in self.completed]
        selected_pending = pending[:self.args.max_runs] if self.args.max_runs else pending
        self.allowed = {j['key'] for j in selected_pending}
        for job in pending:
            if job['key'] not in self.allowed and self.rows[job['key']]['status'] == 'pending':
                self.rows[job['key']]['status'] = 'deferred_max_runs'
        self.save_summaries()

    def cleanup_group(self, group_jobs: list[dict]) -> None:
        if not all(j['key'] in self.completed for j in group_jobs):
            return
        for job in group_jobs:
            if self.archive_stamps[job['key']] != (self.stamp(job['archive']), self.stamp(job['receipt'])):
                raise ValueError('Archive or receipt changed during campaign; cache preserved')
        group = group_jobs[0]['group_work']
        cache = group / 'cache'
        if cache.exists():
            marker = cache / 'owner.json'
            if not marker.exists() or json.loads(marker.read_text()) != {'group_signature': group_jobs[0]['group_data']['group_signature']}:
                raise ValueError(f'Cache ownership mismatch; refusing deletion: {cache}')
            safe_tree(cache, self.root)
            known = {'owner.json', 'ibf', 'reference.bases', 'index_state.json'}
            if set(p.name for p in cache.iterdir()) - known:
                raise ValueError(f'Unexpected files inside group cache; preserving for inspection: {cache}')
            index = cache / 'ibf'
            if index.exists() and any(not p.is_file() or not (
                    p.name in {'manifest.tsv', '.lock'} or re.fullmatch(r'shard-\d+\.ibf', p.name))
                    for p in index.iterdir()):
                raise ValueError(f'Unexpected files inside private IBF cache; preserving: {index}')
            remove_private(cache, self.root)
            self.log.info('[CACHE REMOVED] %s: every requested dataset has a verified archive/receipt', group_jobs[0]['group_key'])

    def execute_job(self, job: dict) -> None:
        work, group = job['work'], job['group_work']
        work.mkdir(parents=True, exist_ok=True)
        safe_tree(work, self.root)
        old_job = work / 'job.json'
        if old_job.exists() and json.loads(old_job.read_text()).get('job_signature') != job['signature']:
            raise ValueError(f'Incompatible incomplete job: {work}; choose a new output root')
        results = work / 'results'
        if results.exists():
            # Preserve each failed attempt's diagnostics while keeping the compatible group cache.
            attempt = 1
            history = work / 'failed_attempts'
            while (history / str(attempt)).exists():
                attempt += 1
            history.mkdir(exist_ok=True)
            results.rename(history / str(attempt))
        results.mkdir()
        cache = group / 'cache'
        cache.mkdir(exist_ok=True)
        compatible_json(cache / 'owner.json', {'group_signature': job['group_data']['group_signature']})
        index_state = cache / 'index_state.json'
        expected_reuse = index_state.exists() and (cache / 'ibf/manifest.tsv').exists()
        (work / 'template.toml').write_text(self.plan['source'])
        generated = generate_config(self.plan, self.args, job)
        (work / 'generated_config.toml').write_text(generated)
        write_json(work / 'job.json', dict(job['job_data'], job_signature=job['signature'],
                                         reference_sha256=job['reference_sha256'],
                                         scientific_settings=tomllib.loads(generated)['v2']))
        (work / 'README_archive.txt').write_text(
            'Archive members use relative paths listed in archive_manifest.json.\n'
            'Generated config/command paths record the original invocation and are not portable replay paths.\n'
            'Input filenames and SHA-256 identify external immutable source files; inputs and IBF shards are not embedded.\n'
            'Recovery is candidate-local verified membership, not independent global-exact conservation or statistical FDR.\n')
        command = [str(self.plan['binary']), '--config', str(work / 'generated_config.toml')]
        started = iso_now()
        write_json(work / 'execution.json', dict(command=command, cwd=str(Path.cwd()), started_utc=started,
                                                binary=self.plan['campaign']['binary'], dataset=job['dataset']))
        write_json(work / 'state.json', {'status': 'running', 'message': 'C++ query execution started'})
        with self.lock:
            self.attempted += 1
            self.rows[job['key']].update(status='running', message='C++ query execution started')
            self.save_summaries()
        self.log.info('[START] %s shared_index=%s', job['key'], cache / 'ibf')
        t0 = time.monotonic()
        try:
            disk_check(self.root, self.args.min_free_disk_gb)
            rc = run_executable(self.plan['binary'], work / 'generated_config.toml', Path.cwd(),
                                results / 'fraglncrna.log', self.stop)
            if rc:
                raise RuntimeError(f'FragIncRNA exited {rc}; work/cache preserved: {work}')
            safe_tree(results, self.root)
            summary, metrics = analyze_results(results, self.plan['identifiers'][job['dataset']], job['k'], job['h'])
            if expected_reuse and metrics['shards_built']:
                raise ValueError('Executable rebuilt a previously complete compatible shard group when queries changed; '
                                 'preserving cache/results for diagnosis')
            manifest_file = regular_file(cache / 'ibf/manifest.tsv')
            shutil.copyfile(manifest_file, work / 'reference_index_manifest.tsv')
            write_json(index_state, {'group_signature': job['group_data']['group_signature'],
                                     'shard_count': metrics['shard_count']})
            elapsed = time.monotonic() - t0
            write_json(work / 'execution.json', dict(command=command, cwd=str(Path.cwd()), started_utc=started,
                       finished_utc=iso_now(), returncode=rc, elapsed_seconds=elapsed,
                       binary=self.plan['campaign']['binary'], dataset=job['dataset']))
            disk_check(self.root, self.args.min_free_disk_gb)
            receipt = create_archive(work, job['archive'], job['signature'], elapsed, summary, metrics)
            write_json(job['receipt'], receipt)
            if json.loads(job['receipt'].read_text()) != receipt:
                raise ValueError('New completion receipt failed read-back validation')
            with self.lock:
                self.completed[job['key']] = receipt
                self.archive_stamps[job['key']] = (self.stamp(job['archive']), self.stamp(job['receipt']))
                self.rows[job['key']].update(status='complete', elapsed_seconds=receipt['elapsed_seconds'], message='verified archive')
                self.save_summaries()
            write_json(work / 'state.json', {'status': 'complete', 'elapsed_seconds': elapsed})
            self.log.info('[COMPLETE] %s built=%d reused=%d archive=%s', job['key'],
                          metrics['shards_built'], metrics['shards_reused'], job['archive'])
            if not self.args.keep_work:
                # Keep unrelated job siblings and earlier failed-attempt logs.
                remove_private(results, self.root)
                for name in ('generated_config.toml', 'template.toml', 'job.json', 'execution.json',
                             'README_archive.txt', 'state.json', 'archive_manifest.json', 'reference_index_manifest.tsv'):
                    no_symlinks(work / name).unlink(missing_ok=True)
                if not any(work.iterdir()):
                    work.rmdir()
        except Exception as error:
            elapsed = time.monotonic() - t0
            state = dict(status='interrupted' if self.stop.is_set() else 'failed', elapsed_seconds=round(elapsed, 3),
                         message=str(error))
            write_json(work / 'state.json', state)
            with self.lock:
                self.failure_count += 1
                self.rows[job['key']].update(state)
                self.save_summaries()
            self.log.error('[PRESERVED] %s: %s', job['key'], error)
            if not self.args.continue_on_error:
                self.stop.set()

    def process_genome(self, genome: Path) -> None:
        jobs = [j for j in self.jobs if j['genome'] == genome]
        groups = dict.fromkeys(j['group_key'] for j in jobs)
        try:
            for group_key in groups:
                if self.stop.is_set():
                    break
                group_jobs = [j for j in jobs if j['group_key'] == group_key]
                group = group_jobs[0]['group_work']
                if all(j['key'] in self.completed for j in group_jobs):
                    self.log.info('[SKIP GROUP] %s: all requested archives verified', group_key)
                    self.cleanup_group(group_jobs)
                    continue
                if not any(j['key'] in self.allowed for j in group_jobs):
                    continue
                group.mkdir(parents=True, exist_ok=True)
                compatible_json(group / 'group.json', group_jobs[0]['group_data'])
                for job in group_jobs:
                    if self.stop.is_set():
                        break
                    if job['key'] in self.completed:
                        self.log.info('[SKIP DATASET] %s: verified completion', job['key'])
                    elif job['key'] in self.allowed:
                        self.execute_job(job)
                self.cleanup_group(group_jobs)
        except Exception as error:
            with self.lock:
                self.failure_count += 1
                for job in group_jobs:
                    if job['key'] not in self.completed and self.rows[job['key']]['status'] not in {'failed', 'interrupted'}:
                        self.rows[job['key']].update(status='group_error', message=str(error))
                self.save_summaries()
            if group.exists():
                write_json(group / 'error.json', {'status': 'group_error', 'message': str(error), 'time': iso_now()})
            self.log.error('[GROUP PRESERVED] %s: %s', genome.name, error)
            if not self.args.continue_on_error:
                self.stop.set()

    def run(self) -> int:
        def interrupted(signum, frame):
            self.interrupted = True
            self.stop.set()
            self.log.warning('[INTERRUPTED] Signal %s; stopping child process groups and preserving work/cache', signum)
        old_handlers = {s: signal.signal(s, interrupted) for s in (signal.SIGINT, signal.SIGTERM)}
        try:
            with ThreadPoolExecutor(max_workers=self.args.parallel_genomes) as executor:
                futures = [executor.submit(self.process_genome, genome) for genome in self.plan['genomes']]
                for future in futures:
                    future.result()
        finally:
            for s, previous in old_handlers.items():
                signal.signal(s, previous)
            with self.lock:
                self.save_summaries()
        self.log.info('[SUMMARY] planned=%d attempted=%d completed=%d failed=%d',
                      len(self.jobs), self.attempted, len(self.completed), self.failure_count)
        return 130 if self.interrupted else 1 if self.failure_count else 0


def main() -> int:
    args = argument_parser().parse_args()
    try:
        plan = preflight(args)
        jobs = build_jobs(plan, args, hash_references=not args.dry_run)
        group_count = len({j['group_key'] for j in jobs})
        print(f'[PLAN] genomes={len(plan["genomes"])} datasets={len(plan["selected"])} '
              f'k={plan["ks"]} hashes={plan["hashes"]} query_executions={len(jobs)} index_groups={group_count} '
              f'threads_per_process={args.threads} parallel_genomes={args.parallel_genomes}', flush=True)
        if args.dry_run:
            for label in plan['selected']:
                print(f'  {label}: {plan["datasets"][label]["query_count"]} queries {plan["datasets"][label]["path"]}')
            print('[DRY-RUN] No output files written; no executable invoked; reference headers probed, not fully hashed/scanned.')
            return 0
        root = plan['root']
        root.mkdir(parents=True, exist_ok=True)
        with campaign_lock(root):
            compatible_json(root / 'manifests/campaign.json', plan['campaign'])
            compatible_json(root / 'manifests/dataset_registry.json', plan['registry'])
            write_json(root / 'manifests/dataset_window_counts.json', {
                label: {key: info[key] for key in ('valid_windows', 'original_windows', 'zero_valid_queries')}
                for label, info in plan['datasets'].items()})
            write_tsv(root / 'manifests/planned_jobs.tsv', [dict(genome=j['genome'].name, dataset=j['dataset'],
                       k=j['k'], hash_functions=j['h'], group=j['group_key'], job_signature=j['signature'],
                       archive_path=str(j['archive'].relative_to(root))) for j in jobs],
                      ['genome', 'dataset', 'k', 'hash_functions', 'group', 'job_signature', 'archive_path'])
            logger = logging.getLogger('null_ibf')
            logger.setLevel(logging.INFO)
            logger.propagate = False
            formatter = logging.Formatter('%(asctime)sZ %(levelname)s %(message)s')
            formatter.converter = time.gmtime
            handlers = [logging.FileHandler(root / 'runner.log'), logging.StreamHandler(sys.stdout)]
            for handler in handlers:
                handler.setFormatter(formatter)
                logger.addHandler(handler)
            try:
                campaign = Campaign(plan, args, jobs, logger)
                campaign.initialize()
                return campaign.run()
            finally:
                for handler in handlers:
                    logger.removeHandler(handler)
                    handler.close()
    except (OSError, ValueError, RuntimeError, KeyError, csv.Error, tarfile.TarError, EOFError) as error:
        print(f'[ERROR] {error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print('[INTERRUPTED] No new completion claimed; existing archives and private work preserved', file=sys.stderr)
        sys.exit(130)
