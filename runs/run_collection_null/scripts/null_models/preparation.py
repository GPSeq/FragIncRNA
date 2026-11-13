"""Reproducible population preparation, stratified selection and constrained matching."""
from __future__ import annotations

import bisect
import json
import logging
import math
import random
import sqlite3
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path

import numpy as np

from biology import (Transcript, composition, coverage, project_repeats, stable_seed,
                     subtract, valid_windows)
from inputs import Genome, fasta_records, read_annotation, read_repeats, sha256

LOG = logging.getLogger('null_models')
CACHE_VERSION = 2


def prepare(config: dict, paths: dict[str, Path], hashes: dict, cache: Path,
            code_hashes: dict) -> tuple[Genome, object, sqlite3.Connection, dict]:
    genome = Genome(paths['genome_fasta'], cache, hashes['genome_fasta'])
    eligible = {genome.normalize(c) for c in config['chromosomes']}
    if None in eligible:
        raise ValueError('Requested eligible chromosomes missing from reference')
    LOG.info('Reading RepeatMasker intervals')
    repeats, repeat_info = read_repeats(paths['repeats_bed'], genome, eligible)
    fingerprint = {'version': CACHE_VERSION, 'inputs': hashes, 'code': code_hashes,
                   'chromosomes': sorted(eligible), 'gencode_version': config['gencode_version'],
                   'lncrna_transcript_types': config['lncrna_transcript_types'],
                   'missing_coding_annotation': config['missing_coding_annotation'],
                   'max_splice_mismatch': config['max_splice_mismatch']}
    receipt = cache / 'prepared.json'
    db_path = cache / 'transcripts.sqlite3'
    if receipt.exists():
        cached = json.loads(receipt.read_text())
        if cached['fingerprint'] != fingerprint:
            raise ValueError(f'Incompatible preparation cache: {cache}; select a new --cache-dir')
        if not db_path.exists() or sha256(db_path) != cached.get('database_sha256'):
            raise ValueError(f'Missing transcript cache {db_path}')
        db = sqlite3.connect(db_path)
        if db.execute('PRAGMA quick_check').fetchone()[0] != 'ok':
            raise ValueError('Damaged transcript cache')
        LOG.info('Reusing verified preparation cache')
        return genome, repeats, db, cached['info']
    db_path.unlink(missing_ok=True)
    db = sqlite3.connect(db_path)
    db.execute('PRAGMA journal_mode=OFF')
    db.execute('PRAGMA synchronous=OFF')
    db.execute('CREATE TABLE records (tid TEXT PRIMARY KEY, kind TEXT, header TEXT, sequence TEXT)')
    db.execute('CREATE TABLE features (tid TEXT PRIMARY KEY, kind TEXT, data TEXT, annotation TEXT)')
    db.execute('CREATE TABLE ineligible (tid TEXT PRIMARY KEY, kind TEXT, reason TEXT)')
    ids: dict[str, set[str]] = {}
    for kind, key in [('lncrna', 'lncrna_fasta'), ('coding', 'coding_fasta')]:
        LOG.info('Streaming %s FASTA into disk-backed transcript store', kind)
        ids[kind] = set()
        for tid, header, sequence in fasta_records(paths[key]):
            if tid in ids.get('lncrna', set()) and kind == 'coding':
                raise ValueError(f'Transcript {tid} occurs in both coding and lncRNA FASTAs')
            ids[kind].add(tid)
            db.execute('INSERT INTO records VALUES (?,?,?,?)', (tid, kind, header, sequence))
        db.commit()
    LOG.info('Parsing comprehensive GENCODE annotation and all-gene exclusion intervals')
    full, genes, full_info = read_annotation(paths['annotation_gtf'], ids['lncrna'] | ids['coding'],
                                           genome, eligible, config['gencode_version'], True)
    LOG.info('Parsing lncRNA exon annotation')
    lncs, _, lnc_info = read_annotation(paths['lncrna_gtf'], ids['lncrna'], genome, eligible,
                                      config['gencode_version'])
    exclusions: Counter = Counter()
    biotypes: Counter = Counter()
    exclusions['reference_excluded_contigs'] = len(set(genome.chromosomes) - eligible)
    exclusions['repeat_excluded_rows'] = sum(v for k, v in repeat_info['counts'].items()
                                              if k.startswith('excluded'))
    for kind, annotation, info in [('lncrna', lncs, lnc_info), ('coding', full, full_info)]:
        absent = ids[kind] - annotation.keys() - set(info['excluded_transcripts'])
        if absent and (kind != 'coding' or config['missing_coding_annotation'] == 'error'):
            raise ValueError(f'{kind}: {len(absent)} FASTA IDs missing from compatible annotation; '
                             f'examples={sorted(absent)[:8]} (check versioned IDs/GTF release)')
        if absent:
            LOG.warning('%s: excluding %d FASTA IDs without supplied comprehensive annotation; '
                        'all IDs will be recorded; no coordinates inferred', kind, len(absent))
            exclusions[f'{kind}_missing_comprehensive_annotation'] = len(absent)
        exclusions[f'{kind}_excluded_chromosome'] = len(ids[kind] & set(info['excluded_transcripts']))
        LOG.info('Computing %s features and validating every eligible spliced sequence', kind)
        cursor = db.execute('SELECT tid,header,sequence FROM records WHERE kind=? ORDER BY tid', (kind,))
        eligible_count = 0
        for tid, header, sequence in cursor:
            if tid not in annotation:
                reason = 'missing_comprehensive_annotation' if tid in absent else 'excluded_chromosome'
                db.execute('INSERT INTO ineligible VALUES (?,?,?)', (tid, kind, reason))
                continue
            transcript = annotation[tid]
            expected = 'lncRNA' if kind == 'lncrna' else 'protein_coding'
            biotypes[f'{kind}:{transcript.gene_biotype}:{transcript.biotype}'] += 1
            if transcript.gene_biotype != expected:
                reason = f'unaccepted_gene_biotype:{transcript.gene_biotype or "missing"}'
                db.execute('INSERT INTO ineligible VALUES (?,?,?)', (tid, kind, reason))
                exclusions[f'{kind}_{reason}'] += 1
                continue
            accepted_types = config['lncrna_transcript_types'] if kind == 'lncrna' else ['protein_coding']
            if transcript.biotype not in accepted_types:
                reason = f'unaccepted_transcript_biotype:{transcript.biotype}'
                db.execute('INSERT INTO ineligible VALUES (?,?,?)', (tid, kind, reason))
                exclusions[f'{kind}_{reason}'] += 1
                continue
            fields = header.split()[0].split('|')
            if len(fields) > 1 and fields[1] != transcript.gene_id:
                raise ValueError(f'{tid}: versioned FASTA/GTF gene IDs disagree')
            if kind == 'lncrna':
                if tid not in full or asdict(transcript) != asdict(full[tid]):
                    raise ValueError(f'{tid}: lncRNA and comprehensive annotations disagree')
            spliced = genome.splice(transcript)
            if len(spliced) != len(sequence):
                raise ValueError(f'{tid}: exon sum {len(spliced)} != FASTA length {len(sequence)}')
            mismatches = sum(a != b for a, b in zip(spliced, sequence.upper()))
            mismatch_fraction = mismatches / len(sequence)
            if mismatch_fraction > config['max_splice_mismatch']:
                raise ValueError(f'{tid}: exon-spliced GRCh38 mismatch {mismatches}/{len(sequence)}; '
                                 'check assembly/annotation/FASTA compatibility')
            if not valid_windows(sequence, 32):
                exclusions[f'{kind}_no_valid_32mer'] += 1
                db.execute('INSERT INTO ineligible VALUES (?,?,?)', (tid, kind, 'no_valid_32mer'))
                continue
            if len(sequence) < 200:
                exclusions[f'{kind}_short_under_200_flagged_not_excluded'] += 1
            gc, n, ambiguous = composition(sequence)
            projected = project_repeats(transcript, repeats)
            repeat_fraction = coverage([i for values in projected.values() for i in values]) / len(sequence)
            feature = dict(transcript_id=tid, gene_id=transcript.gene_id,
                           gene_symbol=transcript.gene_symbol, source_biotype=transcript.biotype,
                           source_gene_biotype=transcript.gene_biotype,
                           length=len(sequence), gc_fraction=gc, n_fraction=n,
                           ambiguous_fraction=ambiguous, repeat_fraction=repeat_fraction,
                           exon_count=len(transcript.exons), chromosome=transcript.chrom,
                           strand=transcript.strand, splice_mismatches=mismatches,
                           quality_flags=';'.join(flag for flag in [
                               'short_under_200' if len(sequence) < 200 else '',
                               'lncRNA_gene_noncanonical_transcript_type'
                               if kind == 'lncrna' and transcript.biotype != 'lncRNA' else ''] if flag),
                           valid_windows_32=valid_windows(sequence, 32))
            db.execute('INSERT INTO features VALUES (?,?,?,?)',
                       (tid, kind, json.dumps(feature, sort_keys=True), json.dumps(asdict(transcript))))
            eligible_count += 1
            if eligible_count % 25000 == 0:
                LOG.info('%s: validated %d transcripts', kind, eligible_count)
                db.commit()
        db.commit()
        exclusions[f'{kind}_eligible'] = eligible_count
    def brief(info: dict) -> dict:
        return {k: v for k, v in info.items() if k != 'excluded_transcripts'}
    result = {'genes': genes, 'annotation': brief(full_info), 'lncrna_annotation': brief(lnc_info),
              'repeats': repeat_info, 'exclusions': dict(exclusions),
              'source_biotype_counts': dict(biotypes),
              'eligible_chromosomes': sorted(eligible),
              'excluded_reference_contigs': sorted(set(genome.chromosomes) - eligible),
              'chromosome_mapping': genome.mapping}
    receipt_temp = cache / 'prepared.json.partial'
    receipt_temp.write_text(json.dumps({'fingerprint': fingerprint, 'info': result,
                                        'database_sha256': sha256(db_path)}, sort_keys=True))
    receipt_temp.replace(receipt)
    LOG.info('Preparation cache complete')
    return genome, repeats, db, result


def pool(db: sqlite3.Connection, kind: str) -> list[dict]:
    return [json.loads(row[0]) for row in db.execute(
        'SELECT data FROM features WHERE kind=? ORDER BY tid', (kind,))]


def fetch_record(db: sqlite3.Connection, tid: str) -> tuple[str, str, Transcript]:
    header, sequence, annotation = db.execute(
        'SELECT r.header,r.sequence,f.annotation FROM records r JOIN features f USING(tid) WHERE tid=?',
        (tid,)).fetchone()
    return header, sequence, Transcript(**json.loads(annotation))


def stratum(feature: dict) -> str:
    axes = [('L', feature['length'], [200, 500, 1000, 2000, 4000, 8000, 16000]),
            ('G', feature['gc_fraction'], [.35, .45, .55, .65]),
            ('R', feature['repeat_fraction'], [.01, .25, .5, .75]),
            ('E', feature['exon_count'], [1, 2, 4, 8])]
    return '_'.join(f'{label}{bisect.bisect_left(edges, value)}' for label, value, edges in axes)


def select_lncrnas(features: list[dict], n: int, seed: int, max_per_gene: int,
                   allow_partial: bool) -> tuple[list[dict], list[dict]]:
    by_gene: dict[str, list[dict]] = defaultdict(list)
    for feature in features:
        by_gene[feature['gene_id']].append(feature)
    representatives = []
    for gene in sorted(by_gene):
        group = sorted(by_gene[gene], key=lambda f: (stable_seed(seed, 'gene', f['transcript_id']),
                                                     f['transcript_id']))
        representatives.extend(group[:max_per_gene] if max_per_gene else group)
    if len(representatives) < n and not allow_partial:
        raise ValueError(f'Only {len(representatives)} eligible gene-limited lncRNAs for requested {n}')
    n = min(n, len(representatives))
    if not n:
        raise ValueError('No eligible lncRNA transcripts')
    strata: dict[str, list[dict]] = defaultdict(list)
    for feature in representatives:
        strata[stratum(feature)].append(feature)
    total = len(representatives)
    allocation = {key: n * len(values) // total for key, values in strata.items()}
    remainder_order = sorted(strata, key=lambda s: (-(n * len(strata[s]) % total),
                                                   stable_seed(seed, 'allocation', s), s))
    for key in remainder_order[:n - sum(allocation.values())]:
        allocation[key] += 1
    selected = []
    diagnostics = []
    original_counts = Counter(stratum(f) for f in features)
    for key in sorted(strata):
        ranked = sorted(strata[key], key=lambda f: (stable_seed(seed, 'stratum', f['transcript_id']),
                                                    f['transcript_id']))
        selected.extend(dict(f, stratum=key) for f in ranked[:allocation[key]])
        diagnostics.append(dict(stratum=key, source_transcripts=original_counts[key],
                                gene_limited_population=len(ranked), selected=allocation[key]))
    return sorted(selected, key=lambda f: f['transcript_id']), diagnostics


def match_coding(parents: list[dict], candidates: list[dict], config: dict
                 ) -> tuple[dict[str, dict], list[dict]]:
    if not candidates:
        raise ValueError('No eligible protein-coding candidates')
    ordered = sorted(candidates, key=lambda f: f['transcript_id'])
    lengths = np.array([f['length'] for f in ordered], dtype=float)
    gc = np.array([f['gc_fraction'] for f in ordered])
    repeats = np.array([f['repeat_fraction'] for f in ordered])
    exons = np.log1p([f['exon_count'] for f in ordered])
    genes = [f['gene_id'] for f in ordered]
    feature_matrix = np.column_stack((np.log(lengths), gc, repeats, exons))
    scales = np.std(feature_matrix, axis=0)
    scales[scales == 0] = 1
    def feasible(parent: dict) -> np.ndarray:
        result = (np.abs(lengths - parent['length']) / parent['length'] <= config['length_tolerance'])
        result &= np.abs(gc - parent['gc_fraction']) <= config['gc_tolerance']
        if config['repeat_tolerance'] is not None:
            result &= np.abs(repeats - parent['repeat_fraction']) <= config['repeat_tolerance']
        return np.flatnonzero(result)
    # Process rare feasible sets first; compute each set again to bound memory.
    priority = sorted(parents, key=lambda p: (len(feasible(p)), -p['length'], p['transcript_id']))
    used_genes: set[str] = set()
    used_tids: set[str] = set()
    matches = {}
    diagnostics = []
    for parent in priority:
        indexes = feasible(parent)
        indexes = np.array([i for i in indexes if ordered[i]['transcript_id'] not in used_tids
                           and (config['coding_isoforms'] or genes[i] not in used_genes)], dtype=int)
        if not len(indexes):
            diagnostics.append(dict(control_type='coding_positive', parent_query_id=parent['transcript_id'],
                                    status='unmatched', reason='no_unused_candidate_within_tolerances'))
            continue
        target = np.array([math.log(parent['length']), parent['gc_fraction'],
                           parent['repeat_fraction'], math.log1p(parent['exon_count'])])
        cost = np.sum(((feature_matrix[indexes] - target) / scales) ** 2, axis=1)
        cost += np.array([0 if ordered[i]['chromosome'] == parent['chromosome'] else .1
                          for i in indexes])
        best = int(np.argmin(cost))  # sorted transcript IDs break exact cost ties.
        match = dict(ordered[indexes[best]], matching_cost=float(cost[best]))
        matches[parent['transcript_id']] = match
        used_tids.add(match['transcript_id'])
        used_genes.add(match['gene_id'])
        diagnostics.append(dict(control_type='coding_positive', parent_query_id=parent['transcript_id'],
                                status='matched', reason='within_tolerances'))
    return matches, diagnostics


def intergenic_universe(genome: Genome, genes: dict, eligible: list[str], flank: int) -> dict:
    universe = {}
    for chrom in eligible:
        length = genome.length(chrom)
        buffered = [(max(0, a - flank), min(length, b + flank)) for a, b in genes.get(chrom, [])]
        forbidden = buffered + [tuple(x) for x in genome.gaps[chrom]]
        universe[chrom] = subtract(length, forbidden)
    return universe


def match_intergenic(parents: list[dict], genome: Genome, repeats: object, genes: dict,
                     eligible: list[str], config: dict) -> tuple[dict[str, dict], list[dict], dict]:
    universe = intergenic_universe(genome, genes, eligible, config['flank'])
    stats = {'eligible_acgt_intergenic_bases': sum(b - a for v in universe.values() for a, b in v),
             'eligible_intervals': sum(map(len, universe.values())), 'parents': {}}
    matches = {}
    diagnostics = []
    for parent in sorted(parents, key=lambda f: (-f['length'], f['transcript_id'])):
        length = parent['length']
        intervals = [(chrom, a, b) for chrom in sorted(universe) for a, b in universe[chrom]
                     if b - a >= length]
        cumulative = []
        count = 0
        for _, a, b in intervals:
            count += b - a - length + 1
            cumulative.append(count)
        rng = random.Random(stable_seed(config['seed'], 'intergenic', parent['transcript_id']))
        best = None
        considered: set[tuple[str, int]] = set()
        acceptable = 0
        for _ in range(config['intergenic_candidates'] if count else 0):
            draw = rng.randrange(count)
            index = bisect.bisect_right(cumulative, draw)
            chrom, a, _ = intervals[index]
            start = a + draw - (cumulative[index - 1] if index else 0)
            if (chrom, start) in considered:
                continue
            considered.add((chrom, start))
            end = start + length
            sequence = genome.fetch(chrom, start, end)
            gc, n, ambiguous = composition(sequence)
            if ambiguous:
                raise ValueError('Internal error: intergenic universe contains ambiguous bases')
            repeat_hits = repeats.overlaps(chrom, start, end)
            fraction = coverage([(a, b) for a, b, _ in repeat_hits]) / length
            if abs(gc - parent['gc_fraction']) > config['gc_tolerance']:
                continue
            if (config['repeat_tolerance'] is not None
                    and abs(fraction - parent['repeat_fraction']) > config['repeat_tolerance']):
                continue
            acceptable += 1
            cost = ((gc - parent['gc_fraction']) / .1) ** 2
            cost += ((fraction - parent['repeat_fraction']) / .25) ** 2
            cost += 0 if chrom == parent['chromosome'] else .1
            candidate = dict(chromosome=chrom, start=start, end=end, strand='+', length=length,
                             gc_fraction=gc, n_fraction=n, ambiguous_fraction=ambiguous,
                             repeat_fraction=fraction, exon_count='', matching_cost=cost,
                             source_biotype='intergenic_DNA', transcript_id='', gene_id='', gene_symbol='')
            if best is None or (cost, chrom, start) < (best['matching_cost'], best['chromosome'], best['start']):
                best = candidate
        stats['parents'][parent['transcript_id']] = dict(available_starts=count,
                                                       sampled_candidates=len(considered),
                                                       acceptable_candidates=acceptable)
        if best is None:
            diagnostics.append(dict(control_type='intergenic_matched', parent_query_id=parent['transcript_id'],
                                    status='unmatched', reason='bounded_candidates_exhausted'))
            continue
        chrom, start, end = best['chromosome'], best['start'], best['end']
        # Remove the selected locus so all later parent intervals are disjoint.
        updated = []
        for a, b in universe[chrom]:
            if b <= start or a >= end:
                updated.append((a, b))
            else:
                if a < start:
                    updated.append((a, start))
                if b > end:
                    updated.append((end, b))
        universe[chrom] = updated
        matches[parent['transcript_id']] = best
        diagnostics.append(dict(control_type='intergenic_matched', parent_query_id=parent['transcript_id'],
                                status='matched', reason='within_tolerances'))
    return matches, diagnostics, stats
