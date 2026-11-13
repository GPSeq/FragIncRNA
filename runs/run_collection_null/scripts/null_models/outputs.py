"""Deterministic streaming artifacts, paired metadata and independent verification."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from biology import composition, coverage, mask_sequence, project_repeats, shuffle_parent, union, valid_windows
from inputs import fasta_records, open_text, sha256
from preparation import fetch_record

META_FIELDS = ['query_id', 'control_type', 'parent_query_id', 'replicate', 'transcript_id', 'gene_id',
               'gene_symbol', 'source_biotype', 'source_gene_biotype', 'length', 'gc_fraction', 'n_fraction', 'ambiguous_fraction',
               'repeat_fraction', 'repeat_fraction_basis', 'exon_count', 'chromosome', 'strand',
               'source_coordinates', 'matching_cost', 'relative_length_error', 'gc_delta', 'repeat_delta',
               'exon_count_delta', 'splice_mismatches', 'quality_flags', 'stratum', 'rng_seed', 'sequence_sha256']
SHUFFLE_FIELDS = META_FIELDS + ['attempts', 'unchanged_attempts', 'duplicate_attempts',
                               'unavoidable_duplicate', 'unchanged']
WINDOW_FIELDS = ['query_id', 'control_type', 'parent_query_id', 'k', 'original_windows',
                 'valid_acgt_windows', 'masking_induced_lost_windows']


@contextmanager
def gzip_text(path: Path) -> Iterator[TextIO]:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('wb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as text:
                yield text


def table_writer(handle: TextIO, fields: list[str]) -> csv.DictWriter:
    writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader()
    return writer


def write_table(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = table_writer(handle, fields)
        writer.writerows(rows)


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def write_fasta(handle: TextIO, header: str, sequence: str) -> None:
    handle.write('>' + header + '\n')
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + '\n')


def emit_outputs(stage: Path, db: object, genome: object, repeats: object, parents: list[dict],
                 coding: dict, intergenic: dict, config: dict) -> tuple[dict, dict]:
    n = len(parents)
    paths = {
        'real_lncrna': f'fasta/real_lncrna_{n}.fa.gz',
        'dinuc_shuffle': f'fasta/dinuc_shuffled_{n * config["shuffles_per_transcript"]}.fa.gz',
        'intergenic_matched': f'fasta/intergenic_matched_{len(intergenic)}.fa.gz',
        'coding_positive': f'fasta/coding_matched_{len(coding)}.fa.gz',
        'repeat_masked_query': f'fasta/real_lncrna_{n}_repeat_masked.fa.gz',
    }
    counts = {kind: 0 for kind in paths}
    shuffle_qc = dict(unavoidable_duplicate=0, unchanged=0, duplicate_attempts=0, unchanged_attempts=0)
    summary = {'counts': counts, 'shuffle_qc': shuffle_qc, 'sequence_policy':
               'Baseline preserves original sequence/case/header; shuffle retains case at each position; '
               'reference controls uppercase; repeats replaced by uppercase N.'}
    matching_rows = []
    with ExitStack() as stack:
        fastas = {kind: stack.enter_context(gzip_text(stage / path)) for kind, path in paths.items()}
        def writer(name: str, fields: list[str]) -> csv.DictWriter:
            return table_writer(stack.enter_context(gzip_text(stage / 'metadata' / name)), fields)
        all_meta = writer('query_metadata.tsv.gz', META_FIELDS)
        selected = writer('selected_lncrnas.tsv.gz', META_FIELDS)
        shuffle_map = writer('shuffled_parent_map.tsv.gz', SHUFFLE_FIELDS)
        coding_map = writer('coding_matches.tsv.gz', META_FIELDS)
        intergenic_map = writer('intergenic_matches.tsv.gz', META_FIELDS)
        mask_map = writer('repeat_mask_stats.tsv.gz', META_FIELDS + ['masked_bases', 'mask_fraction'])
        windows = writer('per_query_valid_windows.tsv.gz', WINDOW_FIELDS)
        class_stats = writer('transcript_repeat_classes.tsv.gz',
                             ['query_id', 'parent_query_id', 'repeat_class', 'covered_bases', 'fraction'])
        bed = stack.enter_context(gzip_text(stage / 'metadata/transcript_repeat_intervals.bed.gz'))
        ineligible = writer('ineligible_transcripts.tsv.gz',
                            ['transcript_id', 'source_kind', 'reason', 'length', 'k',
                             'original_windows', 'valid_acgt_windows'])
        for tid, kind, reason, sequence in db.execute(
                'SELECT i.tid,i.kind,i.reason,r.sequence FROM ineligible i JOIN records r USING(tid) '
                'ORDER BY i.kind,i.tid'):
            for k in sorted(set(config['k_values']) | {15, 20, 32}):
                ineligible.writerow(dict(transcript_id=tid, source_kind=kind, reason=reason,
                                         length=len(sequence), k=k, original_windows=max(0, len(sequence) - k + 1),
                                         valid_acgt_windows=valid_windows(sequence, k)))

        def emit(kind: str, header: str, sequence: str, feature: dict, parent_id: str = '',
                 parent_feature: dict | None = None, baseline: str | None = None,
                 extra: dict | None = None) -> dict:
            query_id = header.split()[0]
            row = dict(feature, query_id=query_id, control_type=kind, parent_query_id=parent_id,
                       replicate='', rng_seed=config['seed'], sequence_sha256=sequence_hash(sequence),
                       repeat_fraction_basis='genomic_annotation_union')
            gc, n_fraction, ambiguous = composition(sequence)
            row.update(length=len(sequence), gc_fraction=gc, n_fraction=n_fraction,
                       ambiguous_fraction=ambiguous)
            if parent_feature is not None:
                row.update(relative_length_error=(len(sequence) - parent_feature['length']) / parent_feature['length'],
                           gc_delta=gc - parent_feature['gc_fraction'],
                           repeat_delta=feature['repeat_fraction'] - parent_feature['repeat_fraction'],
                           exon_count_delta=(feature['exon_count'] - parent_feature['exon_count']
                                             if isinstance(feature['exon_count'], int) else ''))
            if extra:
                row.update(extra)
            write_fasta(fastas[kind], header, sequence)
            all_meta.writerow(row)
            counts[kind] += 1
            for k in config['k_values']:
                valid = valid_windows(sequence, k)
                lost = valid_windows(baseline, k) - valid if baseline is not None else 0
                windows.writerow(dict(query_id=query_id, control_type=kind, parent_query_id=parent_id,
                                      k=k, original_windows=max(0, len(sequence) - k + 1),
                                      valid_acgt_windows=valid, masking_induced_lost_windows=lost))
            return row

        tasks = []
        source_by_tid = {}
        for feature in parents:
            tid = feature['transcript_id']
            header, sequence, transcript = fetch_record(db, tid)
            parent_id = header.split()[0]
            source_coordinates = ','.join(f'{transcript.chrom}:{a}-{b}:{transcript.strand}'
                                          for a, b in transcript.ordered_exons())
            feature = dict(feature, source_coordinates=source_coordinates)
            selected.writerow(emit('real_lncrna', header, sequence, feature))
            source_by_tid[tid] = (feature, parent_id, sequence)
            tasks.append((parent_id, sequence, config['seed'], config['shuffles_per_transcript'],
                          config['shuffle_retries']))
            projected = project_repeats(transcript, repeats)
            intervals = union([i for group in projected.values() for i in group])
            masked = mask_sequence(sequence, intervals)
            mask_header = f'{tid}|repeat_masked_query|parent_query_id={parent_id}'
            row = emit('repeat_masked_query', mask_header, masked, feature, parent_id,
                       baseline=sequence, extra={'masked_bases': coverage(intervals),
                                                 'mask_fraction': coverage(intervals) / len(sequence)})
            mask_map.writerow(row)
            for a, b in intervals:
                bed.write(f'{row["query_id"]}\t{a}\t{b}\trepeat_union\n')
            for label, values in projected.items():
                class_stats.writerow(dict(query_id=row['query_id'], parent_query_id=parent_id,
                                          repeat_class=label, covered_bases=coverage(values),
                                          fraction=coverage(values) / len(sequence)))
            for matches, kind, pair_writer in [(coding, 'coding_positive', coding_map),
                                                (intergenic, 'intergenic_matched', intergenic_map)]:
                if tid not in matches:
                    continue
                match = dict(matches[tid])
                if kind == 'coding_positive':
                    _, control_sequence, control_tx = fetch_record(db, match['transcript_id'])
                    control_header = f'{match["transcript_id"]}|{kind}|parent_query_id={parent_id}'
                    match['source_coordinates'] = ','.join(
                        f'{control_tx.chrom}:{a}-{b}:{control_tx.strand}' for a, b in control_tx.ordered_exons())
                else:
                    chrom, start, end = match['chromosome'], match['start'], match['end']
                    control_sequence = genome.fetch(chrom, start, end)
                    control_header = f'intergenic|{chrom}:{start}-{end}:+|parent_query_id={parent_id}'
                    match['source_coordinates'] = f'{chrom}:{start}-{end}:+'
                pair = emit(kind, control_header, control_sequence, match, parent_id, feature)
                pair_writer.writerow(pair)
                matching_rows.append({key: pair.get(key, '') for key in [
                    'parent_query_id', 'query_id', 'control_type', 'length', 'gc_fraction', 'repeat_fraction',
                    'relative_length_error', 'gc_delta', 'repeat_delta', 'exon_count_delta', 'matching_cost']})
        executor = (stack.enter_context(ProcessPoolExecutor(max_workers=config['threads']))
                    if config['threads'] > 1 else None)
        results = executor.map(shuffle_parent, tasks, chunksize=1) if executor else map(shuffle_parent, tasks)
        for original_feature, shuffles in zip(parents, results):
            tid = original_feature['transcript_id']
            feature, parent_id, _ = source_by_tid[tid]
            for sequence, extra in shuffles:
                replicate = extra['replicate']
                header = f'{tid}|dinuc_shuffle|replicate={replicate}|parent_query_id={parent_id}'
                extra = dict(extra, repeat_fraction='',
                             repeat_fraction_basis='not_applicable_parent_annotation_in_selected_table',
                             source_coordinates='')
                row = emit('dinuc_shuffle', header, sequence, feature, parent_id, extra=extra)
                # A shuffled sequence has no repeat-annotation coordinates; source repeat
                # coverage is retained separately through its parent, never inferred here.
                row['repeat_fraction'] = ''
                row['repeat_fraction_basis'] = 'not_applicable_parent_annotation_in_selected_table'
                shuffle_map.writerow(row)
                for key in shuffle_qc:
                    shuffle_qc[key] += int(extra[key])
    write_table(stage / 'qc/length_gc_repeat_matching.tsv', matching_rows)
    return paths, summary


def verify_outputs(root: Path, fasta_paths: dict, k_values: list[int]) -> dict:
    """Re-open gzip members to CRC-check and compare every sequence to metadata."""
    metadata = {}
    with open_text(root / 'metadata/query_metadata.tsv.gz') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['query_id'] in metadata:
                raise ValueError('Duplicate output metadata query ID')
            metadata[row['query_id']] = row
    windows = {}
    with open_text(root / 'metadata/per_query_valid_windows.tsv.gz') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            key = (row['query_id'], int(row['k']))
            if key in windows:
                raise ValueError('Duplicate valid-window row')
            windows[key] = row
    seen = set()
    counts = {}
    for kind, relative in fasta_paths.items():
        count = 0
        for query_id, _, sequence in fasta_records(root / relative, canonical_ids=False, allow_empty=True):
            row = metadata[query_id]
            if (query_id in seen or row['control_type'] != kind or len(sequence) != int(row['length'])
                    or sequence_hash(sequence) != row['sequence_sha256']):
                raise ValueError(f'Output metadata disagreement for {query_id}')
            for k in k_values:
                w = windows[(query_id, k)]
                if (int(w['valid_acgt_windows']) != valid_windows(sequence, k)
                        or int(w['original_windows']) != max(0, len(sequence) - k + 1)):
                    raise ValueError(f'Output valid-window disagreement for {query_id}, k={k}')
                expected_loss = (int(windows[(row['parent_query_id'], k)]['valid_acgt_windows'])
                                 - int(w['valid_acgt_windows']) if kind == 'repeat_masked_query' else 0)
                if int(w['masking_induced_lost_windows']) != expected_loss:
                    raise ValueError(f'Output masked-window disagreement for {query_id}, k={k}')
            if row['parent_query_id'] and row['parent_query_id'] not in metadata:
                raise ValueError(f'Missing parent {row["parent_query_id"]}')
            seen.add(query_id)
            count += 1
        counts[kind] = count
    if seen != set(metadata) or len(windows) != len(seen) * len(k_values):
        raise ValueError('FASTA and metadata populations disagree')
    linked_tables = {'selected_lncrnas.tsv.gz': 'real_lncrna',
                     'shuffled_parent_map.tsv.gz': 'dinuc_shuffle',
                     'intergenic_matches.tsv.gz': 'intergenic_matched',
                     'coding_matches.tsv.gz': 'coding_positive',
                     'repeat_mask_stats.tsv.gz': 'repeat_masked_query'}
    for filename, kind in linked_tables.items():
        table_ids = set()
        with open_text(root / 'metadata' / filename) as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                query_id = row['query_id']
                if query_id in table_ids or any(row[key] != metadata[query_id][key] for key in META_FIELDS):
                    raise ValueError(f'Linked metadata mismatch in {filename}: {query_id}')
                table_ids.add(query_id)
        if table_ids != {qid for qid, row in metadata.items() if row['control_type'] == kind}:
            raise ValueError(f'Linked metadata population mismatch in {filename}')
    # All other gzip metadata are independently fully read, checking truncation/CRC.
    for path in sorted((root / 'metadata').glob('*.gz')):
        with open_text(path) as handle:
            for _ in handle:
                pass
    return {'verified_counts': counts, 'all_gzip_members_read': True,
            'sequence_hashes_and_valid_windows_verified': True}


def write_output_hashes(root: Path) -> None:
    # The receipt excludes itself; manifest is already final before hashing.
    rows = [dict(path=str(path.relative_to(root)), sha256=sha256(path), size_bytes=path.stat().st_size)
            for path in sorted(root.rglob('*')) if path.is_file()
            and path != root / 'provenance/outputs.sha256.tsv']
    write_table(root / 'provenance/outputs.sha256.tsv', rows)


def verify_receipt(root: Path) -> None:
    receipt = root / 'provenance/outputs.sha256.tsv'
    with receipt.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    expected = {row['path'] for row in rows} | {'provenance/outputs.sha256.tsv'}
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if expected != actual:
        raise ValueError('Resume output file set differs from recorded provenance')
    for row in rows:
        path = root / row['path']
        if path.resolve().is_relative_to(root.resolve()) and sha256(path) == row['sha256']:
            continue
        raise ValueError(f'Resume checksum mismatch: {path}')
