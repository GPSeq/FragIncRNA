"""Verified atomic archives and schema-checked candidate-local result summaries."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import statistics
import tarfile
from pathlib import Path, PurePosixPath

from null_ibf_data import no_symlinks, safe_tree
from run_fraglncrna_grid import iso_now, sha256_file, write_json


def analyze_results(results: Path, identifiers: set[str], k: int, hashes: int) -> tuple[dict, dict]:
    main = results / 'lncrna.tsv'
    if not main.is_file() or not main.stat().st_size:
        raise ValueError(f'Executable succeeded but main result TSV is missing/empty: {main}')
    metrics_path = results / 'metrics.json'
    if not metrics_path.is_file():
        raise ValueError('Executable did not write requested metrics.json; use a compatible FragIncRNA v2 binary')
    metrics = json.loads(metrics_path.read_text())
    for key in ('whole_reference_exact_scans', 'whole_reference_bases_visited'):
        if metrics.get(key, 0) != 0:
            raise ValueError(f'IBF-only campaign detected forbidden global-exact work: {key}')
    for key in ('shards_built', 'shards_reused', 'shard_count'):
        if not isinstance(metrics.get(key), int) or metrics[key] < 0:
            raise ValueError(f'Binary metrics lack a valid {key}; shard reuse cannot be verified')
    if metrics['shards_built'] + metrics['shards_reused'] != metrics['shard_count']:
        raise ValueError('Inconsistent shard build/reuse counts from executable')
    if metrics.get('mode', 'ibf') != 'ibf' or metrics.get('scoring_method', 'per_bin_verified') != 'per_bin_verified':
        raise ValueError('Executable metrics disagree with the IBF candidate-local verified campaign')
    values = []
    statuses = {}
    seen = set()
    valid_total = zero_valid = 0
    with main.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = set(reader.fieldnames or [])
        supported = {'query_id', 'valid_windows', 'per_bin_verified_local_recovery', 'status',
                     'mode', 'scoring_method', 'verification_scope', 'k', 'hash_functions'} <= fields
        for row in reader:
            if 'query_id' in fields:
                qid = row['query_id']
                if qid in seen or qid not in identifiers:
                    raise ValueError(f'Duplicate or unexpected result query ID: {qid}')
                seen.add(qid)
            if not supported:
                continue
            if (row['mode'] != 'ibf' or row['scoring_method'] != 'per_bin_verified'
                    or row['verification_scope'] != 'candidate_intervals'
                    or int(row['k']) != k or int(row['hash_functions']) != hashes):
                raise ValueError('Result TSV contains incompatible method/scope/k/hash values')
            valid = int(row['valid_windows'])
            if valid < 0:
                raise ValueError('Negative valid-window count in result')
            valid_total += valid
            zero_valid += valid == 0
            score = row['per_bin_verified_local_recovery']
            if score not in ('NA', ''):
                value = float(score)
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError('Candidate-local verified recovery must be in [0,1] or NA')
                values.append(value)
            statuses[row['status']] = statuses.get(row['status'], 0) + 1
        if 'query_id' in fields and seen != identifiers:
            raise ValueError('Result TSV query IDs/count differ from the input FASTA')
    if not supported:
        return {'aggregation_status': 'unsupported_schema', 'query_count': len(identifiers)}, metrics
    n = len(identifiers)
    basic = statuses.get('PASS_BASIC', 0)
    strict = statuses.get('PASS_STRICT', 0)
    return dict(aggregation_status='available', metric='candidate_local_verified_recovery',
                query_count=n, numeric_recovery_count=len(values),
                mean_recovery=statistics.fmean(values) if values else '',
                median_recovery=statistics.median(values) if values else '',
                pass_basic_count=basic, pass_basic_fraction=basic / n,
                pass_strict_count=strict, pass_strict_fraction=strict / n,
                at_least_basic_count=basic + strict, at_least_basic_fraction=(basic + strict) / n,
                valid_windows=valid_total, zero_valid_queries=zero_valid), metrics


def archive_members(work: Path) -> list[Path]:
    safe_tree(work, work.parent)  # full campaign containment is checked by the caller
    files = [work / name for name in ('generated_config.toml', 'template.toml', 'job.json',
                                      'execution.json', 'README_archive.txt', 'reference_index_manifest.tsv')]
    files.extend(p for p in (work / 'results').rglob('*') if p.is_file())
    for path in files:
        no_symlinks(path)
        if not path.is_file():
            raise ValueError(f'Missing archive member: {path}')
    return sorted(files, key=lambda p: str(p.relative_to(work)))


def verify_archive(path: Path, expected_signature: str) -> dict:
    no_symlinks(path)
    with gzip.open(path, 'rb') as compressed:
        with tarfile.open(fileobj=compressed, mode='r|') as archive:
            first = archive.next()
            if (first is None or not first.isfile() or first.name != 'archive_manifest.json'
                    or first.size > 16 * 1024 * 1024):
                raise ValueError(f'Invalid archive manifest member: {path}')
            handle = archive.extractfile(first)
            manifest = json.load(handle)
            if (manifest.get('format') != 'fraglncrna-null-ibf-archive-v1'
                    or manifest.get('job_signature') != expected_signature):
                raise ValueError(f'Archive provenance mismatch: {path}')
            seen = set()
            for member in iter(archive.next, None):
                name = member.name
                if (not member.isfile() or name in seen or PurePosixPath(name).is_absolute()
                        or '..' in PurePosixPath(name).parts or '\\' in name
                        or name not in manifest['files']):
                    raise ValueError(f'Unsafe/unexpected/duplicate archive member: {name}')
                expected = manifest['files'][name]
                if member.size != expected['size']:
                    raise ValueError(f'Archive member size mismatch: {name}')
                digest = hashlib.sha256()
                with archive.extractfile(member) as source:
                    for chunk in iter(lambda: source.read(4 * 1024 * 1024), b''):
                        digest.update(chunk)
                if digest.hexdigest() != expected['sha256']:
                    raise ValueError(f'Archive member checksum mismatch: {name}')
                seen.add(name)
            if seen != set(manifest['files']):
                raise ValueError('Archive member list incomplete')
        # Drain trailing tar padding and gzip trailer to verify CRC/truncation.
        for _ in iter(lambda: compressed.read(4 * 1024 * 1024), b''):
            pass
    return manifest


def make_receipt(archive: Path, manifest: dict) -> dict:
    return dict(status='complete', archive=archive.name, archive_sha256=sha256_file(archive),
                archive_bytes=archive.stat().st_size, job_signature=manifest['job_signature'],
                completed_utc=manifest['finished_utc'], elapsed_seconds=manifest['elapsed_seconds'],
                result_summary=manifest['result_summary'], index_counts=manifest['index_counts'])


def create_archive(work: Path, archive: Path, signature: str, elapsed: float,
                   summary: dict, metrics: dict) -> dict:
    files = archive_members(work)
    manifest = dict(format='fraglncrna-null-ibf-archive-v1', job_signature=signature,
                    elapsed_seconds=round(elapsed, 3), finished_utc=iso_now(), result_summary=summary,
                    index_counts={key: metrics[key] for key in ('shards_built', 'shards_reused', 'shard_count')},
                    files={str(p.relative_to(work)): {'sha256': sha256_file(p), 'size': p.stat().st_size}
                           for p in files})
    write_json(work / 'archive_manifest.json', manifest)
    archive.parent.mkdir(parents=True, exist_ok=True)
    partial = archive.with_name(archive.name + '.part')
    no_symlinks(partial)
    # A failure leaves this partial and the original results/cache for inspection.
    with tarfile.open(partial, 'w:gz') as output:
        output.add(work / 'archive_manifest.json', arcname='archive_manifest.json', recursive=False)
        for path in files:
            output.add(path, arcname=str(path.relative_to(work)), recursive=False)
    if verify_archive(partial, signature) != manifest:
        raise ValueError('New archive did not preserve the complete manifest')
    if archive.exists():
        raise ValueError(f'Refusing to replace an existing archive: {archive}')
    partial.replace(archive)
    return make_receipt(archive, manifest)


def verify_completion(archive: Path, receipt: Path, signature: str, repair: bool) -> dict | None:
    no_symlinks(archive)
    no_symlinks(receipt)
    if not archive.exists() and not receipt.exists():
        return None
    if receipt.exists() and not archive.is_file():
        raise ValueError(f'Completion receipt exists without archive: {archive}')
    if not receipt.exists():
        if not repair:
            # Read-only preflight can validate an orphan, but never publishes a receipt.
            manifest = verify_archive(archive, signature)
            return dict(make_receipt(archive, manifest), needs_receipt_repair=True)
        manifest = verify_archive(archive, signature)
        payload = make_receipt(archive, manifest)
        write_json(receipt, payload)
        return payload
    payload = json.loads(receipt.read_text())
    if payload.get('job_signature') != signature:
        raise ValueError(f'Completed archive exists for different inputs/settings: {archive}; choose a new output root')
    if (payload.get('status') != 'complete' or payload.get('archive') != archive.name
            or payload.get('archive_bytes') != archive.stat().st_size
            or payload.get('archive_sha256') != sha256_file(archive)):
        raise ValueError(f'Archive/receipt checksum mismatch: {archive}; preserve and investigate')
    manifest = verify_archive(archive, signature)
    expected = {key: manifest[key] for key in ('elapsed_seconds', 'result_summary', 'index_counts')}
    if any(payload.get(key) != value for key, value in expected.items()):
        raise ValueError(f'Completion receipt content disagrees with archive manifest: {receipt}')
    return payload
