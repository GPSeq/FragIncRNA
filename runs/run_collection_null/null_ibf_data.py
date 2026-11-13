"""Strict named-dataset discovery, streaming validation and filesystem safety."""
from __future__ import annotations

import csv
import gzip
import json
import os
import re
import shutil
import tomllib
from contextlib import contextmanager
from pathlib import Path

from run_fraglncrna_grid import FASTA_SUFFIXES, sha256_file, write_json

LABELS = ('real_lncrna', 'dinuc_shuffled', 'intergenic_matched', 'coding_matched', 'repeat_masked_lncrna')
ALIASES = dict(zip(('real_lncrna', 'dinuc_shuffle', 'intergenic_matched', 'coding_positive',
                    'repeat_masked_query'), LABELS))
PATTERNS = {
    'real_lncrna': r'real_lncrna(?:_\d+)?',
    'dinuc_shuffled': r'dinuc_shuffl(?:ed|e)(?:_\d+)?',
    'intergenic_matched': r'intergenic_matched(?:_\d+)?',
    'coding_matched': r'coding_matched(?:_\d+)?',
    'repeat_masked_lncrna': r'(?:real_lncrna(?:_\d+)?_repeat_masked|repeat_masked_lncrna(?:_\d+)?)',
}


def no_symlinks(path: Path) -> Path:
    """Check before resolve: resolving first would hide traversal through a link."""
    path = Path(os.path.abspath(path.expanduser()))
    for component in (path, *path.parents):
        if component.is_symlink():
            raise ValueError(f'Symlink traversal refused: {component}')
    return path


def regular_file(path: Path) -> Path:
    path = no_symlinks(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing, empty or nonregular input: {path}')
    return path


def fasta_file(path: Path) -> Path:
    path = regular_file(path)
    if not path.name.lower().endswith(FASTA_SUFFIXES):
        raise ValueError(f'Unsupported FASTA extension: {path}')
    return path


@contextmanager
def text_input(path: Path):
    with path.open('rb') as raw:
        compressed = raw.read(2) == b'\x1f\x8b'
    opener = gzip.open if compressed else open
    try:
        with opener(path, 'rt', encoding='utf-8') as handle:
            yield handle
    except (EOFError, gzip.BadGzipFile, UnicodeError, OSError) as error:
        raise ValueError(f'Corrupt/unreadable input {path}: {error}') from error


def scan_fasta(path: Path, k_values: list[int]) -> tuple[dict, set[str]]:
    """Keep only IDs and per-record counters, never whole sequences."""
    identifiers: set[str] = set()
    result = dict(path=str(path), sha256=sha256_file(path), bytes=path.stat().st_size,
                  query_count=0, bases=0, n_bases=0, ambiguous_bases=0, lowercase_acgt_bases=0,
                  valid_windows={str(k): 0 for k in k_values},
                  original_windows={str(k): 0 for k in k_values},
                  zero_valid_queries={str(k): 0 for k in k_values})
    current = None
    length = tail = 0
    valid = {k: 0 for k in k_values}

    def finish() -> None:
        if current is None:
            return
        if not length:
            raise ValueError(f'{path}: empty FASTA record {current}')
        result['query_count'] += 1
        result['bases'] += length
        for k in k_values:
            result['valid_windows'][str(k)] += valid[k]
            result['original_windows'][str(k)] += max(0, length - k + 1)
            result['zero_valid_queries'][str(k)] += valid[k] == 0

    with text_input(path) as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                finish()
                header = line[1:]
                if not header or header[0].isspace():
                    raise ValueError(f'{path}:{number}: empty/malformed FASTA header')
                current = header.split()[0]
                if current in identifiers:
                    raise ValueError(f'{path}: duplicate FASTA ID {current}')
                identifiers.add(current)
                length = tail = 0
                valid = {k: 0 for k in k_values}
            else:
                if current is None:
                    raise ValueError(f'{path}:{number}: sequence before FASTA header')
                sequence = ''.join(line.split())
                if re.search('[^ACGTRYKMSWBDHVNacgtrykmswbdhvn]', sequence):
                    raise ValueError(f'{path}:{number}: invalid DNA sequence')
                upper = sequence.upper()
                result['n_bases'] += upper.count('N')
                result['lowercase_acgt_bases'] += sum(m.end() - m.start() for m in re.finditer('[acgt]+', sequence))
                canonical = 0
                next_tail = 0
                for match in re.finditer('[ACGT]+', upper):
                    size = match.end() - match.start()
                    canonical += size
                    prior = tail if match.start() == 0 else 0
                    for k in k_values:
                        valid[k] += max(0, size + prior - k + 1) - max(0, prior - k + 1)
                    if match.end() == len(upper):
                        next_tail = prior + size
                tail = next_tail
                result['ambiguous_bases'] += len(upper) - canonical
                length += len(upper)
    finish()
    if not identifiers:
        raise ValueError(f'{path}: no FASTA records')
    return result, identifiers


def probe_genome(path: Path) -> None:
    with text_input(path) as handle:
        header = None
        for _ in range(100):
            line = handle.readline(1024 * 1024)
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            if header is None:
                if not line.startswith('>') or not line[1:].strip():
                    raise ValueError(f'Invalid reference FASTA header: {path}')
                header = line
            else:
                if line.startswith('>') or re.search(r'[^ACGTRYKMSWBDHVNacgtrykmswbdhvn\s]', line):
                    raise ValueError(f'Empty/invalid first reference record: {path}')
                return
    raise ValueError(f'Reference has no initial sequence: {path}')


def discover_genomes(directory: Path, selections: list[str] | None) -> list[Path]:
    directory = no_symlinks(directory)
    if not directory.is_dir():
        raise ValueError(f'Genome directory does not exist: {directory}')
    candidates = sorted((p for p in directory.iterdir() if p.name.lower().endswith(FASTA_SUFFIXES)),
                        key=lambda p: p.name)
    if selections:
        chosen = set()
        for selector in selections:
            if Path(selector).is_absolute():
                matches = [Path(selector)]  # explicit absolute paths do not expand wildcards
            else:
                if Path(selector).name != selector:
                    raise ValueError('Genome selectors must be flat filenames/globs or absolute file paths')
                matches = sorted(directory.glob(selector))
            if len(matches) != 1:
                raise ValueError(f'Genome selector {selector!r} must match exactly one file; '
                                 f'matches={[str(p) for p in matches]}')
            chosen.add(fasta_file(matches[0]))
        candidates = sorted(chosen, key=lambda p: (p.name, str(p)))
    genomes = [fasta_file(p) for p in candidates]
    if not genomes:
        raise ValueError(f'No reference genome FASTAs in {directory}')
    if len({p.name for p in genomes}) != len(genomes):
        raise ValueError('Reference filenames must be unique, including absolute selections')
    for genome in genomes:
        probe_genome(genome)
    return genomes


def _relative(path: str, parent: Path) -> Path:
    if not isinstance(path, str) or not path:
        raise ValueError('Manifest paths must be nonempty strings')
    p = Path(path).expanduser()
    return regular_file(p if p.is_absolute() else parent / p)


def discover_datasets(manifest: Path | None, directory: Path | None, selected: list[str],
                      k_values: list[int]) -> tuple[dict, dict, dict]:
    paths = {}
    metadata = {}
    generator = None
    if manifest:
        manifest = regular_file(manifest)
        data = tomllib.loads(manifest.read_text())
        if set(data) - {'datasets', 'metadata'} or not isinstance(data.get('datasets'), dict):
            raise ValueError('Dataset manifest requires [datasets] and optional [metadata] only')
        unknown = set(data['datasets']) - set(LABELS)
        if unknown:
            raise ValueError(f'Unknown dataset labels: {sorted(unknown)}')
        paths = {key: _relative(value, manifest.parent) for key, value in data['datasets'].items()}
        if not isinstance(data.get('metadata', {}), dict):
            raise ValueError('[metadata] must map labels to file paths')
        metadata = {key: _relative(value, manifest.parent) for key, value in data.get('metadata', {}).items()}
    else:
        directory = no_symlinks(directory)
        if not directory.is_dir():
            raise ValueError(f'Datasets directory does not exist: {directory}')
        roots = [directory, directory.parent] if directory.name == 'fasta' else [directory]
        root = next((p for p in roots if (p / 'provenance/run_manifest.json').exists()), None)
        if root:
            source = regular_file(root / 'provenance/run_manifest.json')
            generator = json.loads(source.read_text())
            if generator.get('status') != 'complete' or generator.get('tool') != 'FragIncRNA-null-models':
                raise ValueError(f'Unrecognized/incomplete generator manifest: {source}')
            for key, value in generator.get('fasta_paths', {}).items():
                if key not in ALIASES:
                    raise ValueError(f'Unrecognized generator control type: {key}')
                paths[ALIASES[key]] = _relative(value, root)
            metadata['source_manifest'] = source
            for path in sorted((root / 'metadata').glob('*')):
                if path.is_file():
                    name = path.name.removesuffix('.gz').removesuffix('.tsv').removesuffix('.bed')
                    metadata[name] = regular_file(path)
            if (root / 'provenance/outputs.sha256.tsv').exists():
                metadata['generator_output_receipt'] = regular_file(root / 'provenance/outputs.sha256.tsv')
        else:
            search = directory / 'fasta' if (directory / 'fasta').is_dir() else directory
            for label in LABELS:
                matches = []
                for path in sorted(search.iterdir()):
                    stem = path.name
                    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
                        if stem.lower().endswith(suffix):
                            stem = stem[:-len(suffix)]
                            if re.fullmatch(PATTERNS[label], stem):
                                matches.append(path)
                            break
                if len(matches) > 1:
                    raise ValueError(f'Ambiguous {label} dataset: {matches}; supply --datasets-manifest')
                if matches:
                    paths[label] = fasta_file(matches[0])
    missing = set(selected) - paths.keys()
    if missing:
        raise ValueError(f'Missing required datasets {sorted(missing)}; available={sorted(paths)}')
    if len(set(paths.values())) != len(paths):
        raise ValueError('Distinct logical datasets may not map to the same FASTA path')
    datasets = {}
    ids = {}
    for label in LABELS:
        if label not in paths:
            continue
        info, identifiers = scan_fasta(fasta_file(paths[label]), k_values)
        datasets[label] = info
        ids[label] = identifiers
        if generator:
            source_label = next(k for k, value in ALIASES.items() if value == label)
            if generator.get('counts', {}).get(source_label) != info['query_count']:
                raise ValueError(f'Generator manifest / FASTA record count mismatch: {label}')
        if label == 'repeat_masked_lncrna':
            if info['lowercase_acgt_bases'] and not info['n_bases']:
                raise ValueError('Repeat-masked dataset appears soft-masked only: lowercase ACGT is '
                                 'ordinary sequence with repeat=include/ambiguity=skip. Supply N-masked queries.')
            info['masking_detection'] = 'hard_N_present' if info['n_bases'] else 'no_N_mask_detected'
    meta = {key: {'path': str(path), 'sha256': sha256_file(path), 'bytes': path.stat().st_size}
            for key, path in sorted(metadata.items())}
    if 'generator_output_receipt' in metadata:
        with metadata['generator_output_receipt'].open() as handle:
            receipt = {row['path']: row for row in csv.DictReader(handle, delimiter='\t')}
        root = metadata['source_manifest'].parents[1]
        verified_inputs = list(datasets.values()) + [v for key, v in meta.items() if key != 'generator_output_receipt']
        for info in verified_inputs:
            relative = str(Path(info['path']).relative_to(root))
            if relative not in receipt or receipt[relative]['sha256'] != info['sha256']:
                raise ValueError(f'Generated FASTA/metadata checksum differs from generator receipt: {relative}')
    metadata_table = metadata.get('query_metadata', metadata.get('control_metadata'))
    if metadata_table:
        seen = {label: set() for label in datasets}
        with text_input(metadata_table) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if {'query_id', 'control_type'} <= set(reader.fieldnames or []):
                for row in reader:
                    label = ALIASES.get(row['control_type'], row['control_type'])
                    if label in seen:
                        qid = row['query_id']
                        if qid in seen[label]:
                            raise ValueError(f'Duplicate query metadata ID: {label}/{qid}')
                        seen[label].add(qid)
                if any(seen[label] != ids[label] for label in datasets):
                    raise ValueError('Prepared FASTA and query metadata ID populations disagree')
    return datasets, meta, ids


def safe_tree(path: Path, root: Path) -> None:
    path = no_symlinks(path)
    root = no_symlinks(root)
    if path == root or not path.is_relative_to(root):
        raise ValueError(f'Unsafe managed path outside output root: {path}')
    if path.exists():
        for parent, directories, files in os.walk(path, followlinks=False):
            for name in directories + files:
                if (Path(parent) / name).is_symlink():
                    raise ValueError(f'Symlink in managed tree refused: {Path(parent) / name}')


def remove_private(path: Path, root: Path) -> None:
    safe_tree(path, root)
    if path.exists():
        shutil.rmtree(path)


@contextmanager
def campaign_lock(root: Path):
    path = root / '.runner.lock'
    no_symlinks(path)
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise ValueError(f'Output is locked: {path}; check the recorded PID before removing a stale lock') from error
    try:
        with os.fdopen(fd, 'w') as handle:
            handle.write(f'{os.getpid()}\n')
        yield
    finally:
        path.unlink(missing_ok=True)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    no_symlinks(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + f'.part-{os.getpid()}')
    with part.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    part.replace(path)


def compatible_json(path: Path, value: dict) -> None:
    no_symlinks(path)
    if path.exists() and json.loads(path.read_text()) != value:
        raise ValueError(f'Provenance mismatch: {path}; choose a new output root')
    if not path.exists():
        write_json(path, value)
