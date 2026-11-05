#!/usr/bin/env python3
"""ENA metadata and optional verified paired FASTQ downloads; standard library only.

One table row per run, with sample/experiment IDs for identifying technical runs.
Missing donor, region and strandedness remain unknown; never infer from filenames.
Run only one copy per READ_ROOT. Existing external downloads must finish before alignment.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

SPECIES = {'Homo sapiens': 'human', 'Pan troglodytes': 'chimpanzee',
           'Pan paniscus': 'bonobo', 'Macaca mulatta': 'macaque'}
FIELDS = ('run_accession,sample_accession,experiment_accession,study_accession,'
          'scientific_name,library_layout,library_strategy,library_source,'
          'library_selection,fastq_ftp,fastq_md5,fastq_bytes')
RUN_RE = re.compile(r'(?:SRR|ERR|DRR)\d+')
COLUMNS = ['run_accession', 'species', 'sample_accession', 'experiment_accession',
           'study_accession', 'donor', 'brain_region', 'sex', 'age',
           'library_layout', 'library_strategy', 'library_source', 'library_selection',
           'strandedness', 'metadata_status', 'download_status', 'fastq_1', 'fastq_2',
           'fastq_ftp', 'fastq_md5', 'fastq_bytes', 'sample_attributes_json', 'error']


def fetch(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return response.read()
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def cached(path, url):
    if path.exists():
        return path.read_bytes()
    data = fetch(url)
    part = path.with_suffix(path.suffix + '.part')
    part.write_bytes(data)
    part.replace(path)
    return data


def report(run, folder):
    query = urllib.parse.urlencode(dict(accession=run, result='read_run',
                                        fields=FIELDS, format='tsv'))
    data = cached(folder / f'{run}.tsv', 'https://www.ebi.ac.uk/ena/portal/api/filereport?' + query)
    rows = list(csv.DictReader(io.StringIO(data.decode()), delimiter='\t'))
    if len(rows) != 1 or rows[0]['run_accession'] != run:
        raise ValueError(f'No unique ENA record for {run}')
    return rows[0]


def attributes(sample, folder):
    data = cached(folder / f'{sample}.xml', 'https://www.ebi.ac.uk/ena/browser/api/xml/' + sample)
    root = ET.fromstring(data)
    if root.find('.//SAMPLE') is None:
        raise ValueError(f'No SAMPLE record for {sample}')
    return {x.findtext('TAG', '').strip().lower(): x.findtext('VALUE', '').strip()
            for x in root.findall('.//SAMPLE_ATTRIBUTE')}


def first(attrs, *keys):
    return next((attrs[k] for k in keys if attrs.get(k)), 'unknown')


def md5(path):
    digest = hashlib.md5()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def download(url, dest, checksum, size):
    # Never overwrite an existing file that fails validation: it may be an active
    # external download or a different FASTQ representation (e.g. fasterq-dump).
    if dest.exists():
        if dest.stat().st_size == size and md5(dest) == checksum:
            return
        raise ValueError(f'Existing file differs from ENA: {dest}; inspect before retrying')
    part = dest.with_suffix(dest.suffix + '.part')
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as response, part.open('wb') as out:
                while True:
                    block = response.read(8 * 1024 * 1024)
                    if not block:
                        break
                    out.write(block)
            if part.stat().st_size != size or md5(part) != checksum:
                raise ValueError(f'Size/MD5 mismatch: {dest.name}')
            part.replace(dest)
            return
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def save_table(path, rows):
    part = path.with_suffix('.tsv.part')
    with part.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows[k] for k in sorted(rows))
    part.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--read-root', required=True)
    parser.add_argument('--species', choices=['all', *SPECIES.values()], default='all')
    parser.add_argument('--max-runs', type=int, default=0)
    parser.add_argument('--accessions', type=Path)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    root = Path(args.read_root).resolve()
    folder = root / 'metadata'
    cache = folder / 'ena_cache'
    cache.mkdir(parents=True, exist_ok=True)
    table = folder / 'samples.tsv'
    rows = {}
    if table.exists():
        with table.open() as handle:
            rows = {r['run_accession']: {k: r.get(k, '') for k in COLUMNS}
                    for r in csv.DictReader(handle, delimiter='\t')}
    if args.accessions:
        runs = []
        for line in args.accessions.read_text().splitlines():
            run = line.split('#')[0].strip()
            if not run:
                continue
            if not RUN_RE.fullmatch(run):
                parser.error(f'Expected one run accession per line, got: {run}')
            runs.append(run)
        runs = list(dict.fromkeys(runs))
    else:
        runs = sorted({p.name.split('_')[0] for sp in SPECIES.values()
                       for p in (root / sp).glob('*_[12].fastq.gz')
                       if RUN_RE.fullmatch(p.name.split('_')[0])})
    if not runs:
        parser.error('Provide --accessions, or place accession-named FASTQs in the species folders')
    selected = []
    counts = {}
    failures = 0
    for run in runs:
        row = dict.fromkeys(COLUMNS, '')
        row.update(run_accession=run, donor='unknown', brain_region='unknown',
                   strandedness='unknown', sex='unknown', age='unknown',
                   metadata_status='failed', download_status='not_checked')
        try:
            record = report(run, cache)
            species = SPECIES.get(record['scientific_name'])
            if not species:
                raise ValueError('Unexpected species: ' + record['scientific_name'])
            if args.species != 'all' and species != args.species:
                continue
            if args.max_runs and counts.get(species, 0) >= args.max_runs:
                continue
            counts[species] = counts.get(species, 0) + 1
            row.update({k: v for k, v in record.items() if k in COLUMNS})
            row['species'] = species
            row['metadata_status'] = 'run_record_only'
            try:
                attrs = attributes(record['sample_accession'], cache)
                row.update(donor=first(attrs, 'donor id', 'donor_id', 'donor', 'individual id', 'individual'),
                           brain_region=first(attrs, 'brain region', 'brain_region', 'brain area', 'brain_area'),
                           sex=first(attrs, 'sex'), age=first(attrs, 'age'),
                           strandedness=first(attrs, 'strandedness', 'library strandedness', 'library_strandedness'),
                           sample_attributes_json=json.dumps(attrs, sort_keys=True),
                           metadata_status='sample_record_loaded')
            except Exception as error:
                row['error'] = f'Sample attributes unavailable: {error}'
            pair = [root / species / f'{run}_{i}.fastq.gz' for i in (1, 2)]
            row.update(fastq_1=str(pair[0]), fastq_2=str(pair[1]))
            row['download_status'] = ('present_unverified' if all(p.exists() and p.stat().st_size for p in pair)
                                      else 'partial' if any(p.exists() for p in pair) else 'pending')
            rows[run] = row
            save_table(table, rows)  # Row is visible before FASTQ download starts.
            if record['library_layout'] != 'PAIRED' or record['library_strategy'].upper() != 'RNA-SEQ':
                raise ValueError('This alignment workflow requires paired-end RNA-Seq')
            if args.download:
                urls = record['fastq_ftp'].split(';')
                hashes = record['fastq_md5'].split(';')
                sizes = record['fastq_bytes'].split(';')
                if len(urls) != len(hashes) or len(urls) != len(sizes):
                    raise ValueError('Incomplete ENA file metadata')
                files = {url.rsplit('/', 1)[-1]: (url, checksum, int(size))
                         for url, checksum, size in zip(urls, hashes, sizes) if url}
                if set(files) != {p.name for p in pair}:
                    raise ValueError('ENA files are not exactly a standard mate pair; inspect manually')
                row['download_status'] = 'downloading'
                save_table(table, rows)
                pair[0].parent.mkdir(parents=True, exist_ok=True)
                for dest in pair:
                    url, checksum, size = files[dest.name]
                    url = re.sub(r'^(?:ftp|https?)://', '', url)
                    print(f'Downloading/verifying {species}/{dest.name}', flush=True)
                    download('https://' + url, dest, checksum, size)
                    row['download_status'] = 'mate_1_verified' if dest == pair[0] else 'verified'
                    save_table(table, rows)
            if all(p.exists() and p.stat().st_size for p in pair):
                selected.append((species, run, record['sample_accession'], record['experiment_accession']))
        except Exception as error:
            failures += 1
            row['error'] = '; '.join(filter(None, [row['error'], str(error)]))
            row['download_status'] = 'error'
            print(f'{run}: {error}', file=sys.stderr)
        rows[run] = row
        save_table(table, rows)
    plan = folder / f'alignment_plan_{args.species}.tsv'
    with plan.open('w', newline='') as handle:
        csv.writer(handle, delimiter='\t').writerows(selected)
    print(f'Sample table: {table}', flush=True)
    print('Unknown donor/region/strandedness needs study metadata or protocol verification.', flush=True)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
