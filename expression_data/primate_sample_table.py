#!/usr/bin/env python3
"""Collect ENA sample metadata and download paired RNA-seq FASTQs with progress.

Standard library only. One run accession per line in --accessions.
Keep one downloader per --read-root; do not run alongside another downloader
writing the same files. Downloads are verified by size and MD5 before renaming.
"""
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import csv
import hashlib
import io
import json
import queue
import re
import sys
import threading
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


def human_bytes(value):
    for unit in ('B', 'KB', 'MB', 'GB', 'TB'):
        if value < 1000 or unit == 'TB':
            return f'{value:.1f} {unit}'
        value /= 1000


def duration(seconds):
    if seconds is None:
        return '--'
    seconds = max(0, int(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, seconds = divmod(rest, 60)
    return f'{hours:d}h {minutes:02d}m' if hours else f'{minutes:d}m {seconds:02d}s'


class Progress:
    def __init__(self, total_bytes, total_files, interval=None, stream=None):
        self.stream = stream if stream is not None else sys.stderr
        self.tty = self.stream.isatty()
        self.interval = interval if interval is not None else (0.5 if self.tty else 10)
        self.total = total_bytes
        self.total_files = total_files
        self.done = 0
        self.files = 0
        self.last = float('-inf')
        self.active = False
        self.lock = threading.RLock()
        self.inflight = {}

    def update(self, name, received, size, started, force=False, stage=''):
        with self.lock:
            self.inflight[name] = (received, size, started)
            self._update(name, received, size, started, force, stage)

    def _update(self, name, received, size, started, force, stage):
        now = time.monotonic()
        if not force and now - self.last < self.interval:
            return
        self.last = now
        elapsed = max(now - started, 0.001)
        rate = received / elapsed
        fraction = min(received / size, 1) if size else 0
        filled = int(12 * fraction)
        bar = '=' * filled + '-' * (12 - filled)
        overall = self.done + sum(item[0] for item in self.inflight.values())
        percent = 100 * overall / self.total if self.total else 100
        eta = duration((size - received) / rate if rate else None)
        text = (f'{name} [{bar}] {fraction:6.1%} '
                f'{human_bytes(received)}/{human_bytes(size)} | '
                f'{human_bytes(rate)}/s ETA {eta} | '
                f'ALL {percent:.1f}% {human_bytes(overall)}/{human_bytes(self.total)} '
                f'({self.files}/{self.total_files} verified, {len(self.inflight)} active)')
        if stage:
            text += f' | {stage}'
        if self.tty:
            # Limit width so carriage-return updates cannot wrap into extra rows.
            import shutil
            width = shutil.get_terminal_size((160, 24)).columns
            if len(text) >= width:
                text = (f'{name} {fraction:.0%} | ALL {percent:.0f}% ({len(self.inflight)} active) | '
                        f'{human_bytes(rate)}/s ETA {eta}')
            text = text[:max(20, width - 1)]
            self.stream.write('\r' + text + '\033[K')
            self.stream.flush()
            self.active = True
        else:
            print(text, file=self.stream, flush=True)

    def message(self, text):
        with self.lock:
            if self.active:
                self.stream.write('\n')
                self.active = False
            print(text, file=self.stream, flush=True)

    def failed(self, name):
        with self.lock:
            self.inflight.pop(name, None)

    def verified(self, name, size, existing=False):
        with self.lock:
            self.inflight.pop(name, None)
            self.done += size
            self.files += 1
            self.message(f'{"Already verified" if existing else "Verified"}: {name} | '
                         f'{self.files}/{self.total_files} files | '
                         f'{human_bytes(self.done)}/{human_bytes(self.total)} total')


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
    query = urllib.parse.urlencode(dict(accession=run, result='read_run', fields=FIELDS, format='tsv'))
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


def download(url, dest, checksum, size, progress, stop=None):
    if stop is not None and stop.is_set():
        raise InterruptedError('Download stopped')
    if dest.exists():
        progress.message(f'Checking existing file: {dest.name}')
        if dest.stat().st_size == size and md5(dest) == checksum:
            progress.verified(dest.name, size, existing=True)
            return
        raise ValueError(f'Existing file differs from ENA: {dest}; inspect before retrying')
    part = dest.with_suffix(dest.suffix + '.part')
    for attempt in range(3):
        received = 0
        started = time.monotonic()
        digest = hashlib.md5()
        try:
            progress.update(dest.name, 0, size, started, force=True)
            with urllib.request.urlopen(url, timeout=120) as response, part.open('wb') as out:
                while True:
                    if stop is not None and stop.is_set():
                        raise InterruptedError('Download stopped')
                    block = response.read(1024 * 1024)
                    if not block:
                        break
                    out.write(block)
                    digest.update(block)
                    received += len(block)
                    progress.update(dest.name, received, size, started)
            progress.update(dest.name, received, size, started, force=True, stage='checking MD5')
            if received != size or digest.hexdigest() != checksum:
                raise ValueError(f'Size/MD5 mismatch: {dest.name}')
            part.replace(dest)
            progress.verified(dest.name, size)
            return
        except Exception as error:
            progress.failed(dest.name)
            if isinstance(error, InterruptedError):
                raise
            progress.message(f'{dest.name}: attempt {attempt + 1}/3 failed: {error}')
            if attempt == 2:
                raise
            progress.message('Retrying this file from the beginning; verified files are kept.')
            time.sleep(2 ** attempt)


def save_table(path, rows):
    part = path.with_suffix('.tsv.part')
    with part.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows[k] for k in sorted(rows))
    part.replace(path)


def paired_files(record, pair):
    urls, hashes, sizes = (record[k].split(';') for k in ('fastq_ftp', 'fastq_md5', 'fastq_bytes'))
    if len(urls) != len(hashes) or len(urls) != len(sizes):
        raise ValueError('Incomplete ENA file metadata')
    files = {url.rsplit('/', 1)[-1]: (url, checksum, int(size))
             for url, checksum, size in zip(urls, hashes, sizes) if url}
    if len(urls) != 2 or set(files) != {p.name for p in pair}:
        raise ValueError('ENA files are not exactly a standard mate pair; inspect manually')
    for _, checksum, size in files.values():
        if not re.fullmatch(r'[0-9a-fA-F]{32}', checksum) or size <= 0:
            raise ValueError('Invalid ENA checksum/size metadata')
    return files


def download_run(job, progress, events, stop):
    """Workers download files; only the main thread updates the shared table."""
    row, record, pair, files = job
    run = row['run_accession']
    pair[0].parent.mkdir(parents=True, exist_ok=True)
    events.put((run, 'downloading', ''))
    try:
        for dest in pair:
            url, checksum, size = files[dest.name]
            url = re.sub(r'^(?:ftp|https?)://', '', url)
            download('https://' + url, dest, checksum, size, progress, stop)
            status = 'mate_1_verified' if dest == pair[0] else 'verified'
            events.put((run, status, ''))
        return row['species'], run, record['sample_accession'], record['experiment_accession']
    except Exception as error:
        events.put((run, 'error', str(error)))
        raise


def download_jobs(jobs, workers, progress, rows, table):
    events = queue.Queue()
    selected = []
    failures = 0
    stop = threading.Event()

    def update_table():
        changed = False
        while True:
            try:
                run, status, error = events.get_nowait()
            except queue.Empty:
                break
            row = rows[run]
            row['download_status'] = status
            if error:
                row['error'] = '; '.join(filter(None, [row['error'], error]))
            changed = True
        if changed:
            save_table(table, rows)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        pending = {executor.submit(download_run, job, progress, events, stop): job[0]['run_accession']
                   for job in jobs}
        try:
            while pending:
                completed, _ = wait(pending, timeout=0.2, return_when=FIRST_COMPLETED)
                update_table()
                for future in completed:
                    run = pending.pop(future)
                    try:
                        selected.append(future.result())
                    except Exception as error:
                        failures += 1
                        # Also handles errors before a worker could send its first event.
                        rows[run]['download_status'] = 'error'
                        if str(error) not in rows[run]['error']:
                            rows[run]['error'] = '; '.join(filter(None, [rows[run]['error'], str(error)]))
                        progress.message(f'{run}: {error}')
                if completed:
                    save_table(table, rows)
        except BaseException:
            stop.set()
            for future in pending:
                future.cancel()
            raise
        update_table()
    return sorted(selected), failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--read-root', required=True)
    parser.add_argument('--species', choices=['all', *SPECIES.values()], default='all')
    parser.add_argument('--max-runs', type=int, default=0)
    parser.add_argument('--accessions', type=Path)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--workers', type=int, default=4,
                        help='Concurrent run downloads (default: 4); each worker downloads one mate at a time')
    args = parser.parse_args()
    if args.max_runs < 0:
        parser.error('--max-runs must be nonnegative')
    if args.workers < 1:
        parser.error('--workers must be at least 1')
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
            if run:
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
    jobs, selected, counts = [], [], {}
    failures = 0
    # Resolve sizes first so overall percentage has a known denominator.
    for index, run in enumerate(runs, 1):
        print(f'Metadata {index}/{len(runs)}: {run}', file=sys.stderr, flush=True)
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
            row.update(species=species, metadata_status='run_record_only')
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
            if record['library_layout'] != 'PAIRED' or record['library_strategy'].upper() != 'RNA-SEQ':
                raise ValueError('This alignment workflow requires paired-end RNA-Seq')
            files = paired_files(record, pair) if args.download else {}
            jobs.append((row, record, pair, files))
        except Exception as error:
            failures += 1
            row['error'] = '; '.join(filter(None, [row['error'], str(error)]))
            row['download_status'] = 'error'
            print(f'{run}: {error}', file=sys.stderr)
        rows[run] = row
        save_table(table, rows)
    total = sum(size for _, _, _, files in jobs for _, _, size in files.values())
    progress = Progress(total, 2 * len(jobs))
    if args.download:
        progress.message(f'Workers: {args.workers} (up to {args.workers} FASTQ files downloading at once).')
        progress.message(f'Download plan: {len(jobs)} runs, {2 * len(jobs)} files, {human_bytes(total)} compressed.')
        progress.message('Overall includes all active downloads and verified files. Speed/ETA refer to the displayed file.')
        selected, download_failures = download_jobs(jobs, args.workers, progress, rows, table)
        failures += download_failures
    else:
        for row, record, pair, files in jobs:
            if all(p.exists() and p.stat().st_size for p in pair):
                selected.append((row['species'], row['run_accession'], record['sample_accession'], record['experiment_accession']))
    with (folder / f'alignment_plan_{args.species}.tsv').open('w', newline='') as handle:
        csv.writer(handle, delimiter='\t').writerows(selected)
    progress.message(f'Finished with {failures} run error(s). Sample table: {table}')
    if args.download:
        progress.message(f'Verified {progress.files}/{progress.total_files} files, '
                         f'{human_bytes(progress.done)}/{human_bytes(progress.total)}.')
    print('Unknown donor/region/strandedness needs study metadata or protocol verification.', flush=True)
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
