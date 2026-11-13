"""Strict input readers and bounded-memory, disk-backed human reference access."""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO
from urllib.parse import unquote

from biology import RepeatIndex, Transcript, reverse_complement, union

INPUTS = {
    'lncrna_fasta': 'gencode.v49.lncRNA_transcripts.fa',
    'lncrna_gtf': 'gencode.v49.long_noncoding_RNAs.gtf',
    'annotation_gtf': 'gencode.v49.annotation.gtf',
    'coding_fasta': 'gencode.v49.pc_transcripts.fa',
    'genome_fasta': 'GRCh38.primary_assembly.genome.fa',
    'repeats_bed': 'GRCh38.repeats.bed',
}
ATTRIBUTE_TOKEN = re.compile(
    r'\s*([^;\s=]+)(?:\s*=\s*([^;]*)|\s+("(?:\\.|[^"\\])*"|[^;]*))\s*(?:;|$)')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_input(directory: Path, name: str, explicit: Path | None = None) -> Path:
    if explicit is not None:
        if not explicit.is_file():
            raise ValueError(f'Input override does not exist: {explicit}')
        return explicit.resolve()
    matches = [directory / (name + suffix) for suffix in ('', '.gz', '.bgz', '.bgzf')
               if (directory / (name + suffix)).is_file()]
    if len(matches) != 1:
        candidates = sorted(str(p) for p in directory.glob(name.split('.')[0] + '*'))
        raise ValueError(f'Cannot unambiguously resolve {name}: exact matches={matches}; '
                         f'candidates={candidates}. Supply the explicit input override.')
    return matches[0].resolve()


@contextmanager
def open_text(path: Path) -> Iterator[TextIO]:
    with path.open('rb') as raw:
        magic = raw.read(2)
        raw.seek(0)
        stream = gzip.GzipFile(fileobj=raw) if magic == b'\x1f\x8b' else raw
        with io.TextIOWrapper(stream, encoding='utf-8') as handle:
            try:
                yield handle
            except (EOFError, gzip.BadGzipFile, UnicodeError, OSError) as error:
                raise ValueError(f'Unreadable/corrupt input {path}: {error}') from error


def fasta_records(path: Path, canonical_ids: bool = True, allow_empty: bool = False
                  ) -> Iterator[tuple[str, str, str]]:
    seen: set[str] = set()
    header: str | None = None
    parts: list[str] = []

    def record() -> tuple[str, str, str]:
        assert header is not None
        identifier = header.split()[0]
        if canonical_ids:
            identifier = identifier.split('|')[0]
        if identifier in seen:
            raise ValueError(f'{path}: duplicate FASTA identifier {identifier}')
        seen.add(identifier)
        sequence = ''.join(parts)
        if not sequence or re.search('[^ACGTRYKMSWBDHVNacgtrykmswbdhvn]', sequence):
            raise ValueError(f'{path}: empty or invalid IUPAC DNA in {identifier}')
        fields = header.split()[0].split('|')
        if len(fields) > 6 and fields[6].isdigit() and int(fields[6]) != len(sequence):
            raise ValueError(f'{identifier}: FASTA header length disagrees with sequence')
        return identifier, header, sequence

    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield record()
                header, parts = line[1:], []
                if not header or header[0].isspace():
                    raise ValueError(f'{path}:{number}: malformed FASTA header')
            elif header is None:
                raise ValueError(f'{path}:{number}: sequence before FASTA header')
            else:
                parts.append(''.join(line.split()))
        if header is not None:
            yield record()
        elif not allow_empty:
            raise ValueError(f'{path}: no FASTA records')


def attributes(text: str) -> dict[str, str]:
    """Quoted GTF attributes (including embedded semicolons) or GFF3 key=value."""
    result = {}
    position = 0
    for match in ATTRIBUTE_TOKEN.finditer(text):
        if text[position:match.start()].strip(' ;\t\r\n'):
            raise ValueError(f'Malformed/Unterminated GTF attribute: {text}')
        position = match.end()
        key, gff_value, gtf_value = match.groups()
        if gff_value is not None:
            value = unquote(gff_value.strip())
        else:
            value = gtf_value.strip()
            if value.startswith('"'):
                if not value.endswith('"') or len(value) < 2:
                    raise ValueError(f'Unterminated GTF attribute: {text}')
                value = value[1:-1].replace('\\"', '"').replace('\\\\', '\\')
        if key in result and result[key] != value and key in {'gene_id', 'transcript_id'}:
            raise ValueError(f'Conflicting {key} attributes')
        result[key] = value
    if text[position:].strip(' ;\t\r\n'):
        raise ValueError(f'Malformed/Unterminated GTF attribute: {text}')
    return result


def canonical_chrom(name: str) -> str:
    value = name[3:] if name.startswith('chr') else name
    return 'MT' if value in {'M', 'MT'} else value


def chromosome_mapping(names: list[str]) -> dict[str, str]:
    result = {}
    for name in names:
        key = canonical_chrom(name)
        if key in result:
            raise ValueError(f'Ambiguous chromosome aliases: {result[key]}, {name}')
        result[key] = name
    return result


class Genome:
    """A whitespace-free byte file plus offsets, never a whole-genome string."""

    def __init__(self, source: Path, cache: Path, source_hash: str) -> None:
        cache.mkdir(parents=True, exist_ok=True)
        self.raw_path = cache / 'genome.bases'
        index_path = cache / 'genome.json'
        if index_path.exists() and self.raw_path.exists():
            info = json.loads(index_path.read_text())
            if (info.get('sha256') != source_hash
                    or self.raw_path.stat().st_size != info.get('size')
                    or sha256(self.raw_path) != info.get('cache_sha256')):
                raise ValueError(f'Incompatible/damaged genome cache at {cache}; use a new cache')
        else:
            temporary = cache / 'genome.bases.partial'
            info = {'sha256': source_hash, 'chromosomes': {}, 'gaps': {}}
            try:
                with open_text(source) as handle, temporary.open('wb') as output:
                    chrom = None
                    for line in handle:
                        line = line.strip()
                        if not line:
                            continue
                        if line.startswith('>'):
                            chrom = line[1:].split()[0]
                            if chrom in info['chromosomes']:
                                raise ValueError(f'Duplicate genome chromosome {chrom}')
                            info['chromosomes'][chrom] = [output.tell(), 0]
                            info['gaps'][chrom] = []
                        else:
                            if chrom is None:
                                raise ValueError('Reference sequence before FASTA header')
                            sequence = ''.join(line.split()).upper()
                            if re.search('[^ACGTRYKMSWBDHVN]', sequence):
                                raise ValueError(f'Invalid reference DNA on {chrom}')
                            length = info['chromosomes'][chrom][1]
                            gaps = info['gaps'][chrom]
                            for match in re.finditer('[^ACGT]+', sequence):
                                start, end = length + match.start(), length + match.end()
                                if gaps and gaps[-1][1] == start:
                                    gaps[-1][1] = end
                                else:
                                    gaps.append([start, end])
                            output.write(sequence.encode('ascii'))
                            info['chromosomes'][chrom][1] += len(sequence)
                    info['size'] = output.tell()
                if not info['chromosomes'] or any(n == 0 for _, n in info['chromosomes'].values()):
                    raise ValueError('Empty reference chromosome')
                temporary.replace(self.raw_path)
                info['cache_sha256'] = sha256(self.raw_path)
                index_temp = cache / 'genome.json.partial'
                index_temp.write_text(json.dumps(info, sort_keys=True))
                index_temp.replace(index_path)
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
        self.chromosomes = info['chromosomes']
        self.gaps = info['gaps']
        self.mapping = chromosome_mapping(list(self.chromosomes))
        self.handle = self.raw_path.open('rb')

    def normalize(self, chrom: str) -> str | None:
        return self.mapping.get(canonical_chrom(chrom))

    def length(self, chrom: str) -> int:
        return self.chromosomes[chrom][1]

    def fetch(self, chrom: str, start: int, end: int) -> str:
        if chrom not in self.chromosomes or not 0 <= start < end <= self.length(chrom):
            raise ValueError(f'Invalid reference coordinates {chrom}:{start}-{end}')
        self.handle.seek(self.chromosomes[chrom][0] + start)
        sequence = self.handle.read(end - start).decode('ascii')
        if len(sequence) != end - start:
            raise ValueError('Truncated reference cache')
        return sequence

    def splice(self, transcript: Transcript) -> str:
        parts = [self.fetch(transcript.chrom, a, b) for a, b in transcript.ordered_exons()]
        if transcript.strand == '-':
            parts = [reverse_complement(part) for part in parts]
        return ''.join(parts)

    def close(self) -> None:
        self.handle.close()


def read_annotation(path: Path, wanted: set[str], genome: Genome, eligible: set[str],
                    expected_version: int, all_genes: bool = False
                    ) -> tuple[dict[str, Transcript], dict[str, list[tuple[int, int]]], dict]:
    transcripts: dict[str, Transcript] = {}
    genes: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gene_ids: set[str] = set()
    gene_bounds: dict[str, tuple[str, int, int]] = {}
    transcript_bounds: dict[str, tuple[int, int]] = {}
    exon_only: dict[str, list[tuple[int, int]]] = defaultdict(list)
    counts: Counter = Counter()
    comments = []
    aliases: dict[str, str] = {}
    gff_transcripts: dict[str, str] = {}
    excluded_transcripts: set[str] = set()
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith('#'):
                comments.append(line.strip())
                assembly = re.search(r'GRCh\d+', line)
                if assembly and assembly.group() != 'GRCh38':
                    raise ValueError(f'{path}: incompatible assembly {assembly.group()}')
                version = re.search(r'version (\d+)', line)
                if version and expected_version and int(version[1]) != expected_version:
                    raise ValueError(f'{path}: incompatible annotation version {version[1]}')
                continue
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected 9 GTF/GFF3 columns')
            chrom, _, feature, left, right, _, strand, _, raw = fields
            if feature not in {'gene', 'transcript', 'mRNA', 'lnc_RNA', 'exon'}:
                continue
            a = attributes(raw)
            if feature in {'transcript', 'mRNA', 'lnc_RNA'}:
                candidate_id = a.get('transcript_id', a.get('ID', ''))
                if candidate_id in wanted and (genome.normalize(chrom) not in eligible):
                    excluded_transcripts.add(candidate_id)
            normalized = genome.normalize(chrom)
            if normalized is None:
                if canonical_chrom(chrom) in {str(x) for x in range(1, 23)} | {'X', 'Y'}:
                    raise ValueError(f'{path}: primary chromosome {chrom} missing from reference')
                counts['excluded_unknown_contig_rows'] += 1
                continue
            if normalized in aliases and aliases[normalized] != chrom:
                raise ValueError(f'{path}: mixed chromosome aliases for {normalized}')
            aliases[normalized] = chrom
            if normalized not in eligible:
                counts['excluded_chromosome_rows'] += 1
                continue
            start, end = int(left) - 1, int(right)
            if not 0 <= start < end <= genome.length(normalized):
                raise ValueError(f'{path}:{number}: invalid GTF coordinates {chrom}:{left}-{right}')
            if feature == 'gene':
                gene = a.get('gene_id', a.get('ID', ''))
                if not gene or gene in gene_ids:
                    raise ValueError(f'{path}:{number}: missing/duplicate gene identifier {gene}')
                gene_ids.add(gene)
                gene_bounds[gene] = (normalized, start, end)
                genes[normalized].append((start, end))
                continue
            tid = a.get('transcript_id', '')
            if feature in {'transcript', 'mRNA', 'lnc_RNA'}:
                tid = tid or a.get('ID', '')
                if 'ID' in a:
                    gff_transcripts[a['ID']] = tid
            elif not tid:
                parents = a.get('Parent', '').split(',')
                if len(parents) != 1:
                    raise ValueError('Multi-parent GFF3 exons require explicit transcript expansion')
                tid = gff_transcripts.get(parents[0], parents[0])
            if tid not in wanted:
                continue
            if strand not in {'+', '-'}:
                raise ValueError(f'{tid}: invalid transcript strand {strand}')
            if feature == 'exon':
                exon_only[tid].append((start, end))
                if tid in transcripts and (transcripts[tid].chrom != normalized
                                          or transcripts[tid].strand != strand):
                    raise ValueError(f'{tid}: inconsistent exon chromosome/strand')
            else:
                if tid in transcripts:
                    raise ValueError(f'{path}: duplicate transcript {tid}')
                gene = a.get('gene_id', a.get('Parent', ''))
                biotype = a.get('transcript_type', a.get('transcript_biotype',
                                a.get('gene_type', a.get('gene_biotype', ''))))
                transcripts[tid] = Transcript(tid, gene, a.get('gene_name', a.get('Name', '')),
                                              biotype, normalized, strand,
                                              gene_biotype=a.get('gene_type', a.get('gene_biotype', '')))
                transcript_bounds[tid] = (start, end)
    if all_genes and not gene_ids:
        raise ValueError(f'{path}: no comprehensive gene features; cannot define intergenic DNA')
    for tid, transcript in transcripts.items():
        transcript.exons = exon_only.get(tid, [])
        transcript.ordered_exons()
        start, end = transcript_bounds[tid]
        if any(a < start or b > end for a, b in transcript.exons):
            raise ValueError(f'{tid}: exon outside annotated transcript bounds')
        if transcript.gene_id not in gene_bounds:
            raise ValueError(f'{tid}: missing gene feature {transcript.gene_id}')
        chrom, a, b = gene_bounds[transcript.gene_id]
        if chrom != transcript.chrom or start < a or end > b:
            raise ValueError(f'{tid}: transcript outside annotated gene bounds')
    counts['genes'] = len(gene_ids)
    counts['transcripts'] = len(transcripts)
    return transcripts, {c: union(v) for c, v in genes.items()}, {
        'counts': dict(counts), 'comments': comments, 'chromosome_mapping': aliases,
        'excluded_transcripts': sorted(excluded_transcripts)}


def read_repeats(path: Path, genome: Genome, eligible: set[str]) -> tuple[RepeatIndex, dict]:
    index = RepeatIndex()
    counts: Counter = Counter()
    aliases: dict[str, str] = {}
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f'{path}:{number}: expected BED with at least three columns')
            chrom = genome.normalize(fields[0])
            if chrom is None:
                if canonical_chrom(fields[0]) in {str(x) for x in range(1, 23)} | {'X', 'Y'}:
                    raise ValueError(f'{path}: primary BED chromosome {fields[0]} missing from reference')
                counts['excluded_unknown_contig_rows'] += 1
                continue
            if chrom not in eligible:
                counts['excluded_chromosome_rows'] += 1
                continue
            if chrom in aliases and aliases[chrom] != fields[0]:
                raise ValueError(f'{path}: mixed BED chromosome aliases')
            aliases[chrom] = fields[0]
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end <= genome.length(chrom):
                raise ValueError(f'{path}:{number}: invalid BED coordinates')
            # Supplied BED: chrom/start/end/name/class/family; BED3 is unclassified.
            label = fields[4] if len(fields) >= 5 and not fields[4].isdigit() else 'unclassified'
            index.add(chrom, start, end, label)
            counts['intervals'] += 1
    index.finish()
    if not counts['intervals']:
        raise ValueError(f'{path}: no repeat intervals on eligible chromosomes')
    return index, {'counts': dict(counts), 'chromosome_mapping': aliases,
                   'classes': index.classes}
