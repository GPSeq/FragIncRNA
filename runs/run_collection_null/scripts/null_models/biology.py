"""Coordinate, composition, and randomization primitives (no file side effects)."""
from __future__ import annotations

import bisect
import hashlib
import random
import re
from array import array
from collections import Counter, defaultdict
from dataclasses import dataclass, field

ACGT = re.compile('[ACGT]+')
COMPLEMENT = str.maketrans('ACGTRYKMSWBDHVNacgtrykmswbdhvn',
                           'TGCAYRMKSWVHDBNtgcayrmkswvhdbn')


def stable_seed(seed: int, *namespace: object) -> int:
    payload = '\0'.join(map(str, (seed, *namespace))).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:16], 'big')


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f'Invalid half-open interval [{start}, {end})')
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def coverage(intervals: list[tuple[int, int]]) -> int:
    return sum(b - a for a, b in union(intervals))


def subtract(length: int, forbidden: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result = []
    cursor = 0
    for start, end in union(forbidden):
        if end > length:
            raise ValueError(f'Forbidden interval {start}:{end} exceeds chromosome {length}')
        if start > cursor:
            result.append((cursor, start))
        cursor = end
    if cursor < length:
        result.append((cursor, length))
    return result


def composition(sequence: str) -> tuple[float, float, float]:
    counts = Counter(sequence.upper())
    canonical = sum(counts[x] for x in 'ACGT')
    return ((counts['G'] + counts['C']) / canonical if canonical else 0.0,
            counts['N'] / len(sequence) if sequence else 0.0,
            1 - canonical / len(sequence) if sequence else 1.0)


def valid_windows(sequence: str, k: int) -> int:
    return sum(max(0, m.end() - m.start() - k + 1) for m in ACGT.finditer(sequence.upper()))


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    gene_symbol: str
    biotype: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)
    gene_biotype: str = ''

    def __post_init__(self) -> None:
        self.exons = [tuple(exon) for exon in self.exons]

    def ordered_exons(self) -> list[tuple[int, int]]:
        ordered = sorted(self.exons, reverse=self.strand == '-')
        if not ordered or len(set(ordered)) != len(ordered):
            raise ValueError(f'{self.transcript_id}: missing or duplicate exons')
        ascending = sorted(ordered)
        if any(b > c for (_, b), (c, _) in zip(ascending, ascending[1:])):
            raise ValueError(f'{self.transcript_id}: overlapping exons')
        return ordered


class RepeatIndex:
    """Compact sorted interval arrays, with prefix maxima for nested overlaps."""

    def __init__(self) -> None:
        self.data: dict[str, tuple[array, array, array, array]] = {}
        self.classes: list[str] = []
        self.class_ids: dict[str, int] = {}

    def add(self, chrom: str, start: int, end: int, repeat_class: str) -> None:
        if repeat_class not in self.class_ids:
            self.class_ids[repeat_class] = len(self.classes)
            self.classes.append(repeat_class)
        if chrom not in self.data:
            self.data[chrom] = tuple(array('q') for _ in range(4))
        starts, ends, labels, maxima = self.data[chrom]
        starts.append(start)
        ends.append(end)
        labels.append(self.class_ids[repeat_class])

    def finish(self) -> None:
        for chrom, (starts, ends, labels, maxima) in list(self.data.items()):
            if any(a > b for a, b in zip(starts, starts[1:])):
                order = sorted(range(len(starts)), key=starts.__getitem__)
                starts, ends, labels = (array('q', (values[i] for i in order))
                                        for values in (starts, ends, labels))
            maximum = 0
            for end in ends:
                maximum = max(maximum, end)
                maxima.append(maximum)
            self.data[chrom] = starts, ends, labels, maxima

    def overlaps(self, chrom: str, start: int, end: int) -> list[tuple[int, int, str]]:
        if chrom not in self.data:
            return []
        starts, ends, labels, maxima = self.data[chrom]
        lo = bisect.bisect_right(maxima, start)
        hi = bisect.bisect_left(starts, end)
        return [(max(start, starts[i]), min(end, ends[i]), self.classes[labels[i]])
                for i in range(lo, hi) if ends[i] > start]


def project_repeats(transcript: Transcript, repeats: RepeatIndex) -> dict[str, list[tuple[int, int]]]:
    projected: dict[str, list[tuple[int, int]]] = defaultdict(list)
    offset = 0
    for start, end in transcript.ordered_exons():
        for a, b, label in repeats.overlaps(transcript.chrom, start, end):
            left, right = ((a - start, b - start) if transcript.strand == '+'
                           else (end - b, end - a))
            projected[label].append((offset + left, offset + right))
        offset += end - start
    return {label: union(intervals) for label, intervals in sorted(projected.items())}


def mask_sequence(sequence: str, intervals: list[tuple[int, int]]) -> str:
    masked = list(sequence)
    for start, end in union(intervals):
        if end > len(sequence):
            raise ValueError('Transcript repeat interval exceeds sequence length')
        masked[start:end] = 'N' * (end - start)
    return ''.join(masked)


def dinucleotide_counts(sequence: str) -> Counter:
    return Counter(a + b for a, b in zip(sequence, sequence[1:]))


def shuffle_run(sequence: str, rng: random.Random) -> str:
    """Wilson last-exit tree + randomized Euler trail; edges are dinucleotides.

    The last outgoing edge of every vertex except the original final base is
    constrained to a random tree toward that base. This avoids stranded edges
    and samples Euler trails using the Altschul/Erickson-Wilson construction.
    """
    if len(sequence) < 3:
        return sequence
    outgoing: dict[str, list[str]] = defaultdict(list)
    for a, b in zip(sequence, sequence[1:]):
        outgoing[a].append(b)
    root = sequence[-1]
    in_tree = {root}
    successor: dict[str, str] = {}
    for vertex in sorted(outgoing):
        current = vertex
        while current not in in_tree:
            successor[current] = rng.choice(outgoing[current])
            current = successor[current]
        current = vertex
        while current not in in_tree:
            in_tree.add(current)
            current = successor[current]
    for vertex in sorted(outgoing):
        edges = outgoing[vertex]
        if vertex != root:
            last = successor[vertex]
            edges.remove(last)
            rng.shuffle(edges)
            edges.append(last)
        else:
            rng.shuffle(edges)
    positions = Counter()
    result = [sequence[0]]
    for _ in range(len(sequence) - 1):
        current = result[-1]
        result.append(outgoing[current][positions[current]])
        positions[current] += 1
    return ''.join(result)


def shuffle_sequence(sequence: str, rng: random.Random) -> str:
    """Preserve case at each position, every ambiguous character, and run ends."""
    parts = list(sequence)
    for match in ACGT.finditer(sequence.upper()):
        original = match.group()
        shuffled = shuffle_run(original, rng)
        if (len(shuffled) != len(original) or Counter(shuffled) != Counter(original)
                or dinucleotide_counts(shuffled) != dinucleotide_counts(original)):
            raise ValueError('Internal error: dinucleotide shuffle invariant failed')
        for index, base in enumerate(shuffled, match.start()):
            parts[index] = base.lower() if sequence[index].islower() else base
    result = ''.join(parts)
    if any(a != b for a, b in zip(sequence, result) if a.upper() not in 'ACGT'):
        raise ValueError('Internal error: ambiguous position moved')
    return result


def shuffle_parent(task: tuple[str, str, int, int, int]) -> list[tuple[str, dict]]:
    parent_id, sequence, seed, replicates, retries = task
    seen = {sequence.upper()}
    output = []
    for replicate in range(1, replicates + 1):
        namespace = stable_seed(seed, 'dinuc', parent_id, replicate)
        rng = random.Random(namespace)
        unchanged_attempts = duplicate_attempts = 0
        for attempt in range(1, retries + 1):
            shuffled = shuffle_sequence(sequence, rng)
            unchanged_attempts += shuffled.upper() == sequence.upper()
            duplicate_attempts += shuffled.upper() in seen
            if shuffled.upper() not in seen:
                break
        duplicate = shuffled.upper() in seen
        seen.add(shuffled.upper())
        output.append((shuffled, dict(replicate=replicate, rng_seed=namespace,
                       attempts=attempt, unchanged_attempts=unchanged_attempts,
                       duplicate_attempts=duplicate_attempts, unavoidable_duplicate=duplicate,
                       unchanged=shuffled.upper() == sequence.upper())))
    return output
