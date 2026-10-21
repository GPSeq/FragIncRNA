#!/usr/bin/env python3
"""Uniformly sample FASTA records with a fixed seed; preserve their source order."""
import argparse
import gzip
import json
import pathlib
import random
import resource
import sys
import time


def open_text(path, mode):
    return gzip.open(path, mode + "t") if path.suffix == ".gz" else path.open(mode, encoding="utf-8")


def records(path):
    header = None
    sequence = []
    with open_text(path, "r") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line.rstrip("\r\n")
                sequence = []
            elif header is None:
                raise ValueError("Sequence before the first FASTA header")
            else:
                sequence.append("".join(line.split()))
    if header is not None:
        yield header, "".join(sequence)


def peak_rss_bytes():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--count", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261017)
    parser.add_argument("--metrics", type=pathlib.Path)
    args = parser.parse_args()
    if args.count < 1:
        parser.error("--count must be positive")
    if args.input.resolve() == args.output.resolve():
        parser.error("Input and output paths must differ")
    metrics_path = args.metrics or pathlib.Path(str(args.output) + ".metrics.json")
    if metrics_path.resolve() in (args.input.resolve(), args.output.resolve()):
        parser.error("Metrics path must differ from input and output")

    start = time.perf_counter()
    rng = random.Random(args.seed)
    reservoir = []
    total = 0
    for total, record in enumerate(records(args.input), start=1):
        if len(reservoir) < args.count:
            reservoir.append((total, record))
        else:
            chosen = rng.randrange(total)
            if chosen < args.count:
                reservoir[chosen] = (total, record)
    if total < args.count:
        raise ValueError(f"Requested {args.count} records, found only {total}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "w") as output:
        for _, (header, sequence) in sorted(reservoir):
            output.write(header + "\n")
            for start_pos in range(0, len(sequence), 80):
                output.write(sequence[start_pos:start_pos + 80] + "\n")
            if not sequence:
                output.write("\n")
    metrics_path.parent.mkdir(parents=True, exist_ok=True)
    report = {"source": str(args.input), "output": str(args.output), "total_records": total,
              "sampled_records": args.count, "seed": args.seed,
              "wall_seconds": time.perf_counter() - start,
              "peak_resident_bytes": peak_rss_bytes()}
    metrics_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Sampled {args.count}/{total} FASTA records; elapsed={report['wall_seconds']:.3f}s; "
          f"peak RSS={report['peak_resident_bytes']} bytes", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        raise SystemExit(f"Error: {exc}")
