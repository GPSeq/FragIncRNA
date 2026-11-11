#!/usr/bin/env python3
"""Prepare read-only-input, order-preserving query pilots and one contiguous reference prefix."""
import argparse
import gzip
import hashlib
import json
import pathlib
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=pathlib.Path, required=True)
    parser.add_argument('--queries', type=pathlib.Path, required=True)
    parser.add_argument('--output-dir', type=pathlib.Path, required=True)
    parser.add_argument('--prefix-bases', type=int, default=25000000)
    parser.add_argument('--counts', nargs='+', type=int, default=[100])
    parser.add_argument('--seed', type=int, default=20261021)
    args = parser.parse_args()
    if args.prefix_bases < 1 or any(n < 1 or n > 5000 for n in args.counts):
        parser.error('Positive prefix length and sample counts in 1..5000 required')
    for path in [args.reference,args.queries]:
        if not path.is_file(): parser.error(f'Missing input: {path}')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    prefix = args.output_dir/'reference-prefix.fa'
    destinations = [prefix, prefix.with_suffix('.provenance.json')]
    destinations += [args.output_dir/f'queries{n}.fa' for n in args.counts]
    if any(p.exists() for p in destinations):
        parser.error('Pilot output already exists; choose a new output directory')
    if any(p.resolve() in {args.reference.resolve(),args.queries.resolve()} for p in destinations):
        parser.error('Outputs must differ from inputs')
    opener = gzip.open if args.reference.suffix == '.gz' else open
    count = 0
    with opener(args.reference,'rt') as source, prefix.open('w') as output:
        header = next(source).strip()
        if not header.startswith('>'): raise ValueError('Reference must begin with a FASTA header')
        contig = header[1:].split()[0]
        # Keep a simple original contig identifier; coordinates start at zero.
        output.write('>'+contig+'\n')
        for line in source:
            if line.startswith('>'): break
            chunk = ''.join(line.split())[:args.prefix_bases-count]
            output.write(chunk+'\n'); count += len(chunk)
            if count == args.prefix_bases: break
    metadata = {'source':str(args.reference.resolve()),'contig':contig,'start':0,'end':count,
        'sampling':'one contiguous first-record prefix; original case, Ns and base order retained; not a whole-genome benchmark',
        'prefix_sha256':hashlib.sha256(prefix.read_bytes()).hexdigest()}
    prefix.with_suffix('.provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
    sampler = pathlib.Path(__file__).with_name('sample_fasta.py')
    for n in args.counts:
        subprocess.run([sys.executable,str(sampler),'--input',str(args.queries),'--output',
            str(args.output_dir/f'queries{n}.fa'),'--count',str(n),'--seed',str(args.seed)],check=True)

if __name__ == '__main__': main()
