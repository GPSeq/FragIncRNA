#!/usr/bin/env python3
"""Optional single-reference comparison; never scans a directory or launches a genome/k grid.

Inputs are read-only. Use sample_fasta.py for a seeded order-preserving query pilot.
Every command, worker metrics and output checksum is saved in benchmark.json.
"""
import argparse
import csv
import hashlib
import json
import pathlib
import subprocess
import time


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--program', type=pathlib.Path, required=True)
    parser.add_argument('--reference', type=pathlib.Path, required=True)
    parser.add_argument('--queries', type=pathlib.Path, required=True)
    parser.add_argument('--output-dir', type=pathlib.Path, required=True)
    parser.add_argument('--caps', type=int, nargs='+', default=[250,1000,5000])
    parser.add_argument('--k', type=int, default=20)
    parser.add_argument('--fragment-size', type=int, default=80000)
    parser.add_argument('--threads', type=int, default=15)
    parser.add_argument('--timeout', type=int, default=1800, help='Hard seconds per worker (not a full-grid time limit)')
    args = parser.parse_args()
    for path in [args.reference, args.queries]:
        if not path.is_file(): parser.error(f'Missing FASTA: {path}')
    if any(cap <= 0 for cap in args.caps): parser.error('Caps must be positive')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    if any(args.output_dir.iterdir()):
        parser.error('Benchmark destination must be empty; choose a new output directory')
    report = {'reference':str(args.reference.resolve()),'queries':str(args.queries.resolve()),
              'k':args.k,'fragment_size':args.fragment_size,'runs':[],
              'interpretation':'candidate-local and global-exact recovery are separate quantities; sharding itself does not lower FPR'}
    report_path = args.output_dir/'benchmark.json'
    if report_path.exists(): parser.error('Existing benchmark.json; choose a new output directory')
    common = [str(args.program.resolve()),'--reference',str(args.reference.resolve()),
        '--queries',str(args.queries.resolve()),'--k',str(args.k),'--fragment-size',str(args.fragment_size),
        '--hash-functions','3','--fpr','0.01','--strand','canonical','--repeat','include',
        '--threads',str(args.threads),'--batch-windows','2000000','--progress-every-queries','5000',
        '--per-bin-min-hits','2','--per-bin-min-fraction','0.05','--per-bin-top-k','64',
        '--per-bin-max-locus-gap','50000','--ibf-auto-tune','recommend']
    runs = [('A_single','ibf','per_bin_ibf',None,False),('B_single','ibf','per_bin_verified',None,True),
            ('Truth','exact','global_exact',None,False)]
    for cap in args.caps:
        runs += [(f'C_{cap}','ibf','per_bin_ibf',cap,False),(f'D_{cap}','ibf','per_bin_verified',cap,False),
                 (f'D_{cap}_reuse','ibf','per_bin_verified',cap,True),(f'E_{cap}','ibf','global_ibf',cap,True),
                 (f'Validate_{cap}','validate','per_bin_verified',cap,True)]
    for name,mode,scoring,cap,reuse in runs:
        out = args.output_dir/(name+'.tsv')
        cmd = [*common,'--output',str(out),'--mode',mode,'--scoring-method',scoring,
               '--reference-cache',str(args.output_dir/'reference.bases')]
        if mode == 'exact': cmd += ['--ibf-layout','single']
        elif cap is None:
            cmd += ['--ibf-layout','single','--load-index' if reuse else '--save-index',str(args.output_dir/'single.ibf')]
        else:
            cmd += ['--ibf-layout','sharded','--max-bins-per-ibf',str(cap),
                    '--ibf-shard-dir',str(args.output_dir/f'index_{cap}'),
                    '--ibf-shard-reuse' if reuse else '--no-ibf-shard-reuse']
        started = time.monotonic()
        with (args.output_dir/(name+'.stderr')).open('w') as stderr:
            try:
                p = subprocess.run(cmd,stdout=stderr,stderr=stderr,timeout=args.timeout)
                code = p.returncode
            except subprocess.TimeoutExpired:
                code = 'timeout'
        item = {'name':name,'command':cmd,'returncode':code,'wall_seconds':time.monotonic()-started}
        if code == 0:
            item['metrics'] = json.loads(pathlib.Path(str(out)+'.metrics.json').read_text())
            item['output_sha256'] = hashlib.sha256(out.read_bytes()).hexdigest()
        report['runs'].append(item)
        report_path.write_text(json.dumps(report,indent=2)+'\n')
        print(name,code,item['wall_seconds'],flush=True)
        if code != 0:
            print('Stopping after worker failure; completed evidence and partial indexes are preserved.',flush=True)
            break
    successful = {r['name']:r for r in report['runs'] if r['returncode'] == 0}
    comparisons = {}
    if 'Truth' in successful:
        truth = rows(args.output_dir/'Truth.tsv')
        for name in successful:
            data = rows(args.output_dir/(name+'.tsv'))
            if [r['query_id'] for r in data] != [r['query_id'] for r in truth]:
                raise ValueError('Query identity/order changed across variants')
            confusion = {}
            for observed, exact in zip(data,truth):
                key = exact['status']+' -> '+observed['status']; confusion[key] = confusion.get(key,0)+1
            comparisons[name] = {'classification_confusion_vs_global_exact':confusion,
                'recovery_different_transcripts':sum(a['recovery_fraction'] != b['recovery_fraction'] for a,b in zip(data,truth)),
                'matched_window_sum':sum(int(r['matched_windows']) for r in data)}
    layout_comparisons = {}
    ignored = {'ibf_layout','shard_count','max_bins_per_ibf','global_bin_count','index_manifest','manifest_sha256'}
    for name in successful:
        baseline = 'A_single' if name.startswith('C_') else 'B_single' if name.startswith('D_') else None
        if baseline and baseline in successful:
            left, right = rows(args.output_dir/(baseline+'.tsv')), rows(args.output_dir/(name+'.tsv'))
            diffs = [(a['query_id'],key) for a,b in zip(left,right) for key in a
                     if key not in ignored and a[key] != b.get(key)]
            layout_comparisons[name] = {'baseline':baseline,'scientific_field_differences':len(diffs),
                'first_differences':diffs[:10]}
    report['layout_comparisons'] = layout_comparisons
    report['comparisons'] = comparisons
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    return 0 if all(run['returncode'] == 0 for run in report['runs']) else 1

if __name__ == '__main__': raise SystemExit(main())
