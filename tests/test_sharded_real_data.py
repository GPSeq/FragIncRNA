"""Optional real-data acceptance pilot: 100 queries against a contiguous 25-Mb prefix.

Only registered when LNCRNA_MERS_REAL_DATA_TESTS=ON. Never launches a full grid.
A fresh timestamped work directory preserves every preceding result.
"""
import argparse
import json
import pathlib
import subprocess
import sys
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--program', type=pathlib.Path, required=True)
parser.add_argument('--reference', type=pathlib.Path, required=True)
parser.add_argument('--queries', type=pathlib.Path, required=True)
parser.add_argument('--work-dir', type=pathlib.Path, required=True)
args = parser.parse_args()
if not args.reference.is_file() or not args.queries.is_file():
    print('SKIP: real reference/query FASTA unavailable')
    raise SystemExit(77)
args.work_dir.mkdir(parents=True,exist_ok=True)
root = pathlib.Path(tempfile.mkdtemp(prefix='pilot-',dir=args.work_dir))
repo = pathlib.Path(__file__).resolve().parents[1]
subprocess.run([sys.executable,str(repo/'scripts/prepare_sharded_pilot.py'),
    '--reference',str(args.reference),'--queries',str(args.queries),'--output-dir',str(root/'inputs'),
    '--counts','100','--prefix-bases','25000000'],check=True)
subprocess.run([sys.executable,str(repo/'scripts/benchmark_sharded_ibf.py'),'--program',str(args.program),
    '--reference',str(root/'inputs/reference-prefix.fa'),'--queries',str(root/'inputs/queries100.fa'),
    '--output-dir',str(root/'benchmark'),'--caps','32','250','1000','5000','--threads','15','--timeout','180'],check=True)
report = json.loads((root/'benchmark/benchmark.json').read_text())
assert all(run['returncode']==0 for run in report['runs'])
assert all(value['scientific_field_differences']==0 for value in report['layout_comparisons'].values())
for run in report['runs']:
    if run['name'].startswith('Validate'):
        assert run['metrics']['false_negative_bin_pairs'] == 0
print('Optional bounded real-data evidence:',root)
