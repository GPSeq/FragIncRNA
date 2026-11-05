#!/usr/bin/env python3
"""Collect only completed frozen production trials; preserve failures explicitly."""
import argparse,csv,gzip,hashlib,itertools,json,pathlib,shutil

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def science_diff(a,b):
    count=0
    with a.open(newline='') as x,b.open(newline='') as y:
        for left,right in itertools.zip_longest(csv.DictReader(x,delimiter='\t'),csv.DictReader(y,delimiter='\t')):
            if left is None or right is None:return -1
            count+=sum(left.get(k)!=right.get(k) for k in left.keys()|right.keys() if k not in {'index_manifest','manifest_sha256'})
    return count

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work-dir',type=pathlib.Path,required=True);p.add_argument('--evidence-dir',type=pathlib.Path,required=True);a=p.parse_args()
    root=a.work_dir;dest=a.evidence_dir;dest.mkdir(parents=True,exist_ok=True);records=[];comparisons=[];oracle=[]
    for path in sorted(root.glob('*-t*.tsv.run.json')):
        data=json.loads(path.read_text());name=path.name.removesuffix('.tsv.run.json');label=name.split('-')[0]
        if label not in {'baseline','optimized','final'}:continue
        shutil.copy2(path,dest/path.name);metrics=root/(name+'.tsv.metrics.json')
        if data['returncode']==0:
            data['metrics']=json.loads(metrics.read_text());shutil.copy2(metrics,dest/metrics.name)
            for ext in ['.tsv.log','.tsv.log.jsonl']:
                f=root/(name+ext)
                if f.exists():shutil.copy2(f,dest/f.name)
        data['variant']=label;records.append(data)
    for label in ['optimized','final']:
        for sample in [100,5000]:
            for thread in [1,15]:
                base=root/f'baseline-{sample}-t{thread}.tsv';new=root/f'{label}-{sample}-t{thread}.tsv'
                if not all(any(x['variant']==v and x['sample']==sample and x['threads']==thread and x['returncode']==0 for x in records) for v in ['baseline',label]):continue
                item={'variant':label,'sample':sample,'threads':thread,'scientific_field_differences':science_diff(base,new),'sidecars':{}}
                for suffix in ['.per_bin.tsv','.per_bin.tsv.candidates.tsv']:
                    x=pathlib.Path(str(base)+suffix);y=pathlib.Path(str(new)+suffix)
                    item['sidecars'][suffix]={'baseline_sha256':sha(x),'optimized_sha256':sha(y),'identical':sha(x)==sha(y)}
                comparisons.append(item)
    for label in ['optimized','final']:
        for sample in [100,5000]:
            if not all(any(x['variant']==label and x['sample']==sample and x['threads']==t and x['returncode']==0 for x in records) for t in [1,15]):continue
            left=root/f'{label}-{sample}-t1.tsv';right=root/f'{label}-{sample}-t15.tsv'
            item={'variant':label+'_thread_reproducibility','comparison':label+' t1 vs t15','sample':sample,'threads':'1/15','scientific_field_differences':science_diff(left,right),'main_tsv_byte_identical':sha(left)==sha(right),'sidecars':{}}
            for suffix in ['.per_bin.tsv','.per_bin.tsv.candidates.tsv']:
                x=pathlib.Path(str(left)+suffix);y=pathlib.Path(str(right)+suffix)
                item['sidecars'][suffix]={'baseline_sha256':sha(x),'optimized_sha256':sha(y),'identical':sha(x)==sha(y)}
            comparisons.append(item)
    old=root/'oracle-windows.tsv'
    if old.exists():
        for thread in [1,15]:
            new=root/f'final-100-t{thread}.windows.tsv.gz'
            if not new.exists():continue
            differences=0;rows=0
            with old.open() as f,gzip.open(new,'rt') as g:
                for expected,actual in itertools.zip_longest(csv.DictReader(f,delimiter='\t'),csv.DictReader(g,delimiter='\t')):
                    if expected is None or actual is None:raise RuntimeError('Oracle window count mismatch')
                    rows+=1
                    for k,v in expected.items():
                        if expected['valid']=='0' and k in {'ibf_any_bin_positive','ibf_candidate_positive','exact_verified_in_candidates'}:continue
                        differences+=actual[k]!=v
            oracle.append({'threads':thread,'window_rows':rows,'predicate_differences':differences})
            shutil.copy2(new,dest/new.name)
    summary={'runs':records,'comparisons':comparisons,'frozen_window_oracle':oracle}
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    text=['# Completed measurements','', 'Wall times are non-isolated observations on native ext4 with uncontrolled OS caches. Optimized is the initial 64-KiB-page trial; final adds the exact-key rejection filter and 16-KiB pages. Different optimized versions are not repeated measurements of one fixed executable.','',
          '| Variant | Queries | Threads requested/used | Wall s | User+sys CPU s | Peak RSS MiB | Query loads | SHA checks | Candidate verify s | Verifier read GiB |',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for x in records:
        if x['returncode']!=0:
            text+=['',f"{x['variant']} {x['sample']} queries t{x['threads']}: **{x['returncode']}**, {">=" if x.get('wall_seconds_is_lower_bound') else ""}{x['wall_seconds']:.1f} s; incomplete output excluded."]
            continue
        m=x['metrics'];read=m.get('candidate_reference_userspace_read_bytes');cpu=m.get('user_cpu_seconds',0)+m.get('system_cpu_seconds',0)
        text.append(f"| {x['variant']} | {x['sample']} | {x['threads']}/{m['threads']} | {x['wall_seconds']:.2f} | {cpu:.2f} | {m['peak_resident_bytes']/1048576:.1f} | {m.get('shard_query_loads','NA')} | {m.get('shard_checksum_count','not instrumented')} | {m['candidate_verification_seconds']:.2f} | {read/1073741824:.3f} |" if read is not None else f"| {x['variant']} | {x['sample']} | {x['threads']}/{m['threads']} | {x['wall_seconds']:.2f} | {cpu:.2f} | {m['peak_resident_bytes']/1048576:.1f} | {m.get('shard_query_loads','NA')} | not instrumented | {m['candidate_verification_seconds']:.2f} | not instrumented |")
    text+=['','| Comparison | Queries | Threads | Scientific differences | Candidate/per-bin sidecars identical |','|---|---:|---:|---:|---|']
    for x in comparisons:text.append(f"| {x.get('comparison', 'baseline vs '+x['variant'])} | {x['sample']} | {x['threads']} | {x['scientific_field_differences']} | {all(y['identical'] for y in x['sidecars'].values())} |")
    for x in oracle:text+=['',f"Frozen-source oracle vs final, t{x['threads']}: {x['window_rows']:,} original window rows, {x['predicate_differences']} predicate/index differences."]
    (dest/'measured-results.md').write_text('\n'.join(text)+'\n')
    prefix=root/'prefix-validation/summary.json'
    if prefix.exists():shutil.copy2(prefix,dest/'prefix-validation.json')
    return int(any(x['scientific_field_differences'] or not all(y['identical'] for y in x['sidecars'].values()) for x in comparisons) or any(x['predicate_differences'] for x in oracle))
if __name__=='__main__':raise SystemExit(main())
