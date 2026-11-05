#!/usr/bin/env python3
"""Frozen-parameter A/B trials on an explicit reference and seeded samples.

Never downloads inputs, drops OS caches, touches user outputs, or launches a
multi-genome/k grid. Copies one verified index snapshot into a NEW output root.
"""
import argparse,csv,hashlib,json,os,pathlib,platform,shutil,subprocess,time

def sha(path):
    with pathlib.Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def compare(a,b):
    x,y=rows(a),rows(b)
    if len(x)!=len(y):return [{'row_count':[len(x),len(y)]}]
    excluded={'index_manifest','manifest_sha256'}
    differences=[]
    for i,(left,right) in enumerate(zip(x,y)):
        for key in left.keys()|right.keys():
            if key not in excluded and left.get(key)!=right.get(key):
                differences.append({'row':i,'field':key,'baseline':left.get(key),'optimized':right.get(key)})
                if len(differences)>=100:return differences
    return differences

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['baseline','optimized','reference','queries100','queries5000','index_dir','output_dir']:
        p.add_argument('--'+name.replace('_','-'),type=pathlib.Path,required=True)
    p.add_argument('--reference-cache',type=pathlib.Path)
    p.add_argument('--samples',default='100,5000')
    p.add_argument('--threads',default='1,15')
    p.add_argument('--repetitions',type=int,default=2)
    p.add_argument('--timeout',type=float,default=2400)
    a=p.parse_args();samples=[int(s) for s in a.samples.split(',')];threads=[int(t) for t in a.threads.split(',')]
    if not set(samples)<={100,5000} or not samples or min(threads)<1 or a.repetitions<1 or a.timeout<=0:p.error('Invalid samples/threads/repetitions/timeout')
    for path in [a.baseline,a.optimized,a.reference,a.queries100,a.queries5000]:
        if not path.is_file():p.error('Missing input: '+str(path))
    if a.output_dir.exists():p.error('--output-dir must not exist; benchmark artifacts are immutable')
    a.output_dir.mkdir(parents=True);out=a.output_dir.resolve()
    source=a.index_dir.resolve();manifest=source/'manifest.tsv'
    if manifest.is_symlink() or not manifest.is_file():p.error('A regular authoritative manifest.tsv is required')
    raw=manifest.read_bytes();lines=raw.decode().splitlines()
    payload=b'\n'.join(raw.split(b'\n')[2:])
    if lines[0]!='FRAGLNCRNA-SHARDS1' or hashlib.sha256(payload).hexdigest()!=lines[1]:p.error('Invalid source manifest checksum')
    at=lines.index('shard_id\tglobal_begin\treal_bins\trelative_path\tsha256\tbytes\tbuild_run')
    archives=[line.split('\t') for line in lines[at+1:]];total=sum(int(row[5]) for row in archives)
    if shutil.disk_usage(out).free<total+1024**3:p.error('Insufficient disk for explicit index snapshot plus result headroom')
    print(f'Copying {total:,} archive bytes into {out / "index-snapshot"}; source files remain untouched.',flush=True)
    snapshot=out/'index-snapshot';snapshot.mkdir()
    for row in archives:
        file=row[3]
        if pathlib.Path(file).name!=file or row[4]=='pending':p.error('Unsafe or incomplete source archive set')
        src=source/file
        if src.is_symlink() or not src.is_file():p.error('Unsafe/missing shard: '+file)
        dst=snapshot/file;shutil.copyfile(src,dst)
        if dst.stat().st_size!=int(row[5]) or sha(dst)!=row[4]:p.error('Source changed or corrupt shard: '+file)
    if manifest.read_bytes()!=raw:p.error('Source manifest changed during snapshot')
    (snapshot/'manifest.tsv').write_bytes(raw)
    environment={'platform':platform.platform(),'cpu_count':os.cpu_count(),'baseline_binary_sha256':sha(a.baseline),'optimized_binary_sha256':sha(a.optimized),
        'reference_file_sha256':sha(a.reference),'query100_sha256':sha(a.queries100),'query5000_sha256':sha(a.queries5000),'manifest_sha256':sha(manifest),
        'filesystem_cache':'uncontrolled OS cache; no drop_caches or cold-cache claim','storage':subprocess.check_output(['df','-T',str(out)],text=True),
        'reference':str(a.reference.resolve()),'index_snapshot_bytes':total}
    (out/'environment.json').write_text(json.dumps(environment,indent=2)+'\n')
    records=[];failed=False
    for sample in samples:
        for thread in threads:
            golden=None
            for rep in range(1,a.repetitions+1):
                pair={}
                for label,program in [('baseline',a.baseline),('optimized',a.optimized)]:
                    name=f'{label}-{sample}-t{thread}-r{rep}';tsv=out/(name+'.tsv')
                    command=[str(program.resolve()),'--reference',str(a.reference.resolve()),'--queries',str(getattr(a,'queries'+str(sample)).resolve()),'--output',str(tsv),
                        '--mode','ibf','--scoring-method','per_bin_verified','--k','20','--fragment-size','80000','--hash-functions','3','--fpr','0.01',
                        '--strand','canonical','--ambiguity','skip','--repeat','include','--threads',str(thread),'--batch-windows','500000','--per-bin-min-hits','2',
                        '--per-bin-min-fraction','0.02','--per-bin-top-k','256','--per-bin-max-locus-gap','50000','--basic-threshold','0.75','--strict-threshold','0.85',
                        '--ibf-layout','sharded','--max-bins-per-ibf','1000','--ibf-shard-dir',str(snapshot),'--ibf-shard-storage','disk','--ibf-shard-retention','keep',
                        '--max-ibf-bin-pairs','100000000','--no-ibf-resume','--no-log-to-console']
                    if a.reference_cache:command+=['--reference-cache',str(a.reference_cache.resolve())]
                    start=time.monotonic()
                    with (out/(name+'.console')).open('w') as log:
                        try:code=subprocess.run(command,stdout=log,stderr=log,timeout=a.timeout).returncode
                        except subprocess.TimeoutExpired:code='timeout'
                    result={'name':name,'sample':sample,'threads':thread,'repetition':rep,'command':command,'returncode':code,'wall_seconds':time.monotonic()-start}
                    if code==0:
                        result['metrics']=json.loads(pathlib.Path(str(tsv)+'.metrics.json').read_text());result['output_sha256']=sha(tsv);pair[label]=tsv
                        if golden is None:golden=tsv
                        result['scientific_differences_from_first']=compare(golden,tsv)
                        failed|=bool(result['scientific_differences_from_first'])
                    else:failed=True
                    records.append(result);(out/(name+'.run.json')).write_text(json.dumps(result,indent=2)+'\n')
                    (out/'summary.json').write_text(json.dumps({'environment':environment,'runs':records,'passed':not failed},indent=2)+'\n')
                    print(name,code,result['wall_seconds'],flush=True)
                if len(pair)==2:
                    for suffix in ['.per_bin.tsv','.per_bin.tsv.candidates.tsv']:
                        if pathlib.Path(str(pair['baseline'])+suffix).read_bytes()!=pathlib.Path(str(pair['optimized'])+suffix).read_bytes():failed=True
    (out/'summary.json').write_text(json.dumps({'environment':environment,'runs':records,'passed':not failed},indent=2)+'\n')
    return int(failed)
if __name__=='__main__':raise SystemExit(main())
