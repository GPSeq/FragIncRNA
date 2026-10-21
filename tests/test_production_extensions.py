"""Directory transactions, original-window exports, cache and thread equivalence."""
import csv,gzip,hashlib,json,pathlib,random,subprocess,sys,tempfile,shutil,time
P=pathlib.Path(sys.argv[1]).resolve()
def table(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def run(args,good=True):
    p=subprocess.run([str(P),*map(str,args)],text=True,capture_output=True,timeout=120)
    assert (p.returncode==0)==good,(args,p.stderr)
    return p

def science(path):
    return [{k:v for k,v in row.items() if k not in {'index_manifest','manifest_sha256'}} for row in table(path)]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rev(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with tempfile.TemporaryDirectory(prefix='frag-production-ext-') as tmp:
    r=pathlib.Path(tmp); genomes=r/'genomes';genomes.mkdir()
    rng=random.Random(214);seq=''.join(rng.choice('ACGT') for _ in range(31000))
    seq=seq[:15000]+'ACGT'*150+seq[15600:]
    ref=genomes/'Assembly.full-v1.fna';ref.write_text('>chr1\n'+seq+'\n')
    q=r/'q.fa'; queries={'duplicate':('ACGT'*35)+'N'+('ACGT'*12),'forward':seq[9980:10130],
        'reverse':rev(seq[20900:21030]),'ambiguous':'NN'+seq[:100]+'NN','short':'ACG'}
    q.write_text(''.join('>'+name+'\n'+s+'\n' for name,s in queries.items()))
    common=['--reference',ref,'--queries',q,'--mode','ibf','--scoring-method','per_bin_verified',
        '--fragment-size','10000','--ibf-layout','sharded','--max-bins-per-ibf','1',
        '--per-bin-min-fraction','0.02','--per-bin-top-k','256','--max-ibf-bin-pairs','1','--no-log-to-console']
    for k in [15,18,20,22,24]:
        baseline=None
        for threads in [1,2,4,8,15]:
            out=r/f'k{k}t{threads}.tsv';win=r/f'k{k}t{threads}.windows.tsv.gz';fa=r/f'k{k}t{threads}.fa.gz'
            args=[*common,'--k',k,'--threads',threads,'--output',out,'--ibf-shard-dir',r/f'i{k}',
                '--ibf-shard-cache-mb','1','--window-output',win,'--window-include-sequence','--kmer-fasta-output',fa,'--batch-windows','300']
            run(args); rows=table(win);metrics=json.loads(pathlib.Path(str(out)+'.metrics.json').read_text())
            assert metrics['shard_checksum_count']==metrics['shard_count']
            assert metrics['shard_cache_peak_bytes']<=metrics['ibf_shard_cache_budget_bytes']
            assert metrics['whole_reference_exact_scans']==0
            assert metrics['shard_cache_evictions']>0
            if baseline is None:baseline=(science(out),rows)
            else:assert baseline==(science(out),rows),(k,threads)
            candidates=table(pathlib.Path(str(out)+'.per_bin.tsv.candidates.tsv'))
            for name,s in queries.items():
                wr=[w for w in rows if w['query_id']==name]
                assert [int(w['k_index']) for w in wr]==list(range(max(0,len(s)-k+1)))
                regions=[x for x in candidates if x['query_id']==name]
                for wi,w in enumerate(wr):
                    mer=s[wi:wi+k];valid=set(mer)<=set('ACGT')
                    assert (w['valid']=='1')==valid
                    assert w['query_start']==str(wi) and w['query_end']==str(wi+k)
                    if not valid:
                        assert w['canonical_kmer_id']=='NA' and w['kmer_sequence']=='NA';continue
                    matches=sum(mer in seq[int(c['start']):int(c['sequence_end'])] or rev(mer) in seq[int(c['start']):int(c['sequence_end'])] for c in regions)
                    assert int(w['matching_candidate_regions_count'])==matches,(k,name,wi,w,regions)
                    assert int(w['exact_verified_in_candidates'])==bool(matches)
                    code=0
                    for c in min(mer,rev(mer)):code=(code<<2)|'ACGT'.index(c)
                    assert w['canonical_kmer_id']==f'dna2bit-v1:{k}:{code:x}'
                    assert w['kmer_sequence']==mer
            with gzip.open(fa,'rt') as f:lines=f.read().splitlines()
            assert len({x for x in lines if x.startswith('>')})==len(lines)//2
        # Calibration still enforces the actually materialized pair cap.
        run([*common,'--k',k,'--mode','calibrate','--output',r/f'cal{k}.tsv','--ibf-shard-dir',r/f'i{k}'],False)
    off=r/'off.tsv';offwin=r/'off.windows.tsv.gz'
    run([*common,'--k',24,'--output',off,'--ibf-shard-dir',r/'i24','--window-output',offwin,'--window-include-sequence',
         '--no-candidate-prefilter','--ibf-shard-cache-mb',0,'--reference-page-cache-mb',0])
    assert science(off)==baseline[0] and table(offwin)==baseline[1]
    # Even exact legacy-schema exports stage the main output on a row-cap error.
    run(['--reference',ref,'--queries',q,'--mode','exact','--k',20,'--output',r/'exactlimited.tsv','--window-output',r/'exactlimited.windows.tsv','--max-window-rows',1],False)
    assert not (r/'exactlimited.tsv').exists()
    # Optional local disk copies are explicit, byte bounded and revalidated.
    local=r/'local-cache'; localout=r/'local.tsv'
    localargs=[*common,'--k',20,'--output',localout,'--ibf-shard-dir',r/'i20','--ibf-local-cache-dir',local,'--ibf-local-cache-max-mb',8]
    run(localargs);lm=json.loads(pathlib.Path(str(localout)+'.metrics.json').read_text())
    assert lm['shard_local_cache_copied_bytes']>0 and lm['shard_local_cache_checksum_count']==lm['shard_count']
    run(localargs);assert json.loads(pathlib.Path(str(localout)+'.metrics.json').read_text())['shard_local_cache_copied_bytes']==0
    original_cached_copy=next(local.glob('*/*.ibf'))
    # A second compatible subtree cannot exceed the same global cache quota.
    run([*common,'--k',18,'--output',r/'local18.tsv','--ibf-shard-dir',r/'i18','--ibf-local-cache-dir',local,'--ibf-local-cache-max-mb',6],False)
    run([*localargs,'--ibf-local-cache-max-mb',1],False)
    copy=original_cached_copy;copy.write_bytes(copy.read_bytes()[:-1]);run(localargs,False)
    # Mutations after startup verification must fail even when the object is cached.
    longq=r/'long.fa';longq.write_text(''.join(f'>long{i}\n{seq[:1500]}\n' for i in range(400)))
    for mutation in ['same_size','truncate','swap','symlink']:
        index=r/('mutate-'+mutation);shutil.copytree(r/'i20',index)
        out=r/(mutation+'.tsv');cmd=[str(P),*map(str,[*common,'--queries',longq,'--k',20,'--output',out,'--ibf-shard-dir',index,'--progress-every-queries',1,'--ibf-shard-cache-mb',8])]
        p=subprocess.Popen(cmd,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        log=pathlib.Path(str(out)+'.log');deadline=time.monotonic()+20
        while p.poll() is None and time.monotonic()<deadline:
            if log.exists() and 'bounded input batches' in log.read_text():break
            time.sleep(.002)
        assert p.poll() is None,'query finished before integrity mutation test'
        shard=index/'shard-0.ibf'
        if mutation=='same_size':
            with shard.open('r+b') as f:f.seek(-1,2);b=f.read(1);f.seek(-1,2);f.write(bytes([b[0]^1]))
        elif mutation=='truncate':shard.write_bytes(shard.read_bytes()[:-1])
        elif mutation=='swap':
            replacement=index/'replacement';shutil.copy2(shard,replacement);replacement.replace(shard)
        else:
            original=index/'original';shard.rename(original);shard.symlink_to(original)
        assert p.wait(timeout=30)!=0 and not out.exists()
        assert 'Shard changed while querying' in log.read_text()
    # Scope predicates, query-ID filtering, byte-identical masks across cache off/on.
    all_rows=baseline[1]
    for scope in ['verified','candidate','failed']:
        out=r/(scope+'.tsv');win=r/(scope+'.windows.tsv')
        run([*common,'--k',24,'--output',out,'--ibf-shard-dir',r/'i24','--window-output',win,'--window-scope',scope])
        wanted={'verified':'exact_verified_in_candidates','candidate':'ibf_candidate_positive','failed':'exact_verified_in_candidates'}[scope]
        expected=[x for x in all_rows if x['valid']=='1' and x[wanted]==('0' if scope=='failed' else '1')]
        assert [(x['query_id'],x['k_index']) for x in table(win)]==[(x['query_id'],x['k_index']) for x in expected]
    ids=r/'ids.txt';ids.write_text('duplicate\n')
    run(['--mode','export-kmers','--queries',q,'--k',20,'--window-output',r/'export.tsv.gz','--window-query-ids',ids,'--window-include-sequence'])
    erows=table(r/'export.tsv.gz');assert len(erows)==len(queries['duplicate'])-19
    assert all(x['reference_id']=='NA' and x['ibf_any_bin_positive']=='NA' and x['exact_verified_in_candidates']=='NA' for x in erows)
    # Limit errors cannot publish a completed report or main TSV.
    run([*common,'--k',20,'--output',r/'limited.tsv','--ibf-shard-dir',r/'i20','--window-output',r/'limited.windows.tsv','--max-window-rows',2],False)
    assert not (r/'limited.windows.tsv').exists() and not (r/'limited.tsv').exists()
    # Directory discovery skips files/dirs symlinks, ignores other extensions.
    with gzip.open(genomes/'Second.assembly.fa.gz','wt') as f:f.write('>different\n'+rev(seq)+'\n')
    (genomes/'ignore.txt').write_text('not fasta');(genomes/'skip.fa').symlink_to(ref)
    (genomes/'nested').mkdir();(genomes/'nested'/'ignored.fa').write_text('>x\nAAAA\n')
    outdir=r/'multi'; indexroot=r/'multi-index'
    args=['--reference-dir',genomes,'--output-dir',outdir,'--queries',q,'--mode','ibf','--scoring-method','per_bin_verified','--k','20',
        '--fragment-size','10000','--ibf-layout','sharded','--max-bins-per-ibf','2','--ibf-shard-dir',indexroot,'--no-log-to-console',
        '--window-output','windows.tsv.gz','--window-include-sequence']
    # Internal ledger names are reserved even on case-insensitive filesystems.
    run([*args,'--output-dir',r/'reserved-name','--output-name-template','completed.sha256.tsv'],False)
    run([*args,'--output-dir',r/'reserved-sidecar','--metrics','COMPLETED.sha256.tsv'],False)
    run(args)
    manifest=table(outdir/'run_manifest.tsv');assert [x['reference_id'] for x in manifest]==['Assembly.full-v1','Second.assembly']
    assert all(x['status']=='COMPLETED' for x in manifest)
    outputs=list(outdir.glob('*/k20/results.tsv'));before={str(p):sha(p) for p in outdir.glob('*/k20/*') if p.is_file()}
    run(args);assert before=={str(p):sha(p) for p in outdir.glob('*/k20/*') if p.is_file()}
    assert all(int(x['query_count'])==len(queries) for x in table(outdir/'run_manifest.tsv'))
    for p in outputs:
        metrics=json.loads(pathlib.Path(str(p)+'.metrics.json').read_text())
        assert '.k-staging-' not in json.dumps(metrics)
    # Any completed result tamper is detected, without replacing old data.
    p=outputs[0];p.write_text(p.read_text()+'tamper\n');tampered=sha(p);run(args,False);assert sha(p)==tampered
    # Ambiguous assembly IDs fail before any genome work.
    (genomes/'assembly.full-v1.fa').write_text('>x\nACGT\n')
    run([*args,'--output-dir',r/'collision'],False)
    assert not list((r/'collision').glob('*/k20/results.tsv'))
    (genomes/'assembly.full-v1.fa').unlink()
    # Continue-on-error retains completed genomes and publishes no failed assembly.
    (genomes/'Broken.fa').write_text('ACGT\n')
    failargs=[*args,'--output-dir',r/'continue','--continue-on-reference-error']
    run(failargs,False)
    statuses={x['reference_id']:x['status'] for x in table(r/'continue/run_manifest.tsv')}
    assert statuses=={'Assembly.full-v1':'COMPLETED','Broken':'FAILED','Second.assembly':'COMPLETED'}
    assert not (r/'continue/Broken/k20').exists()
    # Legacy single-reference exact directory path and TOML types.
    (genomes/'Broken.fa').unlink()
    cfg=r/'multi.toml';cfg.write_text(f'[v2]\nreference_dir="{genomes}"\noutput_dir="{r / "exact-multi"}"\nqueries="{q}"\nmode="exact"\nk=20\nwindow_output="windows.tsv.gz"\nreference_recursive=false\nibf_shard_cache_mb=0\n')
    run(['--config',cfg]);assert len(list((r/'exact-multi').glob('*/k20/results.tsv')))==2
    run(['--config',cfg,'--output-dir',r/'recursive','--reference-recursive']);assert len(list((r/'recursive').glob('*/k20/results.tsv')))==3
    (genomes/'space name.fa').write_text('>x\nACGT\n');(genomes/'space_name.fna').write_text('>x\nACGT\n')
    run([*args,'--output-dir',r/'sanitize-collision'],False)
    cfg.write_text(cfg.read_text().replace('reference_recursive=false','reference_recursive="false"'));run(['--config',cfg],False)
print('Production extensions passed')
