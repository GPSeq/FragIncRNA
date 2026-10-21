"""Physical-shard invariance, lifecycle, corruption and CLI contracts."""
import csv
import hashlib
import json
import os
import pathlib
import random
import subprocess
import sys
import tempfile
import fcntl

PROGRAM = pathlib.Path(sys.argv[1]).resolve()
PROVENANCE = {'ibf_layout','shard_count','max_bins_per_ibf','global_bin_count','index_manifest','manifest_sha256'}

def table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def science(rows):
    return [{k:v for k,v in row.items() if k not in PROVENANCE} for row in rows]

with tempfile.TemporaryDirectory(prefix='frag-shards-') as directory:
    root = pathlib.Path(directory)
    rng = random.Random(20261021)
    a = ''.join(rng.choices('ACGT', k=1457))
    b = ''.join(rng.choices('ACGT', k=327))
    ref = root/'ref.fa'; query = root/'query.fa'
    ref.write_text(f'>chrA\n{a[:600]}NNNN{a[604:]}\n>chrB\n{b}\n')
    rc = lambda s: s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    query.write_text(''.join(f'>{name}\n{seq}\n' for name,seq in [
        ('boundary',a[81:130]), ('splice',a[170:235]+a[430:495]),
        ('rc',rc(b[30:120])), ('repeat',(a[300:340]+'N')*3),
        ('absent',''.join(rng.choices('ACGT',k=140))), ('N','N'*40), ('short','AC')]))
    common = ['--reference',str(ref),'--queries',str(query),'--k','20','--fragment-size','100',
              '--fpr','0.15','--per-bin-min-hits','1','--per-bin-min-fraction','0',
              '--per-bin-max-locus-gap','400','--no-log-to-console']

    def run(name, *options, good=True):
        out = root/(name+'.tsv')
        cmd = [str(PROGRAM),*common,'--output',str(out),*map(str,options)]
        p = subprocess.run(cmd,capture_output=True,text=True,timeout=40)
        assert (p.returncode == 0) == good, (cmd,p.stderr)
        if not good: return p
        if '--ibf-dry-run' in options or '--ibf-build-only' in options: return None
        metrics = json.loads(pathlib.Path(str(out)+'.metrics.json').read_text())
        return out, table(out), metrics

    for mode in ['ibf','validate','calibrate']:
        for scoring in ['per_bin_ibf','per_bin_verified','global_ibf']:
            baseline = run('single_'+mode+scoring,'--mode',mode,'--scoring-method',scoring,'--ibf-layout','single')
            for cap in [1,2,3,7,100]:
                name = f'{mode}_{scoring}_{cap}'; index = root/(name+'_index')
                out,rows,m = run(name,'--mode',mode,'--scoring-method',scoring,'--ibf-layout','sharded',
                    '--max-bins-per-ibf',cap,'--ibf-shard-dir',index,'--threads','2')
                assert science(rows) == science(baseline[1]), (mode,scoring,cap)
                assert m['shard_count'] == (m['global_bin_count']+cap-1)//cap
                assert 1 <= m['max_resident_shards'] <= m['shard_count']
                assert m['shard_cache_peak_bytes'] <= m['ibf_shard_cache_budget_bytes']
                assert m['physical_occupied_bits_all_modes'] == baseline[2]['physical_occupied_bits_all_modes']
                for key in ['false_positive_windows','false_negative_windows','false_positive_bin_pairs',
                            'false_negative_bin_pairs','measured_per_bin_fpr','nominal_any_bin_fpr_independent_bins']:
                    assert m[key] == baseline[2][key], (key,mode,scoring,cap)
                if scoring != 'global_ibf' or mode != 'ibf':
                    for suffix in ['.per_bin.tsv','.per_bin.tsv.candidates.tsv']:
                        assert pathlib.Path(str(out)+suffix).read_bytes() == pathlib.Path(str(baseline[0])+suffix).read_bytes()
                if mode == 'calibrate':
                    assert pathlib.Path(str(out)+'.bins.tsv').read_bytes() == pathlib.Path(str(baseline[0])+'.bins.tsv').read_bytes()

    # Global top-K=1 must choose one strongest bin across three shards.
    top_ref = root/'top.fa'; top_query = root/'top_query.fa'
    top_ref.write_text('>contig\n'+a[:300]+'\n')
    top_query.write_text('>strong\n'+a[210:280]+a[20:42]+'\n')
    top_common = ['--reference',top_ref,'--queries',top_query,'--mode','validate',
        '--scoring-method','per_bin_verified','--ibf-layout','sharded','--max-bins-per-ibf','1',
        '--per-bin-top-k','1','--fpr','0.001']
    top, top_rows, _ = run('top',*top_common)
    loci = table(pathlib.Path(str(top)+'.per_bin.tsv.candidates.tsv'))
    assert len(loci) == 1 and loci[0]['source_bin_ids'] == '2'
    assert top_rows[0]['per_bin_candidate_count'] == '1'
    assert top_rows[0]['candidate_false_positive_windows'] == '0'
    # The original boundary window starts at 81, exactly 19 bases before end=100.
    boundary_rows = table(root/'calibrate_per_bin_verified_1.tsv')
    assert boundary_rows[0]['exact_global_recovery'] == '1.000000'
    boundary_loci = table(root/'calibrate_per_bin_verified_1.tsv.per_bin.tsv.candidates.tsv')
    boundary = [r for r in boundary_loci if r['query_id'] == 'boundary']
    assert any('0' in r['source_bin_ids'].split(',') for r in boundary)
    assert all(r['candidate_false_positive_windows'] in ('0','NA') for r in boundary_rows)
    # Strand and repeat policies retain identical logical scores across layouts.
    for strand in ['forward','both','canonical']:
        for repeat in ['include','exclude']:
            ref.write_text(ref.read_text().replace(a[170:190],a[170:190].lower()))
            opts = ['--mode','validate','--scoring-method','per_bin_verified','--strand',strand,'--repeat',repeat,'--per-bin-top-k','1']
            one = run(f'policy_single_{strand}{repeat}',*opts,'--ibf-layout','single')
            many = run(f'policy_shards_{strand}{repeat}',*opts,'--ibf-layout','sharded','--max-bins-per-ibf','2')
            assert science(one[1]) == science(many[1])
    ref.write_text(f'>chrA\n{a[:600]}NNNN{a[604:]}\n>chrB\n{b}\n')

    # Sparse validation must not rescan the reference for nonresident logical
    # pair totals; calibration still bounds the pairs it actually retains.
    poly_ref = root/'poly.fa'; poly_ref.write_text('>poly\n'+'A'*1000+'\n')
    poly_queries = root/'poly_queries.fa'; poly_queries.write_text(''.join(f'>p{i}\n'+ 'A'*20+'\n' for i in range(10)))
    for layout in ['single','sharded']:
        _, poly_rows, poly_metrics = run('sparse_'+layout,'--reference',poly_ref,'--queries',poly_queries,
            '--mode','validate','--scoring-method','per_bin_verified','--ibf-layout',layout,
            '--max-ibf-bin-pairs','15','--max-bins-per-ibf','2')
        assert poly_metrics['whole_reference_exact_scans'] == 1
        assert len(poly_rows) == 10 and all(r['recovery_fraction']=='1.000000' for r in poly_rows)
    index = root/'shared'; options = ['--mode','ibf','--scoring-method','per_bin_verified',
        '--ibf-layout','sharded','--max-bins-per-ibf','3','--ibf-shard-dir',index]
    out, first, m = run('lifecycle',*options)
    files = sorted(index.glob('*.ibf'))
    before = {p.name: (p.stat().st_mtime_ns,hashlib.sha256(p.read_bytes()).hexdigest()) for p in files}
    manifest = (index/'manifest.tsv').read_bytes()
    _,again,m2 = run('reuse',*options,'--threads','4')
    assert first == again and m2['shards_reused'] == m['shard_count'] and m2['shards_built'] == 0
    assert (index/'manifest.tsv').read_bytes() == manifest
    assert {p.name:(p.stat().st_mtime_ns,hashlib.sha256(p.read_bytes()).hexdigest()) for p in files} == before
    _, memory, mm = run('memory',*options,'--ibf-shard-storage','memory')
    assert memory == first and mm['max_resident_shards'] == mm['shard_count']
    assert mm['shard_query_loads'] == 0
    _, built_memory, built_mm = run('built_memory',*options,'--ibf-shard-dir',root/'new-memory-index','--ibf-shard-storage','memory')
    assert science(built_memory) == science(first)
    assert built_mm['shards_built'] == built_mm['shard_count'] and built_mm['shard_query_loads'] == 0

    config = root/'shards.toml'
    config.write_text('[v2]\n'+''.join(f'{key} = {json.dumps(value)}\n' for key,value in {
        'reference':str(ref),'queries':str(query),'output':str(root/'toml.tsv'),
        'mode':'ibf','scoring_method':'per_bin_verified','k':20,'fragment_size':100,'fpr':0.15,
        'per_bin_min_hits':1,'per_bin_min_fraction':0,'per_bin_max_locus_gap':400,
        'ibf_layout':'sharded','max_bins_per_ibf':3,'ibf_shard_dir':str(index),
        'ibf_shard_storage':'disk','ibf_shard_retention':'keep','ibf_shard_reuse':True,'ibf_resume':True,
        'log_to_console':False}.items()))
    p = subprocess.run([str(PROGRAM),'--config',str(config)],capture_output=True,text=True)
    assert p.returncode == 0, p.stderr
    assert table(root/'toml.tsv') == first
    other_ref = root/'other.fa'; other_ref.write_bytes(ref.read_bytes())
    other = run('other',*options,'--reference',other_ref,'--ibf-shard-dir',root/'other-index')[0]
    merge = pathlib.Path(__file__).resolve().parents[1]/'scripts/merge_kmer_v2.py'
    command = [sys.executable,str(merge),'--input',str(out),'--input',str(other),
        '--output',str(root/'merged.tsv'),'--summary',str(root/'merged.json')]
    p = subprocess.run(command,capture_output=True,text=True)
    assert p.returncode == 0, p.stderr
    altered = out.read_text().replace(m['manifest_sha256'],'0'*64)
    out.write_text(altered)
    p = subprocess.run(command,capture_output=True,text=True)
    assert p.returncode != 0 and 'layout provenance' in p.stderr

    built_only = root/'build-only'
    run('build_only',*options,'--ibf-shard-dir',built_only,'--ibf-build-only')
    assert (built_only/'manifest.tsv').exists() and not (root/'build_only.tsv').exists()
    run('query_only',*options,'--ibf-shard-dir',built_only,'--no-ibf-resume')
    next(built_only.glob('*.ibf')).unlink()
    run('no_resume_missing',*options,'--ibf-shard-dir',built_only,'--no-ibf-resume',good=False)
    run('output_collision',*options,'--output',index/'manifest.tsv',good=False)
    for cap, threads in [(65,1),(65,2),(65,4),(129,4)]:
        # >64 local bins exercises safe concurrent insertion into separate words.
        large_ref = root/'large.fa'; large_ref.write_text('>large\n'+a*15+'\n')
        args = ['--reference',large_ref,'--mode','validate','--scoring-method','per_bin_verified','--threads',threads]
        single = run(f'word_single{cap}_{threads}',*args,'--ibf-layout','single')
        sharded = run(f'word_shards{cap}_{threads}',*args,'--ibf-layout','sharded','--max-bins-per-ibf',cap)
        assert science(single[1]) == science(sharded[1])
    files[1].unlink()
    files[2].write_bytes(files[2].read_bytes()[:-1]+b'X')
    _, repaired, mr = run('repair',*options)
    assert science(first) == science(repaired) and mr['shards_built'] == 2
    # Simulate interruption after a pending checkpoint and a partial archive.
    pending = (index/'manifest.tsv').read_bytes().split(b'\n',2)[2]
    lines = pending.decode().splitlines()
    header = lines.index('shard_id\tglobal_begin\treal_bins\trelative_path\tsha256\tbytes\tbuild_run')
    fields = lines[header+1].split('\t')
    fields[4:7] = ['pending','0','']
    lines[header+1] = '\t'.join(fields)
    payload = ('\n'.join(lines)+'\n').encode()
    (index/'manifest.tsv').write_bytes(b'FRAGLNCRNA-SHARDS1\n'+hashlib.sha256(payload).hexdigest().encode()+b'\n'+payload)
    files[0].write_bytes(b'partial archive')
    sentinel = index/'unrelated.partial'; sentinel.write_text('preserve')
    run('pending_no_resume',*options,'--no-ibf-resume',good=False)
    _, pending_rows, pending_metrics = run('pending_resume',*options)
    assert science(pending_rows) == science(first) and pending_metrics['shards_built'] == 1
    assert sentinel.read_text() == 'preserve'
    # Cleanup never deletes reused archives or arbitrary data.
    unrelated = index/'reference.fa'; unrelated.write_text('do not delete')
    run('keep_shared',*options,'--ibf-shard-retention','cleanup_on_success')
    assert all(p.exists() for p in files) and unrelated.read_text() == 'do not delete'
    clean = root/'cleanup'; clean.mkdir(); (clean/'unrelated').write_text('safe')
    run('clean',*options,'--ibf-shard-dir',clean,'--ibf-shard-retention','cleanup_on_success')
    assert not list(clean.glob('*.ibf')) and (clean/'unrelated').read_text() == 'safe'
    assert (clean/'manifest.tsv').exists()
    # Failure after building preserves owned shards for resume.
    broken_query = root/'broken.fa'; broken_query.write_text('>dup\nACGT\n>dup\nAAAA\n')
    failed_index = root/'failed'
    run('failure',*options,'--queries',broken_query,'--ibf-shard-dir',failed_index,
        '--ibf-shard-retention','cleanup_on_success',good=False)
    assert list(failed_index.glob('*.ibf'))
    # Lock is released automatically by the OS; second writer is rejected.
    with (index/'.lock').open('r+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        assert 'locked' in run('locked',*options,good=False).stderr
    for extra in [['--k','15'],['--fragment-size','101'],['--fpr','0.2'],['--strand','forward'],
                  ['--repeat','exclude'],['--mask-target','reference'],['--hash-functions','2']]:
        assert 'Incompatible shard manifest' in run('mismatch',*options,*extra,good=False).stderr
    ref.write_text(ref.read_text().replace(a[:20],rc(a[:20]),1))
    assert 'Incompatible shard manifest' in run('reference_mismatch',*options,good=False).stderr
    ref.write_text(f'>chrA\n{a[:600]}NNNN{a[604:]}\n>chrB\n{b}\n')
    (index/'manifest.tsv').write_bytes(manifest[:-5]+b'X\n')
    assert 'manifest checksum' in run('corrupt_manifest',*options,good=False).stderr
    for value in ['0','-1','18446744073709551616','18446744073709551615']:
        run('bad_cap',*options,'--max-bins-per-ibf',value,good=False)
    assert 'single .ibf' in run('legacy_conflict',*options,'--save-index',root/'old.ibf',good=False).stderr
    assert 'physical RAM' in run('unsafe_physical','--mode','ibf','--ibf-layout','sharded','--fragment-size','1000000000000000','--ibf-dry-run',good=False).stderr
    assert 'exceed' in run('budget',*options,'--ibf-max-index-memory-gb','0.000001',good=False).stderr
    _, _, auto = run('auto','--mode','ibf','--ibf-layout','auto','--max-bins-per-ibf','2')
    assert auto['ibf_layout'] == 'sharded' and auto['scoring_method'] == 'per_bin_verified'
    _, _, exact = run('exact','--mode','exact','--ibf-layout','sharded')
    assert exact['shard_count'] == 0 and exact['ibf_layout'] == 'none'
    assert not (root/'exact.tsv.ibf-shards').exists()
    run('auto_build_without_destination','--mode','ibf','--ibf-layout','auto','--max-bins-per-ibf','1000000','--ibf-build-only',good=False)
    # Explicit any-bin budgeting uses global B, never bins in an individual shard.
    for strand in ['canonical','both']:
        args = ['--mode','ibf','--scoring-method','global_ibf','--genome-fpr','0.01','--ibf-auto-tune','apply','--strand',strand]
        one = run('budget_one_'+strand,*args,'--ibf-layout','single')
        many = run('budget_many_'+strand,*args,'--ibf-layout','sharded','--max-bins-per-ibf','2')
        assert science(one[1]) == science(many[1])
        assert one[2]['configured_per_bin_fpr'] == many[2]['configured_per_bin_fpr']
    _, _, explicit_unverified = run('explicit_unverified','--mode','ibf','--ibf-layout','sharded','--no-per-bin-exact-verify')
    assert explicit_unverified['scoring_method'] == 'per_bin_ibf'
    _, _, explicit_global = run('explicit_global','--mode','ibf','--ibf-layout','sharded','--no-fdr-control-per-bin')
    assert explicit_global['scoring_method'] == 'global_ibf'
    _, _, requested = run('requested_target','--mode','ibf','--ibf-layout','sharded','--genome-fpr','0.01','--ibf-auto-tune','recommend')
    assert requested['global_error_budget_policy'] == 'requested_only'
    assert requested['configured_per_bin_fpr'] == 0.15
    _, _, default_validation = run('default_validation','--mode','validate','--ibf-layout','sharded')
    assert default_validation['candidate_exact_verification'] is True
    run('dry','--mode','ibf','--ibf-layout','sharded','--ibf-dry-run')
    assert not (root/'dry.tsv.ibf-shards').exists()
print('sharded IBF physical-layout invariance and lifecycle checks passed')
