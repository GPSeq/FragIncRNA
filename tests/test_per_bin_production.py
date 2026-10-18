"""Independent scientific and operational oracles for schema-3 production."""
import csv
import gzip
import hashlib
import json
import pathlib
import random
import signal
import subprocess
import sys
import tempfile
import time

PROGRAM = pathlib.Path(sys.argv[1]).resolve()
ROOT = pathlib.Path(__file__).resolve().parents[1]


def rows(path):
    with path.open() as handle:
        result = list(csv.DictReader(handle, delimiter="\t"))
    assert all(None not in row and None not in row.values() for row in result), result
    return {row["query_id"]: row for row in result}


def fasta(path, records):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


def oracle(seq, refs, k, strand="canonical"):
    rc = lambda text: text.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
    seq = seq.upper()
    valid = [seq[i:i+k] for i in range(max(0, len(seq)-k+1)) if set(seq[i:i+k]) <= set("ACGT")]
    return len(valid), sum(any(w in r.upper() or (strand != "forward" and rc(w) in r.upper()) for r in refs) for w in valid)


with tempfile.TemporaryDirectory(prefix="fraglncrna-production-") as directory:
    root = pathlib.Path(directory)
    rng = random.Random(20261018)
    reference = "".join(rng.choices("ACGT", k=30*512))
    absent = "".join(rng.choices("ACGT", k=200))
    local = reference[2660:2760]
    ref = root / "ref.fa"
    query = root / "q.fa"
    records = [("absent", absent), ("local", local), ("N", "N"*40), ("short", "AC")]
    fasta(ref, [("chrA",reference)])
    fasta(query,records)
    common = ["--reference",ref,"--queries",query,"--k",15,"--fragment-size",512,
              "--fpr",0.5,"--per-bin-min-fraction",0.8,"--per-bin-max-locus-gap",0]

    def run(name, *options, good=True, base=None):
        out = root / (name + ".tsv")
        cmd = [str(PROGRAM), *map(str,common if base is None else base),"--output",str(out),*map(str,options)]
        process = subprocess.run(cmd,capture_output=True,text=True,timeout=20)
        assert (process.returncode == 0) == good, process.stderr
        if good:
            return out, rows(out), json.loads(pathlib.Path(str(out)+".metrics.json").read_text()), process
        return out, process

    index = root / "same.ibf"
    cache = root / "cache.bases"
    out, prediction, metrics, process = run("production", "--mode","ibf","--fdr-control-per-bin",
        "--save-index",index,"--reference-cache",cache)
    assert prediction["absent"]["status"] == "NO_MATCH"
    assert float(prediction["absent"]["global_ibf_recovery"]) > .8
    assert prediction["absent"]["per_bin_candidate_count"] == "0"
    assert prediction["local"]["status"] == "PASS_STRICT"
    assert prediction["local"]["recovery_fraction"] == "1.000000"
    assert prediction["N"]["recovery_fraction"] == "NA"
    assert prediction["short"]["status"] == "NO_VALID_KMERS"
    assert metrics["whole_reference_exact_scans"] == 0
    assert metrics["tested_absent_windows"] is None and metrics["exact_positive_windows"] is None
    assert 0 < metrics["candidate_reference_bytes_read"] < len(reference)
    assert metrics["candidate_reference_bytes_read"] == metrics["candidate_reference_bases_fetched"]
    assert all(row["scoring_method"] == "per_bin_verified" and row["verification_scope"] == "candidate_intervals"
               and row["exact_global_recovery"] == "NA" for row in prediction.values())
    assert all(row["schema_version"] == "3" for row in prediction.values())
    tune = json.loads(pathlib.Path(str(out)+".ibf_tuning.json").read_text())
    assert not tune["changed"] and tune["effective_per_bin_fpr"] == .5
    assert len(tune["hash_sensitivity"]) == 5 and len(tune["fragment_alternatives"]) == 3
    assert "nominal any-bin false-positive probability" in process.stderr
    assert "recommendation=" in process.stderr and "B=30" in process.stderr

    logs = [json.loads(line) for line in pathlib.Path(str(out)+".log.jsonl").read_text().splitlines()]
    assert logs[0]["event"] == "stage_started" and logs[-1]["event"] == "run_completed"
    assert any(item["phase"] == "verify" and item["event"] == "stage_completed" for item in logs)
    assert not any(item["phase"] == "global-exact" for item in logs)
    assert len({item["run_id"] for item in logs}) == 1
    assert metrics["run_id"] == logs[0]["run_id"]
    assert all(item["elapsed_seconds"] >= 0 for item in logs)

    vout, validation, vm, _ = run("validate","--mode","validate","--scoring-method","per_bin_verified",
        "--load-index",index,"--reference-cache",cache,"--threads",4)
    assert vm["reference_cache_reused"] and vm["whole_reference_exact_scans"] == 1
    assert vm["measured_per_bin_fpr"] is not None
    assert vm["whole_reference_exact_bases_visited"] == len(reference)
    for name,seq in records:
        n, exact = oracle(seq,[reference],15)
        row = validation[name]
        assert int(row["valid_windows"]) == n
        assert row["exact_global_status"] == ("NO_VALID_KMERS" if not n else "NO_MATCH" if not exact else "PASS_STRICT")
        assert row["per_bin_verified_local_recovery"] == prediction[name]["per_bin_verified_local_recovery"]
        assert row["candidate_false_positive_windows"] == "0"
        assert row["candidate_false_negative_windows"] == "0"
        if n:
            c = int(row["per_bin_candidate_window_support"])
            assert int(row["per_bin_ibf_false_positive_windows"])-int(row["per_bin_ibf_false_negative_windows"]) == c-exact
    _, _, calibration, _ = run("calibration", "--mode", "calibrate", "--scoring-method", "per_bin_verified",
        "--load-index", index, "--reference-cache", cache)
    for key in ("tested_absent_bin_pairs", "false_positive_bin_pairs", "false_negative_bin_pairs"):
        assert vm[key] == calibration[key], (key, vm[key], calibration[key])

    # Several thread counts and repeated runs compare all scientific TSV bytes.
    for threads in (1,2,4,15,4):
        other, _, m, _ = run(f"threads{threads}","--mode","ibf","--scoring-method","per_bin_verified",
            "--load-index",index,"--reference-cache",cache,"--threads",threads,"--log-level","off")
        assert other.read_bytes() == out.read_bytes()
        for suffix in (".per_bin.tsv", ".per_bin.tsv.candidates.tsv"):
            assert pathlib.Path(str(other)+suffix).read_bytes() == pathlib.Path(str(out)+suffix).read_bytes()
        assert m["threads"] <= threads and m["requested_threads"] == threads

    # Gzip source follows documented one-time normalized cache strategy.
    compressed = root / "ref.fa.gz"
    with gzip.open(compressed,"wt") as handle: handle.write(ref.read_text())
    gzbase = list(common); gzbase[1] = compressed
    _, gr, gm, _ = run("gzip","--mode","ibf","--scoring-method","per_bin_verified",
        "--reference-cache",cache,base=gzbase)
    assert gm["reference_cache_reused"]
    assert gr["local"]["matched_windows"] == prediction["local"]["matched_windows"]

    # Corrupted reusable cache is never silently replaced.
    data = cache.read_bytes(); cache.write_bytes(b"X"+data[1:])
    bad, failure = run("bad_cache","--mode","ibf","--scoring-method","per_bin_verified",
        "--reference-cache",cache,good=False)
    assert "checksum mismatch" in failure.stderr
    assert cache.read_bytes()[0:1] == b"X"
    footer = [json.loads(line) for line in pathlib.Path(str(bad)+".log.jsonl").read_text().splitlines()]
    assert footer[-1]["event"] == "run_failed"
    cache.write_bytes(data)

    # Unsafe global objective: recommend retains settings; apply meets bound;
    # impossible memory requests fail before allocation and preserve report.
    _, failure = run("unsafe","--mode","ibf","--scoring-method","global_ibf",
        "--ibf-fail-on-unsafe-settings",good=False)
    assert "Unsafe global_ibf" in failure.stderr
    gout,_,_,_ = run("tuned","--mode","ibf","--scoring-method","global_ibf","--ibf-auto-tune","apply")
    gt = json.loads(pathlib.Path(str(gout)+".ibf_tuning.json").read_text())
    assert gt["changed"] and gt["union_bound"] <= .01000000001
    _, failure = run("budget","--mode","ibf","--scoring-method","per_bin_verified",
        "--ibf-max-index-memory-gb",1e-9,good=False)
    assert "exceed --ibf-max-index-memory-gb" in failure.stderr

    # A feasible chosen h=5 design must survive unrepresentable h=1
    # alternatives; tiny probabilities must not be serialized as zero.
    tiny = root/"tiny.tsv"
    tiny_command = [str(PROGRAM),"--reference",str(ROOT/"tests/fixtures/kmer_v2/reference.fa"),
        "--queries",str(query),"--output",str(tiny),"--mode","ibf","--scoring-method","global_ibf",
        "--k","4","--fragment-size","4","--hash-functions","5","--fpr","1e-30",
        "--ibf-auto-tune","off","--ibf-dry-run","--ibf-max-index-memory-gb","1"]
    tiny_process = subprocess.run(tiny_command,capture_output=True,text=True,timeout=20)
    assert tiny_process.returncode == 0, tiny_process.stderr
    tiny_report = json.loads(pathlib.Path(str(tiny)+".ibf_tuning.json").read_text())
    assert tiny_report["effective_per_bin_fpr"] == 1e-30
    assert tiny_report["hash_sensitivity"][0]["same_p_bin_bytes"] is None
    assert "overflow" in tiny_report["hash_sensitivity"][0]["same_p_bin_error"]
    assert tiny_report["hash_sensitivity"][4]["same_p_bin_bytes"] > 0
    actual_tiny = root/"actual_tiny.tsv"
    actual_command = [argument for argument in tiny_command if argument != "--ibf-dry-run"]
    actual_command[actual_command.index(str(tiny))] = str(actual_tiny)
    actual_command[actual_command.index("1e-30")] = "1e-18"
    actual_process = subprocess.run(actual_command,capture_output=True,text=True,timeout=20)
    assert actual_process.returncode == 0, actual_process.stderr
    actual_metrics = json.loads(pathlib.Path(str(actual_tiny)+".metrics.json").read_text())
    assert actual_metrics["configured_per_bin_fpr"] == 1e-18
    assert all(float(row["configured_per_bin_fpr"]) == 1e-18 and float(row["index_fpr_per_bin"]) == 1e-18
               for row in rows(actual_tiny).values())

    # Many original windows with one repeated key must count fully without
    # a reference-occurrence x query-occurrence Cartesian-product scan.
    poly_ref = root/"poly.fa"; poly_query = root/"poly_q.fa"
    fasta(poly_ref,[("poly", "A"*320000)])
    fasta(poly_query,[("poly_query", "A"*250000)])
    poly_base = ["--reference",poly_ref,"--queries",poly_query,"--k",15,"--fragment-size",80000,
                 "--mode","ibf","--scoring-method","per_bin_verified","--threads",4]
    _, poly_rows, poly_metrics, _ = run("poly",base=poly_base)
    assert poly_rows["poly_query"]["matched_windows"] == str(250000-15+1)
    assert poly_rows["poly_query"]["recovery_fraction"] == "1.000000"
    assert poly_metrics["candidate_repeated_reference_keys_skipped"] > 100000

    # Adjacent normalized reads can have different query owners. A shared
    # motif found in the first locus must still be verified in the second.
    boundary_rng = random.Random(91)
    motif = "".join(boundary_rng.choices("ACGT",k=40))
    left = "".join(boundary_rng.choices("ACGT",k=80))+motif
    right = "".join(boundary_rng.choices("ACGT",k=80))+motif
    boundary_ref = root/"boundary.fa"; boundary_query = root/"boundary_q.fa"
    fasta(boundary_ref,[("chr", "N"*100+left+"N"*(512-220)+"N"*100+right+"N"*(1024-732))])
    fasta(boundary_query,[("left",left),("right",right)])
    boundary_base = ["--reference",boundary_ref,"--queries",boundary_query,"--k",15,"--fragment-size",512,
        "--mode","validate","--scoring-method","per_bin_verified","--per-bin-min-fraction",.8,
        "--per-bin-max-locus-gap",0,"--threads",4]
    _, boundary_rows, _, _ = run("owner_boundaries",base=boundary_base)
    assert all(row["recovery_fraction"] == "1.000000" and row["candidate_false_negative_windows"] == "0"
               for row in boundary_rows.values())

    # Split exons, cross-contig and top-K losses checked with string truth.
    exon1 = "".join(rng.choices("ACGT",k=80)); exon2 = "".join(rng.choices("ACGT",k=80))
    sr = root / "spliced.fa"; sq = root / "spliced_q.fa"
    sequence = "G"*100 + exon1 + "N"*1000 + exon2 + "C"*100
    fasta(sr,[("chrS",sequence),("empty", "N"*200)])
    fasta(sq,[("splice",exon1+exon2),("repeat",exon1+"NNN"+exon1)])
    spbase = ["--reference",sr,"--queries",sq,"--k",15,"--fragment-size",100,"--mode","validate",
        "--scoring-method","per_bin_verified","--per-bin-min-fraction",.2,"--per-bin-max-locus-gap",2000]
    _, sp, _, _ = run("splice",base=spbase)
    n, e = oracle(exon1+exon2,[sequence,"N"*200],15)
    assert abs(float(sp["splice"]["exact_global_recovery"])-e/n) < 1e-6
    assert sp["splice"]["candidate_false_negative_windows"] == "0"
    assert sp["repeat"]["per_bin_verified_local_recovery"] == "1.000000"
    _, limited, _, _ = run("topk","--per-bin-top-k",1,base=spbase)
    assert int(limited["splice"]["candidate_false_negative_windows"]) > 0
    assert float(limited["splice"]["candidate_window_recall"]) < 1
    assert float(limited["splice"]["per_bin_verified_local_recovery"]) < float(limited["splice"]["exact_global_recovery"])

    for strand in ("canonical","forward","both"):
        # Mask first exon using BED, N ambiguity retained, original contig end.
        bed = root / "repeat.bed"; bed.write_text("chrS\t100\t180\tREPEAT\n")
        _, masked, _, _ = run("masked_"+strand,"--strand",strand,"--repeat","exclude","--mask-target","reference",
            "--reference-repeat-bed",bed,base=spbase)
        assert int(masked["splice"]["candidate_false_positive_windows"]) == 0
        assert 0 <= float(masked["splice"]["per_bin_verified_local_recovery"]) <= float(masked["splice"]["exact_global_recovery"])

    # Exceptions in a parallel query worker propagate with terminal log footer.
    _, failure = run("worker_failure","--threads",4,"--mode","ibf","--scoring-method","per_bin_verified",
        "--max-ibf-bin-pairs",1,good=False)
    assert "IBF-positive bin pairs exceed" in failure.stderr

    config = root / "production.toml"
    config.write_text(f'[v2]\nreference="{ref}"\nqueries="{query}"\noutput="{root / "config.tsv"}"\n'
        'mode="ibf"\nscoring_method="per_bin_verified"\nk=15\nfragment_size=512\nfpr=0.5\n'
        'fdr_control_per_bin=true\nper_bin_exact_verify=true\nper_bin_min_fraction=0.8\nper_bin_max_locus_gap=0\nthreads=4\n')
    process = subprocess.run([str(PROGRAM),"--config",str(config)],capture_output=True,text=True)
    assert process.returncode == 0, process.stderr
    assert (root / "config.tsv").read_bytes() == out.read_bytes()
    cm = json.loads((root / "config.tsv.metrics.json").read_text())
    assert cm["config_sha256"] == hashlib.sha256(config.read_bytes()).hexdigest()

    # Startup failure paths retain structured terminal events.
    for label, changes in (("missing_index", ["--load-index",root/"missing.ibf"]),
                           ("missing_reference", ["--reference",root/"missing.fa"])):
        failed, result = run(label,"--mode","ibf","--scoring-method","per_bin_verified",*changes,good=False)
        events = [json.loads(line) for line in pathlib.Path(str(failed)+".log.jsonl").read_text().splitlines()]
        assert events[-1]["event"] == "run_failed"
    malformed = root/"malformed.fa"; malformed.write_text("ACGT\n")
    failed,result = run("malformed","--reference",malformed,"--mode","ibf","--scoring-method","per_bin_verified",good=False)
    assert "Sequence before FASTA header" in result.stderr
    assert not failed.exists()  # main output was never finalized

    # Graceful signal handling writes an error summary outside the handler.
    long_ref = root/"long.fa"
    fasta(long_ref,[("long", "ACGT"*1000000)])
    cancelled = root/"cancelled.tsv"
    command = [str(PROGRAM),"--reference",str(long_ref),"--queries",str(query),"--output",str(cancelled),
               "--mode","exact","--scoring-method","global_exact"]
    child = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    logpath = pathlib.Path(str(cancelled)+".log.jsonl")
    deadline = time.monotonic()+5
    while not logpath.exists() and child.poll() is None and time.monotonic()<deadline:
        time.sleep(.005)
    assert logpath.exists()
    child.send_signal(signal.SIGTERM)
    _,stderr = child.communicate(timeout=10)
    assert child.returncode != 0 and "interrupted by signal" in stderr
    events = [json.loads(line) for line in logpath.read_text().splitlines()]
    assert events[-1]["event"] == "run_failed"

print("production scientific, interval-access, logging, tuning, and determinism checks passed")
