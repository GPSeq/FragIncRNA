"""Independent black-box checks for optional per-bin coherence screening."""

import csv
import json
import pathlib
import random
import subprocess
import sys
import tempfile


PROGRAM = pathlib.Path(sys.argv[1])
DNA = "ACGT"


def fasta(path, records):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(root, name, reference, queries, *options, good=True):
    output = root / f"{name}.tsv"
    command = [str(PROGRAM), "--reference", str(reference), "--queries", str(queries),
               "--output", str(output), *options]
    process = subprocess.run(command, capture_output=True, text=True)
    if good:
        assert process.returncode == 0, process.stderr
        paired = {row["query_id"]: row for row in table(pathlib.Path(str(output) + ".per_bin.tsv"))}
        candidates = table(pathlib.Path(str(output) + ".per_bin.tsv.candidates.tsv"))
        return output, paired, candidates, json.loads(pathlib.Path(str(output) + ".metrics.json").read_text()), process
    assert process.returncode != 0
    return process


with tempfile.TemporaryDirectory(prefix="fraglncrna-per-bin-") as directory:
    root = pathlib.Path(directory)
    rng = random.Random(20261018)
    random_reference = "".join(rng.choices(DNA, k=30 * 512))
    absent = "".join(rng.choices(DNA, k=200))
    localized = random_reference[5 * 512 + 100:5 * 512 + 200]
    ref = root / "random_ref.fa"
    q = root / "random_queries.fa"
    fasta(ref, [("chrA", random_reference)])
    fasta(q, [("absent", absent), ("localized", localized), ("all_N", "N" * 40)])
    index = root / "high.ibf"
    options = ("--k", "15", "--fragment-size", "512", "--hash-functions", "3",
               "--fpr", "0.5", "--mode", "validate", "--fdr-control-per-bin",
               "--per-bin-min-hits", "2", "--per-bin-min-fraction", "0.8",
               "--per-bin-max-locus-gap", "0")
    output, rows, candidates, metrics, process = run(root, "scatter", ref, q, *options,
                                                     "--save-index", str(index))
    assert rows["absent"]["exact_global_matched_windows"] == "0"
    assert float(rows["absent"]["ibf_union_recovery"]) > 0.8
    assert rows["absent"]["union_status"] == "PASS_STRICT"
    assert rows["absent"]["exact_status"] == "NO_MATCH"
    assert rows["absent"]["candidate_bin_count"] == "0"
    assert rows["absent"]["per_bin_status"] == "NO_CANDIDATE"
    assert rows["localized"]["per_bin_status"] == "VERIFIED_CANDIDATE"
    assert rows["localized"]["best_locus_exact_recovery"] == "1.000000"
    assert rows["all_N"]["per_bin_status"] == "NO_VALID_KMERS"
    assert rows["all_N"]["max_bin_support"] == "NA"
    assert "nominal any-bin false-positive probability" in process.stderr
    assert metrics["false_positive_windows"] > 0

    _, diagnostic, _, loaded_metrics, warning = run(
        root, "unverified", ref, q, "--k", "15", "--fragment-size", "512",
        "--hash-functions", "3", "--fpr", "0.5", "--mode", "ibf",
        "--fdr-control-per-bin", "--no-per-bin-exact-verify",
        "--per-bin-min-hits", "2", "--per-bin-min-fraction", "0.8",
        "--per-bin-max-locus-gap", "0", "--load-index", str(index))
    assert diagnostic["absent"]["per_bin_status"] == "NO_CANDIDATE"
    assert diagnostic["absent"]["exact_global_recovery"] == "NA"
    assert diagnostic["localized"]["per_bin_status"] == "CANDIDATE_UNVERIFIED"
    assert loaded_metrics["index_provenance"] == "loaded"
    assert loaded_metrics["index_file_sha256"] == metrics["index_file_sha256"]
    assert "unverified IBF evidence" in warning.stderr

    exon1 = "".join(rng.choices(DNA, k=80))
    exon2 = "".join(rng.choices(DNA, k=80))
    splice_ref = ("".join(rng.choices(DNA, k=100)) + exon1 +
                  "".join(rng.choices(DNA, k=800)) + exon2 +
                  "".join(rng.choices(DNA, k=100)))
    splice_query = exon1 + exon2
    splice_r = root / "splice_ref.fa"
    splice_q = root / "splice_queries.fa"
    fasta(splice_r, [("chrS", splice_ref)])
    fasta(splice_q, [("splice", splice_query), ("reverse", exon1.translate(
        str.maketrans("ACGT", "TGCA"))[::-1]), ("N_gap", exon1[:30] + "N" + exon2[:30])])
    expected = root / "expected.bed"
    expected.write_text("chrS\t100\t1100\tsplice\tgene\n")
    splice_options = ("--k", "15", "--fragment-size", "100", "--mode", "validate",
                      "--fdr-control-per-bin", "--per-bin-min-hits", "2",
                      "--per-bin-min-fraction", "0.20", "--per-bin-max-locus-gap", "1000")
    first, rows, loci, _, _ = run(root, "splice", splice_r, splice_q, *splice_options,
                                  "--expected-loci-bed", str(expected))
    spliced = rows["splice"]
    assert int(spliced["candidate_bin_count"]) >= 2
    assert spliced["candidate_locus_count"] == "1"
    assert spliced["best_locus_exact_matched_windows"] == spliced["exact_global_matched_windows"]
    assert spliced["candidate_union_exact_matched_windows"] == spliced["exact_global_matched_windows"]
    assert 0 < float(spliced["max_bin_support"]) < float(spliced["best_locus_exact_recovery"])
    assert rows["reverse"]["exact_global_matched_windows"] != "0"
    assert rows["N_gap"]["valid_windows"] == "32"  # two 30-base blocks, k=15
    assert pathlib.Path(str(first) + ".expected_loci.tsv").exists()

    second, _, _, _, _ = run(root, "splice_repeat", splice_r, splice_q, *splice_options,
                              "--expected-loci-bed", str(expected))
    assert pathlib.Path(str(first) + ".per_bin.tsv").read_bytes() == pathlib.Path(
        str(second) + ".per_bin.tsv").read_bytes()
    assert pathlib.Path(str(first) + ".per_bin.tsv.candidates.tsv").read_bytes() == pathlib.Path(
        str(second) + ".per_bin.tsv.candidates.tsv").read_bytes()

    _, limited, _, _, _ = run(root, "top_one", splice_r, splice_q, *splice_options,
                               "--per-bin-top-k", "1")
    assert int(limited["splice"]["top_k_dropped_bins"]) >= 1
    assert int(limited["splice"]["global_exact_windows_outside_candidates"]) > 0
    assert int(limited["splice"]["candidate_union_exact_matched_windows"]) < int(
        limited["splice"]["exact_global_matched_windows"])

    split_r = root / "split_ref.fa"
    split_q = root / "split_query.fa"
    fasta(split_r, [("chr1", "A" * 100 + exon1 + "C" * 100),
                    ("chr2", "T" * 100 + exon2 + "G" * 100)])
    fasta(split_q, [("split", splice_query)])
    _, split, split_loci, _, _ = run(root, "split_contigs", split_r, split_q,
                                     *splice_options)
    assert int(split["split"]["candidate_locus_count"]) >= 2
    assert int(split["split"]["best_locus_exact_matched_windows"]) < int(
        split["split"]["candidate_union_exact_matched_windows"])
    assert len({row["contig"] for row in split_loci}) == 2

    masked_r = root / "masked_ref.fa"
    fasta(masked_r, [("chrS", splice_ref[:100] + exon1.lower() + splice_ref[180:])])
    _, masked, _, _, _ = run(root, "masked", masked_r, splice_q, *splice_options,
                              "--repeat", "exclude", "--mask-target", "reference")
    assert int(masked["splice"]["exact_global_matched_windows"]) < int(
        spliced["exact_global_matched_windows"])
    forward = run(root, "forward", splice_r, splice_q, *splice_options, "--strand", "forward")
    assert forward[1]["reverse"]["exact_global_matched_windows"] == "0"

    for extra, message in ((('--mode', 'exact'), 'requires an IBF mode'),
                           (('--per-bin-min-fraction', '1.1'), 'Invalid numeric option'),
                           (('--per-bin-top-k', '0'), 'Invalid numeric option'),
                           (('--per-bin-top-k', '1025'), 'Invalid numeric option'),
                           (('--per-bin-min-hits', '-1'), 'Invalid value'),
                           (('--hash-functions', '22'), 'Invalid numeric option')):
        failed = run(root, "bad", splice_r, splice_q, *splice_options, *extra, good=False)
        assert message in failed.stderr, failed.stderr

print("per-bin IBF synthetic checks passed")
