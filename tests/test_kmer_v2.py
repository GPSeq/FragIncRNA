"""Deterministic black-box checks for the version 2 matching path."""
import csv
import gzip
import json
import pathlib
import random
import subprocess
import sys
import tempfile


PROGRAM = pathlib.Path(sys.argv[1])
FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "kmer_v2"


def run(tmp, label, *, reference=None, queries=None, k=4, fragment_size=4, extra=()):
    output = tmp / (label + ".tsv")
    command = [str(PROGRAM), "--reference", str(reference or FIXTURE / "reference.fa"),
               "--queries", str(queries or FIXTURE / "queries.fa"),
               "--output", str(output), "--k", str(k), "--fragment-size", str(fragment_size),
               *extra]
    subprocess.run(command, check=True, capture_output=True, text=True)
    with output.open(newline="") as handle:
        rows = {r["query_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    return rows, json.loads(pathlib.Path(str(output) + ".metrics.json").read_text()), output


with tempfile.TemporaryDirectory(prefix="fraglncrna-v2-test-") as directory:
    tmp = pathlib.Path(directory)
    rows, metrics, first = run(tmp, "validate")
    assert first.read_text() == (FIXTURE / "expected_validate.tsv").read_text()
    for key, expected in json.loads((FIXTURE / "expected_metrics_fields.json").read_text()).items():
        assert metrics[key] == expected, (key, metrics[key], expected)
    assert rows["reverse_only"]["exact_recovery"] == "1.000000"
    assert rows["reverse_only"]["ibf_recovery"] == "1.000000"
    assert rows["contig_join"]["exact_recovery"] == "0.000000"
    assert rows["n_join"]["exact_recovery"] == "0.000000"
    assert rows["duplicated"]["matched_windows"] == "2"
    assert rows["duplicated"]["distinct_matched_kmers"] == "1"
    assert rows["palindrome"]["matched_windows"] == "1"
    assert rows["all_n"]["status"] == "NO_VALID_KMERS"
    assert rows["all_n"]["recovery_fraction"] == "NA"
    assert rows["short"]["total_windows"] == "0"
    assert rows["short"]["status"] == "NO_VALID_KMERS"
    assert metrics["bin_count"] == 6
    assert metrics["query_count"] == 8
    assert metrics["peak_resident_bytes"] > 0
    assert metrics["index_construction_seconds"] >= 0
    assert metrics["exact_scan_seconds"] >= 0

    again, _, second = run(tmp, "repeat")
    assert first.read_text() == second.read_text()
    assert rows == again

    forward, _, _ = run(tmp, "forward", extra=("--strand", "forward"))
    assert forward["reverse_only"]["exact_recovery"] == "0.000000"
    both, _, _ = run(tmp, "both", extra=("--strand", "both"))
    assert both["reverse_only"]["exact_recovery"] == "1.000000"
    assert both["palindrome"]["matched_windows"] == "1"

    exact, _, _ = run(tmp, "exact", extra=("--mode", "exact"))
    assert exact["duplicated"]["matched_windows"] == "2"
    assert exact["duplicated"]["distinct_matched_kmers"] == "1"
    assert exact["duplicated"]["nonrepeat_recovery"] == "0.400000"
    exact_local, exact_metrics, _ = run(tmp, "exact_local",
                                        extra=("--mode", "exact", "--local-window", "4"))
    assert exact_local["reverse_only"]["best_region_recovery"] == "1.000000"
    assert exact_metrics["measured_genome_wide_fpr"] is None
    ibf_only, ibf_metrics, _ = run(tmp, "ibf_only", extra=("--mode", "ibf"))
    assert ibf_only["reverse_only"]["exact_recovery"] == "NA"
    assert ibf_metrics["false_positive_windows"] is None
    small_batches, _, _ = run(tmp, "small_batches", extra=("--batch-windows", "1"))
    assert small_batches == rows

    excluded, _, _ = run(tmp, "exclude", extra=("--repeat", "exclude"))
    assert excluded["softmasked"]["status"] == "NO_VALID_KMERS"
    assert excluded["softmasked"]["repeat_fraction"] == "1.000000"
    assert excluded["duplicated"]["matched_windows"] == "2"
    annotated, _, annotate_output = run(tmp, "annotate", extra=("--repeat", "annotate"))
    assert annotated["softmasked"]["repeat_matched_windows"] == "1"
    assert annotated["softmasked"]["valid_nonrepeat_windows"] == "0"
    with pathlib.Path(str(annotate_output) + ".repeat_classes.tsv").open(newline="") as handle:
        softmask_classes = list(csv.DictReader(handle, delimiter="\t"))
    assert any(row["query_id"] == "softmasked" and row["repeat_class"] == "SOFTMASK" and
               row["matched_windows"] == "1" for row in softmask_classes)

    ambiguity_query = tmp / "ambiguity.fa"
    ambiguity_query.write_text(">mixed\naAgCRacgt\n")
    ambiguous, _, _ = run(tmp, "ambiguity", queries=ambiguity_query)
    assert ambiguous["mixed"]["total_windows"] == "6"
    assert ambiguous["mixed"]["valid_windows"] == "2"
    assert ambiguous["mixed"]["invalid_windows"] == "4"
    empty_query = tmp / "empty.fa"
    empty_query.write_text(">empty\n>next\nAAGC\n")
    empty_rows, _, _ = run(tmp, "empty", queries=empty_query)
    assert empty_rows["empty"]["status"] == "NO_VALID_KMERS"
    assert empty_rows["next"]["exact_recovery"] == "1.000000"

    repeat_ref = tmp / "repeat.fa"
    repeat_ref.write_text(">masked\naagc\n")
    included, _, _ = run(tmp, "repeat_include", reference=repeat_ref)
    excluded_ref, _, _ = run(tmp, "repeat_exclude", reference=repeat_ref,
                             extra=("--repeat", "exclude"))
    assert included["duplicated"]["exact_recovery"] == "0.400000"
    assert excluded_ref["duplicated"]["exact_recovery"] == "0.000000"

    with gzip.open(tmp / "reference.fa.gz", "wt") as handle:
        handle.write((FIXTURE / "reference.fa").read_text())
    gz_rows, gz_metrics, _ = run(tmp, "gzip", reference=tmp / "reference.fa.gz")
    assert gz_rows["duplicated"]["exact_recovery"] == "0.400000"
    assert gz_metrics["reference_sequence_fnv1a64"] == metrics["reference_sequence_fnv1a64"]
    with gzip.open(tmp / "queries.fa.gz", "wt") as handle:
        handle.write((FIXTURE / "queries.fa").read_text())
    gz_queries, gz_query_metrics, _ = run(tmp, "gzip_queries", queries=tmp / "queries.fa.gz")
    assert gz_queries == rows
    assert gz_query_metrics["query_sequence_fnv1a64"] == metrics["query_sequence_fnv1a64"]

    spanning_ref = tmp / "spanning.fa"
    spanning_query = tmp / "spanning_query.fa"
    spanning_ref.write_text(">chr\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
    for k in (4, 15, 18, 20, 22, 24, 32):
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        spanning_query.write_text(">boundary\n" + sequence[3:3 + k] + "\n")
        k_rows, _, _ = run(tmp, "k" + str(k), reference=spanning_ref,
                           queries=spanning_query, k=k, fragment_size=k)
        assert k_rows["boundary"]["ibf_recovery"] == "1.000000", k
        assert k_rows["boundary"]["exact_recovery"] == "1.000000", k

    threshold_rows, _, _ = run(tmp, "threshold", extra=("--basic-threshold", "0.4",
                                                      "--strict-threshold", "0.5"))
    assert threshold_rows["duplicated"]["status"] == "PASS_BASIC"
    table = tmp / "thresholds.tsv"
    table.write_text("reference_id\tk\tbasic_threshold\tstrict_threshold\ttraining_control_set\n"
                     "reference\t4\t0.4\t0.5\tshuffled_training\n")
    table_rows, table_metrics, _ = run(tmp, "threshold_table",
                                       extra=("--threshold-table", str(table),
                                              "--evaluation-control-set", "lncrna_evaluation"))
    assert table_rows["duplicated"]["status"] == "PASS_BASIC"
    assert table_rows["duplicated"]["threshold_source"] == "table"
    assert table_metrics["training_control_set"] == "shuffled_training"
    assert table_metrics["evaluation_control_set"] == "lncrna_evaluation"
    same_controls = subprocess.run([str(PROGRAM), "--reference", str(FIXTURE / "reference.fa"),
                                    "--queries", str(FIXTURE / "queries.fa"),
                                    "--output", str(tmp / "same_controls.tsv"), "--k", "4",
                                    "--fragment-size", "4", "--threshold-table", str(table),
                                    "--evaluation-control-set", "shuffled_training"],
                                   capture_output=True, text=True)
    assert same_controls.returncode != 0 and "must be nonempty and different" in same_controls.stderr

    rng = random.Random(4217)
    calibration_ref = tmp / "calibration_ref.fa"
    calibration_query = tmp / "calibration_query.fa"
    calibration_ref.write_text(">r\n" + "".join(rng.choices("ACGT", k=2000)) + "\n")
    calibration_query.write_text(">absent_control\n" +
                                 "".join(rng.choices("ACGT", k=2000)) + "\n")
    calibrated, calibration_metrics, _ = run(tmp, "calibrate", reference=calibration_ref,
                                              queries=calibration_query, k=15,
                                              fragment_size=2000,
                                              extra=("--mode", "calibrate", "--fpr", "0.5"))
    assert calibrated["absent_control"]["exact_recovery"] == "0.000000"
    assert calibrated["absent_control"]["recovery_fraction"] == "0.000000"
    assert calibrated["absent_control"]["matched_windows"] == "0"
    assert calibrated["absent_control"]["status"] == "NO_MATCH"
    assert float(calibrated["absent_control"]["ibf_recovery"]) > 0
    assert calibration_metrics["primary_match_source"] == "exact"
    assert calibration_metrics["false_positive_windows"] > 0
    assert calibration_metrics["false_negative_windows"] == 0
    assert abs(calibration_metrics["measured_genome_wide_fpr"] -
               calibration_metrics["false_positive_windows"] /
               calibration_metrics["tested_absent_windows"]) < 1e-6
    with pathlib.Path(str(tmp / "calibrate.tsv") + ".bins.tsv").open(newline="") as handle:
        calibration_bins = list(csv.DictReader(handle, delimiter="\t"))
    assert len(calibration_bins) == 1
    assert calibration_bins[0]["tested_absent_query_windows"] == "1986"
    assert calibration_bins[0]["false_positive_bin_pairs"] == str(
        calibration_metrics["false_positive_bin_pairs"])
    assert calibration_metrics["measured_per_bin_fpr"] == calibration_metrics["measured_genome_wide_fpr"]

    bin_ref = tmp / "bin_ref.fa"
    bin_query = tmp / "bin_query.fa"
    bin_ref.write_text(">chr1\nAAAACCCC\n>chr2\nGGGGTTTT\n")
    bin_query.write_text(">q\nAAAANNNNACAC\n")
    _, bin_metrics, bin_output = run(tmp, "calibrate_bins", reference=bin_ref,
                                      queries=bin_query, extra=("--mode", "calibrate", "--strand", "forward"))
    with pathlib.Path(str(bin_output) + ".bins.tsv").open(newline="") as handle:
        bins = list(csv.DictReader(handle, delimiter="\t"))
    assert [int(row["tested_absent_query_windows"]) for row in bins] == [1, 2, 2, 2]
    assert bin_metrics["tested_absent_bin_pairs"] == 7
    assert bin_metrics["false_negative_bin_pairs"] == 0
    assert sum(int(row["false_positive_bin_pairs"]) for row in bins) == bin_metrics["false_positive_bin_pairs"]
    ibf_pair_cap = subprocess.run([str(PROGRAM), "--reference", str(FIXTURE / "reference.fa"),
                                   "--queries", str(FIXTURE / "queries.fa"),
                                   "--output", str(tmp / "ibf_pair_cap.tsv"), "--k", "4",
                                   "--fragment-size", "4", "--mode", "calibrate",
                                   "--max-ibf-bin-pairs", "1"], capture_output=True, text=True)
    assert ibf_pair_cap.returncode != 0 and "exceed --max-ibf-bin-pairs" in ibf_pair_cap.stderr
    exact_pair_cap = subprocess.run([str(PROGRAM), "--reference", str(FIXTURE / "reference.fa"),
                                     "--queries", str(FIXTURE / "queries.fa"),
                                     "--output", str(tmp / "exact_pair_cap.tsv"), "--k", "4",
                                     "--fragment-size", "4", "--mode", "calibrate",
                                     "--max-calibration-pairs", "1"], capture_output=True, text=True)
    assert exact_pair_cap.returncode != 0 and "exceed --max-calibration-pairs" in exact_pair_cap.stderr

    repeat_bed_ref = tmp / "repeat_ref.bed"
    repeat_bed_query = tmp / "repeat_query.bed"
    repeat_bed_ref.write_text("r\t0\t1\tLINE\n")
    repeat_bed_query.write_text("q\t4\t5\tSINE\nq\t5\t6\tDNA\n")
    bed_ref = tmp / "bed_ref.fa"
    bed_query = tmp / "bed_query.fa"
    bed_ref.write_text(">r\nAAAACCCC\n")
    bed_query.write_text(">q\nAAAACCCC\n")
    bed_options = ("--strand", "forward", "--reference-repeat-bed", str(repeat_bed_ref),
                   "--query-repeat-bed", str(repeat_bed_query))
    annotated_bed, _, annotated_path = run(tmp, "annotate_bed", reference=bed_ref,
                                            queries=bed_query, extra=(*bed_options, "--repeat", "annotate"))
    assert annotated_bed["q"]["repeat_fraction"] == "0.800000"
    assert annotated_bed["q"]["valid_nonrepeat_windows"] == "1"
    with pathlib.Path(str(annotated_path) + ".repeat_classes.tsv").open(newline="") as handle:
        classes = {row["repeat_class"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert classes["SINE"]["valid_windows"] == "1"
    assert classes["MULTIPLE_CLASSES"]["valid_windows"] == "3"
    reference_only, _, _ = run(tmp, "exclude_reference_bed", reference=bed_ref, queries=bed_query,
                                extra=(*bed_options, "--repeat", "exclude", "--mask-target", "reference"))
    assert reference_only["q"]["valid_windows"] == "5"
    assert reference_only["q"]["exact_recovery"] == "0.800000"
    query_only, _, _ = run(tmp, "exclude_query_bed", reference=bed_ref, queries=bed_query,
                            extra=(*bed_options, "--repeat", "exclude", "--mask-target", "query"))
    assert query_only["q"]["valid_windows"] == "1"
    assert query_only["q"]["exact_recovery"] == "1.000000"
    both_masked, _, _ = run(tmp, "exclude_both_bed", reference=bed_ref, queries=bed_query,
                             extra=(*bed_options, "--repeat", "exclude"))
    assert both_masked["q"]["valid_windows"] == "1"
    assert both_masked["q"]["exact_recovery"] == "0.000000"
    invalid_bed = tmp / "invalid.bed"
    invalid_bed.write_text("wrong_contig\t0\t1\tLINE\n")
    bad_contig = subprocess.run([str(PROGRAM), "--reference", str(bed_ref),
                                 "--queries", str(bed_query), "--output", str(tmp / "bad_bed.tsv"),
                                 "--k", "4", "--fragment-size", "4",
                                 "--reference-repeat-bed", str(invalid_bed)], capture_output=True, text=True)
    assert bad_contig.returncode != 0 and "missing from reference" in bad_contig.stderr
    too_long = subprocess.run([str(PROGRAM), "--reference", str(bed_ref),
                               "--queries", str(bed_query), "--output", str(tmp / "long.tsv"),
                               "--k", "4", "--fragment-size", "4", "--max-query-bases", "4"],
                              capture_output=True, text=True)
    assert too_long.returncode != 0 and "exceeds --max-query-bases" in too_long.stderr

    positions_path = tmp / "positions.tsv.gz"
    localized, _, _ = run(tmp, "localized", extra=("--local-window", "4",
                                                  "--positions", str(positions_path)))
    assert localized["reverse_only"]["candidate_contig"] == "ctg1"
    assert localized["reverse_only"]["candidate_start"] == "0"
    assert localized["reverse_only"]["candidate_end"] == "4"
    assert localized["reverse_only"]["multimapped_query_windows"] == "1"
    with gzip.open(positions_path, "rt", newline="") as handle:
        positions = list(csv.DictReader(handle, delimiter="\t"))
    reverse_positions = [p for p in positions if p["query_id"] == "reverse_only"]
    assert [(p["contig"], p["start"], p["end"], p["strand"]) for p in reverse_positions] == [
        ("ctg1", "0", "4", "-"), ("ctg3", "0", "4", "-")]
    assert [p["strand"] for p in positions if p["query_id"] == "palindrome"] == ["."]

    spread_ref = tmp / "spread_ref.fa"
    spread_query = tmp / "spread_query.fa"
    spread_ref.write_text(">c1\nAAAACCCC\n>c2\nGGGGTTTT\n")
    spread_query.write_text(">concentrated\nAAAACCCC\n>dispersed\nAAAAGGGG\n")
    spread, _, _ = run(tmp, "spread", reference=spread_ref, queries=spread_query,
                       extra=("--strand", "forward", "--local-window", "10"))
    assert spread["concentrated"]["global_recovery"] == "1.000000"
    assert spread["concentrated"]["best_region_recovery"] == "1.000000"
    assert spread["dispersed"]["global_recovery"] == "0.400000"
    assert spread["dispersed"]["best_region_recovery"] == "0.200000"
    assert spread["dispersed"]["local_concentration"] == "0.500000"
    assert spread["dispersed"]["off_region_fraction"] == "0.500000"
    assert spread["dispersed"]["number_of_supported_regions"] == "2"
    boundary_ref = tmp / "boundary_ref.fa"
    boundary_query = tmp / "boundary_query.fa"
    boundary_ref.write_text(">chr\nAAAACCCC\n")
    boundary_query.write_text(">boundary\nAAACCC\n")
    boundary, _, _ = run(tmp, "boundary_local", reference=boundary_ref,
                          queries=boundary_query, k=3, extra=("--strand", "forward", "--local-window", "4"))
    assert boundary["boundary"]["best_region_recovery"] == "0.750000"
    assert boundary["boundary"]["sliding_best_region_recovery"] == "1.000000"
    assert boundary["boundary"]["boundary_gain_windows"] == "1"

    splice_ref = tmp / "splice_ref.fa"
    splice_query = tmp / "splice_query.fa"
    splice_blocks = tmp / "splice_blocks.tsv"
    splice_ref.write_text(">chr\nAAAANNNNCCCCNNNNGGGG\n")
    splice_query.write_text(">plus\nAAAACCCC\n>minus\nAAAACCCC\n")
    splice_blocks.write_text("query_id\tquery_start\tquery_end\tcontig\tref_start\tref_end\tstrand\n"
                             "plus\t0\t4\tchr\t0\t4\t+\n"
                             "plus\t4\t8\tchr\t8\t12\t+\n"
                             "minus\t0\t4\tchr\t0\t4\t+\n"
                             "minus\t4\t8\tchr\t16\t20\t-\n")
    spliced, spliced_metrics, _ = run(tmp, "splice", reference=splice_ref, queries=splice_query,
                                       extra=("--exon-blocks", str(splice_blocks), "--local-window", "4"))
    for row in spliced.values():
        assert row["exact_recovery"] == "0.400000"
        assert row["best_region_recovery"] == "0.200000"
        assert row["exon_eligible_windows"] == "2"
        assert row["exon_matched_windows"] == "2"
        assert row["exon_aware_recovery"] == "1.000000"
    assert spliced_metrics["exon_eligible_windows"] == 4
    assert spliced_metrics["exon_aware_recovery"] == 1.0

    overflow = subprocess.run([str(PROGRAM), "--reference", str(FIXTURE / "reference.fa"),
                               "--queries", str(FIXTURE / "queries.fa"),
                               "--output", str(tmp / "overflow.tsv"), "--k", "4",
                               "--fragment-size", "4", "--local-window", "4",
                               "--max-position-rows", "1"], capture_output=True, text=True)
    assert overflow.returncode != 0 and "exceeds --max-position-rows" in overflow.stderr

    collision = subprocess.run([str(PROGRAM), "--reference", str(FIXTURE / "reference.fa"),
                                "--queries", str(FIXTURE / "queries.fa"),
                                "--output", str(FIXTURE / "queries.fa")],
                               capture_output=True, text=True)
    assert collision.returncode != 0 and "paths must be distinct" in collision.stderr

    error = subprocess.run([str(PROGRAM), "--reference", str(FIXTURE / "reference.fa"),
                            "--queries", str(FIXTURE / "queries.fa"),
                            "--output", str(tmp / "error.tsv"), "--k", "4",
                            "--fragment-size", "4", "--ambiguity", "error"],
                           capture_output=True, text=True)
    assert error.returncode != 0 and "Ambiguous base" in error.stderr

    class_ref = tmp / "class_ref.fa"
    class_query = tmp / "class_query.fa"
    class_bed = tmp / "class_ref.bed"
    class_ref.write_text(">chr\naaaaCCCCAAAA\n")
    class_query.write_text(">q\nAAAA\n")
    class_bed.write_text("chr\t0\t1\tLINE\nchr\t8\t9\tSINE\n")
    fragment_path = tmp / "fragments.tsv"
    hit_class_path = tmp / "reference_hits.tsv"
    class_rows, class_metrics, _ = run(
        tmp, "reference_classes", reference=class_ref, queries=class_query,
        extra=("--mode", "exact", "--strand", "forward", "--repeat", "annotate",
               "--reference-repeat-bed", str(class_bed),
               "--reference-hit-classes", str(hit_class_path), "--fragments", str(fragment_path)))
    assert class_rows["q"]["matched_windows"] == "1"
    with fragment_path.open(newline="") as handle:
        fragments = list(csv.DictReader(handle, delimiter="\t"))
    assert [(r["contig"], r["start"], r["end"], r["sequence_end"]) for r in fragments] == [
        ("chr", "0", "4", "7"), ("chr", "4", "8", "11"), ("chr", "8", "12", "12")]
    assert all(r["inserted_reference_windows"] == "NA" for r in fragments)
    with hit_class_path.open(newline="") as handle:
        hit_classes = {r["reference_repeat_class"]: r for r in csv.DictReader(handle, delimiter="\t")}
    assert set(hit_classes) == {"SINE", "MULTIPLE_CLASSES"}
    assert all(r["matched_query_windows"] == "1" and r["exact_occurrences"] == "1" and
               r["fraction_of_exact_matched_windows"] == "1.000000" for r in hit_classes.values())
    assert class_metrics["fragment_output_seconds"] >= 0
    assert class_metrics["peak_rss_after_fragment_output_bytes"] > 0

    ibf_fragments = tmp / "ibf_fragments.tsv"
    ibf_rows, ibf_metrics, _ = run(
        tmp, "ibf_target", reference=class_ref, queries=class_query,
        extra=("--mode", "ibf", "--strand", "forward", "--fragments", str(ibf_fragments)))
    assert ibf_rows["q"]["matched_windows"] == "1"
    assert abs(ibf_metrics["configured_per_bin_fpr"] -
               ibf_metrics["target_genome_wide_fpr"] / 3) < 1e-12
    assert ibf_metrics["observed_occupancy_union_estimate"] <= ibf_metrics["target_genome_wide_fpr"]
    assert ibf_metrics["occupancy_guard_seconds"] >= 0
    with ibf_fragments.open(newline="") as handle:
        indexed_fragments = list(csv.DictReader(handle, delimiter="\t"))
    assert indexed_fragments[0]["inserted_reference_windows"] == "4"
    assert indexed_fragments[0]["distinct_reference_kmers"] == "4"
    assert indexed_fragments[2]["inserted_reference_windows"] == "1"

    unsafe_index = tmp / "unsafe.ibf"
    run(tmp, "unsafe_index", reference=class_ref, queries=class_query,
        extra=("--mode", "validate", "--strand", "forward", "--save-index", str(unsafe_index)))
    unsafe_load = subprocess.run(
        [str(PROGRAM), "--reference", str(class_ref), "--queries", str(class_query),
         "--output", str(tmp / "unsafe_load.tsv"), "--k", "4", "--fragment-size", "4",
         "--mode", "ibf", "--strand", "forward", "--load-index", str(unsafe_index)],
        capture_output=True, text=True)
    assert unsafe_load.returncode != 0 and "configured_fpr" in unsafe_load.stderr

print("kmer v2 integration checks passed")
