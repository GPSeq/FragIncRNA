#!/usr/bin/env python3
"""Summarize measured FragIncRNA per-bin comparison runs without rerunning them."""

import argparse
import csv
import json
from pathlib import Path


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def value(item):
    return "NA" if item is None else item


def ratio(numerator, denominator):
    return "NA" if denominator == 0 else f"{numerator / denominator:.6f}"


def union_status(item):
    if int(item["valid_windows"]) == 0:
        return "NO_VALID_KMERS"
    recovery = float(item["ibf_recovery"])
    if recovery == 0:
        return "NO_MATCH"
    if recovery >= float(item["strict_threshold"]):
        return "PASS_STRICT"
    if recovery >= float(item["basic_threshold"]):
        return "PASS_BASIC"
    return "LOW_SUPPORT"


def paired_rows(root):
    exact_path = root / "T_exact_100.tsv"
    if not exact_path.exists():
        return []
    exact = {row["query_id"]: row for row in read_tsv(exact_path)}
    data = {}
    for name in ("A0_high_union_100", "A1_budgeted_union_100"):
        path = root / f"{name}.tsv"
        if path.exists():
            data[name] = {row["query_id"]: row for row in read_tsv(path)}
    for name in ("B0_high_coherence_unverified_100", "B1_high_coherence_verified_100",
                 "B2_budgeted_coherence_verified_100"):
        path = root / f"{name}.tsv.per_bin.tsv"
        if path.exists():
            data[name] = {row["query_id"]: row for row in read_tsv(path)}
    rows = []
    for query_id, truth in exact.items():
        row = {"query_id": query_id, "k": truth["k"],
               "exact_global_recovery": truth["recovery_fraction"],
               "exact_global_matched_windows": truth["matched_windows"],
               "exact_status": truth["status"]}
        for name, short in (("A0_high_union_100", "A0"), ("A1_budgeted_union_100", "A1")):
            item = data.get(name, {}).get(query_id)
            row[f"{short}_ibf_union_recovery"] = item["ibf_recovery"] if item else "NA"
            row[f"{short}_union_status"] = union_status(item) if item else "NA"
        for name, short in (("B0_high_coherence_unverified_100", "B0"),
                            ("B1_high_coherence_verified_100", "B1"),
                            ("B2_budgeted_coherence_verified_100", "B2")):
            item = data.get(name, {}).get(query_id)
            for source, suffix in (("max_bin_support", "max_bin_support"),
                                   ("candidate_bin_count", "candidate_bin_count"),
                                   ("candidate_locus_count", "candidate_locus_count"),
                                   ("best_locus_exact_recovery", "best_locus_exact_recovery"),
                                   ("global_exact_windows_outside_candidates", "exact_windows_missed"),
                                   ("per_bin_status", "per_bin_status")):
                row[f"{short}_{suffix}"] = item[source] if item else "NA"
        rows.append(row)
    return rows


def summarize(root):
    rows = []
    truth_metrics_path = root / "T_exact_100.tsv.metrics.json"
    truth_metrics = json.loads(truth_metrics_path.read_text()) if truth_metrics_path.exists() else None
    for command_path in sorted(root.glob("*.command.json")):
        name = command_path.name.removesuffix(".command.json")
        if name.startswith("occupancy_"):
            continue
        metrics_path = root / f"{name}.tsv.metrics.json"
        if not metrics_path.exists():
            continue
        command = json.loads(command_path.read_text())
        metrics = json.loads(metrics_path.read_text())
        comparable_truth = (truth_metrics and metrics["k"] == truth_metrics["k"] and
                            metrics["reference_sequence_sha256"] == truth_metrics["reference_sequence_sha256"] and
                            metrics["query_sequence_sha256"] == truth_metrics["query_sequence_sha256"])
        truth_source = truth_metrics if metrics["mode"] == "ibf" and comparable_truth else metrics
        paired_path = root / f"{name}.tsv.per_bin.tsv"
        paired = read_tsv(paired_path) if paired_path.exists() else []
        truth_positive = predicted_positive = tp = fp = fn = tn = 0
        if paired and metrics["mode"] != "ibf":
            for item in paired:
                truth = int(item["exact_global_matched_windows"]) > 0
                predicted = item["per_bin_status"] == "VERIFIED_CANDIDATE"
                truth_positive += truth
                predicted_positive += predicted
                tp += truth and predicted
                fp += not truth and predicted
                fn += truth and not predicted
                tn += not truth and not predicted
        wall = command["wall_seconds"]
        row = {
            "run": name, "k": metrics["k"], "queries": metrics["query_count"],
            "valid_windows": metrics["valid_query_windows"],
            "exact_present_windows": truth_source["exact_positive_windows"] if truth_source is not metrics else
                                     (metrics["exact_positive_windows"] if metrics["mode"] != "ibf" else "NA"),
            "exact_absent_windows": truth_source["tested_absent_windows"] if truth_source is not metrics else
                                    (metrics["tested_absent_windows"] if metrics["mode"] != "ibf" else "NA"),
            "window_union_fp": value(metrics["false_positive_windows"]),
            "window_union_fpr": value(metrics["measured_genome_wide_fpr"]),
            "window_union_fn": value(metrics["false_negative_windows"]),
            "candidate_queries": metrics.get("per_bin_queries_with_candidates", 0),
            "candidate_bins": metrics.get("per_bin_candidate_bins", 0),
            "candidate_loci": metrics.get("per_bin_candidate_loci", 0),
            "verified_loci": metrics.get("per_bin_verified_candidate_loci", 0),
            "rejected_loci": metrics.get("per_bin_rejected_candidate_loci", 0),
            "rejected_locus_fraction": ratio(metrics.get("per_bin_rejected_candidate_loci", 0),
                                             metrics.get("per_bin_candidate_loci", 0))
                                       if paired and metrics["mode"] != "ibf" else "NA",
            "missed_exact_queries": metrics.get("per_bin_missed_global_exact_queries", 0)
                                    if paired and metrics["mode"] != "ibf" else "NA",
            "exact_windows_outside_candidates": metrics.get("per_bin_global_exact_windows_outside_candidates", 0)
                                                if paired and metrics["mode"] != "ibf" else "NA",
            "top_k_dropped_bins": metrics.get("per_bin_top_k_dropped_bins", 0),
            "technical_endpoint": "any global exact k-mer versus any verified candidate",
            "TP": tp if paired and metrics["mode"] != "ibf" else "NA",
            "FP": fp if paired and metrics["mode"] != "ibf" else "NA",
            "FN": fn if paired and metrics["mode"] != "ibf" else "NA",
            "TN": tn if paired and metrics["mode"] != "ibf" else "NA",
            "technical_precision": ratio(tp, tp + fp) if paired and metrics["mode"] != "ibf" else "NA",
            "technical_recall": ratio(tp, tp + fn) if paired and metrics["mode"] != "ibf" else "NA",
            "technical_specificity": ratio(tn, tn + fp) if paired and metrics["mode"] != "ibf" else "NA",
            "wall_seconds": f"{wall:.3f}", "peak_rss_bytes": command["peak_child_resident_bytes"],
            "queries_per_second": f"{metrics['query_count'] / wall:.3f}" if wall else "NA",
            "valid_windows_per_second": f"{metrics['valid_query_windows'] / wall:.3f}" if wall else "NA",
            "index_load_seconds_internal": metrics["index_load_seconds"],
            "index_build_seconds_internal": metrics["index_construction_seconds"],
            "index_save_seconds_internal": metrics["index_save_seconds"],
            "query_seconds_internal": metrics["query_seconds"],
            "exact_scan_seconds_internal": metrics["exact_scan_seconds"],
            "index_bytes": metrics["disk_index_bytes"], "bin_count": metrics["bin_count"],
            "bin_bits": metrics["bin_bits"],
            "effective_per_bin_fpr": metrics["configured_per_bin_fpr"],
            "target_whole_reference_fpr": value(metrics["target_genome_wide_fpr"]),
            "index_sha256": metrics["index_file_sha256"],
            "reference_semantic_sha256": metrics["reference_sequence_sha256"],
            "query_semantic_sha256": metrics["query_sequence_sha256"],
        }
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root
    summary = summarize(root)
    if summary:
        write_tsv(root / "summary.tsv", summary, list(summary[0]))
    paired = paired_rows(root)
    if paired:
        write_tsv(root / "paired_100.tsv", paired, list(paired[0]))
    grid = [row for row in summary if row["run"].startswith("grid_high_") or
            row["run"] == "B1_high_coherence_verified_100"]
    if grid:
        write_tsv(root / "threshold_grid.tsv", grid, list(grid[0]))
    occupancy = []
    for label in ("high", "low"):
        bins_path = root / f"occupancy_{label}.bins.tsv"
        metrics_path = root / f"occupancy_{label}.tsv.metrics.json"
        command_path = root / f"occupancy_{label}.command.json"
        if not all(path.exists() for path in (bins_path, metrics_path, command_path)):
            continue
        bins = read_tsv(bins_path)
        metrics = json.loads(metrics_path.read_text())
        command = json.loads(command_path.read_text())
        fractions = sorted(float(item["bit_occupancy_fraction"]) for item in bins)
        occupied = sum(int(item["occupied_bits"]) for item in bins)
        possible = sum(int(item["bin_bits"]) for item in bins)
        occupancy.append({"index": label, "bins": len(bins), "occupied_bits": occupied,
                          "total_bin_bits": possible,
                          "aggregate_occupancy_fraction": ratio(occupied, possible),
                          "median_bin_occupancy_fraction": f"{fractions[len(fractions) // 2]:.6f}",
                          "max_bin_occupancy_fraction": f"{fractions[-1]:.6f}",
                          "occupancy_scan_seconds_internal": metrics["bit_occupancy_measurement_seconds"],
                          "run_wall_seconds": f"{command['wall_seconds']:.3f}",
                          "peak_rss_bytes": command["peak_child_resident_bytes"],
                          "index_bytes": metrics["disk_index_bytes"],
                          "index_sha256": metrics["index_file_sha256"]})
    if occupancy:
        write_tsv(root / "occupancy_summary.tsv", occupancy, list(occupancy[0]))
    manifest = []
    for row in summary:
        name = row["run"]
        metrics = json.loads((root / f"{name}.tsv.metrics.json").read_text())
        manifest.append({"id": name, "command_json": f"{name}.command.json",
                         "result_tsv": f"{name}.tsv", "metrics_json": f"{name}.tsv.metrics.json",
                         "per_bin_tsv": metrics.get("per_bin_output") or None,
                         "index_sha256": row["index_sha256"],
                         "reference_semantic_sha256": row["reference_semantic_sha256"],
                         "query_semantic_sha256": row["query_semantic_sha256"]})
    (root / "aggregate.json").write_text(json.dumps(
        {"schema_version": 1, "runs": len(summary), "paired_queries": len(paired),
         "source": "measured command and v2 metrics JSON; no inferred biological truth",
         "method_ids": {"A0": "historical global union", "A1": "budgeted global union",
                        "B0": "historical per-bin unverified", "B1": "historical per-bin exact verified",
                        "B2": "budgeted per-bin exact verified", "T": "full-reference exact"},
         "run_manifest": manifest}, indent=2) + "\n")
    print(f"Summarized {len(summary)} measured runs and {len(paired)} paired query IDs")


if __name__ == "__main__":
    main()
