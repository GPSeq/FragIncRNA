#!/usr/bin/env python3
"""Fit reference/k recovery thresholds on one control set and audit a separate set."""

import argparse
import csv
import json
import math
import pathlib
import resource
import sys
import time


def peak_rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else rss * 1024)


def read_controls(path):
    grouped = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"reference_id", "k", "query_id", "recovery_fraction", "schema_version"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"Missing schema 2 control columns in {path}")
        for line_no, row in enumerate(reader, 2):
            if row["schema_version"] != "2":
                raise ValueError(f"Expected schema 2 on {path}:{line_no}")
            if not row["query_id"] or not row["reference_id"]:
                raise ValueError(f"Empty identifier on {path}:{line_no}")
            if row["recovery_fraction"] == "NA":
                continue  # A sequence without valid k-mers cannot be classified by a fraction.
            value = float(row["recovery_fraction"])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Invalid recovery fraction on {path}:{line_no}")
            key = (row["reference_id"], int(row["k"]))
            if key[1] < 1:
                raise ValueError(f"Invalid k on {path}:{line_no}")
            group = grouped.setdefault(key, {})
            if row["query_id"] in group:
                raise ValueError(f"Duplicate control ID for {key} on {path}:{line_no}")
            group[row["query_id"]] = value
    if not grouped:
        raise ValueError(f"No classifiable control rows in {path}")
    return grouped


def fit_threshold(values, target_rate):
    scores = sorted(values, reverse=True)
    allowed = math.floor(target_rate * len(scores))
    if allowed == len(scores):
        return 0.0
    threshold = math.nextafter(scores[allowed], math.inf)
    if threshold > 1:
        raise ValueError("Cannot achieve target null pass rate: control recovery reaches 1")
    return threshold


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training", required=True, type=pathlib.Path)
    parser.add_argument("--evaluation", required=True, type=pathlib.Path)
    parser.add_argument("--training-label", required=True)
    parser.add_argument("--evaluation-label", required=True)
    parser.add_argument("--basic-null-rate", type=float, default=0.05)
    parser.add_argument("--strict-null-rate", type=float, default=0.01)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--summary", required=True, type=pathlib.Path)
    args = parser.parse_args()
    if (not args.training_label or not args.evaluation_label or
            args.training_label == args.evaluation_label):
        parser.error("Training and evaluation labels must be nonempty and different")
    if not (0 <= args.strict_null_rate <= args.basic_null_rate < 1):
        parser.error("Require 0 <= strict-null-rate <= basic-null-rate < 1")
    if args.output.resolve() in {args.training.resolve(), args.evaluation.resolve(), args.summary.resolve()} or \
       args.summary.resolve() in {args.training.resolve(), args.evaluation.resolve()}:
        parser.error("Input and output paths must be distinct")
    training = read_controls(args.training)
    evaluation = read_controls(args.evaluation)
    if training.keys() != evaluation.keys():
        raise ValueError("Training and evaluation must cover the same reference/k groups")
    threshold_rows = []
    audit_rows = []
    for reference_id, k in sorted(training):
        key = (reference_id, k)
        overlap = training[key].keys() & evaluation[key].keys()
        if overlap:
            raise ValueError(f"Training/evaluation query IDs overlap for {key}: {len(overlap)} IDs")
        basic = fit_threshold(training[key].values(), args.basic_null_rate)
        strict = fit_threshold(training[key].values(), args.strict_null_rate)
        threshold_rows.append((reference_id, k, basic, strict, args.training_label))
        audit_rows.append({
            "reference_id": reference_id,
            "k": k,
            "training_count": len(training[key]),
            "evaluation_count": len(evaluation[key]),
            "basic_threshold": basic,
            "strict_threshold": strict,
            "training_basic_pass_count": sum(x >= basic for x in training[key].values()),
            "training_strict_pass_count": sum(x >= strict for x in training[key].values()),
            "evaluation_basic_pass_count": sum(x >= basic for x in evaluation[key].values()),
            "evaluation_strict_pass_count": sum(x >= strict for x in evaluation[key].values()),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["reference_id", "k", "basic_threshold", "strict_threshold",
                         "training_control_set"])
        for reference_id, k, basic, strict, label in threshold_rows:
            writer.writerow([reference_id, k, format(basic, ".17g"), format(strict, ".17g"), label])
    summary = {
        "schema_version": 1,
        "training": str(args.training),
        "evaluation": str(args.evaluation),
        "training_label": args.training_label,
        "evaluation_label": args.evaluation_label,
        "target_basic_null_pass_rate": args.basic_null_rate,
        "target_strict_null_pass_rate": args.strict_null_rate,
        "groups": audit_rows,
        "wall_seconds": time.monotonic() - started,
        "peak_resident_bytes": peak_rss_bytes(),
    }
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {len(audit_rows)} threshold groups; wall={summary['wall_seconds']:.3f}s "
          f"peak RSS={summary['peak_resident_bytes']}B")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as error:
        raise SystemExit(str(error)) from error
