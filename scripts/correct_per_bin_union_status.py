#!/usr/bin/env python3
"""Repair the union_status label in measured preview per-bin TSVs.

Runs made before the corrected writer accidentally put exact primary status in
union_status for validate/calibrate mode. All other fields were written
correctly. This script preserves a backup, changes only that column, and
records hashes so the measured raw output remains auditable.
"""

import argparse
import csv
import hashlib
import json
import pathlib
import shutil


def digest(path):
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            sha.update(block)
    return sha.hexdigest()


def classify(matched, valid, basic, strict):
    if valid == 0:
        return "NO_VALID_KMERS"
    if matched == 0:
        return "NO_MATCH"
    if matched / valid >= strict:
        return "PASS_STRICT"
    if matched / valid >= basic:
        return "PASS_BASIC"
    return "LOW_SUPPORT"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=pathlib.Path)
    parser.add_argument("--backup-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    args.backup_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for metrics_path in sorted(args.root.glob("*.tsv.metrics.json")):
        metrics = json.loads(metrics_path.read_text())
        if not metrics.get("fdr_control_per_bin") or metrics["mode"] == "ibf":
            continue
        source = pathlib.Path(metrics["per_bin_output"])
        if not source.is_absolute():
            source = pathlib.Path.cwd() / source
        if not source.exists():
            continue
        backup = args.backup_dir / source.name
        if backup.exists():
            raise RuntimeError(f"Backup already exists: {backup}")
        original_sha = digest(source)
        shutil.copyfile(source, backup)
        temporary = source.with_name(source.name + ".union_status.partial")
        changed = rows = 0
        with source.open(newline="") as input_handle, temporary.open("w", newline="") as output_handle:
            reader = csv.DictReader(input_handle, delimiter="\t")
            writer = csv.DictWriter(output_handle, fieldnames=reader.fieldnames,
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                rows += 1
                correct = classify(int(row["ibf_union_matched_windows"]), int(row["valid_windows"]),
                                   metrics["basic_threshold"], metrics["strict_threshold"])
                changed += row["union_status"] != correct
                row["union_status"] = correct
                writer.writerow(row)
        temporary.replace(source)
        records.append({"path": str(source), "backup": str(backup), "rows": rows,
                        "changed_union_status_rows": changed, "raw_sha256": original_sha,
                        "corrected_sha256": digest(source)})
    output = args.root / "union_status_corrections.json"
    output.write_text(json.dumps({"schema_version": 1, "records": records}, indent=2) + "\n")
    print(f"Checked {len(records)} verified per-bin TSVs; changed "
          f"{sum(r['changed_union_status_rows'] for r in records)} union-status cells")


if __name__ == "__main__":
    main()
