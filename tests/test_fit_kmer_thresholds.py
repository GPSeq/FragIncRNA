"""Check holdout separation and deterministic null threshold fitting."""

import csv
import json
import pathlib
import subprocess
import sys
import tempfile


SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "fit_kmer_thresholds.py"


def write_scores(path, scores):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["schema_version", "reference_id", "k", "query_id", "recovery_fraction"])
        for name, score in scores:
            writer.writerow(["2", "reference", "4", name, score])


with tempfile.TemporaryDirectory(prefix="fraglncrna-fit-test-") as directory:
    tmp = pathlib.Path(directory)
    training = tmp / "train.tsv"
    evaluation = tmp / "eval.tsv"
    thresholds = tmp / "thresholds.tsv"
    summary = tmp / "summary.json"
    write_scores(training, [("train0", "0.0"), ("train1", "0.1"),
                            ("train2", "0.2"), ("train3", "0.3")])
    write_scores(evaluation, [("eval0", "0.05"), ("eval1", "0.35")])
    command = [sys.executable, str(SCRIPT), "--training", str(training),
               "--evaluation", str(evaluation), "--training-label", "train_set",
               "--evaluation-label", "held_out", "--basic-null-rate", "0.25",
               "--strict-null-rate", "0", "--output", str(thresholds),
               "--summary", str(summary)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    with thresholds.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert float(rows[0]["basic_threshold"]) > 0.2
    assert float(rows[0]["basic_threshold"]) < 0.3
    assert float(rows[0]["strict_threshold"]) > 0.3
    assert rows[0]["training_control_set"] == "train_set"
    audit = json.loads(summary.read_text())
    assert audit["groups"][0]["training_basic_pass_count"] == 1
    assert audit["groups"][0]["training_strict_pass_count"] == 0
    assert audit["groups"][0]["evaluation_basic_pass_count"] == 1
    assert audit["wall_seconds"] > 0
    assert audit["peak_resident_bytes"] > 0
    again = tmp / "again.tsv"
    subprocess.run([*command[:-4], "--output", str(again), "--summary", str(tmp / "again.json")],
                   check=True, capture_output=True, text=True)
    assert thresholds.read_text() == again.read_text()
    write_scores(evaluation, [("train0", "0.05")])
    overlap = subprocess.run(command, capture_output=True, text=True)
    assert overlap.returncode != 0 and "overlap" in overlap.stderr

print("threshold fitting checks passed")
