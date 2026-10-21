"""Check streaming multi-reference merge, settings checks, and resource report."""

import csv
import json
import pathlib
import subprocess
import sys
import tempfile


PROGRAM = pathlib.Path(sys.argv[1])
MERGER = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "merge_kmer_v2.py"


with tempfile.TemporaryDirectory(prefix="fraglncrna-merge-test-") as directory:
    tmp = pathlib.Path(directory)
    queries = tmp / "queries.fa"
    queries.write_text(">matched\nAAAA\n>absent\nACAC\n")
    references = []
    for name, sequence in (("ref_a", "AAAA"), ("ref_b", "CCCC")):
        reference = tmp / f"{name}.fa"
        reference.write_text(f">chr\n{sequence}\n")
        output = tmp / f"{name}.tsv"
        subprocess.run([str(PROGRAM), "--reference", str(reference), "--queries", str(queries),
                        "--output", str(output), "--k", "4", "--fragment-size", "4"],
                       check=True, capture_output=True, text=True)
        references.append(output)
    merged = tmp / "merged.tsv.gz"
    summary = tmp / "summary.json"
    command = [sys.executable, str(MERGER), "--input", str(references[1]), "--input", str(references[0]),
               "--output", str(merged), "--summary", str(summary)]
    subprocess.run(command, check=True, capture_output=True, text=True)
    import gzip
    with gzip.open(merged, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["query_id"], row["reference_id"]) for row in rows] == [
        ("matched", "ref_a"), ("matched", "ref_b"), ("absent", "ref_a"), ("absent", "ref_b")]
    assert [row["matched_windows"] for row in rows] == ["1", "0", "0", "0"]
    audit = json.loads(summary.read_text())
    assert audit["query_count"] == 2 and audit["reference_count"] == 2
    assert audit["primary_match_source"] == "exact"
    assert audit["references"]["ref_a"]["statuses"]["PASS_STRICT"] == 1
    assert audit["references"]["ref_b"]["statuses"]["NO_MATCH"] == 2
    assert audit["wall_seconds"] > 0 and audit["peak_resident_bytes"] > 0
    second = tmp / "second.tsv.gz"
    subprocess.run([*command[:-4], "--output", str(second), "--summary", str(tmp / "second.json")],
                   check=True, capture_output=True, text=True)
    with gzip.open(second, "rt") as handle:
        assert rows == list(csv.DictReader(handle, delimiter="\t"))
    assert second.read_bytes() == merged.read_bytes()
    second_manifest = pathlib.Path(str(references[1]) + ".metrics.json")
    original_manifest = second_manifest.read_text()
    changed = json.loads(original_manifest)
    changed["primary_match_source"] = "ibf"
    second_manifest.write_text(json.dumps(changed))
    mixed_source = subprocess.run(command, capture_output=True, text=True)
    assert mixed_source.returncode != 0 and (
        "Primary match sources differ" in mixed_source.stderr or
        "Missing IBF-wide FPR sizing provenance" in mixed_source.stderr)
    changed.pop("primary_match_source")
    second_manifest.write_text(json.dumps(changed))
    missing_source = subprocess.run(command, capture_output=True, text=True)
    assert missing_source.returncode != 0 and "primary_match_source" in missing_source.stderr
    second_manifest.write_text(original_manifest)
    altered_queries = tmp / "altered.fa"
    altered_queries.write_text(">matched\nAAAT\n>absent\nACAC\n")
    subprocess.run([str(PROGRAM), "--reference", str(tmp / "ref_b.fa"), "--queries", str(altered_queries),
                    "--output", str(references[1]), "--k", "4", "--fragment-size", "4"],
                   check=True, capture_output=True, text=True)
    mismatch = subprocess.run(command, capture_output=True, text=True)
    assert mismatch.returncode != 0 and "hashes differ" in mismatch.stderr

    ibf_outputs = []
    for name, sequence in (("ibf_a", "AAAA"), ("ibf_b", "CCCCCCCC")):
        reference = tmp / f"{name}.fa"
        reference.write_text(f">chr\n{sequence}\n")
        output = tmp / f"{name}.tsv"
        subprocess.run([str(PROGRAM), "--reference", str(reference), "--queries", str(queries),
                        "--output", str(output), "--k", "4", "--fragment-size", "4", "--mode", "ibf"],
                       check=True, capture_output=True, text=True)
        ibf_outputs.append(output)
    ibf_summary = tmp / "ibf_summary.json"
    subprocess.run([sys.executable, str(MERGER), "--input", str(ibf_outputs[0]),
                    "--input", str(ibf_outputs[1]), "--output", str(tmp / "ibf_merged.tsv"),
                    "--summary", str(ibf_summary)], check=True, capture_output=True, text=True)
    ibf_a = json.loads(pathlib.Path(str(ibf_outputs[0]) + ".metrics.json").read_text())
    ibf_b = json.loads(pathlib.Path(str(ibf_outputs[1]) + ".metrics.json").read_text())
    assert ibf_a["configured_per_bin_fpr"] != ibf_b["configured_per_bin_fpr"]
    assert json.loads(ibf_summary.read_text())["ibf_target_genome_wide_fpr"] == 0.01

print("multi-reference merge checks passed")
