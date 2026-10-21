"""Demonstrate and preserve the historical bin-summed score semantics."""
import csv
import pathlib
import subprocess
import sys
import tempfile


with tempfile.TemporaryDirectory(prefix="fraglncrna-legacy-test-") as directory:
    root = pathlib.Path(directory)
    refs = root / "refs"
    refs.mkdir()
    (refs / "reference.fa").write_text(">chr\nACGTACGTAC\n")
    (root / "queries.fa").write_text(">q\nACGTAC\n")
    out = root / "out"
    config = root / "config.toml"
    config.write_text(f"""[general]
index_method = "ibf"
ref_dir = "{refs}"
query_file = "{root / 'queries.fa'}"
output_dir = "{out}"
output_file = "results.tsv"
log_file = "run.log"
fragment_size = 5
kmer_size = 3
hit_threshold = 1
threads = 1
single_results_writer = false
[ibf]
hash_functions = 3
fpr = 0.01
multi_indexing = false
""")
    subprocess.run([sys.argv[1], str(config)], check=True, capture_output=True, text=True)
    with (out / "results_reference.tsv").open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert int(row["reference_count"]) > 4
    assert float(row["reference_pct"]) > 1

print(f"legacy score reproduced: count={row['reference_count']}, pct={row['reference_pct']}")
