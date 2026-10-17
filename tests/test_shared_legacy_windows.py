"""Flat legacy and negative_set executables use the shared position-safe windows."""

import csv
import pathlib
import subprocess
import sys
import tempfile


LEGACY = pathlib.Path(sys.argv[1])
NEGATIVE = pathlib.Path(sys.argv[2])


with tempfile.TemporaryDirectory(prefix="fraglncrna-shared-windows-") as directory:
    root = pathlib.Path(directory)
    refs = root / "refs"
    refs.mkdir()
    (refs / "reference.fa").write_text(">chr\n" + "A" * 32 + "N" + "C" * 32 + "\n")
    (root / "queries.fa").write_text(">reverse\n" + "T" * 32 + "\n"
                                     ">cross_n\n" + "A" * 31 + "C\n"
                                     ">ambiguous\n" + "A" * 15 + "N" + "A" * 16 + "\n")
    legacy_config = root / "legacy.toml"
    legacy_config.write_text(f"""[general]
index_method = "ibf"
ref_dir = "{refs}"
query_file = "{root / 'queries.fa'}"
output_dir = "{root / 'legacy_out'}"
fragment_size = 33
kmer_size = 32
hit_threshold = 1
threads = 1
single_results_writer = false
[ibf]
hash_functions = 3
fpr = 0.000001
multi_indexing = false
""")
    subprocess.run([str(LEGACY), str(legacy_config)], check=True, capture_output=True, text=True)
    with (root / "legacy_out" / "results_reference.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert int(rows[0]["reference_count"]) > 0  # canonical reverse complement
    assert rows[1]["reference_count"] == "0"  # no join across N
    assert rows[2]["reference_count"] == "0"  # ambiguous query retains its position
    with (root / "legacy_out" / "reference_kmers.tsv").open(newline="") as handle:
        hits = list(csv.DictReader(handle, delimiter="\t"))
    assert hits[1]["kmer_indices"] == "" and hits[2]["kmer_indices"] == ""

    genome = root / "genome.fa"
    genome.write_text(">chr\nAAAAACCCCCNNNNNGGGGGTTTTT\n")
    negative_config = root / "negative.toml"
    negative_config.write_text(f"""[general]
genome_file = "{genome}"
output_dir = "{root / 'negative_out'}"
kmer_size = 5
fragment_size = 6
chromosome_count = 1
sequence_count = 3
min_length = 6
max_length = 6
seed = 42
[ibf]
hash_functions = 3
fpr = 0.01
""")
    result = subprocess.run([str(NEGATIVE), str(negative_config)], check=True,
                            capture_output=True, text=True)
    assert "overlap=4" in (root / "negative_out" / "negative_set.log").read_text()
    with (root / "negative_out" / "negative_set.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3, result.stdout
    assert all(row["total_kmers"] == "2" and 0 <= float(row["match_fraction"]) <= 1
               for row in rows)

print("shared flat window checks passed")
