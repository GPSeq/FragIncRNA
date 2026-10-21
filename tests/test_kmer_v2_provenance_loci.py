"""Index round trips, SHA-256 provenance, and exact locus verification."""

import csv
import gzip
import hashlib
import json
import pathlib
import subprocess
import sys
import tempfile


PROGRAM = pathlib.Path(sys.argv[1])


def run(tmp, label, reference, queries, *extra):
    output = tmp / f"{label}.tsv"
    command = [str(PROGRAM), "--reference", str(reference), "--queries", str(queries),
               "--output", str(output), "--k", "4", "--fragment-size", "4", *extra]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise AssertionError((command, result.stderr))
    with output.open(newline="") as handle:
        rows = {row["query_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    return rows, json.loads(pathlib.Path(str(output) + ".metrics.json").read_text()), output


def rejected(tmp, label, reference, queries, *extra, contains):
    result = subprocess.run([str(PROGRAM), "--reference", str(reference), "--queries", str(queries),
                             "--output", str(tmp / f"{label}.tsv"), "--k", "4",
                             "--fragment-size", "4", *extra], capture_output=True, text=True)
    assert result.returncode != 0 and contains in result.stderr, result.stderr


def semantic_bytes(records):
    payload = bytearray(b"FragIncRNA-FASTA-semantic-v1")
    for identifier, sequence in records:
        encoded_id = identifier.encode()
        payload += b"\xff" + len(encoded_id).to_bytes(8, "little") + encoded_id
        payload += sequence.encode() + b"\xfe" + len(sequence).to_bytes(8, "little")
    return bytes(payload)


with tempfile.TemporaryDirectory(prefix="fraglncrna-provenance-test-") as directory:
    tmp = pathlib.Path(directory)
    reference = tmp / "reference.fa"
    queries = tmp / "queries.fa"
    reference.write_text(">chr1\nAAAACCCCNNNNGGGGTTTT\n>chr2\nAAAACCCC\n")
    queries.write_text(">concentrated\nAAAACCCC\n>dispersed\nAAAAGGGG\n>unannotated\nTTTT\n")
    index = tmp / "reference.ibf"
    built, built_metrics, built_output = run(tmp, "built", reference, queries,
                                              "--strand", "forward", "--save-index", str(index))
    loaded, loaded_metrics, loaded_output = run(tmp, "loaded", reference, queries,
                                                 "--strand", "forward", "--load-index", str(index))
    assert built_output.read_bytes() == loaded_output.read_bytes()
    assert built == loaded
    assert built_metrics["index_provenance"] == "built"
    assert loaded_metrics["index_provenance"] == "loaded"
    assert loaded_metrics["index_construction_seconds"] == 0
    assert loaded_metrics["index_load_seconds"] >= 0
    assert built_metrics["disk_index_bytes"] == index.stat().st_size > 0
    assert loaded_metrics["disk_index_bytes"] == index.stat().st_size
    assert built_metrics["index_file_sha256"] == hashlib.sha256(index.read_bytes()).hexdigest()
    metadata = json.loads(pathlib.Path(str(index) + ".metadata.json").read_text())
    assert metadata["index_format_version"] == 3
    assert metadata["fragment_overlap"] == 3
    assert metadata["contigs"] == [{"id": "chr1", "length": 20}, {"id": "chr2", "length": 8}]
    reference_semantic = semantic_bytes([("chr1", "AAAACCCCNNNNGGGGTTTT"), ("chr2", "AAAACCCC")])
    query_semantic = semantic_bytes([("concentrated", "AAAACCCC"), ("dispersed", "AAAAGGGG"),
                                     ("unannotated", "TTTT")])
    assert built_metrics["reference_sequence_sha256"] == hashlib.sha256(reference_semantic).hexdigest()
    assert built_metrics["query_sequence_sha256"] == hashlib.sha256(query_semantic).hexdigest()
    assert metadata["reference_sequence_sha256"] == built_metrics["reference_sequence_sha256"]
    with gzip.open(tmp / "reference.fa.gz", "wt") as handle:
        handle.write(reference.read_text())
    _, compressed_metrics, _ = run(tmp, "compressed_load", tmp / "reference.fa.gz", queries,
                                    "--strand", "forward", "--load-index", str(index))
    assert compressed_metrics["reference_sequence_sha256"] == built_metrics["reference_sequence_sha256"]
    rejected(tmp, "wrong_strand", reference, queries, "--load-index", str(index), contains="strand")
    rejected(tmp, "wrong_k", reference, queries, "--k", "3", "--strand", "forward",
             "--load-index", str(index), contains="metadata: k")
    rejected(tmp, "wrong_fpr", reference, queries, "--strand", "forward", "--fpr", "0.001",
             "--load-index", str(index), contains="configured_fpr")
    changed_reference = tmp / "changed.fa"
    changed_reference.write_text(reference.read_text().replace("GGGG", "GAGG"))
    rejected(tmp, "changed_reference", changed_reference, queries, "--strand", "forward",
             "--load-index", str(index), contains="reference_sha256")
    corrupt = tmp / "corrupt.ibf"
    corrupt.write_bytes(b"?" + index.read_bytes()[1:])
    rejected(tmp, "corrupt", reference, queries, "--strand", "forward", "--load-index", str(corrupt),
             contains="bad magic")
    damaged_payload = tmp / "damaged_payload.ibf"
    payload = bytearray(index.read_bytes())
    payload[-1] ^= 1
    damaged_payload.write_bytes(payload)
    rejected(tmp, "damaged_payload", reference, queries, "--strand", "forward",
             "--load-index", str(damaged_payload), contains="payload checksum mismatch")
    repeat_bed = tmp / "reference_repeats.bed"
    repeat_bed.write_text("chr1\t0\t1\tLINE\n")
    masked_index = tmp / "masked.ibf"
    run(tmp, "masked_build", reference, queries, "--strand", "forward", "--repeat", "exclude",
        "--mask-target", "reference", "--reference-repeat-bed", str(repeat_bed),
        "--save-index", str(masked_index))
    run(tmp, "masked_load", reference, queries, "--strand", "forward", "--repeat", "exclude",
        "--mask-target", "reference", "--reference-repeat-bed", str(repeat_bed),
        "--load-index", str(masked_index))
    repeat_bed.write_text("chr1\t1\t2\tLINE\n")
    rejected(tmp, "changed_bed", reference, queries, "--strand", "forward", "--repeat", "exclude",
             "--mask-target", "reference", "--reference-repeat-bed", str(repeat_bed),
             "--load-index", str(masked_index), contains="reference_repeat_bed_sha256")

    expected_bed = tmp / "expected.bed"
    expected_bed.write_text("chr1\t0\t8\tconcentrated\tmain\n"
                            "chr2\t0\t8\tconcentrated\talternate\n"
                            "chr1\t12\t16\tdispersed\ttarget\n")
    loci_path = tmp / "loci_candidates.tsv"
    rows, metrics, output = run(tmp, "loci", reference, queries, "--strand", "forward",
                                "--load-index", str(index), "--discover-loci", str(loci_path),
                                "--locus-gap", "4", "--expected-loci-bed", str(expected_bed))
    assert rows["concentrated"]["number_of_discovered_loci"] == "2"
    assert rows["concentrated"]["best_locus_matched_windows"] == "5"
    assert rows["concentrated"]["best_locus_recovery"] == "1.000000"
    assert rows["dispersed"]["number_of_discovered_loci"] == "3"
    assert rows["dispersed"]["expected_locus_matched_windows"] == "1"
    assert rows["dispersed"]["off_expected_locus_windows"] == "1"
    assert rows["unannotated"]["expected_locus_recovery"] == "NA"
    assert metrics["expected_locus_matched_windows"] == 6
    assert metrics["expected_loci_bed_sha256"] == hashlib.sha256(expected_bed.read_bytes()).hexdigest()
    with loci_path.open(newline="") as handle:
        loci = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["contig"], row["start"], row["end"], row["matched_windows"])
            for row in loci if row["query_id"] == "concentrated"] == [
                ("chr1", "0", "8", "5"), ("chr2", "0", "8", "5")]
    with pathlib.Path(str(output) + ".expected_loci.tsv").open(newline="") as handle:
        expected_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["label"], row["matched_windows"]) for row in expected_rows] == [
        ("main", "5"), ("alternate", "5"), ("target", "1")]
    second_loci = tmp / "loci_again_candidates.tsv"
    second_rows, _, second_output = run(tmp, "loci_again", reference, queries, "--strand", "forward",
                                         "--load-index", str(index), "--discover-loci", str(second_loci),
                                         "--locus-gap", "4", "--expected-loci-bed", str(expected_bed))
    assert rows == second_rows and output.read_bytes() == second_output.read_bytes()
    assert loci_path.read_bytes() == second_loci.read_bytes()
    rejected(tmp, "locus_cap", reference, queries, "--strand", "forward", "--discover-loci",
             str(tmp / "capped_candidates.tsv"), "--locus-gap", "4", "--max-locus-rows", "1",
             contains="exceed --max-locus-rows")
    boundary_bed = tmp / "boundary.bed"
    boundary_bed.write_text("chr1\t12\t15\tdispersed\tshort\n")
    edge_rows, _, _ = run(tmp, "edge", reference, queries, "--strand", "forward",
                          "--expected-loci-bed", str(boundary_bed))
    assert edge_rows["dispersed"]["expected_locus_matched_windows"] == "0"
    assert edge_rows["dispersed"]["off_expected_locus_windows"] == "2"
    overlapping = tmp / "overlap.bed"
    overlapping.write_text("chr1\t0\t8\tconcentrated\nchr1\t7\t9\tconcentrated\n")
    rejected(tmp, "overlap", reference, queries, "--expected-loci-bed", str(overlapping),
             contains="Overlapping expected loci")

    repeated_reference = tmp / "repeated.fa"
    repeated_reference.write_text(">repeats\nAAAAAA\n>empty\nNNNNNN\n")
    repeated_query = tmp / "repeated_query.fa"
    repeated_query.write_text(">q\nAAAA\n")
    repeated_index = tmp / "repeated.ibf"
    _, repeated_metrics, repeated_output = run(tmp, "repeated_build", repeated_reference, repeated_query,
                                                "--mode", "calibrate", "--save-index", str(repeated_index))
    _, loaded_repeat_metrics, loaded_repeat_output = run(tmp, "repeated_load", repeated_reference,
                                                          repeated_query, "--mode", "calibrate",
                                                          "--load-index", str(repeated_index))
    with pathlib.Path(str(repeated_output) + ".bins.tsv").open(newline="") as handle:
        bin_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert pathlib.Path(str(repeated_output) + ".bins.tsv").read_bytes() == pathlib.Path(
        str(loaded_repeat_output) + ".bins.tsv").read_bytes()
    assert len(bin_rows) == 4
    assert bin_rows[0]["inserted_reference_windows"] == "3"
    assert bin_rows[0]["distinct_reference_kmers"] == "1"
    assert 0 < int(bin_rows[0]["occupied_bits"]) <= 3
    assert all(row["distinct_reference_kmers"] == "0" and row["occupied_bits"] == "0"
               for row in bin_rows[1:])
    assert repeated_metrics["distinct_bin_kmers_sum"] == 1
    assert repeated_metrics["occupied_bin_bits_sum"] == loaded_repeat_metrics["occupied_bin_bits_sum"]
    assert repeated_metrics["bit_occupancy_measurement_seconds"] >= 0
    assert repeated_metrics["peak_rss_after_bit_occupancy_bytes"] > 0
    assert json.loads(pathlib.Path(str(repeated_index) + ".metadata.json").read_text())["distinct_bin_kmers_sum"] == 1

print("index provenance and locus checks passed")
