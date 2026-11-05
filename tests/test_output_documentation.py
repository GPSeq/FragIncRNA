#!/usr/bin/env python3
"""Exercise all v2 report writers and check README coverage of emitted names.

Run explicitly: python3 tests/test_output_documentation.py build-production/lncrna_mers_v2
This is a documentation contract check, not a biological benchmark.
"""
import csv
import gzip
import json
import pathlib
import re
import subprocess
import sys
import tempfile

repo = pathlib.Path(__file__).resolve().parents[1]
program = pathlib.Path(sys.argv[1]).resolve()
readme = (repo / "README.md").read_text()
documented = set(re.findall(r"`([^`\n]+)`", readme))
seen_columns, seen_metrics = set(), set()

def check_name(name, context):
    assert name in documented, f"README lacks exact field {name!r} from {context}"

with tempfile.TemporaryDirectory(prefix="frag-output-docs-") as directory:
    root = pathlib.Path(directory)
    ref = root / "reference.fa"
    queries = root / "queries.fa"
    ref.write_text(">chrA\naACGTACGTNNACGTTACGATGCGTACCTAGGCT\n>chrB\nTTTTCCCCGGGGAAAATGCATGCA\n")
    queries.write_text(">q1\nACGTACGT\n>q2\nTTTANNNGC\n>q3\nNNNN\n>q4\nAC\n")
    query_bed, ref_bed = root / "query.bed", root / "reference.bed"
    query_bed.write_text("q1\t0\t2\tSINE\n")
    ref_bed.write_text("chrA\t2\t5\tLINE\n")
    expected, exons = root / "expected.bed", root / "exons.tsv"
    expected.write_text("chrA\t0\t9\tq1\tknown\nchrB\t0\t10\tq2\tnegative\n")
    exons.write_text("query_id\tquery_start\tquery_end\tcontig\tref_start\tref_end\tstrand\n"
                     "q1\t0\t8\tchrA\t1\t9\t+\n")
    output = root / "all.tsv"
    subprocess.run([str(program),"--reference",str(ref),"--queries",str(queries),
        "--output",str(output),"--mode","calibrate","--scoring-method","per_bin_verified",
        "--k","4","--fragment-size","8","--threads","4","--per-bin-min-hits","1",
        "--per-bin-min-fraction","0","--repeat","annotate","--query-repeat-bed",str(query_bed),
        "--reference-repeat-bed",str(ref_bed),"--positions",str(root/"positions.tsv.gz"),
        "--fragments",str(root/"fragments.tsv"),"--reference-hit-classes",str(root/"classes.tsv"),
        "--discover-loci",str(root/"loci.tsv"),"--expected-loci-bed",str(expected),
        "--exon-blocks",str(exons),"--local-window","8","--save-index",str(root/"index.ibf"),
        "--reference-cache",str(root/"reference.bases"),"--no-log-to-console"],check=True)
    tables = sorted(root.glob("*.tsv"))
    tables.remove(exons)  # supplied annotation, not generated output
    tables.append(root/"positions.tsv.gz")
    for path in tables:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path,"rt") as handle:
            reader = csv.DictReader(handle,delimiter="\t")
            assert reader.fieldnames, path
            for name in reader.fieldnames:
                check_name(name,path.name)
                seen_columns.add(name)
            for row in reader:
                assert None not in row, f"Extra fields in {path.name}"
                assert all(value is not None for value in row.values()), f"Missing field in {path.name}"
    for suffix in (".metrics.json",".ibf_tuning.json"):
        values = json.loads(pathlib.Path(str(output)+suffix).read_text())
        for name in values:
            check_name(name,suffix)
            seen_metrics.add(name)
    for name in json.loads((root/"index.ibf.metadata.json").read_text()):
        check_name(name,"index metadata")
    for event in map(json.loads,pathlib.Path(str(output)+".log.jsonl").read_text().splitlines()):
        for name in event:
            check_name(name,"JSONL")
    with output.open() as handle:
        rows = {row["query_id"]:row for row in csv.DictReader(handle,delimiter="\t")}
    assert rows["q3"]["status"] == rows["q4"]["status"] == "NO_VALID_KMERS"
    assert rows["q1"]["exon_eligible_windows"] == rows["q1"]["exon_matched_windows"] == "5"
    assert rows["q1"]["verification_scope"] == "candidate_intervals"
    assert rows["q1"]["exact_global_recovery"] == "1.000000"
print(f"README covers {len(seen_columns)} emitted TSV column names and {len(seen_metrics)} metrics/tuning keys; all optional v2 report writers exercised")
