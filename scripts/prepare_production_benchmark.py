#!/usr/bin/env python3
"""Inventory supplied inputs and extract a bounded reference; never edit inputs."""
import argparse
import gzip
import hashlib
import json
import pathlib
import platform
import re
import subprocess
import time


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def inventory(path, extract=None, limit=0, details=True):
    sha = digest(path)
    opener = gzip.open if path.suffix == ".gz" else open
    contigs = []
    selected = 0
    output = extract.open("x") if extract else None
    seen = set()
    chunks = []
    chunk_bases = 0
    def flush_stats():
        nonlocal chunk_bases
        if not chunks:
            return
        seq = "".join(chunks)
        if re.search(r"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", seq):
            raise ValueError("Invalid FASTA alphabet")
        upper = seq.upper()
        contigs[-1]["ambiguous"] += len(seq) - sum(upper.count(c) for c in "ACGT")
        contigs[-1]["lowercase"] += len(seq) - len(seq.translate(str.maketrans("", "", "acgtryswkmbdhvn")))
        chunks.clear(); chunk_bases = 0
    try:
        with opener(path, "rt") as handle:
            for raw in handle:
                if raw.startswith(">"):
                    flush_stats()
                    name = raw[1:].split()[0]
                    if name in seen:
                        raise ValueError(f"Duplicate ID {name}")
                    seen.add(name)
                    contigs.append({"id": name, "length": 0, "ambiguous": 0, "lowercase": 0})
                    if len(contigs) == 1 and output:
                        output.write(">" + name + "\n")
                elif raw.strip():
                    if not contigs:
                        raise ValueError("FASTA sequence before header")
                    seq = "".join(raw.split())
                    contigs[-1]["length"] += len(seq)
                    chunks.append(seq); chunk_bases += len(seq)
                    if chunk_bases >= 1 << 20:
                        flush_stats()
                    if len(contigs) == 1 and output:
                        part = seq if not limit else seq[:max(0, limit-selected)]
                        if part:
                            output.write(part+"\n"); selected += len(part)
        flush_stats()
    finally:
        if output:
            output.close()
    return {"path": str(path.resolve()), "file_bytes": path.stat().st_size, "raw_sha256": sha,
            "type": "gzip FASTA" if path.suffix == ".gz" else "plain FASTA",
            "record_count": len(contigs), "total_bases": sum(c["length"] for c in contigs),
            "contigs": contigs if details else contigs[:3], "contig_details_complete": details,
            "ambiguous_bases": sum(c["ambiguous"] for c in contigs),
            "lowercase_bases": sum(c["lowercase"] for c in contigs), "extracted_bases": selected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=pathlib.Path, default=pathlib.Path("rerun_test"))
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--max-reference-bases", type=int, default=0,
                        help="0 extracts entire first contig; otherwise take its original-coordinate prefix")
    parser.add_argument("--reuse-extracted", action="store_true", help="verify and reuse an existing first-contig extraction")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    destination = args.output / "reference_first_contig.fa"
    staging = args.output / "reference_first_contig.fa.partial"
    if destination.exists() and not args.reuse_extracted:
        raise ValueError(f"Refusing to overwrite benchmark reference {destination}")
    genomes = sorted(p for p in (args.source / "genomes").iterdir() if p.is_file())
    if not genomes:
        raise ValueError("No supplied genomes found")
    result = {"created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": platform.platform(),
              "git_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
              "query_release": "GENCODE v49", "selected_assembly": str(genomes[0]),
              "reference_scope": "first_contig" if args.max_reference_bases == 0 else "first_contig_prefix",
              "reference_base_cap": args.max_reference_bases,
              "genomes": [], "queries": inventory(args.source / "gencode.v49.lncRNA_transcripts.fa", details=False)}
    for i, path in enumerate(genomes):
        result["genomes"].append(inventory(path, staging if i == 0 and not args.reuse_extracted else None, args.max_reference_bases))
    if args.reuse_extracted:
        extracted = inventory(destination)
        original = result["genomes"][0]["contigs"][0]
        expected = min(original["length"], args.max_reference_bases) if args.max_reference_bases else original["length"]
        if extracted["record_count"] != 1 or extracted["contigs"][0]["id"] != original["id"] or extracted["total_bases"] != expected:
            raise ValueError("Existing extraction is incomplete or mismatched; use a new output directory")
        result["genomes"][0]["extracted_bases"] = extracted["total_bases"]
    else:
        original = result["genomes"][0]["contigs"][0]
        extracted = inventory(staging)
        expected = min(original["length"], args.max_reference_bases) if args.max_reference_bases else original["length"]
        if extracted["record_count"] != 1 or extracted["total_bases"] != expected:
            raise ValueError("Extraction failed length validation; partial output retained")
        staging.rename(destination)
    result["benchmark_reference"] = {"path": str(destination), "raw_sha256": digest(destination),
                                     "original_contig_id_preserved": True}
    (args.output / "input_inventory.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"query_records": result["queries"]["record_count"],
                      "assembly_contigs": result["genomes"][0]["record_count"],
                      "selected_contig": result["genomes"][0]["contigs"][0],
                      "extracted_bases": result["genomes"][0]["extracted_bases"]}, indent=2))


if __name__ == "__main__":
    main()
