#!/usr/bin/env python3
"""Stream schema 2 per-reference TSVs into one query-major comparison and summary."""

import argparse
import contextlib
import csv
import gzip
import io
import json
import pathlib
import resource
import sys
import time


COMPARABLE = ("k", "strand_mode", "repeat_policy", "ambiguity_policy", "mode",
              "fragment_size", "hash_functions")
STATUSES = ("NO_VALID_KMERS", "NO_MATCH", "LOW_SUPPORT", "PASS_BASIC", "PASS_STRICT")


@contextlib.contextmanager
def open_text(path, mode):
    if path.suffix != ".gz":
        with path.open(mode, newline="") as handle:
            yield handle
    elif mode == "r":
        with gzip.open(path, "rt", newline="") as handle:
            yield handle
    else:
        with path.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, newline="") as handle:
                    yield handle


def peak_rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else rss * 1024)


def next_row(reader, path, expected_index):
    row = next(reader, None)
    if row is None:
        return None
    if row.get("schema_version") not in ("2", "3") or int(row["query_index"]) != expected_index:
        raise ValueError(f"Unexpected schema or query order in {path} at index {expected_index}")
    return row


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=pathlib.Path, action="append", required=True,
                        help="One schema 2 TSV per reference; repeat this option")
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--summary", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if len(args.input) < 2:
        parser.error("At least two reference TSVs are required")
    inputs = {path.resolve() for path in args.input}
    manifests = {pathlib.Path(str(path) + ".metrics.json").resolve() for path in args.input}
    if len(inputs) != len(args.input) or args.output.resolve() in inputs | manifests or \
       args.summary.resolve() in inputs | manifests or args.output.resolve() == args.summary.resolve():
        parser.error("Inputs and outputs must be distinct")
    with contextlib.ExitStack() as stack:
        bundles = []
        expected_header = None
        common_query_hash = None
        common_settings = None
        common_primary_source = None
        common_ibf_target = None
        common_scoring = None
        common_schema = None
        for path in args.input:
            manifest_path = pathlib.Path(str(path) + ".metrics.json")
            manifest = json.loads(manifest_path.read_text())
            primary_source = manifest.get("primary_match_source")
            if primary_source not in ("ibf", "exact", "ibf_candidates", "candidate_exact"):
                raise ValueError(f"Missing or invalid primary_match_source in {manifest_path}")
            if common_primary_source is None:
                common_primary_source = primary_source
            elif primary_source != common_primary_source:
                raise ValueError("Primary match sources differ across reference runs")
            scoring = manifest.get("scoring_method", "global_ibf" if primary_source == "ibf" else "global_exact")
            schema = manifest.get("schema_version", 2)
            if common_scoring is None:
                common_scoring, common_schema = scoring, schema
            elif scoring != common_scoring or schema != common_schema:
                raise ValueError("Scoring methods or schema versions differ across reference runs")
            query_hash = manifest["query_sequence_sha256"]
            if common_query_hash is None:
                common_query_hash = query_hash
            elif query_hash != common_query_hash:
                raise ValueError("Query FASTA hashes differ across reference runs")
            handle = stack.enter_context(open_text(path, "r"))
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "reference_id" not in reader.fieldnames:
                raise ValueError(f"Invalid schema 2 TSV header in {path}")
            if expected_header is None:
                expected_header = reader.fieldnames
            elif reader.fieldnames != expected_header:
                raise ValueError("Schema 2 TSV headers differ")
            first = next_row(reader, path, 0)
            if first is None:
                raise ValueError(f"No query rows in {path}")
            if first["schema_version"] != str(schema):
                raise ValueError("TSV and manifest schema versions differ")
            if schema == 3 and (first.get("scoring_method") != scoring or first.get("primary_match_source") != primary_source):
                raise ValueError("TSV and manifest scoring provenance differ")
            settings = tuple(first[field] for field in COMPARABLE)
            if primary_source == "ibf" and manifest.get("target_genome_wide_fpr") is not None:
                ibf_target = (manifest.get("target_genome_wide_fpr"),
                              manifest.get("requested_per_bin_fpr_cap"))
                if None in ibf_target:
                    raise ValueError(f"Missing IBF-wide FPR sizing provenance in {manifest_path}")
                if common_ibf_target is None:
                    common_ibf_target = ibf_target
                elif ibf_target != common_ibf_target:
                    raise ValueError("IBF-wide FPR targets or per-bin caps differ across reference runs")
            else:
                settings += (first["configured_per_bin_fpr"],)
            if schema == 3:
                settings += tuple(manifest.get(key) for key in ("per_bin_min_hits", "per_bin_min_fraction",
                    "per_bin_top_k", "per_bin_max_locus_gap", "per_bin_exact_verify", "max_locus_span", "mask_target"))
            if common_settings is None:
                common_settings = settings
            elif settings != common_settings:
                raise ValueError("Index/query settings differ across reference runs")
            bundles.append({"path": path, "reader": reader, "row": first, "reference_id": first["reference_id"],
                            "manifest": manifest,
                            "local_settings": tuple(first[field] for field in COMPARABLE) +
                                              (first["configured_per_bin_fpr"],)})
        bundles.sort(key=lambda bundle: bundle["reference_id"])
        if len({b["reference_id"] for b in bundles}) != len(bundles):
            raise ValueError("Duplicate reference_id across inputs")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        output = stack.enter_context(open_text(args.output, "w"))
        writer = csv.DictWriter(output, fieldnames=expected_header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        aggregates = {b["reference_id"]: {"query_count": 0, "valid_windows": 0,
                                           "matched_windows": 0,
                                           "false_positive_windows": 0, "false_negative_windows": 0,
                                           "statuses": {status: 0 for status in STATUSES}}
                      for b in bundles}
        query_count = 0
        while True:
            present = [b["row"] is not None for b in bundles]
            if not any(present):
                break
            if not all(present):
                raise ValueError("Reference runs have different numbers of query rows")
            query_id = bundles[0]["row"]["query_id"]
            for bundle in bundles:
                row = bundle["row"]
                if row["query_id"] != query_id or row["reference_id"] != bundle["reference_id"]:
                    raise ValueError(f"Query/reference identity mismatch at index {query_count}")
                manifest = bundle["manifest"]
                if row["schema_version"] != str(manifest.get("schema_version",2)):
                    raise ValueError("TSV and manifest schema versions differ")
                if row["schema_version"] == "3" and (row.get("scoring_method") != manifest["scoring_method"] or
                        row.get("primary_match_source") != manifest["primary_match_source"]):
                    raise ValueError("TSV and manifest scoring provenance differ")
                if (tuple(row[field] for field in COMPARABLE) +
                        (row["configured_per_bin_fpr"],)) != bundle["local_settings"]:
                    raise ValueError(f"Settings change within {bundle['path']}")
                writer.writerow(row)
                aggregate = aggregates[bundle["reference_id"]]
                aggregate["query_count"] += 1
                aggregate["valid_windows"] += int(row["valid_windows"])
                aggregate["matched_windows"] += int(row["matched_windows"])
                for key in ("false_positive_windows", "false_negative_windows"):
                    if row[key] != "NA":
                        aggregate[key] += int(row[key])
                if row["status"] not in aggregate["statuses"]:
                    raise ValueError(f"Unknown status in {bundle['path']}: {row['status']}")
                aggregate["statuses"][row["status"]] += 1
                bundle["row"] = next_row(bundle["reader"], bundle["path"], query_count + 1)
            query_count += 1
    for aggregate in aggregates.values():
        valid = aggregate["valid_windows"]
        aggregate["window_weighted_recovery"] = aggregate["matched_windows"] / valid if valid else None
    summary = {"schema_version": 1, "query_count": query_count, "reference_count": len(bundles),
               "query_sequence_sha256": common_query_hash,
               "primary_match_source": common_primary_source,
               "ibf_target_genome_wide_fpr": common_ibf_target[0] if common_ibf_target else None,
               "ibf_requested_per_bin_fpr_cap": common_ibf_target[1] if common_ibf_target else None,
               "inputs": [str(b["path"]) for b in bundles], "references": aggregates,
               "wall_seconds": time.monotonic() - started, "peak_resident_bytes": peak_rss_bytes()}
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Merged {query_count} queries x {len(bundles)} references; "
          f"wall={summary['wall_seconds']:.3f}s peak RSS={summary['peak_resident_bytes']}B")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(str(error)) from error
