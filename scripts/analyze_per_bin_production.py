#!/usr/bin/env python3
"""Compare actual IBF predictions with independent whole-reference k-mer truth."""
import argparse
import csv
import hashlib
import json
import pathlib
import statistics


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ratio(n, d):
    return n/d if d else None


def confusion(predictions, truth):
    tp = sum(p and t for p,t in zip(predictions,truth))
    fp = sum(p and not t for p,t in zip(predictions,truth))
    fn = sum(not p and t for p,t in zip(predictions,truth))
    tn = sum(not p and not t for p,t in zip(predictions,truth))
    return dict(tp=tp,fp=fp,fn=fn,tn=tn,sensitivity=ratio(tp,tp+fn),specificity=ratio(tn,tn+fp),
                precision=ratio(tp,tp+fp),f1=ratio(2*tp,2*tp+fp+fn))


def analyze(path):
    records = table(path)
    paired_path = pathlib.Path(str(path)+".per_bin.tsv")
    if not paired_path.exists() or not records or records[0].get("exact_global_recovery", "NA") == "NA":
        return None
    paired = {r["query_id"]:r for r in table(paired_path)}
    usable = [r for r in records if int(r["valid_windows"]) > 0]
    truth = [int(paired[r["query_id"]]["exact_global_matched_windows"]) for r in usable]
    valid = [int(r["valid_windows"]) for r in usable]
    methods = {}
    for method, count_column, status_column in (
        ("global_ibf","ibf_union_matched_windows","global_ibf_status"),
        ("per_bin_ibf","per_bin_candidate_window_support","per_bin_ibf_status"),
        ("per_bin_verified","candidate_union_exact_matched_windows","per_bin_verified_status"),
        ("global_exact","exact_global_matched_windows","exact_global_status")):
        counts = [int(r[count_column]) if count_column in r else int(paired[r["query_id"]][count_column]) for r in usable]
        errors = [abs(p/n-t/n) for p,t,n in zip(counts,truth,valid)]
        statuses = [r[status_column] for r in usable]
        truth_statuses = [r["exact_global_status"] for r in usable]
        classes = {}
        for label,threshold_key in (("basic","basic_threshold"),("strict","strict_threshold")):
            classes[label] = confusion([p/n >= float(r[threshold_key]) for p,n,r in zip(counts,valid,usable)],
                                      [t/n >= float(r[threshold_key]) for t,n,r in zip(truth,valid,usable)])
        matrix = {}
        for t,p in zip(truth_statuses,statuses):
            matrix[t+" -> "+p] = matrix.get(t+" -> "+p,0)+1
        strata = {}
        for name,lo,hi in (("lt500",0,500),("500to2000",500,2000),("ge2000",2000,10**12)):
            chosen = [i for i,n in enumerate(valid) if lo <= n < hi]
            strata[name] = {"queries":len(chosen),"mean_absolute_recovery_error":
                            statistics.mean(errors[i] for i in chosen) if chosen else None}
        methods[method] = dict(query_count=len(usable),matched_windows=sum(counts),
            weighted_recovery=ratio(sum(counts),sum(valid)),
            any_match_confusion=confusion([p>0 for p in counts],[t>0 for t in truth]),
            threshold_confusions=classes,class_confusion=matrix,
            classification_agreement=ratio(sum(p==t for p,t in zip(statuses,truth_statuses)),len(usable)),
            mean_absolute_recovery_error=statistics.mean(errors) if errors else None,
            median_absolute_recovery_error=statistics.median(errors) if errors else None,
            p95_absolute_recovery_error=sorted(errors)[max(0,int(.95*len(errors))-1)] if errors else None,
            length_strata_by_valid_windows=strata)
        columns = {"global_ibf":("false_positive_windows","false_negative_windows"),
                   "per_bin_ibf":("per_bin_ibf_false_positive_windows","per_bin_ibf_false_negative_windows"),
                   "per_bin_verified":("candidate_false_positive_windows","candidate_false_negative_windows")}
        if method in columns:
            fp_col,fn_col = columns[method]
            available = all(r.get(fp_col,"NA") != "NA" and r.get(fn_col,"NA") != "NA" for r in usable)
            fp = sum(int(r[fp_col]) for r in usable) if available else None
            fn = sum(int(r[fn_col]) for r in usable) if available else None
            methods[method]["false_positive_windows"] = fp
            methods[method]["false_negative_windows"] = fn
            methods[method]["window_fpr"] = ratio(fp,sum(valid)-sum(truth)) if fp is not None else None
            methods[method]["window_recall"] = ratio(sum(truth)-fn,sum(truth)) if fn is not None else None
    exact = sum(truth)
    verified = methods["per_bin_verified"]["matched_windows"]
    scope = "complete first contig of supplied assembly; see input_inventory.json" if (path.parent/"input_inventory.json").exists() else "prefix extraction attempt; see INPUT_SCOPE.md; excluded from final benchmark"
    if (path.parent/"benchmark_scope.json").exists():
        scope = json.loads((path.parent/"benchmark_scope.json").read_text())["reference_scope"]
    return dict(path=str(path),reference_scope=scope,
        k=int(records[0]["k"]),
        biological_truth=False,query_count=len(records),zero_valid_queries=len(records)-len(usable),
        valid_windows=sum(valid),exact_positive_windows=exact,candidate_window_recall=ratio(verified,exact),
        candidate_false_negative_windows=sum(int(r["candidate_false_negative_windows"]) for r in usable),
        candidate_false_positive_windows=sum(int(r["candidate_false_positive_windows"]) for r in usable),
        global_ibf_false_positive_windows=sum(int(r["false_positive_windows"]) for r in usable),
        top_k_dropped_bins=sum(int(r["top_k_dropped_bins"]) for r in paired.values()),methods=methods)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=pathlib.Path)
    args = parser.parse_args()
    reports = {}
    timings = []
    for path in sorted(args.directory.glob("*.tsv")):
        if not pathlib.Path(str(path)+".metrics.json").exists():
            continue
        metrics = json.loads(pathlib.Path(str(path)+".metrics.json").read_text())
        command_path = path.with_suffix(".command.json")
        if not command_path.exists():
            continue
        command = json.loads(command_path.read_text())
        if command["exit_code"] != 0:
            continue
        result = analyze(path)
        if result:
            reports[path.name] = result
        timings.append(dict(run=path.stem,queries=metrics["query_count"],k=metrics["k"],method=metrics["scoring_method"],
            threads=metrics["threads"],wall_seconds=command["wall_seconds"],peak_rss_bytes=metrics["peak_resident_bytes"],
            index_build_seconds=metrics["index_construction_seconds"],index_load_seconds=metrics["index_load_seconds"],
            index_bytes=metrics["disk_index_bytes"],query_seconds=metrics["query_seconds"],
            candidate_selection_seconds=metrics["candidate_selection_seconds"],
            verification_seconds=metrics["candidate_verification_seconds"],exact_seconds=metrics["exact_scan_seconds"],
            candidate_bases=metrics["candidate_reference_bases_fetched"],global_exact_scans=metrics["whole_reference_exact_scans"],
            requested_per_bin_fpr=metrics["requested_per_bin_fpr_cap"],effective_per_bin_fpr=metrics["configured_per_bin_fpr"]))
    (args.directory/"method_comparison.json").write_text(json.dumps(reports,indent=2)+"\n")
    comparison_rows = []
    for run, report in reports.items():
        for method, values in report["methods"].items():
            row = dict(run=run,k=report["k"],queries=report["query_count"],method=method,
                valid_windows=report["valid_windows"],exact_positive_windows=report["exact_positive_windows"],
                matched_windows=values["matched_windows"],
                false_positive_windows=values.get("false_positive_windows"),
                false_negative_windows=values.get("false_negative_windows"),
                window_recall=values.get("window_recall"),window_fpr=values.get("window_fpr"),
                classification_agreement=values["classification_agreement"],
                mean_absolute_recovery_error=values["mean_absolute_recovery_error"],
                median_absolute_recovery_error=values["median_absolute_recovery_error"],
                p95_absolute_recovery_error=values["p95_absolute_recovery_error"])
            for endpoint, metrics in [("any_match",values["any_match_confusion"]),
                    *values["threshold_confusions"].items()]:
                row.update({endpoint+"_"+key:value for key,value in metrics.items()})
            comparison_rows.append({key:"NA" if value is None else value for key,value in row.items()})
    if comparison_rows:
        with (args.directory/"method_comparison.tsv").open("w") as handle:
            writer = csv.DictWriter(handle,fieldnames=list(comparison_rows[0]),delimiter="\t",lineterminator="\n")
            writer.writeheader(); writer.writerows(comparison_rows)
    if timings:
        with (args.directory/"timings.tsv").open("w") as handle:
            writer = csv.DictWriter(handle,fieldnames=list(timings[0]),delimiter="\t",lineterminator="\n")
            writer.writeheader(); writer.writerows(timings)
    checksums = {str(path.relative_to(args.directory)):hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in sorted(args.directory.iterdir()) if path.is_file() and path.suffix in (".tsv",".json",".toml",".log",".jsonl",".fa")
                 and not path.name.startswith("reference") and path.name != "checksums.json"}
    (args.directory/"checksums.json").write_text(json.dumps(checksums,indent=2)+"\n")
    print(f"Compared {len(reports)} validation tables; recorded {len(timings)} measured runs")


if __name__ == "__main__":
    main()
