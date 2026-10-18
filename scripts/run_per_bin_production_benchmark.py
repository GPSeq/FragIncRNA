#!/usr/bin/env python3
"""Bounded real-data benchmark; reuses indexes, records every command, no downloads."""
import argparse
import csv
import hashlib
import json
import pathlib
import random
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--program",type=pathlib.Path,required=True)
    parser.add_argument("--directory",type=pathlib.Path,required=True)
    parser.add_argument("--queries-directory",type=pathlib.Path)
    parser.add_argument("--assembly",type=pathlib.Path,help="Existing supplied genome; required only for the bounded whole-assembly stage")
    parser.add_argument("--stage",choices=("smoke","scale","optimized","truth","extras","sweep","whole"),required=True)
    parser.add_argument("--timeout-seconds",type=int,default=300)
    parser.add_argument("--total-budget-seconds",type=int,default=1800)
    parser.add_argument("--max-index-memory-gb",type=float,default=4)
    args = parser.parse_args()
    start = time.monotonic()
    root = args.directory.resolve()
    root.mkdir(parents=True,exist_ok=True)
    queries_directory = args.queries_directory.resolve() if args.queries_directory else root
    reference = root/"reference_first_contig.fa"
    index = root/"k20_nominal.ibf"
    base = [str(args.program.resolve()),"--reference",str(reference),"--k","20","--fragment-size","80000",
        "--strand","canonical","--hash-functions","3","--fpr","0.01","--batch-windows","1000000",
        "--max-ibf-bin-pairs","500000000","--ibf-max-index-memory-gb",str(args.max_index_memory_gb),
        "--reference-cache",str(root/"reference.bases"),"--no-log-to-console"]
    query_counts = {}
    def run(name,count,options):
        if time.monotonic()-start > args.total_budget_seconds:
            raise RuntimeError("Total benchmark time budget exhausted; completed records retained")
        command_path = root/(name+".command.json")
        if command_path.exists():
            prior = json.loads(command_path.read_text())
            if prior["exit_code"] == 0 and (root/(name+".tsv.metrics.json")).exists():
                print("Reusing completed "+name,flush=True); return
            raise RuntimeError("Prior failed run exists: choose a new run name/directory")
        timeout = min(args.timeout_seconds,max(1,int(args.total_budget_seconds-(time.monotonic()-start))))
        command = [*base,"--queries",str(queries_directory/f"queries_{count}.fa"),"--output",str(root/(name+".tsv")),*options]
        query_path = pathlib.Path(command[max(i for i,value in enumerate(command) if value == "--queries")+1])
        if query_path not in query_counts:
            with query_path.open() as handle:
                query_counts[query_path] = sum(line.startswith(">") for line in handle)
        actual_queries = query_counts[query_path]
        if not 0 < actual_queries <= count:
            raise RuntimeError(f"Refusing {actual_queries} query records: this stage has a query count cap of {count}")
        print("Running "+name,flush=True)
        sources = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(pathlib.Path("src").glob("*")) if path.is_file()}
        (root/(name+".config.json")).write_text(json.dumps({"command":command,"timeout_seconds":timeout,
            "actual_query_records":actual_queries,
            "executable_sha256":hashlib.sha256(args.program.read_bytes()).hexdigest(),"source_sha256":sources},indent=2)+"\n")
        process = subprocess.run([sys.executable,"scripts/measure_command.py","--output",str(command_path),
            "--","timeout","--signal=TERM","--kill-after=10",str(timeout),*command],check=False)
        if process.returncode:
            raise RuntimeError(f"{name} exited {process.returncode}; inspect its persistent logs")
    if args.stage == "whole":
        if args.assembly is None:
            parser.error("--assembly is required for --stage whole; it runs exactly 100 queries at k=20")
        assembly = args.assembly.resolve()
        checksum = hashlib.sha256()
        with assembly.open("rb") as handle:
            for block in iter(lambda:handle.read(1 << 20),b""):
                checksum.update(block)
        (root/"benchmark_scope.json").write_text(json.dumps({"reference_scope":"entire supplied assembly",
            "reference":str(assembly),"raw_sha256":checksum.hexdigest(),"query_count_cap":100,"k":20},indent=2)+"\n")
        run("whole_validate100",100,["--reference",str(assembly),"--mode","validate","--scoring-method","per_bin_verified",
            "--threads","4","--save-index",str(index)])
    elif args.stage == "smoke":
        run("cold_validate100",100,["--mode","validate","--scoring-method","per_bin_verified","--threads","1",
                                   "--save-index",str(index)])
        run("production100_t1",100,["--mode","ibf","--scoring-method","per_bin_verified","--threads","1","--load-index",str(index)])
        run("exact100",100,["--mode","exact","--scoring-method","global_exact","--threads","1"])
        run("approx100",100,["--mode","ibf","--scoring-method","per_bin_ibf","--threads","1","--load-index",str(index)])
    elif args.stage == "scale":
        for threads in (1,2,4,8,15):
            run(f"production5000_t{threads}",5000,["--mode","ibf","--scoring-method","per_bin_verified",
                "--threads",str(threads),"--load-index",str(index)])
        run("validate5000_t4",5000,["--mode","validate","--scoring-method","per_bin_verified","--threads","4","--load-index",str(index)])
        # Source-independent truth invocation, no IBF input.
        run("exact5000",5000,["--mode","exact","--scoring-method","global_exact","--threads","1"])
        hashes = {}
        for threads in (1,2,4,8,15):
            name = f"production5000_t{threads}.tsv"
            hashes[str(threads)] = {suffix: hashlib.sha256((root/(name+suffix)).read_bytes()).hexdigest()
                                    for suffix in ("",".per_bin.tsv",".per_bin.tsv.candidates.tsv")}
        assert len({json.dumps(value,sort_keys=True) for value in hashes.values()}) == 1
        (root/"serial_parallel_checksums.json").write_text(json.dumps(hashes,indent=2)+"\n")
    elif args.stage == "optimized":
        checksums = {}
        for threads in (1,2,4,8,15):
            name = f"optimized5000_t{threads}"
            run(name,5000,["--mode","ibf","--scoring-method","per_bin_verified","--threads",str(threads),"--load-index",str(index)])
            checksums[str(threads)] = {suffix:hashlib.sha256((root/(name+".tsv"+suffix)).read_bytes()).hexdigest()
                                      for suffix in ("",".per_bin.tsv",".per_bin.tsv.candidates.tsv")}
            for suffix in (".per_bin.tsv",".per_bin.tsv.candidates.tsv"):
                assert (root/(name+".tsv"+suffix)).read_bytes() == (root/("production5000_t1.tsv"+suffix)).read_bytes()
        assert len({json.dumps(value,sort_keys=True) for value in checksums.values()}) == 1
        with (root/"optimized5000_t1.tsv").open() as handle:
            current = list(csv.DictReader(handle,delimiter="\t"))
        with (root/"production5000_t1.tsv").open() as handle:
            previous = list(csv.DictReader(handle,delimiter="\t"))
        fields = ("query_id","valid_windows","matched_windows","recovery_fraction","status","per_bin_verified_local_recovery")
        assert len(current) == len(previous) == 5000
        assert all(all(a[field] == b[field] for field in fields) for a,b in zip(current,previous))
        (root/"optimized_serial_parallel_checksums.json").write_text(json.dumps(checksums,indent=2)+"\n")
    elif args.stage == "truth":
        run("validate5000_t4_retry600",5000,["--mode","validate","--scoring-method","per_bin_verified","--threads","4","--load-index",str(index)])
        run("exact5000",5000,["--mode","exact","--scoring-method","global_exact","--threads","1"])
        hashes = {}
        for threads in (1,2,4,8,15):
            name = f"production5000_t{threads}.tsv"
            hashes[str(threads)] = {suffix: hashlib.sha256((root/(name+suffix)).read_bytes()).hexdigest()
                                    for suffix in ("",".per_bin.tsv",".per_bin.tsv.candidates.tsv")}
        assert len({json.dumps(value,sort_keys=True) for value in hashes.values()}) == 1
        (root/"serial_parallel_checksums.json").write_text(json.dumps(hashes,indent=2)+"\n")
    elif args.stage == "extras":
        run("budgeted_validate100",100,["--mode","validate","--scoring-method","global_ibf","--threads","4",
            "--genome-fpr","0.01","--hash-functions","5","--save-index",str(root/"k20_budgeted.ibf")])
        run("calibrate100",100,["--mode","calibrate","--scoring-method","global_exact","--threads","4","--load-index",str(index)])
        for repeat in range(1,4):
            levels = ["off","info","debug"]
            random.Random(20261018+repeat).shuffle(levels)
            for level in levels:
                run(f"logging100_{level}_r{repeat}",100,["--mode","ibf","--scoring-method","per_bin_verified","--threads","4",
                    "--load-index",str(index),"--log-level",level])
        for k in (15,18,22,24):
            run(f"k{k}_validate100",100,["--mode","validate","--scoring-method","per_bin_verified","--threads","4",
                "--k",str(k),"--save-index",str(root/f"k{k}_nominal.ibf")])
    else:
        records = []
        for line in (queries_directory/"queries_100.fa").read_text().splitlines():
            if line.startswith(">"):
                records.append([line, ""])
            else:
                records[-1][1] += line.strip()
        def gene(header):
            fields = header[1:].split()[0].split("|")
            return fields[1] if len(fields)>1 else fields[0]
        genes = sorted({gene(header) for header,_ in records})
        random.Random(20261018).shuffle(genes)
        development = set(genes[:len(genes)//2])
        split = {"seed":20261018,"grouping":"GENCODE gene ID (not independently annotated gene family)",
                 "biological_family_independence_established":False}
        for label in ("development","heldout"):
            selected = [(h,s) for h,s in records if (gene(h) in development) == (label == "development")]
            path = root/(label+".fa")
            path.write_text("".join(h+"\n"+s+"\n" for h,s in selected))
            split[label] = {"query_ids":[h[1:].split()[0] for h,_ in selected],"gene_ids":sorted({gene(h) for h,_ in selected})}
        assert not set(split["development"]["gene_ids"]) & set(split["heldout"]["gene_ids"])
        (root/"threshold_split.json").write_text(json.dumps(split,indent=2)+"\n")
        grid = [("baseline",2,.05,64,50000),("hits1",1,.05,64,50000),
                ("fraction01",2,.01,64,50000),("top4",2,.05,4,50000),
                ("gap0",2,.05,64,0),("fraction10",2,.10,64,50000)]
        scores = []
        def settings(values,subset):
            label,hits,fraction,top,gap = values
            return ["--mode","validate","--scoring-method","per_bin_verified","--threads","4",
                "--load-index",str(index),"--queries",str(root/(subset+".fa")),
                "--per-bin-min-hits",str(hits),"--per-bin-min-fraction",str(fraction),
                "--per-bin-top-k",str(top),"--per-bin-max-locus-gap",str(gap)]
        for config in grid:
            name = "development_"+config[0]
            run(name,100,settings(config,"development"))
            with (root/(name+".tsv.per_bin.tsv")).open() as handle:
                values = list(csv.DictReader(handle,delimiter="\t"))
            exact = sum(int(row["exact_global_matched_windows"]) for row in values)
            local = sum(int(row["candidate_union_exact_matched_windows"]) for row in values)
            candidates = sum(int(row["candidate_locus_count"]) for row in values)
            scores.append({"configuration":config,"window_recall":local/exact if exact else None,"candidates":candidates})
        chosen = max(scores,key=lambda value: (value["window_recall"] or 0,-value["candidates"]))
        (root/"threshold_selection.json").write_text(json.dumps({"objective":"maximize technical global-exact window recall on DEVELOPMENT only, tie-break fewer candidate loci",
            "development_results":scores,"selected":chosen},indent=2)+"\n")
        run("heldout_selected",100,settings(chosen["configuration"],"heldout"))
        run("heldout_baseline",100,settings(grid[0],"heldout"))


if __name__ == "__main__":
    main()
