#!/usr/bin/env python3
"""Run whole-reference exact k-mer matching for genome FASTAs x k=15..32.

Python >= 3.11, standard library only. Uses lncrna_mers_v2 mode=exact,
scoring_method=global_exact. Each successful genome/k job is archived with
checksums and provenance. Rerun the same command to resume safely.
Exact matching tests k-mer membership anywhere in the reference; it is not
full-transcript alignment. Run from the FragIncRNA repository root.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime as dt
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import time
import threading
import tomllib

FASTA_SUFFIXES = (
    ".fa", ".fasta", ".fna", ".fas",
    ".fa.gz", ".fasta.gz", ".fna.gz", ".fas.gz",
)
# Optional reports are redirected into each job to avoid shared output paths.
REPORT_OUTPUTS = (
    "positions", "discover_loci", "expected_loci_output", "repeat_classes",
    "reference_hit_classes", "window_output", "kmer_fasta_output",
)


def job_overrides(source: str, genome: Path, queries: Path, results: Path,
                  k: int, threads: int) -> dict[str, object]:
    config = tomllib.loads(source)["v2"]
    overrides = {
        "reference": str(genome), "queries": str(queries),
        "output": str(results / "lncrna.tsv"), "metrics": "",
        "k": k, "threads": threads, "reference_cache": "", "log_file": "",
        "save_index": "", "load_index": "",
    }
    for key in REPORT_OUTPUTS:
        if config.get(key):
            suffix = ".fa" if key == "kmer_fasta_output" else ".tsv"
            if str(config[key]).endswith(".gz"):
                suffix += ".gz"
            overrides[key] = str(results / (key + suffix))
    return overrides


def iso_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def slug_for(genome: Path) -> str:
    name = genome.name
    for suffix in sorted(FASTA_SUFFIXES, key=len, reverse=True):
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
            break
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._-")[:72] or "genome"
    return f"{slug}-{hashlib.sha256(genome.name.encode()).hexdigest()[:10]}"


def find_genomes(directory: Path) -> list[Path]:
    items = []
    for item in sorted(directory.iterdir(), key=lambda p: p.name):
        if item.is_symlink():
            print(f"[SKIP] symlink: {item}", file=sys.stderr)
        elif item.is_file() and item.name.lower().endswith(FASTA_SUFFIXES):
            items.append(item.resolve())
    return items


def toml_literal(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    if isinstance(value, str):
        # JSON quoted strings are valid TOML basic strings for these paths.
        return json.dumps(value, ensure_ascii=False)
    raise TypeError(f"Unsupported TOML override: {type(value).__name__}")


def render_config(template: str, overrides: dict[str, object]) -> str:
    """Rewrite only [v2], preserving all other template settings/comments."""
    if set(tomllib.loads(template)) != {"v2"}:
        raise ValueError("Template must contain exactly one [v2] section")
    lines = (template.rstrip() + "\n").splitlines(keepends=True)
    current = ""
    found: set[str] = set()
    result = []
    insert_at = None
    for line in lines:
        section_match = re.match(r"^\s*\[([^]]+)\]\s*(?:#.*)?$", line)
        if section_match:
            if current == "v2" and insert_at is None:
                insert_at = len(result)
            current = section_match.group(1).strip()
        match = re.match(r"^(\s*)([A-Za-z_][\w-]*)\s*=", line)
        if current == "v2" and match and match.group(2) in overrides:
            key = match.group(2)
            if key in found:
                raise ValueError(f"Duplicate {key} key in TOML template")
            found.add(key)
            result.append(f"{match.group(1)}{key} = {toml_literal(overrides[key])}\n")
        else:
            result.append(line)
    if insert_at is None:
        insert_at = len(result)
    missing = [f"{key} = {toml_literal(val)}\n" for key, val in overrides.items() if key not in found]
    result[insert_at:insert_at] = missing
    rendered = "".join(result)
    config = tomllib.loads(rendered)["v2"]
    for key, val in overrides.items():
        if config.get(key) != val:
            raise ValueError(f"TOML round-trip failed: {key}")
    if config.get("mode") != "exact" or config.get("scoring_method") != "global_exact":
        raise ValueError('Template must set mode="exact" and scoring_method="global_exact"')
    if config.get("fdr_control_per_bin", False) or config.get("per_bin_output"):
        raise ValueError("Exact mode cannot use per-bin IBF candidate selection")
    if "genome_fpr" in config or config.get("calibration_bins"):
        raise ValueError("Exact mode cannot use genome_fpr or calibration_bins")
    if config.get("reference_dir") or config.get("output_dir"):
        raise ValueError("Template must use single-reference mode; the runner manages the genome grid")
    if config.get("ibf_build_only", False) or config.get("ibf_dry_run", False):
        raise ValueError("Exact jobs cannot use IBF build-only or dry-run options")
    return rendered


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + f".part-{os.getpid()}")
    try:
        part.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        part.replace(path)
    finally:
        part.unlink(missing_ok=True)


def list_archivable_files(work_dir: Path) -> list[Path]:
    roots = [work_dir / "generated_config.toml", work_dir / "job.json"]
    results = work_dir / "results"
    for entry in results.rglob("*"):
        if entry.is_symlink():
            raise RuntimeError(f"Refusing to archive symlink: {entry}")
        if entry.is_file():
            roots.append(entry)
    for item in roots:
        if not item.is_file() or item.is_symlink():
            raise RuntimeError(f"Missing or unsafe archive member: {item}")
    return sorted(roots, key=lambda p: str(p.relative_to(work_dir)))


def archive_and_verify(work_dir: Path, archive: Path, job: dict, elapsed: float) -> dict:
    """Archive complete outputs; re-read every member for integrity."""
    files = list_archivable_files(work_dir)
    manifest = {
        "format": "fraglncrna-grid-archive-v1",
        "job_signature": job["job_signature"],
        "finished_utc": iso_now(),
        "elapsed_seconds": round(elapsed, 3),
        "files": {
            str(p.relative_to(work_dir)): {"size": p.stat().st_size, "sha256": sha256_file(p)}
            for p in files
        },
    }
    write_json(work_dir / "archive_manifest.json", manifest)
    files.append(work_dir / "archive_manifest.json")
    archive.parent.mkdir(parents=True, exist_ok=True)
    part = archive.with_name(archive.name + f".part-{os.getpid()}")
    try:
        with tarfile.open(part, mode="w:gz", compresslevel=6) as tf:
            for p in files:
                tf.add(p, arcname=str(p.relative_to(work_dir)), recursive=False)
        with tarfile.open(part, mode="r:gz") as tf:
            members = {member.name: member for member in tf.getmembers()}
            expected = set(manifest["files"]) | {"archive_manifest.json"}
            if set(members) != expected:
                raise RuntimeError("Archive member list mismatch")
            for name, meta in manifest["files"].items():
                member = members[name]
                if not member.isfile() or member.size != meta["size"]:
                    raise RuntimeError(f"Archive member size/type mismatch: {name}")
                h = hashlib.sha256()
                handle = tf.extractfile(member)
                if handle is None:
                    raise RuntimeError(f"Cannot extract {name}")
                with handle:
                    for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                        h.update(chunk)
                if h.hexdigest() != meta["sha256"]:
                    raise RuntimeError(f"Archive checksum mismatch: {name}")
            archived_manifest = tf.extractfile(members["archive_manifest.json"])
            if archived_manifest is None or json.load(archived_manifest) != manifest:
                raise RuntimeError("Archive provenance manifest mismatch")
        part.replace(archive)
    finally:
        part.unlink(missing_ok=True)
    receipt = {
        "status": "complete",
        "archive": archive.name,
        "archive_sha256": sha256_file(archive),
        "archive_bytes": archive.stat().st_size,
        "job_signature": job["job_signature"],
        "completed_utc": iso_now(),
        "elapsed_seconds": round(elapsed, 3),
    }
    return receipt


def completion_receipt_valid(archive: Path, receipt_file: Path, job_signature: str) -> bool:
    if not archive.is_file() and not receipt_file.exists():
        return False
    if not archive.is_file() or not receipt_file.is_file():
        raise RuntimeError(f"Incomplete archive/receipt pair; investigate before rerun: {archive}")
    payload = json.loads(receipt_file.read_text(encoding="utf-8"))
    if payload.get("job_signature") != job_signature:
        raise RuntimeError(f"Completed archive exists for different inputs/settings: {archive}")
    if (payload.get("status") != "complete" or payload.get("archive") != archive.name
            or payload.get("archive_sha256") != sha256_file(archive)):
        raise RuntimeError(f"Existing archive failed integrity check: {archive}")
    # The archive's uncompressed members were already checked before receipt publication.
    return True


def run_executable(binary: Path, config: Path, cwd: Path, log_file: Path,
                   stop: threading.Event) -> int:
    # Write directly to the private executable log: concurrent output must not mix.
    with log_file.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            [str(binary), "--config", str(config)], cwd=cwd,
            stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
        )
        while True:
            try:
                return process.wait(timeout=0.25)
            except subprocess.TimeoutExpired:
                if stop.is_set():
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        process.wait()
                    raise RuntimeError("Run cancelled; work/logs preserved")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--genomes-dir", type=Path, required=True, help="Flat directory containing .fa/.fna/.fasta[.gz]")
    parser.add_argument("--queries", type=Path, required=True, help="lncRNA query FASTA")
    parser.add_argument("--output-dir", type=Path, required=True, help="All archives and transient job files go here")
    parser.add_argument("--threads", type=int, required=True, help="Threads requested per run (whole-reference exact scans are serial)")
    parser.add_argument("--parallel-genomes", type=int, default=1, help="Maximum genomes processed simultaneously; --threads applies to EACH genome")
    parser.add_argument("--binary", type=Path, default=Path("build/lncrna_mers_v2"))
    parser.add_argument("--template", type=Path, default=Path(__file__).with_name("default_exact.toml"))
    parser.add_argument("--k-min", type=int, default=15)
    parser.add_argument("--k-max", type=int, default=32)
    parser.add_argument("--dry-run", action="store_true", help="Print plan; do not run or modify outputs")
    parser.add_argument("--max-runs", type=int, default=0, help="Stop after N attempted NEW runs (0: unlimited)")
    parser.add_argument("--keep-work", action="store_true", help="Keep generated configs/results work dir after archival")
    parser.add_argument("--continue-on-error", action="store_true", help="Proceed to next job on failure; preserves failed job folder")
    args = parser.parse_args()

    if args.parallel_genomes <= 0:
        parser.error("--parallel-genomes must be positive")
    if not 1 <= args.threads <= 4096:
        parser.error("--threads must be in [1,4096]")
    if args.k_min < 1 or args.k_max > 32 or args.k_min > args.k_max:
        parser.error("Expected 1 <= --k-min <= --k-max <= 32")
    if args.max_runs < 0:
        parser.error("--max-runs must be >= 0")

    repo_cwd = Path.cwd().resolve()
    genome_dir = args.genomes_dir.expanduser().resolve(strict=True)
    queries = args.queries.expanduser().resolve(strict=True)
    output_root = args.output_dir.expanduser().resolve()
    binary = args.binary.expanduser().resolve(strict=True)
    template_file = args.template.expanduser().resolve(strict=True)
    if not genome_dir.is_dir() or not queries.is_file() or not binary.is_file():
        parser.error("Missing genome directory, query FASTA, or executable")
    if not os.access(binary, os.X_OK):
        parser.error(f"Executable is not executable: {binary}")
    if not template_file.is_file():
        parser.error(f"TOML template does not exist: {template_file}")
    if genome_dir == output_root or output_root in genome_dir.parents:
        parser.error("Output root cannot be an ancestor of the genome input directory")
    if queries == output_root or output_root in queries.parents:
        parser.error("Query input cannot be inside output root")
    genomes = find_genomes(genome_dir)
    if not genomes:
        parser.error("No supported FASTA genome files found")
    source = template_file.read_text(encoding="utf-8")
    tomllib.loads(source)
    # Catch incompatible template values before expensive input hashing.
    check = render_config(source, job_overrides(
        source, genomes[0], queries, output_root / "_work" / "test" / "results",
        args.k_min, args.threads))
    del check
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    requested_threads = args.threads * min(args.parallel_genomes, len(genomes))
    if requested_threads > available:
        print(f"[WARN] Requested up to {requested_threads} total threads; only {available} CPUs in current affinity", file=sys.stderr)
    n_total = len(genomes) * (args.k_max - args.k_min + 1)
    print(f"[PLAN] genomes={len(genomes)} k={args.k_min}..{args.k_max} mode=exact "
          f"runs={n_total} threads_per_genome={args.threads} parallel_genomes={args.parallel_genomes} output={output_root}", flush=True)
    if args.dry_run:
        for g in genomes:
            print(f"  {g.name} => {slug_for(g)} ({(args.k_max - args.k_min + 1)} runs)")
        print("[DRY-RUN] No files created; no genomes queried.")
        return 0

    output_root.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("fraglncrna.grid.exact")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    formatter = logging.Formatter("%(asctime)sZ %(levelname)s %(message)s")
    formatter.converter = time.gmtime
    file_handler = logging.FileHandler(output_root / "runner.log", encoding="utf-8")
    console_handler = logging.StreamHandler(sys.stdout)
    for handler in (file_handler, console_handler):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    try:
        return run_grid(args, genomes, source, binary, queries, template_file,
                        output_root, repo_cwd, n_total, logger)
    finally:
        for handler in (file_handler, console_handler):
            logger.removeHandler(handler)
            handler.close()


def run_grid(args, genomes, source, binary, queries, template_file,
             output_root, repo_cwd, n_total, logger) -> int:
    logger.info("[PLAN] genomes=%d parallel_genomes=%d threads_per_genome=%d total_jobs=%d k=%d..%d mode=exact output=%s",
                len(genomes), args.parallel_genomes, args.threads, n_total,
                args.k_min, args.k_max, output_root)
    for genome in genomes:
        logger.info("[FOUND GENOME] %s => %s", genome, slug_for(genome))
    logger.info("[PROVENANCE] Hashing executable and queries")
    template_sha = hashlib.sha256(source.encode("utf-8")).hexdigest()
    executable_sha = sha256_file(binary)
    queries_sha = sha256_file(queries)
    # Changes to optional annotation/threshold inputs must invalidate resume too.
    auxiliary_inputs = {}
    template_config = tomllib.loads(source)["v2"]
    for key in ("fragments", "reference_repeat_bed", "query_repeat_bed",
                "threshold_table", "expected_loci_bed", "exon_blocks", "window_query_ids"):
        if template_config.get(key):
            path = Path(template_config[key]).resolve(strict=True)
            auxiliary_inputs[key] = {"path": str(path), "sha256": sha256_file(path)}
    attempted = successes = failures = skipped = 0
    state_lock = threading.Lock()
    stop = threading.Event()
    running = {}

    def log_running():
        # Caller holds state_lock, so each snapshot is coherent.
        logger.info("[RUNNING] %d genome(s) simultaneously; genomes=%s",
                    len(running), ", ".join(f"{name} ({job})" for name, job in sorted(running.items())) or "none")

    def run_genome(genome):
        nonlocal attempted, successes, failures, skipped
        if stop.is_set():
            return
        genome_slug = slug_for(genome)
        existing = any((output_root / folder / genome_slug).exists() for folder in ("archives", "_work"))
        logger.info("[%s GENOME] %s; checking each requested exact k completion",
                    "RESUME" if existing else "NEW", genome.name)
        try:
            process_genome(genome, genome_slug)
        except Exception:
            with state_lock:
                failures += 1
                if not args.continue_on_error:
                    stop.set()
            logger.exception("[GENOME ERROR] %s", genome.name)
        finally:
            with state_lock:
                running.pop(genome.name, None)
                log_running()
            logger.info("[GENOME END] %s", genome.name)

    def process_genome(genome, genome_slug):
        nonlocal attempted, successes, failures, skipped
        logger.info("[GENOME] %s: hashing reference for provenance", genome.name)
        genome_sha = sha256_file(genome)
        for k in range(args.k_min, args.k_max + 1):
            label = f"{genome_slug}/exact/k{k}"
            work_dir = output_root / "_work" / genome_slug / "exact" / f"k{k}"
            results_dir = work_dir / "results"
            archive = output_root / "archives" / genome_slug / "exact" / f"k{k}.tar.gz"
            receipt_path = archive.with_name(archive.name + ".complete.json")
            provenance = {
                "format": "fraglncrna-grid-exact-job-v1",
                "genome_name": genome.name,
                "reference": str(genome), "reference_sha256": genome_sha,
                "queries": str(queries), "queries_sha256": queries_sha,
                "executable": str(binary), "executable_sha256": executable_sha,
                "template": str(template_file), "template_sha256": template_sha,
                "auxiliary_inputs": auxiliary_inputs,
                "mode": "exact", "scoring_method": "global_exact", "k": k, "threads": args.threads,
            }
            job_signature = hashlib.sha256(
                json.dumps(provenance, sort_keys=True).encode("utf-8")
            ).hexdigest()
            job = {**provenance, "job_signature": job_signature}
            if stop.is_set():
                return
            try:
                if completion_receipt_valid(archive, receipt_path, job_signature):
                    logger.info("[SKIP COMPLETE] %s", label)
                    with state_lock:
                        skipped += 1
                    continue
                with state_lock:
                    if stop.is_set():
                        return
                    if args.max_runs and attempted >= args.max_runs:
                        logger.info("[LIMIT] --max-runs=%d reached; no new jobs for %s", args.max_runs, genome.name)
                        return
                    attempted += 1
                    attempt = attempted
                    running[genome.name] = f"exact/k{k}"
                    log_running()
                logger.info("[START %d] %s; executable_log=%s", attempt, label, results_dir / "fraglncrna.log")
                if work_dir.exists():
                    prior_job = work_dir / "job.json"
                    if not prior_job.exists() or json.loads(prior_job.read_text(encoding="utf-8")).get("job_signature") != job_signature:
                        raise RuntimeError(f"Work dir exists for different/unknown job; refusing overwrite: {work_dir}")
                    logger.info("[RETRY] %s: rerunning partial exact results", label)
                    # Remove only this failed job's partial results.
                    if results_dir.exists():
                        shutil.rmtree(results_dir)
                results_dir.mkdir(parents=True, exist_ok=True)
                overrides = job_overrides(source, genome, queries, results_dir, k, args.threads)
                generated = render_config(source, overrides)
                (work_dir / "generated_config.toml").write_text(generated, encoding="utf-8")
                write_json(work_dir / "job.json", job)
                t0 = time.monotonic()
                rc = run_executable(binary, work_dir / "generated_config.toml", repo_cwd,
                                    results_dir / "fraglncrna.log", stop)
                elapsed = time.monotonic() - t0
                if rc != 0:
                    raise RuntimeError(f"FragIncRNA exited {rc}; preserved logs/results at {work_dir}")
                output_tsv = results_dir / "lncrna.tsv"
                if not output_tsv.is_file() or output_tsv.stat().st_size == 0:
                    raise RuntimeError(f"FragIncRNA exited 0 but main TSV missing/empty: {output_tsv}")
                receipt = archive_and_verify(work_dir, archive, job, elapsed)
                write_json(receipt_path, receipt)
                if not args.keep_work:
                    shutil.rmtree(work_dir)
                with state_lock:
                    successes += 1
                logger.info("[COMPLETE] %s time=%.1fs archive=%s size=%.1f MiB",
                            label, elapsed, archive, receipt['archive_bytes']/1024**2)
            except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
                with state_lock:
                    failures += 1
                    if not args.continue_on_error:
                        stop.set()
                logger.error("[ERROR] %s: %s", label, exc)
                logger.error("[PRESERVED] %s (inspect or rerun; never use partial outputs)", work_dir)
                if not args.continue_on_error:
                    return
            finally:
                with state_lock:
                    if genome.name in running:
                        running.pop(genome.name)
                        log_running()

    executor = ThreadPoolExecutor(max_workers=args.parallel_genomes)
    try:
        futures = [executor.submit(run_genome, genome) for genome in genomes]
        for future in futures:
            future.result()
    except KeyboardInterrupt:
        stop.set()
        logger.warning("[INTERRUPTED] Stopping active executables; completed archives and work retained")
        raise
    finally:
        executor.shutdown(wait=True, cancel_futures=True)
    logger.info("[SUMMARY] total=%d attempted=%d success=%d skipped=%d failed=%d",
                n_total, attempted, successes, skipped, failures)
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Completed archives preserved; active work folder retained", file=sys.stderr)
        sys.exit(130)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        sys.exit(1)

