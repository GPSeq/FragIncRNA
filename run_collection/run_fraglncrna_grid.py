#!/usr/bin/env python3
"""Run FragIncRNA v2 sequentially for genome FASTAs × k=15..32 × hashes=3,4.

Python >= 3.11, standard library only. One independent, provenance-checked
.tar.gz per successful run. IBF shard caches live INSIDE the private per-job
_work directory; the runner removes them only after success and archive
verification. Failed jobs retain their caches and logs for safe retry.

Run from your FragIncRNA repository root. Use --dry-run first.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tarfile
import time
import tomllib

FASTA_SUFFIXES = (
    ".fa", ".fasta", ".fna", ".fas",
    ".fa.gz", ".fasta.gz", ".fna.gz", ".fas.gz",
)
OVERRIDE_KEYS = {
    "reference", "queries", "output", "metrics", "threads", "k",
    "hash_functions", "ibf_shard_dir", "reference_cache", "log_file",
    "ibf_tuning_report", "ibf_build_only", "ibf_dry_run",
    "save_index", "load_index",
}


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
        return json.dumps(value, ensure_ascii=True)
    raise TypeError(f"Unsupported TOML override: {type(value).__name__}")


def render_config(template: str, overrides: dict[str, object]) -> str:
    """Rewrite only [v2], preserving all other template settings/comments."""
    if "v2" not in tomllib.loads(template):
        raise ValueError("Template must contain [v2]")
    lines = template.splitlines(keepends=True)
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
    if config.get("mode") != "ibf" or config.get("ibf_layout") != "sharded":
        raise ValueError('Template must set mode="ibf" and ibf_layout="sharded"')
    if config.get("scoring_method") != "per_bin_verified":
        raise ValueError('Template must use scoring_method="per_bin_verified"')
    if config.get("ibf_shard_storage") != "disk":
        raise ValueError('Template must set ibf_shard_storage="disk"')
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
    """Archive complete outputs (not IBFs); re-read every member for integrity."""
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


def run_executable(binary: Path, config: Path, cwd: Path, log_file: Path) -> int:
    argv = [str(binary), "--config", str(config)]
    print("$ " + " ".join(argv), flush=True)
    with log_file.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            argv, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", bufsize=1,
            start_new_session=True,
        )
        try:
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                print(line, end="", flush=True)
            return process.wait()
        except KeyboardInterrupt:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            process.wait()
            raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--genomes-dir", type=Path, required=True, help="Flat directory containing .fa/.fna/.fasta[.gz]")
    parser.add_argument("--queries", type=Path, required=True, help="lncRNA query FASTA")
    parser.add_argument("--output-dir", type=Path, required=True, help="All archives and transient IBF caches go here")
    parser.add_argument("--threads", type=int, required=True, help="OpenMP threads per run")
    parser.add_argument("--binary", type=Path, default=Path("build/lncrna_mers_v2"))
    parser.add_argument("--template", type=Path, default=Path(__file__).with_name("default.toml"))
    parser.add_argument("--k-min", type=int, default=15)
    parser.add_argument("--k-max", type=int, default=32)
    parser.add_argument("--hash-functions", nargs="+", type=int, default=[3, 4], metavar="H")
    parser.add_argument("--dry-run", action="store_true", help="Print plan; do not run or modify outputs")
    parser.add_argument("--max-runs", type=int, default=0, help="Stop after N attempted NEW runs (0: unlimited)")
    parser.add_argument("--keep-work", action="store_true", help="Keep generated configs/results work dir after archival (IBF cache still removed)")
    parser.add_argument("--continue-on-error", action="store_true", help="Proceed to next job on failure; preserves failed job folder")
    args = parser.parse_args()

    if args.threads <= 0:
        parser.error("--threads must be positive")
    if args.k_min < 1 or args.k_max > 32 or args.k_min > args.k_max:
        parser.error("Expected 1 <= --k-min <= --k-max <= 32")
    if len(set(args.hash_functions)) != len(args.hash_functions) or any(h < 1 or h > 5 for h in args.hash_functions):
        parser.error("Hash counts must be distinct and in [1,5]")
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
    if genome_dir == output_root or genome_dir in output_root.parents:
        parser.error("Output root cannot be an ancestor of the genome input directory")
    if queries == output_root or output_root in queries.parents:
        parser.error("Query input cannot be inside output root")
    genomes = find_genomes(genome_dir)
    if not genomes:
        parser.error("No supported FASTA genome files found")
    source = template_file.read_text(encoding="utf-8")
    tomllib.loads(source)
    # Catch incompatible template values before expensive input hashing.
    check = render_config(source, {
        "reference": str(genomes[0]), "queries": str(queries), "output": str(output_root / "_work" / "test.tsv"),
        "metrics": "", "threads": args.threads, "k": args.k_min, "hash_functions": args.hash_functions[0],
        "ibf_shard_dir": str(output_root / "_work" / "index"), "reference_cache": "", "log_file": "",
        "ibf_tuning_report": "", "ibf_build_only": False, "ibf_dry_run": False,
        "save_index": "", "load_index": "",
    })
    del check
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    if args.threads > available:
        print(f"[WARN] Requested {args.threads} threads; only {available} CPUs in current affinity", file=sys.stderr)
    n_total = len(genomes) * (args.k_max - args.k_min + 1) * len(args.hash_functions)
    print(f"[PLAN] genomes={len(genomes)} k={args.k_min}..{args.k_max} hashes={args.hash_functions} "
          f"runs={n_total} threads={args.threads} output={output_root}", flush=True)
    if args.dry_run:
        for g in genomes:
            print(f"  {g.name} => {slug_for(g)} ({(args.k_max - args.k_min + 1)*len(args.hash_functions)} runs)")
        print("[DRY-RUN] No files created; no genomes queried.")
        return 0

    output_root.mkdir(parents=True, exist_ok=True)
    template_sha = hashlib.sha256(source.encode("utf-8")).hexdigest()
    executable_sha = sha256_file(binary)
    queries_sha = sha256_file(queries)
    attempted = 0
    successes = 0
    failures = 0
    skipped = 0

    for genome in genomes:
        genome_slug = slug_for(genome)
        print(f"[GENOME] {genome.name}: hashing reference for provenance", flush=True)
        genome_sha = sha256_file(genome)
        for hashes in args.hash_functions:
            for k in range(args.k_min, args.k_max + 1):
                label = f"{genome_slug}/hash{hashes}/k{k}"
                work_dir = output_root / "_work" / genome_slug / f"hash{hashes}" / f"k{k}"
                results_dir = work_dir / "results"
                cache_dir = work_dir / "cache" / "ibf"
                archive = output_root / "archives" / genome_slug / f"hash{hashes}" / f"k{k}.tar.gz"
                receipt_path = archive.with_name(archive.name + ".complete.json")
                provenance = {
                    "format": "fraglncrna-grid-job-v1",
                    "genome_name": genome.name,
                    "reference": str(genome), "reference_sha256": genome_sha,
                    "queries": str(queries), "queries_sha256": queries_sha,
                    "executable": str(binary), "executable_sha256": executable_sha,
                    "template": str(template_file), "template_sha256": template_sha,
                    "k": k, "hash_functions": hashes, "threads": args.threads,
                }
                job_signature = hashlib.sha256(
                    json.dumps(provenance, sort_keys=True).encode("utf-8")
                ).hexdigest()
                job = {**provenance, "job_signature": job_signature}
                if completion_receipt_valid(archive, receipt_path, job_signature):
                    print(f"[SKIP COMPLETE] {label}", flush=True)
                    skipped += 1
                    continue
                if args.max_runs and attempted >= args.max_runs:
                    print(f"[STOP] --max-runs={args.max_runs} reached", flush=True)
                    print(f"[SUMMARY] success={successes} skipped={skipped} failed={failures}")
                    return 1 if failures else 0
                attempted += 1
                print(f"\n[START {attempted}] {label} at {iso_now()}", flush=True)
                try:
                    if work_dir.exists():
                        prior_job = work_dir / "job.json"
                        if not prior_job.exists() or json.loads(prior_job.read_text(encoding="utf-8")).get("job_signature") != job_signature:
                            raise RuntimeError(f"Work dir exists for different/unknown job; refusing overwrite: {work_dir}")
                        # Reuse only the compatible shard cache. Remove any partial result files.
                        if results_dir.exists():
                            shutil.rmtree(results_dir)
                    results_dir.mkdir(parents=True, exist_ok=True)
                    cache_dir.mkdir(parents=True, exist_ok=True)
                    overrides = {
                        "reference": str(genome), "queries": str(queries),
                        "output": str(results_dir / "lncrna.tsv"), "metrics": "",
                        "k": k, "hash_functions": hashes, "threads": args.threads,
                        "ibf_shard_dir": str(cache_dir), "reference_cache": "", "log_file": "",
                        "ibf_tuning_report": "", "ibf_build_only": False, "ibf_dry_run": False,
                        "save_index": "", "load_index": "",
                    }
                    generated = render_config(source, overrides)
                    (work_dir / "generated_config.toml").write_text(generated, encoding="utf-8")
                    write_json(work_dir / "job.json", job)
                    t0 = time.monotonic()
                    rc = run_executable(binary, work_dir / "generated_config.toml", repo_cwd,
                                        results_dir / "fraglncrna.log")
                    elapsed = time.monotonic() - t0
                    if rc != 0:
                        raise RuntimeError(f"FragIncRNA exited {rc}; preserved logs/cache at {work_dir}")
                    output_tsv = results_dir / "lncrna.tsv"
                    if not output_tsv.is_file() or output_tsv.stat().st_size == 0:
                        raise RuntimeError(f"FragIncRNA exited 0 but main TSV missing/empty: {output_tsv}")
                    # Remove only the scratch IBFs created inside THIS JOB'S folder.
                    # This frees ~GBs before the archive compression phase.
                    cache_base = work_dir / "cache"
                    if cache_base.exists():
                        shutil.rmtree(cache_base)
                    receipt = archive_and_verify(work_dir, archive, job, elapsed)
                    write_json(receipt_path, receipt)
                    if not args.keep_work:
                        shutil.rmtree(work_dir)
                    successes += 1
                    print(f"[COMPLETE] {label} time={elapsed:.1f}s archive={archive} "
                          f"size={receipt['archive_bytes']/1024**2:.1f} MiB; IBFs removed", flush=True)
                except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
                    failures += 1
                    print(f"[ERROR] {label}: {exc}", file=sys.stderr, flush=True)
                    print(f"[PRESERVED] {work_dir} (inspect or rerun; never use partial outputs)", file=sys.stderr)
                    if not args.continue_on_error:
                        print(f"[SUMMARY] success={successes} skipped={skipped} failed={failures}")
                        return 1

    print(f"\n[SUMMARY] total={n_total} success={successes} skipped={skipped} failed={failures}")
    return 1 if failures else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[INTERRUPTED] Completed archives preserved; active work folder retained", file=sys.stderr)
        sys.exit(130)
