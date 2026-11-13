# FragIncRNA v2 sharded IBF grid runner

Requires Python 3.11+ (standard library only) and a working `build/lncrna_mers_v2` executable. Run commands **from the FragIncRNA repository root**. Place `default.toml` next to `run_fraglncrna_grid.py` (they are already packaged together).

## Run one dry-run and one scientific pilot

```bash
python /path/to/run_fraglncrna_grid.py \
  --genomes-dir /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/genomes \
  --queries /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/gencode.v49.lncRNA_transcripts.fa \
  --output-dir /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/production_grid \
  --threads 15 --dry-run

python /path/to/run_fraglncrna_grid.py \
  --genomes-dir /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/genomes \
  --queries /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/gencode.v49.lncRNA_transcripts.fa \
  --output-dir /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/production_grid \
  --threads 15 --max-runs 1
```

The pilot runs **only one k=15, 3-hash experiment**, against the first genome in filename order, and saves the archive. Rerun the same command **without** `--max-runs 1` to continue all remaining jobs:

```bash
python /path/to/run_fraglncrna_grid.py \
  --genomes-dir /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/genomes \
  --queries /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/gencode.v49.lncRNA_transcripts.fa \
  --output-dir /mnt/d/primates_bmc/review/FragIncRNA/rerun_test/production_grid \
  --threads 15
```

```bash
# Bash runs
python3 run_collection/run_fraglncrna_grid.py \
  --genomes-dir /tmp/lutfia95/FragIncRNA/data/genomes\
  --queries /tmp/lutfia95/FragIncRNA/data/lncrna/gencode.v49.lncRNA_transcripts.fa \
  --output-dir /tmp/lutfia95/FragIncRNA/review/FragIncRNA/review_results \
  --threads 45 --dry-run
```

Optional: `--binary /absolute/path/to/lncrna_mers_v2` and `--template /absolute/path/to/default.toml`.

## Scheduling

Use `--parallel-genomes N` to process up to N genomes at the same time (default: 1). Genomes enter the queue in lexicographic filename order. Each genome runs one job at a time: `hash_functions = 3` across k=15..32, then `hash_functions = 4` across k=15..32. When a genome finishes, its worker takes the next genome. One genome = 36 jobs. Eleven genomes = 396 jobs.

For example, run three genomes simultaneously, with 15 threads **per genome** (up to 45 threads total):

```bash
python3 run_collection/run_fraglncrna_grid.py \
  --genomes-dir /path/to/genomes \
  --queries /path/to/lncrna.fa \
  --output-dir /path/to/production_grid \
  --threads 15 --parallel-genomes 3
```

Rerun the same command to resume. Existing genome folders are inspected job by job: verified completed hash/k archives are skipped, incomplete jobs are retried, and missing combinations are run in order. A completely processed genome is skipped so its worker can take another genome. This also fills earlier gaps rather than assuming everything before the highest k is complete. You can change `--parallel-genomes` when resuming; inputs, template, executable and per-genome thread settings must still match the recorded provenance.

`--max-runs N` limits new job attempts across **all workers together**, excluding completed jobs that are skipped. With concurrency, the particular genomes receiving those attempts depend on worker timing. Without `--continue-on-error`, a failure stops new work and terminates active executables, preserving their work directories for retry. Ctrl+C also stops active executables and preserves work.

## Runner logging

`OUTPUT_DIR/runner.log` is a separate, timestamped UTC orchestration log, appended on every invocation and also displayed on the console. It records:

- Discovered genome filenames and their output folder names.
- New or existing genome folders, completed hash/k skips, and partial-job retries.
- Job starts with their executable log paths.
- Current active genome counts and names, including each genome's hash/k job.
- Completions, archive paths, failures, interruptions, and final attempted/success/skipped/failed counts.

Each job's executable stdout/stderr goes to its own `_work/GENOME/hashH/kK/results/fraglncrna.log`, which is included in its successful archive. It is separate from `runner.log`, so output from parallel executables does not mix. Failed-job logs remain in the work folder. Monitor orchestration with:

```bash
tail -f /path/to/production_grid/runner.log
```

## Storage, safety and provenance

- Results/workspace share **one output root**, e.g. `output/_work/` and `output/archives/`.
- Each run has an isolated index directory under its private `_work/.../cache/ibf/`.
- A successful CLI exit and a nonempty main TSV are required before archive creation.
- The runner removes **only its own per-job cache** after the CLI succeeds; it never deletes user indexes outside `_work/`.
- Only results, log, sidecars, job/provenance and actual generated configuration are included in the `.tar.gz`. The runner reopens the archive and verifies every member against its SHA-256, publishes an external SHA-256 completion receipt, and removes the work directory.
- On an unsuccessful job, its private files remain. On retry, the script keeps its compatible IBF cache, removes only the partial results folder, and reruns that job. A different executable/template/inputs for an existing failed job require a new output root.
- On restart, completed `.tar.gz` archives are checked against the receipts and skipped. A missing/mismatched receipt causes a failure instead of silent overwriting.
- Genome SHA-256 is calculated once per genome for provenance; query FASTA and executable SHA-256 are calculated once per invocation.
- `--continue-on-error` can run later jobs even if one fails, but failed-work IBFs stay on disk and may consume space.
- This runner uses **tar.gz**, because gzip alone cannot package a directory into one file.
- Existing live jobs and caches outside the chosen `--output-dir` are never touched.

## Checking an archive

```bash
tar -tzf /path/to/output/archives/GENOME/hash3/k15.tar.gz | head -30
# To extract later:
mkdir -p extracted
 tar -xzf /path/to/output/archives/GENOME/hash3/k15.tar.gz -C extracted
```

## Notes

- The full-scale 396-run experiment can take a long time. Measure one full run first and ensure sufficient memory and free disk space for all concurrently active IBF sets and the final archives.
- k=32 is the current maximum supported by the supplied configuration documentation; verify boundary behavior on a small dataset before launching the full grid.
- `fdr_control_per_bin=true` is historical nomenclature, **not** measured statistical FDR. `per_bin_verified` is candidate-local, **not** global-exact recovery.
- No separate k-mer FASTA is enabled by this template. Per-window output is also not enabled; configure optional reports separately if needed.
- The template defaults to `repeat=include`, `fpr=0.01`, `fragment_size=80000`, `per_bin_min_fraction=0.02`, `per_bin_top_k=256`, and 75%/85% transcript thresholds.
- Prefer running heavy workloads on native Linux storage where possible instead of WSL `/mnt/d/`, especially for repeated shard reads.


python run_fraglncrna_grid.py \
  --genomes-dir "/mnt/d/primates_bmc/review/FragIncRNA/rerun_test/genomes" \
  --queries "/mnt/d/primates_bmc/review/FragIncRNA/rerun_test/gencode.v49.lncRNA_transcripts.fa" \
  --output-dir "/mnt/d/primates_bmc/review/FragIncRNA/rerun_test/production_grid" \
  --threads 15