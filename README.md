# FragIncRNA

The selectable production workflow is now available as schema 3. Use
`--mode ibf --scoring-method per_bin_verified` for IBF candidates followed by
exact verification within indexed candidate intervals. The main recovery and
status then use candidate-local verified support. It is a lower bound on
whole-reference exact recovery; it is not the same biological quantity.
See the [production audit](docs/PER_BIN_PRODUCTION_AUDIT.md),
[current method](docs/PER_BIN_IBF_METHOD.md),
[measured benchmark](docs/PER_BIN_PRODUCTION_BENCHMARK.md), and
[acceptance status](docs/PER_BIN_PRODUCTION_STATUS.md).

Ready-to-run fixture examples, after the measured build below:

```bash
python3 scripts/measure_command.py --output /tmp/frag-configure.json -- cmake -S . -B build-production -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLNCRNA_MERS_ENABLE_HIBF=OFF
python3 scripts/measure_command.py --output /tmp/frag-build.json -- cmake --build build-production -j 4
python3 scripts/measure_command.py --output /tmp/frag-tests.json -- ctest --test-dir build-production --output-on-failure
python3 scripts/measure_command.py --output /tmp/frag-production.json -- build-production/lncrna_mers_v2 --config configs/per_bin_production.toml
python3 scripts/measure_command.py --output /tmp/frag-validate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_production.toml --mode validate --output /tmp/frag-validate.tsv
python3 scripts/measure_command.py --output /tmp/frag-calibrate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_production.toml --mode calibrate --output /tmp/frag-calibrate.tsv
```

The fixture uses k=4/4-base bins and a one-hit threshold solely for executable
checks. For research, set the supplied v49 query/reference paths, k=20,
fragment_size=80000 and per_bin_min_hits=2; the 5%/top-64 rule is heuristic.
Run one reference per invocation. The [bounded runner](scripts/run_per_bin_production_benchmark.py)
processes 100 then 5,000 queries and never launches the full manuscript grid.
The data-preparation script stages and validates an extraction before exposing it.

After preparing the supplied v49 subsets and complete first contig, these use
the same real-data settings as the benchmark (outputs are new example files):

```bash
python3 scripts/measure_command.py --output /tmp/v49-production.json -- build-production/lncrna_mers_v2 --config configs/per_bin_v49_bounded.toml --load-index results/per_bin_production_validation/complete_chr1/k20_nominal.ibf --reference-cache results/per_bin_production_validation/complete_chr1/reference.bases --output /tmp/v49-production.tsv
python3 scripts/measure_command.py --output /tmp/v49-validate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_v49_bounded.toml --mode validate --load-index results/per_bin_production_validation/complete_chr1/k20_nominal.ibf --reference-cache results/per_bin_production_validation/complete_chr1/reference.bases --output /tmp/v49-validate.tsv
python3 scripts/measure_command.py --output /tmp/v49-calibrate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_v49_bounded.toml --mode calibrate --load-index results/per_bin_production_validation/complete_chr1/k20_nominal.ibf --reference-cache results/per_bin_production_validation/complete_chr1/reference.bases --output /tmp/v49-calibrate.tsv
```

Omit `--load-index` on first use to build instead. The indexed reference here
is one complete macaque contig, not all 22 contigs. Validation/calibration scan
the full supplied reference independently; production does not. **Do not use
the default candidate filter as a drop-in global-recovery replacement:** on
the 5,000-query sample it missed 36.34% of exact-positive windows and 3,005 of
4,482 queries with any exact support. Biological accuracy is not established.

Choose `global_ibf` for historical any-bin recovery, `per_bin_ibf` for
approximate support in selected bins, `per_bin_verified` for exact membership
in selected intervals, or `global_exact` for independent full-reference
membership. Validation reports all four, even if exact remains primary.
Old unspecified exact/global invocations retain schema 2; per-bin IBF
production and explicit scoring use schema 3 with method/scope/provenance.
The merger checks source, scoring, schema, thresholds and candidate settings.

`--ibf-auto-tune off|recommend|apply` diagnoses settings and writes
`<output>.ibf_tuning.json`. Explicit off/recommend retains settings; apply
can budget a global-IBF objective using `--genome-fpr`. Per-bin methods retain
their design settings and recommend measured candidate-recall validation.
Unspecified legacy global-IBF invocations preserve their previous automatic
whole-reference budget; this compatibility rule is explicitly logged.
`--ibf-max-index-memory-gb` preflights raw index bytes; it is not a total
process RSS cap. Dry-run uses `--ibf-dry-run`. Loaded archives are checked and
never resized/rebuilt silently. Fragment alternatives are estimates, never
automatic scientific-setting changes. Per-bin Bloom FPR, approximate union
FPR, transcript classification errors and biological FDR are distinct.

`--threads 15` uses OpenMP for bounded query and interval work; actual workers
are capped by available resources and recorded. Default is 1. Explicit threads
wins over OMP_NUM_THREADS. Index building/global exact scans remain serial.
Without OpenMP the binary reports serial fallback. A normalized original-base
cache enables random access for plain/ordinary gzip FASTA; use
`--reference-cache FILE` for reuse. Original-reference inspection and cache
checksum cost remain in total time. The candidate-only byte counter excludes
this separately reported provenance/preparation I/O.

Persistent UTC text logs and matching JSONL default to `<output>.log` and
`<output>.log.jsonl`. Levels are trace/debug/info/warning/error/off;
`--no-log-to-console` disables console copies. Progress uses
`--progress-interval-seconds` and bounded `--progress-every-queries` batches.
Debug adds per-query scoring events; `--log-per-query`/`--log-per-bin` explicitly
enable expensive detail. Lines flush immediately. Signals/worker exceptions
produce terminal error events where graceful termination is possible.
Schema-3 main TSV/metrics finalize by rename; failed sidecars/partial files
must not be treated as successful runs. SIGKILL cannot produce a footer.

The older schema-2 documentation below is retained for compatibility/history;
its full-reference candidate verifier and IBF-only restrictions are superseded.
[![CMake tests](https://github.com/GPSeq/FragIncRNA/actions/workflows/cmake-single-platform.yml/badge.svg?branch=reporting)](https://github.com/GPSeq/FragIncRNA/actions/workflows/cmake-single-platform.yml)
`FragIncRNA` builds per-reference fragment indexes from FASTA references, queries lncRNA sequences against them, and writes match summaries plus unique matching k-mers. The index backend can be either a flat `ibf` or a hierarchical `hibf`.

## Corrected k-mer analysis (schema 2 preview)

`lncrna_mers_v2` is the corrected analysis path. It reports per-window recovery in [0,1], default canonical strand matching, independent exact validation, and wall-time/peak-memory measurements. Its default `validate` mode now classifies from exact verified windows and reports Bloom-filter recovery separately. The TOML driven `lncrna_mers [config.toml]` retains its original output columns, but its flat-IBF path now shares v2's position-preserving, ambiguity-safe, canonical k-mer encoder and uses `k-1` fragment overlap. `negative_set` shares those flat-IBF rules too. Historical flat-IBF result values can therefore change when rerun. The legacy `*_pct` still sums bin hits and can exceed 1; do not treat it as a recovery fraction. The v2 path supports flat IBF for one reference per invocation. See [methods](docs/KMER_METHODS.md), [audit](docs/KMER_AUDIT.md), [resource measurements](docs/KMER_BENCHMARKS.md), [real chr1 run](docs/KMER_REAL_DATA_FIRST_RUN.md), and [implementation status](docs/KMER_IMPLEMENTATION_STATUS.md) before using it for biological conclusions.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLNCRNA_MERS_ENABLE_HIBF=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j 4
ctest --test-dir build --output-on-failure
./build/lncrna_mers_v2 --reference tests/fixtures/kmer_v2/reference.fa --queries tests/fixtures/kmer_v2/queries.fa --output /tmp/fraglncrna-v2.tsv --k 4 --fragment-size 4 --mode validate
./build/lncrna_mers_v2 --config configs/example.toml
```

The CLI example writes `/tmp/fraglncrna-v2.tsv` and its metrics JSON. The fully commented [TOML example](configs/example.toml) covers every v2 option; `--config` can be followed by CLI overrides, and the config file's SHA-256 is recorded in the run JSON. Legacy `lncrna_mers` continues to use its existing `configs/config_*.toml` layout. Run `./build/lncrna_mers_v2 --help` for options. Plain and `.gz` FASTA are accepted. The output is one row per query and reference, with stable FASTA IDs and explicit zero-hit rows. `--mode exact` skips IBF construction; use it for the supplied chr1 setting. IBF-only mode now sizes the per-bin filter for a default whole-reference FPR target of 0.01 (`--genome-fpr`) and checks observed bit occupancy before classifying. `--mode calibrate` measures genome-wide per-window and per-bin false positives against exact truth and writes a `.bins.tsv` report with distinct k-mer counts and actual bit occupancy. `--save-index reference.ibf` stores the flat IBF with checked compatibility metadata; `--load-index reference.ibf` reuses it. The metrics JSON contains semantic FASTA SHA-256 hashes, auxiliary input hashes, archive size, wall time, and peak RSS. `--fragments fragments.tsv` writes explicit reference bin and overlap coordinates. `--local-window 8000 --positions /tmp/verified-hits.tsv.gz` adds exact positions, fixed-region support, and sliding-window sensitivity with a bounded occurrence cap. `--discover-loci candidates.tsv` groups exact hits into ranked locus candidates; `--expected-loci-bed expected.bed` writes a separate expected-locus report and per-query recovery columns. `--exon-blocks blocks.tsv` reports eligible-window recovery at supplied splice-alignment coordinates. `--reference-repeat-bed`, `--query-repeat-bed`, and `--mask-target` control repeat annotation and exclusion; `--repeat annotate` writes query class counts and `--reference-hit-classes hits.tsv` writes exact reference-side class attribution. A threshold table keyed by reference/k can be fitted on separate training and evaluation controls with `scripts/fit_kmer_thresholds.py`; `scripts/merge_kmer_v2.py` streams verified per-reference runs into a query-major comparison. Exact validation rescans the reference for each query batch, so benchmark a subset before planning a full-genome run.

For the first real-data check, use the [measured 5,000-transcript runbook](docs/KMER_REAL_DATA_FIRST_RUN.md). `scripts/sample_fasta.py` makes a seeded subset and writes its own time and peak-memory JSON.

### Optional per-bin coherence screen

`--fdr-control-per-bin` adds candidate fragment scores and exact candidate-region
verification without changing the existing main TSV columns. The requested
flag name does **not** imply statistical FDR control: it measures how many
original query windows hit each IBF bin, groups nearby candidate bins, and
checks selected regions against the original reference. The paired per-query
table defaults to `<output>.per_bin.tsv`; candidate regions go to
`<output>.per_bin.tsv.candidates.tsv`. Use `--per-bin-min-hits`,
`--per-bin-min-fraction`, `--per-bin-top-k`, and
`--per-bin-max-locus-gap` to set the heuristic screen. The same Boolean and
numeric options are available in `[v2]` TOML; see
[the complete example](configs/example.toml) and [method definition](docs/PER_BIN_IBF_METHOD.md).
For an unverified historical-index diagnostic, add `--mode ibf
--no-per-bin-exact-verify`; the separate candidate status must not be read as
global recovery. The [comparison report](docs/PER_BIN_IBF_BENCHMARK.md) records
measured chr1 outcomes and resource costs.

## Requirements

- CMake 3.20+
- A C++20 compiler

SeqAn3 is fetched automatically by CMake during configure. This project currently pins SeqAn3 `3.4.1`.
HIBF support is optional and currently requires a newer compiler toolchain than the one used for the older flat IBF-only setup. In practice, use GCC 12+ and CMake 3.25+ for HIBF builds.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

If sdsl-lite is causing any errors from seqan3 build, use:

```bash
git clone https://github.com/xxsds/sdsl-lite.git
cd sdsl-lite
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build -j
cmake --install build
```

This produces the executable `build/lncrna_mers`.

## HIBF Toolchain Notes

If CMake prints:

```text
HIBF support disabled: GCC 11.4.0 is too old. HIBF requires GCC >= 12.
```

then your current compiler is too old for the HIBF library. The flat `ibf` backend can still be built, but `hibf` needs GCC 12 or newer.

On Ubuntu or WSL Ubuntu, install GCC 12:

```bash
sudo apt update
sudo apt install gcc-12 g++-12
```

Check that the compiler is available:

```bash
which g++-12
g++-12 --version
```

If CMake reports:

```text
The CMAKE_CXX_COMPILER:
  g++-12
is not a full path and was not found in the PATH.
```

then either `g++-12` is not installed, or it is not visible on `PATH`. In that case, use the full compiler paths:

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-12 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-12
cmake --build build -j
```

You can also set the environment variables instead:

```bash
export CC=/usr/bin/gcc-12
export CXX=/usr/bin/g++-12
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
```

If your CMake version is very new and configure fails while processing the fetched SeqAn3 project, add:

```bash
-DCMAKE_POLICY_VERSION_MINIMUM=3.5
```

Example:

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j
```

After building with a supported compiler, enable HIBF in [config.toml](config.toml):

```toml
[general]
index_method = "hibf"
```

## Run

The program now reads its settings from a TOML file instead of command-line flags.

1. Edit [config.toml](config.toml).
2. Run:

```bash
./build/lncrna_mers
```

By default it loads `./config.toml`. You can also pass a different TOML file path:

```bash
./build/lncrna_mers /path/to/config.toml
```

Example `config.toml`:

```toml
[general]
index_method = "ibf"

ref_dir = "/path/to/references"
query_file = "/path/to/queries.fa"

output_dir = "./ibf_results"
output_file = "results.tsv"
log_file = "ibf_run.log"

fragment_size = 8000
kmer_size = 15
hit_threshold = 13
threads = 8

store_fragments = false
store_index = false
cleanup_index = false
single_results_writer = false

[ibf]
hash_functions = 3
fpr = 0.01
multi_indexing = false

[hibf]
hash_functions = 2
maximum_fpr = 0.05
relaxed_fpr = 0.30
threads = 8
sketch_bits = 12
tmax = 0
empty_bin_fraction = 0.0
alpha = 1.2
max_rearrangement_ratio = 0.5
disable_estimate_union = false
disable_rearrangement = false
```

## Config Keys

The TOML file is split into three sections:

- `[general]`
- `[ibf]`
- `[hibf]`

General keys:

- `general.index_method`: `"ibf"` or `"hibf"`
- `general.ref_dir`: directory containing `.fa`, `.fna`, `.fasta`, and gzipped variants
- `general.query_file`: FASTA file with query sequences
- `general.fragment_size`: fragment length used for reference splitting; must be at least `4`
- `general.kmer_size`: k-mer length used for indexing and querying; valid range `1..32`
- `general.hit_threshold`: minimum total k-mer hits required for a passing match
- `general.threads`: number of references to process in parallel; use `0` to auto-detect hardware concurrency
- `general.output_dir`: directory for logs, results, stored indexes, and optional fragment FASTA files
- `general.output_file`: combined TSV file name when `single_results_writer = true`
- `general.log_file`: log file name written under `output_dir`
- `general.store_fragments`: when `true`, writes `<reference>_fragments.fasta`
- `general.store_index`: when `true`, stores one serialized index file per reference
- `general.cleanup_index`: when `true`, removes the stored index file after processing
- `general.single_results_writer`: `true` writes one combined TSV; `false` writes one `results_<reference>.tsv` per reference

IBF-specific keys:

- `ibf.hash_functions`: number of IBF hash functions; valid range `1..32`
- `ibf.fpr`: target false positive rate per bin; must be in `(0, 0.5]`
- `ibf.multi_indexing`: when `true`, builds one flat IBF per reference fragment, splits each fragment into near-equal sub-bins, stores the `.ibf` files in `output_dir` while searching, and removes them after that reference is processed

HIBF-specific keys:

- `hibf.hash_functions`: number of HIBF hash functions; valid range `1..5`
- `hibf.maximum_fpr`: maximum false positive rate; must be in `(0, 1)`
- `hibf.relaxed_fpr`: relaxed FPR for merged bins; must be in `(0, 1)` and `>= maximum_fpr`
- `hibf.threads`: internal HIBF build threads; `0` falls back to `general.threads`
- `hibf.sketch_bits`: HyperLogLog sketch bits; valid range `5..32`
- `hibf.tmax`: maximum technical bins for each IBF level; `0` lets HIBF choose
- `hibf.empty_bin_fraction`: reserved empty-bin fraction; valid range `[0.0, 1.0)`
- `hibf.alpha`: HIBF alpha parameter; must be positive
- `hibf.max_rearrangement_ratio`: valid range `[0.0, 1.0]`
- `hibf.disable_estimate_union`: disables union estimation during layout
- `hibf.disable_rearrangement`: disables rearrangement during layout

## Output

Combined mode writes a TSV with these columns per reference:

- `*_count`: total IBF hits across bins for the query/reference pair
- `*_ibf_unique_kmer`: number of distinct query k-mers that
  1. have IBF count exactly `1` for that reference
  2. are counted once per distinct k-mer, even if they appear multiple times in the query
- `*_pass`: `1` if `count >= hit_threshold`, else `0`
- `*_pct`: `count / number_of_query_kmers`

Unique kmers are also written under `unique_mers/<reference>.tsv`, with
`query_index` in the first column using the 0-based order from the query FASTA file.
All matching query k-mers are additionally written to `<output_dir>/<reference>_kmers.tsv`
with three columns: `query_index`, slash-separated 1-based `kmer_indices`, and the
matching slash-separated `kmer_counts` for those indices.
The main results TSV also writes `query_index` instead of the full FASTA query name.

## Tests

Configure and build as usual, then run:

```bash
ctest --test-dir build --output-on-failure
```

To run only the unit test binary directly:

```bash
./build/lncrna_mers_tests
```
