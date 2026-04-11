# FragIncRNA

`FragIncRNA` builds per-reference fragment indexes from FASTA references, queries lncRNA sequences against them, and writes match summaries plus unique matching k-mers. The index backend can be either a flat `ibf` or a hierarchical `hibf`.

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

After building with a supported compiler, enable HIBF in [config.toml](/mnt/c/Users/yynk1/Desktop/FragIncRNA/config.toml:1):

```toml
[general]
index_method = "hibf"
```

## Run

The program now reads its settings from a TOML file instead of command-line flags.

1. Edit [config.toml](/mnt/c/Users/yynk1/Desktop/FragIncRNA/config.toml).
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
