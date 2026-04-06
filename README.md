# FragIncRNA

`FragIncRNA` builds per-reference interleaved Bloom filters (IBFs) from FASTA references, queries lncRNA sequences against them, and writes match summaries plus unique matching k-mers.

## Requirements

- CMake 3.20+
- A C++20 compiler

SeqAn3 is fetched automatically by CMake during configure.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This produces the executable `build/lncrna_mers`.

## Run

The program now reads its settings from a TOML file instead of command-line flags.

1. Edit [`config.toml`](/Users/ahmadlutfi/Downloads/FragIncRNA/FragIncRNA/config.toml).
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
ref_dir = "/path/to/references"
query_file = "/path/to/queries.fa"

output_dir = "./ibf_results"
output_file = "results.tsv"
log_file = "ibf_run.log"

fragment_size = 8000
kmer_size = 15
hash_functions = 3
fpr = 0.01
hit_threshold = 13
threads = 8

store_fragments = false
store_ibf = true
cleanup_ibf = true
single_results_writer = true
```

## Config Keys

The TOML file uses these argument names:

- `ref_dir`: directory containing `.fa`, `.fna`, `.fasta`, and gzipped variants
- `query_file`: FASTA file with query sequences
- `fragment_size`: fragment length used for reference splitting; must be at least `4`
- `kmer_size`: k-mer length used for indexing and querying; valid range `1..32`
- `hash_functions`: number of IBF hash functions; valid range `1..32`
- `fpr`: target false positive rate per bin; must be in `(0, 0.5]`
- `hit_threshold`: minimum total k-mer hits required for a passing match
- `threads`: number of references to process in parallel; use `0` to auto-detect hardware concurrency
- `output_dir`: directory for logs, results, IBFs, and optional fragment FASTA files
- `output_file`: combined TSV file name when `single_results_writer = true`
- `log_file`: log file name written under `output_dir`
- `store_fragments`: when `true`, writes `<reference>_fragments.fasta`
- `store_ibf`: when `true`, writes one `<reference>.ibf` file per reference
- `cleanup_ibf`: when `true`, removes stored IBFs after processing
- `single_results_writer`: `true` writes one combined TSV; `false` writes one `results_<reference>.tsv` per reference

## Output

Combined mode writes a TSV with these columns per reference:

- `*_count`: total IBF hits across bins for the query/reference pair
- `*_ibf_unique_kmer`: number of distinct query k-mers that
  1. have IBF count exactly `1` for that reference
  2. are counted once per distinct k-mer, even if they appear multiple times in the query
- `*_pass`: `1` if `count >= hit_threshold`, else `0`
- `*_pct`: `count / number_of_query_kmers`

IBF unique kmers are also written under `unique_mers/<reference>.tsv`, with
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
