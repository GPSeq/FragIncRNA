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

If sdsl-lite is causing any errors from seqan3 build, use: 

```bash
git clone https://github.com/xxsds/sdsl-lite.git
cd sdsl-lite
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build -j
cmake --install build
```

This produces the executable `build/lncrna_mers`.

## Run

Example:

```bash
./build/lncrna_mers \
  --ref-dir /path/to/references \
  --query-file /path/to/queries.fa \
  --fragment-size 8000 \
  --threads 8 \
  --kmer-size 15 \
  --hash-functions 3 \
  --hit-threshold 13 \
  --output-dir ./ibf_results \
  --output-file results.tsv \
  --log-file ibf_run.log \
  --store-ibf \
  --cleanup-ibf
```

Important options:

- `--ref-dir`: directory containing `.fa`, `.fna`, `.fasta` and gzipped variants
- `--query-file`: FASTA file with query sequences
- `--fragment-size`: fragment length used for reference splitting; must be at least `4`
- `--threads`, `-j`: number of references to process in parallel; defaults to hardware concurrency
- `--kmer-size`: k-mer length used for indexing and querying
- `--hash-functions`: number of IBF hash functions
- `--fpr`: target false positive rate per bin, default `0.01`
- `--hit-threshold`: minimum total k-mer hits required for a passing match
- `--output-dir`: directory for logs, results, IBFs, and optional fragment FASTA files
- `--store-fragments`: writes `<reference>_fragments.fasta`
- `--store-ibf`: writes one `<reference>.ibf` file per reference
- `--cleanup-ibf`: removes stored IBFs after processing
- `--single-results-writer false`: writes one `results_<reference>.tsv` per reference instead of one combined TSV

## Output

Combined mode writes a TSV with these columns per reference:

- `*_count`: total IBF hits across bins for the query/reference pair
- `*_unique_kmers`: number of distinct matching k-mers
- `*_pass`: `1` if `count >= --hit-threshold`, else `0`
- `*_pct`: `count / number_of_query_kmers`

Unique matching k-mers are also written under `unique_mers/<reference>.tsv`.

## Tests

Configure and build as usual, then run:

```bash
ctest --test-dir build --output-on-failure
```

To run only the unit test binary directly:

```bash
./build/lncrna_mers_tests
```
