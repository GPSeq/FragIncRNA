#!/usr/bin/env bash

set -euo pipefail

# Edit these paths before running.
REF_DIR="/mnt/d/primates_bmc/genomes/"
QUERY_FILE="/mnt/d/primates_bmc/lncrna/gencode.v49.lncRNA_transcripts.fa"
OUTPUT_DIR="/mnt/d/primates_bmc/kmers_results/20"

BINARY="./build/lncrna_mers"

# Choose "ibf" or "hibf".
INDEX_METHOD="ibf"

THREADS=1
KMER_SIZE=20
SINGLE_RESULTS_WRITER=false

FRAGMENT_SIZE=80000
HIT_THRESHOLD=15
STORE_FRAGMENTS=false
STORE_INDEX=false
CLEANUP_INDEX=true
LOG_FILE="ibf_run.log"
OUTPUT_FILE="results.tsv"

IBF_HASH_FUNCTIONS=3
IBF_FPR=0.01

HIBF_HASH_FUNCTIONS=2
HIBF_MAXIMUM_FPR=0.05
HIBF_RELAXED_FPR=0.30
HIBF_THREADS=1
HIBF_SKETCH_BITS=12
HIBF_TMAX=0
HIBF_EMPTY_BIN_FRACTION=0.0
HIBF_ALPHA=1.2
HIBF_MAX_REARRANGEMENT_RATIO=0.5
HIBF_DISABLE_ESTIMATE_UNION=false
HIBF_DISABLE_REARRANGEMENT=false

if [[ "${REF_DIR}" == "/path/to/reference_genomes" ]] || \
   [[ "${QUERY_FILE}" == "/path/to/queries.fa" ]] || \
   [[ "${OUTPUT_DIR}" == "/path/to/output_dir" ]]; then
    echo "Set REF_DIR, QUERY_FILE, and OUTPUT_DIR at the top of $0 before running."
    exit 1
fi

if [[ ! -x "${BINARY}" ]]; then
    echo "Binary not found or not executable: ${BINARY}"
    echo "Build first, for example:"
    echo "  cmake -S . -B build -DCMAKE_BUILD_TYPE=Release"
    echo "  cmake --build build -j"
    exit 1
fi

if [[ ! -d "${REF_DIR}" ]]; then
    echo "Reference directory does not exist: ${REF_DIR}"
    exit 1
fi

if [[ ! -f "${QUERY_FILE}" ]]; then
    echo "Query FASTA does not exist: ${QUERY_FILE}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

CONFIG_PATH="${OUTPUT_DIR}/config_k15_per_reference.toml"

cat > "${CONFIG_PATH}" <<EOF
[general]
index_method = "${INDEX_METHOD}"
ref_dir = "${REF_DIR}"
query_file = "${QUERY_FILE}"

output_dir = "${OUTPUT_DIR}"
output_file = "${OUTPUT_FILE}"
log_file = "${LOG_FILE}"

fragment_size = ${FRAGMENT_SIZE}
kmer_size = ${KMER_SIZE}
hit_threshold = ${HIT_THRESHOLD}
threads = ${THREADS}

store_fragments = ${STORE_FRAGMENTS}
store_index = ${STORE_INDEX}
cleanup_index = ${CLEANUP_INDEX}
single_results_writer = ${SINGLE_RESULTS_WRITER}

[ibf]
hash_functions = ${IBF_HASH_FUNCTIONS}
fpr = ${IBF_FPR}

[hibf]
hash_functions = ${HIBF_HASH_FUNCTIONS}
maximum_fpr = ${HIBF_MAXIMUM_FPR}
relaxed_fpr = ${HIBF_RELAXED_FPR}
threads = ${HIBF_THREADS}
sketch_bits = ${HIBF_SKETCH_BITS}
tmax = ${HIBF_TMAX}
empty_bin_fraction = ${HIBF_EMPTY_BIN_FRACTION}
alpha = ${HIBF_ALPHA}
max_rearrangement_ratio = ${HIBF_MAX_REARRANGEMENT_RATIO}
disable_estimate_union = ${HIBF_DISABLE_ESTIMATE_UNION}
disable_rearrangement = ${HIBF_DISABLE_REARRANGEMENT}
EOF

echo "Running lncrna_mers with:"
echo "  ref_dir: ${REF_DIR}"
echo "  query_file: ${QUERY_FILE}"
echo "  output_dir: ${OUTPUT_DIR}"
echo "  index_method: ${INDEX_METHOD}"
echo "  kmer_size: ${KMER_SIZE}"
echo "  threads: ${THREADS}"
echo "  store_index: ${STORE_INDEX}"
echo "  single_results_writer: ${SINGLE_RESULTS_WRITER}"
echo

"${BINARY}" "${CONFIG_PATH}"

echo
echo "Finished."
echo "Per-reference result files are in: ${OUTPUT_DIR}"
echo "Unique k-mers are in: ${OUTPUT_DIR}/unique_mers"
