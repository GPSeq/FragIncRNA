#!/usr/bin/env bash

set -euo pipefail

# Edit these paths before running.
REF_DIR="/mnt/d/primates_bmc/genomes/"
QUERY_FILE="/mnt/d/primates_bmc/lncrna/gencode.v49.lncRNA_transcripts.fa"
OUTPUT_DIR="/mnt/d/primates_bmc/output_kmers"

BINARY="./build/lncrna_mers"

THREADS=1
KMER_SIZE=15
STORE_IBF=false
SINGLE_RESULTS_WRITER=false

FRAGMENT_SIZE=80000
HASH_FUNCTIONS=3
FPR=0.01
HIT_THRESHOLD=15
STORE_FRAGMENTS=false
CLEANUP_IBF=true
LOG_FILE="ibf_run.log"
OUTPUT_FILE="results.tsv"

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
ref_dir = "${REF_DIR}"
query_file = "${QUERY_FILE}"

output_dir = "${OUTPUT_DIR}"
output_file = "${OUTPUT_FILE}"
log_file = "${LOG_FILE}"

fragment_size = ${FRAGMENT_SIZE}
kmer_size = ${KMER_SIZE}
hash_functions = ${HASH_FUNCTIONS}
fpr = ${FPR}
hit_threshold = ${HIT_THRESHOLD}
threads = ${THREADS}

store_fragments = ${STORE_FRAGMENTS}
store_ibf = ${STORE_IBF}
cleanup_ibf = ${CLEANUP_IBF}
single_results_writer = ${SINGLE_RESULTS_WRITER}
EOF

echo "Running lncrna_mers with:"
echo "  ref_dir: ${REF_DIR}"
echo "  query_file: ${QUERY_FILE}"
echo "  output_dir: ${OUTPUT_DIR}"
echo "  kmer_size: ${KMER_SIZE}"
echo "  threads: ${THREADS}"
echo "  store_ibf: ${STORE_IBF}"
echo "  single_results_writer: ${SINGLE_RESULTS_WRITER}"
echo

"${BINARY}" "${CONFIG_PATH}"

echo
echo "Finished."
echo "Per-reference result files are in: ${OUTPUT_DIR}"
echo "Unique k-mers are in: ${OUTPUT_DIR}/unique_mers"
