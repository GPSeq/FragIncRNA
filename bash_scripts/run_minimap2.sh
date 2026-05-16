#!/usr/bin/env bash
set -euo pipefail

MINIMAP2="./minimap2-2.30_x64-linux/minimap2"
INPUT_DIR="./data/all_primates"
OUTPUT_DIR="./output"
QUERY_FA="./gencode.v49.lncRNA_transcripts.fa"

THREADS_MINIMAP=30
THREADS_SAMTOOLS=8
KMER_SIZE=15

mkdir -p "${OUTPUT_DIR}"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${OUTPUT_DIR}/run_all_primates_${TIMESTAMP}.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log "Starting run"
log "Log file: ${LOG_FILE}"

if [[ ! -x "${MINIMAP2}" ]]; then
    log "ERROR: minimap2 not found or not executable: ${MINIMAP2}"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    log "ERROR: samtools not found in PATH"
    exit 1
fi

if [[ ! -f "${QUERY_FA}" ]]; then
    log "ERROR: query FASTA not found: ${QUERY_FA}"
    exit 1
fi

shopt -s nullglob
ref_files=(
    "${INPUT_DIR}"/*.fasta.gz
    "${INPUT_DIR}"/*.fa.gz
    "${INPUT_DIR}"/*.fna.gz
    "${INPUT_DIR}"/*.fas.gz
)

if [[ ${#ref_files[@]} -eq 0 ]]; then
    log "ERROR: no supported gzipped FASTA files found in ${INPUT_DIR}"
    exit 1
fi

count_fasta=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.fasta.gz" | wc -l)
count_fa=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.fa.gz" | wc -l)
count_fna=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.fna.gz" | wc -l)
count_fas=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.fas.gz" | wc -l)
count_all_gz=$(find "${INPUT_DIR}" -maxdepth 1 -type f -name "*.gz" | wc -l)

log "Found ${#ref_files[@]} supported reference file(s)"
log "All .gz files: ${count_all_gz}"
log "*.fasta.gz:  ${count_fasta}"
log "*.fa.gz:     ${count_fa}"
log "*.fna.gz:    ${count_fna}"
log "*.fas.gz:    ${count_fas}"
log "THREADS_MINIMAP=${THREADS_MINIMAP}"
log "THREADS_SAMTOOLS=${THREADS_SAMTOOLS}"
log "KMER_SIZE=${KMER_SIZE}"

for ref in "${ref_files[@]}"; do
    ref_name="$(basename "${ref}")"
    ref_base="${ref_name%.fasta.gz}"
    ref_base="${ref_base%.fa.gz}"
    ref_base="${ref_base%.fna.gz}"
    ref_base="${ref_base%.fas.gz}"

    prefix="${ref_base%%_*}"
    if [[ "${prefix}" == "${ref_base}" ]]; then
        prefix="${ref_base}"
    fi

    mmi="${OUTPUT_DIR}/${ref_base}.mmi"
    sam="${OUTPUT_DIR}/human_lncRNA_vs_${prefix}.sam"
    bam="${OUTPUT_DIR}/human_lncRNA_vs_${prefix}.sorted.bam"
    bai="${bam}.bai"

    log "============================================================"
    log "Processing reference: ${ref}"
    log "Reference basename : ${ref_base}"
    log "Output prefix      : ${prefix}"
    log "MMI                : ${mmi}"
    log "SAM                : ${sam}"
    log "BAM                : ${bam}"
    log "BAI                : ${bai}"
    log "============================================================"

    log "[1/6] Building minimap2 index"
    "${MINIMAP2}" -k "${KMER_SIZE}" -d "${mmi}" "${ref}"

    log "[2/6] Aligning query transcripts"
    "${MINIMAP2}" \
        -k "${KMER_SIZE}" \
        -ax splice:hq \
        --secondary=no \
        -t "${THREADS_MINIMAP}" \
        "${mmi}" \
        "${QUERY_FA}" \
        > "${sam}"

    log "[3/6] Converting SAM to sorted BAM"
    samtools view -@ "${THREADS_SAMTOOLS}" -bS "${sam}" | \
        samtools sort -@ "${THREADS_SAMTOOLS}" -o "${bam}"

    log "[4/6] Indexing BAM"
    samtools index "${bam}"

    log "[5/6] Removing temporary SAM"
    rm -f "${sam}"

    log "[6/6] Removing temporary minimap2 index"
    rm -f "${mmi}"

    log "Finished reference: ${ref_name}"
done

log "All references processed successfully"
