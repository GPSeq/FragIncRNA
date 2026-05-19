#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Summarize k-mer matches per lncRNA transcript across genome result files.

For each transcript in the input FASTA, this script reports:
  - transcript_name
  - gene_name
  - gene_id
  - transcript_id
  - transcript_length
  - total_kmers, computed as transcript_length - k + 1
  - for each genome: matched unique k-mer count and matched/total ratio

The script uses *_kmers.tsv files with columns:
  query_index    kmer_indices    kmer_counts

query_index is interpreted as the 0-based transcript index in the FASTA.
matched unique k-mers are counted from the slash-separated kmer_indices field,
but only when the corresponding slash-separated kmer_counts value exists and is
greater than zero. kmer_counts is otherwise not summed, because it is the number
of reference hits per k-mer, not the number of distinct query k-mers detected.

Usage:
  summarize_kmer_matches_by_transcript.sh [options]

Options:
  -f FASTA       Input lncRNA FASTA
                 default: /mnt/d/primates_bmc/lncrna/gencode.v49.lncRNA_transcripts.fa
  -k INT         k-mer size
                 default: 15
  -i DIR         Directory containing *_kmers.tsv files
                 default: /mnt/d/primates_bmc/kmers_results/output_kmers
  -o FILE        Output TSV file
                 default: DIR/lncRNA_kmer_match_summary_by_transcript.tsv
  -w DIR         Directory for per-genome intermediate summaries
                 default: output file basename + _per_genome
  -h             Show this help

Example:
  bash summarize_kmer_matches_by_transcript.sh \
    -f /mnt/d/primates_bmc/lncrna/gencode.v49.lncRNA_transcripts.fa \
    -i /mnt/d/primates_bmc/kmers_results/output_kmers \
    -o /mnt/d/primates_bmc/kmers_results/output_kmers/lncRNA_kmer_match_summary_by_transcript.tsv
USAGE
}

FASTA="/mnt/d/primates_bmc/lncrna/gencode.v49.lncRNA_transcripts.fa"
KMER_SIZE=15
KMER_DIR="/mnt/d/primates_bmc/kmers_results/output_kmers"
OUT_FILE=""
WORK_DIR=""

while getopts ":f:k:i:o:w:h" opt; do
    case "${opt}" in
        f) FASTA="${OPTARG}" ;;
        k) KMER_SIZE="${OPTARG}" ;;
        i) KMER_DIR="${OPTARG}" ;;
        o) OUT_FILE="${OPTARG}" ;;
        w) WORK_DIR="${OPTARG}" ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: missing argument for -${OPTARG}" >&2; usage >&2; exit 2 ;;
        \?) echo "ERROR: unknown option -${OPTARG}" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "${OUT_FILE}" ]]; then
    OUT_FILE="${KMER_DIR}/lncRNA_kmer_match_summary_by_transcript.tsv"
fi

if [[ -z "${WORK_DIR}" ]]; then
    OUT_BASE="$(basename "${OUT_FILE}")"
    OUT_BASE="${OUT_BASE%.tsv}"
    WORK_DIR="$(dirname "${OUT_FILE}")/${OUT_BASE}_per_genome"
fi

if [[ ! -s "${FASTA}" ]]; then
    echo "ERROR: FASTA not found or empty: ${FASTA}" >&2
    exit 1
fi

if [[ ! -d "${KMER_DIR}" ]]; then
    echo "ERROR: k-mer result directory not found: ${KMER_DIR}" >&2
    exit 1
fi

if ! [[ "${KMER_SIZE}" =~ ^[0-9]+$ ]] || [[ "${KMER_SIZE}" -lt 1 ]]; then
    echo "ERROR: k-mer size must be a positive integer: ${KMER_SIZE}" >&2
    exit 1
fi

shopt -s nullglob
KMER_FILES=("${KMER_DIR}"/*_kmers.tsv)
shopt -u nullglob

if [[ ${#KMER_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no *_kmers.tsv files found in ${KMER_DIR}" >&2
    exit 1
fi

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/kmer_summary.XXXXXX")"
cleanup() {
    rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

META_TSV="${TMP_DIR}/transcript_metadata.tsv"
META_NO_INDEX="${TMP_DIR}/transcript_metadata_no_query_index.tsv"
PASTE_LIST="${TMP_DIR}/paste_files.txt"
HEADER_FILE="${TMP_DIR}/header.tsv"
BODY_FILE="${TMP_DIR}/body.tsv"
mkdir -p "${WORK_DIR}" "$(dirname "${OUT_FILE}")"

echo "Parsing FASTA metadata: ${FASTA}" >&2
grep "^>" "${FASTA}" | awk -v k="${KMER_SIZE}" '
    BEGIN {
        FS = OFS = "\t"
        print "query_index", "transcript_name", "gene_name", "gene_id", \
              "transcript_id", "transcript_length", "total_kmers"
        query_idx = -1
    }

    {
        query_idx++
        header = substr($0, 2)
        sub(/[[:space:]].*$/, "", header)
        n = split(header, parts, "|")
        len = 0
        for (i = n; i >= 1; i--) {
            if (parts[i] ~ /^[0-9]+$/) {
                len = parts[i] + 0
                break
            }
        }
        total = (len >= k) ? len - k + 1 : 0
        transcript_id = (n >= 1) ? parts[1] : header
        gene_id = (n >= 2) ? parts[2] : ""
        transcript_name = (n >= 5) ? parts[5] : transcript_id
        gene_name = (n >= 6) ? parts[6] : gene_id
        print query_idx, transcript_name, gene_name, gene_id, transcript_id, len, total
    }
' > "${META_TSV}"

awk 'BEGIN{FS=OFS="\t"} NR==1{print "transcript_name","gene_name","gene_id","transcript_id","transcript_length","total_kmers"; next} {print $2,$3,$4,$5,$6,$7}' \
    "${META_TSV}" > "${META_NO_INDEX}"

total_files=${#KMER_FILES[@]}
printf '%s\n' "${META_NO_INDEX}" > "${PASTE_LIST}"

echo "Processing ${total_files} k-mer result file(s) one by one." >&2
file_number=0
for kmer_file in "${KMER_FILES[@]}"; do
    file_number=$((file_number + 1))
    genome="$(basename "${kmer_file}")"
    genome="${genome%_kmers.tsv}"
    per_genome="${WORK_DIR}/${genome}.summary.tsv"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${file_number}/${total_files}] ${genome}" >&2

    awk \
        -v genome="${genome}" \
        -v progress_label="[${file_number}/${total_files}] ${genome}" \
        -v progress_every=25000 '
        BEGIN {
            FS = OFS = "\t"
        }

        NR == FNR {
            if (FNR == 1) {
                next
            }
            total[$1] = $7 + 0
            order[++n] = $1
            next
        }

        FNR == 1 {
            next
        }

        {
            query = $1
            indices = $2
            counts = $3
            if (indices == "" || indices == "." || indices == "NA" || \
                counts == "" || counts == "." || counts == "NA") {
                matched[query] = 0
            } else {
                n_indices = split(indices, index_parts, "/")
                n_counts = split(counts, count_parts, "/")
                valid = 0
                for (j = 1; j <= n_indices && j <= n_counts; j++) {
                    if (count_parts[j] + 0 > 0) {
                        valid++
                    }
                }
                matched[query] = valid
            }
            processed = FNR - 1
            if (processed % progress_every == 0) {
                printf("  %s: processed %d/%d transcript rows\n", progress_label, processed, n) > "/dev/stderr"
            }
        }

        END {
            printf("  %s: processed %d/%d transcript rows\n", progress_label, FNR - 1, n) > "/dev/stderr"
            print genome "_matched_kmers", genome "_matched_kmer_ratio"
            for (i = 1; i <= n; i++) {
                query = order[i]
                count = (query in matched) ? matched[query] : 0
                ratio = (total[query] > 0) ? count / total[query] : 0
                print count, sprintf("%.6f", ratio)
            }
        }
    ' "${META_TSV}" "${kmer_file}" > "${per_genome}"

    printf '%s\n' "${per_genome}" >> "${PASTE_LIST}"
done

echo "Combining per-genome summaries: ${OUT_FILE}" >&2
paste $(cat "${PASTE_LIST}") > "${TMP_DIR}/combined_with_header.tsv"
mv "${TMP_DIR}/combined_with_header.tsv" "${OUT_FILE}"

echo "Done." >&2
echo "Output: ${OUT_FILE}" >&2
echo "Per-genome intermediate summaries: ${WORK_DIR}" >&2
