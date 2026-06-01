#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<EOF
Usage:
  $0 <genomes_folder> <lncrna_fasta> <k_size> [gene_name] [output_folder]

Example:
  $0 genomes lncrna/gencode.v49.lncRNA_transcripts.fa 35
  $0 genomes lncrna/gencode.v49.lncRNA_transcripts.fa 35 ZNF667-AS1 results

Output:
  For each genome FASTA in <genomes_folder>, writes:

    <output_folder>/<gene_name>_<genome_basename>_k<k_size>_kmer_locations.tsv
    <output_folder>/<gene_name>_<genome_basename>_k<k_size>_kmer_locations.summary.tsv

Input genome extensions searched:
  *.fna
  *.fa
  *.fasta
  *.fna.gz
  *.fa.gz
  *.fasta.gz

Notes:
  - Only transcripts from the selected gene are used.
  - N/n bases are removed from both transcript and reference sequences.
  - reference_hit_count counts overlapping exact matches.
  - All reference records are concatenated into one N-removed sequence.
  - Positions are 1-based within the concatenated N-removed reference sequence.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 3 || $# -gt 5 ]]; then
    usage
    exit 0
fi

GENOMES_DIR="$1"
LNCRNA_FASTA="$2"
K_SIZE="$3"
GENE_NAME="${4:-ZNF667-AS1}"
OUTPUT_DIR="${5:-.}"

if [[ ! -d "${GENOMES_DIR}" ]]; then
    echo "ERROR: genomes folder not found: ${GENOMES_DIR}" >&2
    exit 1
fi

if [[ ! -f "${LNCRNA_FASTA}" ]]; then
    echo "ERROR: lncRNA FASTA not found: ${LNCRNA_FASTA}" >&2
    exit 1
fi

if ! [[ "${K_SIZE}" =~ ^[0-9]+$ ]] || [[ "${K_SIZE}" -lt 1 ]]; then
    echo "ERROR: k_size must be a positive integer." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

transcripts_tsv="${tmpdir}/${GENE_NAME}_transcripts.tsv"
kmers_tsv="${tmpdir}/${GENE_NAME}_kmers_k${K_SIZE}.tsv"

echo "Extracting ${GENE_NAME} transcripts from ${LNCRNA_FASTA}..." >&2

zcat -f "${LNCRNA_FASTA}" | awk -F'|' -v gene="${GENE_NAME}" '
    /^>/ {
        if (keep && seq != "") {
            print transcript_name "\t" transcript_id "\t" seq
        }

        header = substr($0, 2)
        n = split(header, f, "|")

        transcript_id = f[1]
        transcript_name = f[5]
        gene_name = f[6]

        keep = (gene_name == gene)
        seq = ""
        next
    }

    keep {
        line = toupper($0)
        gsub(/[[:space:]]/, "", line)
        gsub(/N/, "", line)
        seq = seq line
    }

    END {
        if (keep && seq != "") {
            print transcript_name "\t" transcript_id "\t" seq
        }
    }
' > "${transcripts_tsv}"

if [[ ! -s "${transcripts_tsv}" ]]; then
    echo "ERROR: no ${GENE_NAME} transcripts found in ${LNCRNA_FASTA}" >&2
    exit 1
fi

echo "Generating unique ${K_SIZE}-mers per transcript..." >&2

awk -F'\t' -v OFS='\t' -v k="${K_SIZE}" '
    {
        transcript_name = $1
        transcript_id = $2
        seq = $3

        if (length(seq) < k) {
            next
        }

        for (i = 1; i <= length(seq) - k + 1; i++) {
            kmer = substr(seq, i, k)
            if (kmer ~ /^[ACGT]+$/) {
                key = transcript_name SUBSEP transcript_id SUBSEP kmer
                count[key]++
            }
        }
    }

    END {
        for (key in count) {
            split(key, f, SUBSEP)
            print f[1], f[2], k, f[3], count[key]
        }
    }
' "${transcripts_tsv}" | sort -k1,1 -k4,4 > "${kmers_tsv}"

if [[ ! -s "${kmers_tsv}" ]]; then
    echo "ERROR: no valid ${K_SIZE}-mers generated for ${GENE_NAME}" >&2
    exit 1
fi

shopt -s nullglob

genome_files=(
    "${GENOMES_DIR}"/*.fna
    "${GENOMES_DIR}"/*.fa
    "${GENOMES_DIR}"/*.fasta
    "${GENOMES_DIR}"/*.fna.gz
    "${GENOMES_DIR}"/*.fa.gz
    "${GENOMES_DIR}"/*.fasta.gz
)

if [[ "${#genome_files[@]}" -eq 0 ]]; then
    echo "ERROR: no genome FASTA files found in ${GENOMES_DIR}" >&2
    exit 1
fi

echo "Found ${#genome_files[@]} genome files." >&2

for REFERENCE_FASTA in "${genome_files[@]}"; do
    ref_base="$(basename "${REFERENCE_FASTA}")"

    ref_base="${ref_base%.gz}"
    ref_base="${ref_base%.fna}"
    ref_base="${ref_base%.fa}"
    ref_base="${ref_base%.fasta}"

    safe_gene="${GENE_NAME//-/_}"

    OUTPUT="${OUTPUT_DIR}/${safe_gene}_${ref_base}_k${K_SIZE}_kmer_locations.tsv"
    SUMMARY_OUTPUT="${OUTPUT%.tsv}.summary.tsv"

    echo "============================================================" >&2
    echo "Processing genome: ${REFERENCE_FASTA}" >&2
    echo "Output: ${OUTPUT}" >&2
    echo "Summary: ${SUMMARY_OUTPUT}" >&2

    awk -F'\t' -v OFS='\t' -v summary_out="${SUMMARY_OUTPUT}" '
        NR == FNR {
            transcript_name[++n_kmers] = $1
            transcript_id[n_kmers] = $2
            k_size[n_kmers] = $3
            kmer[n_kmers] = $4
            transcript_count[n_kmers] = $5

            k = k_size[n_kmers]

            kmer_rows[kmer[n_kmers]] = kmer_rows[kmer[n_kmers]] \
                (kmer_rows[kmer[n_kmers]] == "" ? "" : ",") n_kmers

            tx_key = transcript_name[n_kmers] SUBSEP transcript_id[n_kmers]

            if (!(tx_key in tx_seen)) {
                tx_seen[tx_key] = ++n_tx
                tx_name[n_tx] = transcript_name[n_kmers]
                tx_id[n_tx] = transcript_id[n_kmers]
            }

            tx_total_windows[tx_key] += transcript_count[n_kmers]
            tx_unique_kmers[tx_key]++

            next
        }

        function search_reference(    pos, ref_kmer, row_count, rows, j, idx) {
            if (seq == "") {
                return
            }

            for (pos = 1; pos <= length(seq) - k + 1; pos++) {
                ref_kmer = substr(seq, pos, k)

                if (ref_kmer in kmer_rows) {
                    row_count = split(kmer_rows[ref_kmer], rows, ",")

                    for (j = 1; j <= row_count; j++) {
                        idx = rows[j]
                        ref_count[idx]++
                        locations[idx] = locations[idx] \
                            (locations[idx] == "" ? "" : ";") "whole_reference:" pos
                    }
                }
            }
        }

        FNR == 1 {
            seq = ""
        }

        /^>/ {
            next
        }

        {
            line = toupper($0)
            gsub(/[[:space:]]/, "", line)
            gsub(/N/, "", line)
            seq = seq line
        }

        END {
            search_reference()

            print "transcript_name", "transcript_id", "k_size", \
                "total_kmer_windows", "matched_kmer_windows", "matched_kmer_ratio", \
                "unique_kmers", "unique_kmers_found", "unique_kmer_ratio", \
                "reference_hit_sum" > summary_out

            print "transcript_name", "transcript_id", "k_size", "kmer", \
                "transcript_kmer_count", "reference_hit_count", "reference_hits"

            for (i = 1; i <= n_kmers; i++) {
                tx_key = transcript_name[i] SUBSEP transcript_id[i]

                if (ref_count[i] > 0) {
                    tx_matched_windows[tx_key] += transcript_count[i]
                    tx_unique_found[tx_key]++
                }

                tx_ref_hit_sum[tx_key] += ref_count[i]

                print transcript_name[i], transcript_id[i], k_size[i], kmer[i], \
                    transcript_count[i], ref_count[i] + 0, \
                    (locations[i] == "" ? "." : locations[i])
            }

            for (j = 1; j <= n_tx; j++) {
                tx_key = tx_name[j] SUBSEP tx_id[j]

                total_windows = tx_total_windows[tx_key] + 0
                matched_windows = tx_matched_windows[tx_key] + 0
                unique_total = tx_unique_kmers[tx_key] + 0
                unique_found = tx_unique_found[tx_key] + 0

                print tx_name[j], tx_id[j], k, \
                    total_windows, matched_windows, \
                    (total_windows ? matched_windows / total_windows : 0), \
                    unique_total, unique_found, \
                    (unique_total ? unique_found / unique_total : 0), \
                    tx_ref_hit_sum[tx_key] + 0 > summary_out
            }
        }
    ' "${kmers_tsv}" <(zcat -f "${REFERENCE_FASTA}") > "${OUTPUT}"

    echo "Wrote ${OUTPUT}" >&2
    echo "Wrote ${SUMMARY_OUTPUT}" >&2
done

echo "============================================================" >&2
echo "Done." >&2
