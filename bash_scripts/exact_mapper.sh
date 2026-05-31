#!/usr/bin/env bash

#./exact_mapper.sh genomes/Mmul10_GCF_003339765.1_Mmul_10_genomic.fna.gz 35 ZNF667_AS1_Mmul10_k35_kmer_locations.tsv
set -euo pipefail

LNCRNA_FASTA="lncrna/gencode.v49.lncRNA_transcripts.fa"
GENE_NAME="ZNF667-AS1"

usage() {
    cat <<EOF
Usage:
  $0 <reference.fna[.gz]> <k_size> [output.tsv]

Example:
  $0 genomes/Mmul10_GCF_003339765.1_Mmul_10_genomic.fna.gz 35 ZNF667_AS1_Mmul10_k35_kmer_locations.tsv

Output columns:
  transcript_name
  transcript_id
  k_size
  kmer
  transcript_kmer_count
  reference_hit_count
  reference_hits

The script also writes a transcript-level summary beside the detailed output:
  <output>.summary.tsv

Notes:
  - The lncRNA FASTA is fixed inside this script:
    ${LNCRNA_FASTA}
  - Only transcripts from gene ${GENE_NAME} are used.
  - N/n bases are removed from both transcript and reference sequences before
    k-mer generation/search, matching the FragIncRNA dna4 behavior used for
    k-mer sizes > 27.
  - reference_hit_count counts overlapping exact matches.
  - All reference records are concatenated into one N-removed sequence before
    searching, so no contig/chromosome assignment is made.
  - reference_hits is formatted as whole_reference:start;whole_reference:start
  - Positions are 1-based within the concatenated N-removed reference sequence,
    not original genome coordinates.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" || $# -lt 2 || $# -gt 3 ]]; then
    usage
    exit 0
fi

REFERENCE_FASTA="$1"
K_SIZE="$2"
OUTPUT="${3:-}"

if [[ ! -f "${LNCRNA_FASTA}" ]]; then
    echo "ERROR: lncRNA FASTA not found: ${LNCRNA_FASTA}" >&2
    exit 1
fi

if [[ ! -f "${REFERENCE_FASTA}" ]]; then
    echo "ERROR: reference FASTA not found: ${REFERENCE_FASTA}" >&2
    exit 1
fi

if ! [[ "${K_SIZE}" =~ ^[0-9]+$ ]] || [[ "${K_SIZE}" -lt 1 ]]; then
    echo "ERROR: k_size must be a positive integer." >&2
    exit 1
fi

if [[ -z "${OUTPUT}" ]]; then
    ref_base="$(basename "${REFERENCE_FASTA}")"
    ref_base="${ref_base%.gz}"
    ref_base="${ref_base%.fna}"
    ref_base="${ref_base%.fa}"
    ref_base="${ref_base%.fasta}"
    OUTPUT="${GENE_NAME}_${ref_base}_k${K_SIZE}_kmer_locations.tsv"
fi

if [[ "${OUTPUT}" == *.tsv ]]; then
    SUMMARY_OUTPUT="${OUTPUT%.tsv}.summary.tsv"
else
    SUMMARY_OUTPUT="${OUTPUT}.summary.tsv"
fi

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

echo "Searching ${REFERENCE_FASTA} as one concatenated N-removed reference; this can take time for large genomes..." >&2
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
