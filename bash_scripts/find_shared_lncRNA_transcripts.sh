#!/usr/bin/env bash
set -euo pipefail

STATUS_MATRIX="${1:-bam_comparison/qc/lncrna_transcript_status_matrix.tsv}"
OUT_DIR="${2:-bam_comparison/qc/shared_lncRNA_transcripts}"

mkdir -p "${OUT_DIR}"

STRICT_OUT="${OUT_DIR}/shared_lncRNA_transcripts_strict_all_genomes.tsv"
BASIC_OUT="${OUT_DIR}/shared_lncRNA_transcripts_basic_all_genomes.tsv"
PASS_COUNTS_OUT="${OUT_DIR}/lncRNA_transcript_pass_counts.tsv"
SUMMARY_OUT="${OUT_DIR}/shared_lncRNA_transcripts_summary.tsv"

if [[ ! -f "${STATUS_MATRIX}" ]]; then
    echo "ERROR: status matrix not found: ${STATUS_MATRIX}" >&2
    exit 1
fi

awk -F '\t' '
    BEGIN { OFS = FS }

    NR == 1 {
        print > strict_out
        print > basic_out
        print $0, "strict_pass_count", "basic_pass_count", "genome_count" > pass_counts_out
        genome_count = NF - 4
        next
    }

    {
        strict_count = 0
        basic_count = 0
        strict_all = 1
        basic_all = 1

        for (i = 5; i <= NF; i++) {
            if ($i == "PASS_STRICT") {
                strict_count++
                basic_count++
            } else {
                strict_all = 0
                if ($i == "PASS_BASIC") {
                    basic_count++
                } else {
                    basic_all = 0
                }
            }
        }

        print $0, strict_count, basic_count, genome_count > pass_counts_out

        total++
        if (strict_all) {
            strict_total++
            print > strict_out
        }
        if (basic_all) {
            basic_total++
            print > basic_out
        }
    }

    END {
        print "metric", "count", "percent_of_transcripts" > summary_out
        print "total_transcripts", total, "100.00" > summary_out
        print "strict_shared_all_genomes", strict_total + 0, sprintf("%.2f", 100 * (strict_total + 0) / total) > summary_out
        print "basic_shared_all_genomes", basic_total + 0, sprintf("%.2f", 100 * (basic_total + 0) / total) > summary_out
    }
' \
    strict_out="${STRICT_OUT}" \
    basic_out="${BASIC_OUT}" \
    pass_counts_out="${PASS_COUNTS_OUT}" \
    summary_out="${SUMMARY_OUT}" \
    "${STATUS_MATRIX}"

echo "Wrote ${STRICT_OUT}"
echo "Wrote ${BASIC_OUT}"
echo "Wrote ${PASS_COUNTS_OUT}"
echo "Wrote ${SUMMARY_OUT}"

#find_shared_lncRNA_transcripts.sh
#   It generated:

#   bam_comparison/qc/shared_lncRNA_transcripts/shared_lncRNA_transcripts_strict_all_genomes.tsv
#   bam_comparison/qc/shared_lncRNA_transcripts/shared_lncRNA_transcripts_basic_all_genomes.tsv
#   bam_comparison/qc/shared_lncRNA_transcripts/lncRNA_transcript_pass_counts.tsv
#   bam_comparison/qc/shared_lncRNA_transcripts/shared_lncRNA_transcripts_summary.tsv