#!/usr/bin/env bash
set -euo pipefail

BAM_LIST="${1:-bam_list.txt}"
OUT_DIR="${2:-bam_comparison/qc}"
MIN_COV="${MIN_COV:-80}"
MIN_ID="${MIN_ID:-80}"
MIN_MAPQ="${MIN_MAPQ:-10}"
STRICT_MIN_ID="${STRICT_MIN_ID:-90}"
STRICT_MIN_MAPQ="${STRICT_MIN_MAPQ:-30}"
LARGE_CONTIG_BP="${LARGE_CONTIG_BP:-10000000}"

mkdir -p "${OUT_DIR}"

SUMMARY="${OUT_DIR}/lncrna_alignment_sample_qc.tsv"
TMP_SUMMARY="${SUMMARY}.tmp"

if [[ ! -f "${BAM_LIST}" ]]; then
    echo "ERROR: BAM list not found: ${BAM_LIST}" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

{
printf '%s\n' \
    "sample	bam	total_transcripts	primary_mapped	primary_mapped_pct	unmapped	unmapped_pct	supplementary	supplementary_per_100_primary	mq0	mq0_pct	mapq_ge10_pct	mapq_ge30_pct	mean_mapq	mean_query_coverage_pct	mean_identity_pct	mean_softclip_pct	pass_cov80_id80_mapq10	pass_cov80_id80_mapq10_pct	pass_cov80_id90_mapq30	pass_cov80_id90_mapq30_pct	primary_with_sa_tag	primary_with_sa_tag_pct	large_contig_bp_threshold	large_contigs	mapped_on_large_contigs_pct"
while IFS= read -r bam || [[ -n "${bam}" ]]; do
    [[ -z "${bam}" ]] && continue
    [[ "${bam}" =~ ^[[:space:]]*# ]] && continue

    if [[ ! -f "${bam}" ]]; then
        echo "WARNING: skipping missing BAM: ${bam}" >&2
        continue
    fi

    sample="$(basename "${bam}" .sorted.bam)"
    echo "QC ${sample}" >&2

    samtools view -h "${bam}" | awk \
        -v sample="${sample}" \
        -v bam="${bam}" \
        -v min_cov="${MIN_COV}" \
        -v min_id="${MIN_ID}" \
        -v min_mapq="${MIN_MAPQ}" \
        -v strict_min_id="${STRICT_MIN_ID}" \
        -v strict_min_mapq="${STRICT_MIN_MAPQ}" \
        -v large_bp="${LARGE_CONTIG_BP}" '
        function has_flag(flag, bit) {
            return int(flag / bit) % 2
        }

        function query_len_from_name(name, seq, parts, n) {
            n = split(name, parts, "|")
            if (n > 1 && parts[n - 1] ~ /^[0-9]+$/) {
                return parts[n - 1] + 0
            }
            if (seq != "*") {
                return length(seq)
            }
            return 0
        }

        function parse_cigar(cigar, token, n, op) {
            cigar_query = 0
            cigar_aligned_query = 0
            cigar_softclip = 0
            while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
                token = substr(cigar, RSTART, RLENGTH)
                n = substr(token, 1, length(token) - 1) + 0
                op = substr(token, length(token), 1)
                if (op ~ /[MIS=X]/) {
                    cigar_query += n
                }
                if (op ~ /[MI=X]/) {
                    cigar_aligned_query += n
                }
                if (op == "S") {
                    cigar_softclip += n
                }
                cigar = substr(cigar, RSTART + RLENGTH)
            }
        }

        function pct(n, d) {
            return d ? (100.0 * n / d) : 0
        }

        BEGIN {
            FS = OFS = "\t"
        }

        /^@SQ/ {
            sn = ""
            ln = 0
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^SN:/) {
                    sn = substr($i, 4)
                } else if ($i ~ /^LN:/) {
                    ln = substr($i, 4) + 0
                }
            }
            if (sn != "" && ln >= large_bp) {
                large_contig[sn] = 1
                large_contig_count++
            }
            next
        }

        /^@/ {
            next
        }

        {
            flag = $2 + 0
            is_unmapped = has_flag(flag, 4)
            is_secondary = has_flag(flag, 256)
            is_supplementary = has_flag(flag, 2048)

            if (is_supplementary) {
                supplementary++
            }

            if (is_secondary || is_supplementary) {
                next
            }

            total_primary++
            if (is_unmapped) {
                unmapped++
                next
            }

            primary_mapped++
            mapq = $5 + 0
            mapq_sum += mapq
            if (mapq == 0) {
                mq0++
            }
            if (mapq >= 10) {
                mapq_ge10++
            }
            if (mapq >= 30) {
                mapq_ge30++
            }
            if (($3 in large_contig)) {
                mapped_large++
            }

            qlen = query_len_from_name($1, $10)
            parse_cigar($6)
            cov = qlen ? pct(cigar_aligned_query, qlen) : 0
            soft = qlen ? pct(cigar_softclip, qlen) : 0
            cov_sum += cov
            soft_sum += soft

            identity = ""
            has_sa = 0
            for (i = 12; i <= NF; i++) {
                if ($i ~ /^de:f:/) {
                    split($i, de_parts, ":")
                    identity = 100.0 * (1.0 - de_parts[3])
                } else if ($i ~ /^SA:Z:/) {
                    has_sa = 1
                }
            }
            if (identity != "") {
                identity_sum += identity
                identity_count++
            }
            if (has_sa) {
                primary_with_sa++
            }

            if (cov >= min_cov && identity != "" && identity >= min_id && mapq >= min_mapq) {
                pass_basic++
            }
            if (cov >= min_cov && identity != "" && identity >= strict_min_id && mapq >= strict_min_mapq) {
                pass_strict++
            }
        }

        END {
            print sample, bam, total_primary, primary_mapped, sprintf("%.2f", pct(primary_mapped, total_primary)), \
                unmapped, sprintf("%.2f", pct(unmapped, total_primary)), \
                supplementary, sprintf("%.2f", 100.0 * supplementary / (primary_mapped ? primary_mapped : 1)), \
                mq0, sprintf("%.2f", pct(mq0, primary_mapped)), \
                sprintf("%.2f", pct(mapq_ge10, primary_mapped)), sprintf("%.2f", pct(mapq_ge30, primary_mapped)), \
                sprintf("%.2f", primary_mapped ? mapq_sum / primary_mapped : 0), \
                sprintf("%.2f", primary_mapped ? cov_sum / primary_mapped : 0), \
                sprintf("%.2f", identity_count ? identity_sum / identity_count : 0), \
                sprintf("%.2f", primary_mapped ? soft_sum / primary_mapped : 0), \
                pass_basic, sprintf("%.2f", pct(pass_basic, total_primary)), \
                pass_strict, sprintf("%.2f", pct(pass_strict, total_primary)), \
                primary_with_sa, sprintf("%.2f", pct(primary_with_sa, primary_mapped)), \
                large_bp, large_contig_count, sprintf("%.2f", pct(mapped_large, primary_mapped))
        }
    '
done < "${BAM_LIST}"
} > "${TMP_SUMMARY}"

mv "${TMP_SUMMARY}" "${SUMMARY}"
echo "Wrote ${SUMMARY}"
#which gave the results file lncrna_alignment_sample_qc.tsv 