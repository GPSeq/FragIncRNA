#!/usr/bin/env bash
set -euo pipefail

BAM_LIST="${1:-bam_list.txt}"
OUT_DIR="${2:-bam_comparison/qc}"
MIN_COV="${MIN_COV:-80}"
MIN_ID="${MIN_ID:-80}"
MIN_MAPQ="${MIN_MAPQ:-10}"
STRICT_MIN_ID="${STRICT_MIN_ID:-90}"
STRICT_MIN_MAPQ="${STRICT_MIN_MAPQ:-30}"
EXPORT_SUPPLEMENTARY="${EXPORT_SUPPLEMENTARY:-0}"

mkdir -p "${OUT_DIR}"

PRIMARY_TSV="${OUT_DIR}/lncrna_transcript_alignment_qc.tsv"
STATUS_MATRIX="${OUT_DIR}/lncrna_transcript_status_matrix.tsv"
COVERAGE_MATRIX="${OUT_DIR}/lncrna_transcript_coverage_matrix.tsv"
IDENTITY_MATRIX="${OUT_DIR}/lncrna_transcript_identity_matrix.tsv"
SUPP_TSV="${OUT_DIR}/lncrna_supplementary_alignment_segments.tsv"

TMP_PRIMARY="${PRIMARY_TSV}.tmp"
TMP_STATUS="${STATUS_MATRIX}.tmp"
TMP_COVERAGE="${COVERAGE_MATRIX}.tmp"
TMP_IDENTITY="${IDENTITY_MATRIX}.tmp"
TMP_SUPP="${SUPP_TSV}.tmp"

if [[ ! -f "${BAM_LIST}" ]]; then
    echo "ERROR: BAM list not found: ${BAM_LIST}" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

sample_order=""
while IFS= read -r bam || [[ -n "${bam}" ]]; do
    [[ -z "${bam}" ]] && continue
    [[ "${bam}" =~ ^[[:space:]]*# ]] && continue
    sample="$(basename "${bam}" .sorted.bam)"
    sample_order="${sample_order}${sample_order:+,}${sample}"
done < "${BAM_LIST}"

write_header() {
    printf '%s\n' \
        "sample	qname	transcript_id	gene_id	transcript_name	gene_name	query_length	flag	alignment_type	mapped	reference	position	strand	mapq	cigar	aligned_query_bases	query_coverage_pct	softclip_bases	softclip_pct	identity_pct	divergence_pct	nm	as_score	has_sa_tag	sa_tag	pass_basic	pass_strict	fail_reasons"
}

{
write_header
if [[ "${EXPORT_SUPPLEMENTARY}" == "1" ]]; then
    write_header > "${TMP_SUPP}"
fi

while IFS= read -r bam || [[ -n "${bam}" ]]; do
    [[ -z "${bam}" ]] && continue
    [[ "${bam}" =~ ^[[:space:]]*# ]] && continue

    if [[ ! -f "${bam}" ]]; then
        echo "WARNING: skipping missing BAM: ${bam}" >&2
        continue
    fi

    sample="$(basename "${bam}" .sorted.bam)"
    echo "Exporting ${sample}" >&2

    samtools view "${bam}" | awk \
        -v sample="${sample}" \
        -v min_cov="${MIN_COV}" \
        -v min_id="${MIN_ID}" \
        -v min_mapq="${MIN_MAPQ}" \
        -v strict_min_id="${STRICT_MIN_ID}" \
        -v strict_min_mapq="${STRICT_MIN_MAPQ}" \
        -v export_supp="${EXPORT_SUPPLEMENTARY}" \
        -v supp_out="${TMP_SUPP}" '
        function has_flag(flag, bit) {
            return int(flag / bit) % 2
        }

        function parse_name(name, parts, n, i) {
            transcript_id = gene_id = transcript_name = gene_name = ""
            query_length = 0
            n = split(name, parts, "|")
            if (n >= 1) transcript_id = parts[1]
            if (n >= 2) gene_id = parts[2]
            if (n >= 5) transcript_name = parts[5]
            if (n >= 6) gene_name = parts[6]
            for (i = n; i >= 1; i--) {
                if (parts[i] ~ /^[0-9]+$/) {
                    query_length = parts[i] + 0
                    break
                }
            }
        }

        function parse_cigar(cigar, token, n, op) {
            aligned_query = 0
            softclip = 0
            if (cigar == "*") {
                return
            }
            while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
                token = substr(cigar, RSTART, RLENGTH)
                n = substr(token, 1, length(token) - 1) + 0
                op = substr(token, length(token), 1)
                if (op ~ /[MI=X]/) {
                    aligned_query += n
                }
                if (op == "S") {
                    softclip += n
                }
                cigar = substr(cigar, RSTART + RLENGTH)
            }
        }

        function pct(n, d) {
            return d ? (100.0 * n / d) : 0
        }

        function tag_value(prefix, value, i) {
            value = ""
            for (i = 12; i <= NF; i++) {
                if (index($i, prefix) == 1) {
                    value = substr($i, length(prefix) + 1)
                    break
                }
            }
            return value
        }

        function make_fail_reasons(mapped, cov, identity, mapq, reasons) {
            reasons = ""
            if (!mapped) {
                return "unmapped"
            }
            if (cov < min_cov) {
                reasons = reasons "low_coverage,"
            }
            if (identity == "" || identity < min_id) {
                reasons = reasons "low_identity,"
            }
            if (mapq < min_mapq) {
                reasons = reasons "low_mapq,"
            }
            sub(/,$/, "", reasons)
            return reasons
        }

        function emit(out_file) {
            print sample, $1, transcript_id, gene_id, transcript_name, gene_name, query_length, \
                flag, alignment_type, mapped, ref, pos, strand, mapq, cigar, aligned_query, \
                sprintf("%.2f", coverage), softclip, sprintf("%.2f", softclip_pct), \
                identity_out, divergence_out, nm, as_score, has_sa, sa_tag, \
                pass_basic, pass_strict, fail_reasons >> out_file
        }

        BEGIN {
            FS = OFS = "\t"
            primary_out = "/dev/stdout"
        }

        {
            flag = $2 + 0
            secondary = has_flag(flag, 256)
            supplementary = has_flag(flag, 2048)
            if (secondary) {
                next
            }

            parse_name($1)

            mapped = !has_flag(flag, 4)
            alignment_type = supplementary ? "supplementary" : "primary"
            ref = mapped ? $3 : ""
            pos = mapped ? $4 : ""
            strand = has_flag(flag, 16) ? "-" : "+"
            mapq = mapped ? ($5 + 0) : 0
            cigar = mapped ? $6 : ""
            aligned_query = 0
            softclip = 0
            coverage = 0
            softclip_pct = 0
            identity = ""
            divergence = ""
            identity_out = ""
            divergence_out = ""
            nm = tag_value("NM:i:")
            as_score = tag_value("AS:i:")
            sa_tag = tag_value("SA:Z:")
            has_sa = sa_tag == "" ? 0 : 1

            if (mapped) {
                parse_cigar(cigar)
                coverage = pct(aligned_query, query_length)
                softclip_pct = pct(softclip, query_length)
                de = tag_value("de:f:")
                if (de != "") {
                    divergence = 100.0 * de
                    identity = 100.0 * (1.0 - de)
                    divergence_out = sprintf("%.2f", divergence)
                    identity_out = sprintf("%.2f", identity)
                }
            }

            pass_basic = mapped && coverage >= min_cov && identity != "" && identity >= min_id && mapq >= min_mapq ? 1 : 0
            pass_strict = mapped && coverage >= min_cov && identity != "" && identity >= strict_min_id && mapq >= strict_min_mapq ? 1 : 0
            fail_reasons = make_fail_reasons(mapped, coverage, identity, mapq)

            if (supplementary) {
                if (export_supp == "1") {
                    emit(supp_out)
                }
                next
            }
            emit(primary_out)
        }
    '
done < "${BAM_LIST}"
} > "${TMP_PRIMARY}"

mv "${TMP_PRIMARY}" "${PRIMARY_TSV}"
if [[ "${EXPORT_SUPPLEMENTARY}" == "1" ]]; then
    mv "${TMP_SUPP}" "${SUPP_TSV}"
fi

awk -F '\t' -v samples="${sample_order}" '
    BEGIN {
        OFS = "\t"
        n = split(samples, sample_names, ",")
        header = "transcript_id\tgene_id\ttranscript_name\tgene_name"
        for (i = 1; i <= n; i++) {
            header = header OFS sample_names[i]
        }
        print header
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) {
            col[$i] = i
        }
        next
    }

    {
        sample = $col["sample"]
        transcript = $col["transcript_id"]
        if (!(transcript in seen)) {
            seen[transcript] = 1
            order[++count] = transcript
            gene_id[transcript] = $col["gene_id"]
            transcript_name[transcript] = $col["transcript_name"]
            gene_name[transcript] = $col["gene_name"]
        }
        if ($col["mapped"] == 0) {
            value = "UNMAPPED"
        } else if ($col["pass_strict"] == 1) {
            value = "PASS_STRICT"
        } else if ($col["pass_basic"] == 1) {
            value = "PASS_BASIC"
        } else {
            value = "LOW_QC"
        }
        status[transcript, sample] = value
    }

    END {
        for (i = 1; i <= count; i++) {
            transcript = order[i]
            line = transcript OFS gene_id[transcript] OFS transcript_name[transcript] OFS gene_name[transcript]
            for (j = 1; j <= n; j++) {
                sample = sample_names[j]
                line = line OFS ((transcript, sample) in status ? status[transcript, sample] : "MISSING")
            }
            print line
        }
    }
' "${PRIMARY_TSV}" > "${TMP_STATUS}"

awk -F '\t' -v samples="${sample_order}" -v value_col="query_coverage_pct" '
    BEGIN {
        OFS = "\t"
        n = split(samples, sample_names, ",")
        header = "transcript_id\tgene_id\ttranscript_name\tgene_name"
        for (i = 1; i <= n; i++) header = header OFS sample_names[i]
        print header
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        next
    }

    {
        sample = $col["sample"]
        transcript = $col["transcript_id"]
        if (!(transcript in seen)) {
            seen[transcript] = 1
            order[++count] = transcript
            gene_id[transcript] = $col["gene_id"]
            transcript_name[transcript] = $col["transcript_name"]
            gene_name[transcript] = $col["gene_name"]
        }
        value[transcript, sample] = $col[value_col]
    }

    END {
        for (i = 1; i <= count; i++) {
            transcript = order[i]
            line = transcript OFS gene_id[transcript] OFS transcript_name[transcript] OFS gene_name[transcript]
            for (j = 1; j <= n; j++) {
                sample = sample_names[j]
                line = line OFS ((transcript, sample) in value ? value[transcript, sample] : "")
            }
            print line
        }
    }
' "${PRIMARY_TSV}" > "${TMP_COVERAGE}"

awk -F '\t' -v samples="${sample_order}" -v value_col="identity_pct" '
    BEGIN {
        OFS = "\t"
        n = split(samples, sample_names, ",")
        header = "transcript_id\tgene_id\ttranscript_name\tgene_name"
        for (i = 1; i <= n; i++) header = header OFS sample_names[i]
        print header
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        next
    }

    {
        sample = $col["sample"]
        transcript = $col["transcript_id"]
        if (!(transcript in seen)) {
            seen[transcript] = 1
            order[++count] = transcript
            gene_id[transcript] = $col["gene_id"]
            transcript_name[transcript] = $col["transcript_name"]
            gene_name[transcript] = $col["gene_name"]
        }
        value[transcript, sample] = $col[value_col]
    }

    END {
        for (i = 1; i <= count; i++) {
            transcript = order[i]
            line = transcript OFS gene_id[transcript] OFS transcript_name[transcript] OFS gene_name[transcript]
            for (j = 1; j <= n; j++) {
                sample = sample_names[j]
                line = line OFS ((transcript, sample) in value ? value[transcript, sample] : "")
            }
            print line
        }
    }
' "${PRIMARY_TSV}" > "${TMP_IDENTITY}"

mv "${TMP_STATUS}" "${STATUS_MATRIX}"
mv "${TMP_COVERAGE}" "${COVERAGE_MATRIX}"
mv "${TMP_IDENTITY}" "${IDENTITY_MATRIX}"

echo "Wrote ${PRIMARY_TSV}"
echo "Wrote ${STATUS_MATRIX}"
echo "Wrote ${COVERAGE_MATRIX}"
echo "Wrote ${IDENTITY_MATRIX}"
if [[ "${EXPORT_SUPPLEMENTARY}" == "1" ]]; then
    echo "Wrote ${SUPP_TSV}"
fi
