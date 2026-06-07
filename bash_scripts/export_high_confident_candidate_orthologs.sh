#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Export high-confidence candidate orthologs from transcript alignment QC.

Definition:
  HIGH_CONFIDENT = query coverage >= 90%
                   alignment identity >= 90%
                   MAPQ >= 30

These are sequence-based candidate orthologs. Synteny and reciprocal-best-hit
checks are still required before calling them confirmed orthologs.

Usage:
  bash bash_scripts/export_high_confident_candidate_orthologs.sh [options]

Options:
  -i FILE   Input lncrna_transcript_alignment_qc.tsv
            default: ../bam_comparison/qc/lncrna_transcript_alignment_qc.tsv
  -o DIR    Output directory
            default: alignment_results/high_confident_candidate_orthologs
  -c FLOAT  Minimum query coverage percentage
            default: 90
  -p FLOAT  Minimum identity percentage
            default: 90
  -q INT    Minimum MAPQ
            default: 30
  -h        Show this help

Outputs:
  high_confident_candidate_ortholog_hits.tsv
      One row per transcript/genome passing all high-confidence thresholds.

  high_confident_candidate_ortholog_status_matrix.tsv
      Per-transcript statuses across genomes:
      HIGH_CONFIDENT, PASS_STRICT, PASS_BASIC, LOW_QC, UNMAPPED, or MISSING.

  high_confident_candidate_ortholog_pass_counts.tsv
      High-confidence, strict, and basic pass counts per transcript.

  high_confident_candidate_orthologs_all_genomes.tsv
      Transcripts classified HIGH_CONFIDENT in every input genome.

  high_confident_candidate_ortholog_summary_by_genome.tsv
      Per-genome transcript and unique-gene counts.

  high_confident_candidate_ortholog_summary.tsv
      Overall thresholds and all-genome counts.
EOF
}

INPUT="../bam_comparison/qc/lncrna_transcript_alignment_qc.tsv"
OUT_DIR="alignment_results/high_confident_candidate_orthologs"
MIN_COVERAGE="90"
MIN_IDENTITY="90"
MIN_MAPQ="30"

while getopts ":i:o:c:p:q:h" opt; do
    case "${opt}" in
        i) INPUT="${OPTARG}" ;;
        o) OUT_DIR="${OPTARG}" ;;
        c) MIN_COVERAGE="${OPTARG}" ;;
        p) MIN_IDENTITY="${OPTARG}" ;;
        q) MIN_MAPQ="${OPTARG}" ;;
        h) usage; exit 0 ;;
        :) echo "ERROR: option -${OPTARG} requires an argument" >&2; usage >&2; exit 2 ;;
        \?) echo "ERROR: unknown option -${OPTARG}" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: input alignment QC table not found or empty: ${INPUT}" >&2
    exit 1
fi

is_number() {
    [[ "$1" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]]
}

if ! is_number "${MIN_COVERAGE}"; then
    echo "ERROR: -c must be a non-negative number" >&2
    exit 2
fi
if ! is_number "${MIN_IDENTITY}"; then
    echo "ERROR: -p must be a non-negative number" >&2
    exit 2
fi
if ! [[ "${MIN_MAPQ}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: -q must be a non-negative integer" >&2
    exit 2
fi

mkdir -p "${OUT_DIR}"

HITS="${OUT_DIR}/high_confident_candidate_ortholog_hits.tsv"
STATUS_MATRIX="${OUT_DIR}/high_confident_candidate_ortholog_status_matrix.tsv"
PASS_COUNTS="${OUT_DIR}/high_confident_candidate_ortholog_pass_counts.tsv"
ALL_GENOMES="${OUT_DIR}/high_confident_candidate_orthologs_all_genomes.tsv"
BY_GENOME="${OUT_DIR}/high_confident_candidate_ortholog_summary_by_genome.tsv"
SUMMARY="${OUT_DIR}/high_confident_candidate_ortholog_summary.tsv"

TMP_PREFIX="${OUT_DIR}/.high_confident_candidate_orthologs.$$"
TMP_HITS="${TMP_PREFIX}.hits"
TMP_STATUS="${TMP_PREFIX}.status"
TMP_COUNTS="${TMP_PREFIX}.counts"
TMP_ALL="${TMP_PREFIX}.all"
TMP_BY_GENOME="${TMP_PREFIX}.by_genome"
TMP_SUMMARY="${TMP_PREFIX}.summary"

cleanup() {
    rm -f \
        "${TMP_HITS}" \
        "${TMP_STATUS}" \
        "${TMP_COUNTS}" \
        "${TMP_ALL}" \
        "${TMP_BY_GENOME}" \
        "${TMP_SUMMARY}"
}
trap cleanup EXIT

echo "Classifying candidate ortholog alignments from ${INPUT}" >&2

awk -F '\t' \
    -v min_coverage="${MIN_COVERAGE}" \
    -v min_identity="${MIN_IDENTITY}" \
    -v min_mapq="${MIN_MAPQ}" \
    -v hits_out="${TMP_HITS}" \
    -v status_out="${TMP_STATUS}" \
    -v counts_out="${TMP_COUNTS}" \
    -v all_out="${TMP_ALL}" \
    -v genome_out="${TMP_BY_GENOME}" \
    -v summary_out="${TMP_SUMMARY}" '
    BEGIN {
        OFS = "\t"
    }

    NR == 1 {
        for (i = 1; i <= NF; i++) {
            col[$i] = i
        }
        required_text = "sample transcript_id gene_id transcript_name gene_name mapped reference position strand mapq cigar query_coverage_pct identity_pct pass_basic pass_strict"
        required_count = split(required_text, required, " ")
        for (i = 1; i <= required_count; i++) {
            if (!(required[i] in col)) {
                print "ERROR: missing required column: " required[i] > "/dev/stderr"
                exit 2
            }
        }

        print \
            "sample", "transcript_id", "gene_id", "transcript_name", "gene_name", \
            "reference", "position", "strand", "mapq", "cigar", \
            "query_coverage_pct", "identity_pct", "nm", "as_score", \
            "has_sa_tag", "candidate_class" > hits_out
        next
    }

    {
        sample = $col["sample"]
        transcript = $col["transcript_id"]
        gene = $col["gene_id"]
        mapped = $col["mapped"] + 0
        coverage = $col["query_coverage_pct"] + 0
        identity_text = $col["identity_pct"]
        identity = identity_text == "" ? -1 : identity_text + 0
        mapq = $col["mapq"] + 0
        pass_basic = $col["pass_basic"] + 0
        pass_strict = $col["pass_strict"] + 0
        high_confident = mapped && coverage >= min_coverage && \
            identity >= min_identity && mapq >= min_mapq

        if (!(sample in sample_seen)) {
            sample_seen[sample] = 1
            samples[++sample_count] = sample
        }
        if (!(transcript in transcript_seen)) {
            transcript_seen[transcript] = 1
            transcripts[++transcript_count] = transcript
            gene_id[transcript] = gene
            transcript_name[transcript] = $col["transcript_name"]
            gene_name[transcript] = $col["gene_name"]
        }

        if (high_confident) {
            status = "HIGH_CONFIDENT"
            high_count[transcript]++
            genome_high_count[sample]++
            genome_high_gene[sample, gene] = 1
            print \
                sample, transcript, gene, $col["transcript_name"], \
                $col["gene_name"], $col["reference"], $col["position"], \
                $col["strand"], mapq, $col["cigar"], \
                sprintf("%.2f", coverage), sprintf("%.2f", identity), \
                (("nm" in col) ? $col["nm"] : ""), \
                (("as_score" in col) ? $col["as_score"] : ""), \
                (("has_sa_tag" in col) ? $col["has_sa_tag"] : ""), \
                "candidate_ortholog" > hits_out
        } else if (pass_strict) {
            status = "PASS_STRICT"
        } else if (pass_basic) {
            status = "PASS_BASIC"
        } else if (!mapped) {
            status = "UNMAPPED"
        } else {
            status = "LOW_QC"
        }

        if (pass_strict) {
            strict_count[transcript]++
        }
        if (pass_basic) {
            basic_count[transcript]++
        }
        status_by_sample[transcript, sample] = status
        genome_total[sample]++
    }

    END {
        if (NR == 1) {
            print "ERROR: input contains a header but no data rows" > "/dev/stderr"
            exit 2
        }

        matrix_header = "transcript_id" OFS "gene_id" OFS "transcript_name" OFS "gene_name"
        for (j = 1; j <= sample_count; j++) {
            matrix_header = matrix_header OFS samples[j]
        }
        print matrix_header > status_out

        print \
            "transcript_id", "gene_id", "transcript_name", "gene_name", \
            "high_confident_pass_count", "strict_pass_count", \
            "basic_pass_count", "genome_count", \
            "high_confident_all_genomes" > counts_out
        print \
            "transcript_id", "gene_id", "transcript_name", "gene_name", \
            "high_confident_pass_count", "genome_count" > all_out

        all_genome_transcripts = 0
        for (i = 1; i <= transcript_count; i++) {
            transcript = transcripts[i]
            matrix_line = transcript OFS gene_id[transcript] OFS \
                transcript_name[transcript] OFS gene_name[transcript]
            for (j = 1; j <= sample_count; j++) {
                sample = samples[j]
                key = transcript SUBSEP sample
                matrix_line = matrix_line OFS \
                    ((key in status_by_sample) ? status_by_sample[key] : "MISSING")
            }
            print matrix_line > status_out

            high = high_count[transcript] + 0
            strict = strict_count[transcript] + 0
            basic = basic_count[transcript] + 0
            all_high = high == sample_count ? 1 : 0
            print \
                transcript, gene_id[transcript], transcript_name[transcript], \
                gene_name[transcript], high, strict, basic, sample_count, \
                all_high > counts_out
            if (all_high) {
                all_genome_transcripts++
                all_genome_gene[gene_id[transcript]] = 1
                print \
                    transcript, gene_id[transcript], transcript_name[transcript], \
                    gene_name[transcript], high, sample_count > all_out
            }
        }

        print \
            "genome", "total_transcripts", "high_confident_transcripts", \
            "high_confident_percent", "high_confident_unique_genes" > genome_out
        for (j = 1; j <= sample_count; j++) {
            sample = samples[j]
            unique_genes = 0
            for (key in genome_high_gene) {
                split(key, parts, SUBSEP)
                if (parts[1] == sample) {
                    unique_genes++
                }
            }
            percent = genome_total[sample] ? 100.0 * \
                genome_high_count[sample] / genome_total[sample] : 0
            print \
                sample, genome_total[sample] + 0, genome_high_count[sample] + 0, \
                sprintf("%.2f", percent), unique_genes > genome_out
        }

        all_genome_genes = 0
        for (gene in all_genome_gene) {
            all_genome_genes++
        }
        print "metric", "value" > summary_out
        print "minimum_query_coverage_pct", min_coverage > summary_out
        print "minimum_identity_pct", min_identity > summary_out
        print "minimum_mapq", min_mapq > summary_out
        print "genome_count", sample_count > summary_out
        print "total_transcripts", transcript_count > summary_out
        print "high_confident_all_genomes_transcripts", all_genome_transcripts > summary_out
        print "high_confident_all_genomes_unique_genes", all_genome_genes > summary_out
    }
' "${INPUT}"

mv "${TMP_HITS}" "${HITS}"
mv "${TMP_STATUS}" "${STATUS_MATRIX}"
mv "${TMP_COUNTS}" "${PASS_COUNTS}"
mv "${TMP_ALL}" "${ALL_GENOMES}"
mv "${TMP_BY_GENOME}" "${BY_GENOME}"
mv "${TMP_SUMMARY}" "${SUMMARY}"

echo "Wrote ${HITS}"
echo "Wrote ${STATUS_MATRIX}"
echo "Wrote ${PASS_COUNTS}"
echo "Wrote ${ALL_GENOMES}"
echo "Wrote ${BY_GENOME}"
echo "Wrote ${SUMMARY}"
