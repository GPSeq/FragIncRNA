#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Check gzipped genome FASTA files for N/n bases.

Usage:
  bash check_genome_N_content.sh [GENOME_DIR] [OUT_TSV]

Defaults:
  GENOME_DIR = ../../genomes
  OUT_TSV    = genome_N_content.tsv

The script scans *.gz files in GENOME_DIR, ignores FASTA headers, and reports
total sequence bases, N/n bases, and N percentage for each file.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

GENOME_DIR="${1:-../../genomes}"
OUT_TSV="${2:-genome_N_content.tsv}"

if [[ ! -d "${GENOME_DIR}" ]]; then
    echo "ERROR: genome directory not found: ${GENOME_DIR}" >&2
    exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip not found in PATH" >&2
    exit 1
fi

shopt -s nullglob
GENOMES=("${GENOME_DIR}"/*.gz)
shopt -u nullglob

if [[ ${#GENOMES[@]} -eq 0 ]]; then
    echo "ERROR: no .gz files found in ${GENOME_DIR}" >&2
    exit 1
fi

printf "file\ttotal_bases\tn_bases\tn_percent\tcontains_N\n" > "${OUT_TSV}"

for genome in "${GENOMES[@]}"; do
    echo "Checking ${genome}" >&2
    gzip -cd "${genome}" | awk -v file="$(basename "${genome}")" '
        BEGIN {
            total = 0
            n_count = 0
        }

        /^>/ {
            next
        }

        {
            gsub(/[ \t\r\n]/, "")
            total += length($0)
            n_count += gsub(/[Nn]/, "", $0)
        }

        END {
            pct = total ? (100.0 * n_count / total) : 0
            print file, total, n_count, sprintf("%.6f", pct), (n_count > 0 ? "yes" : "no")
        }
    ' OFS='\t' >> "${OUT_TSV}"
done

echo "Wrote ${OUT_TSV}" >&2
