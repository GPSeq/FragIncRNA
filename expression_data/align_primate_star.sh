#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
READ_ROOT="${READ_ROOT:-$PWD/brain_rnaseq}"
REF_ROOT="${REF_ROOT:-$PWD/genomes}"
OUT_ROOT="${OUT_ROOT:-$PWD/star_results}"
THREADS="${THREADS:-16}"
SORT_RAM="${SORT_RAM:-8000000000}"
SPECIES="${1:-all}"
MAX_RUNS="${MAX_RUNS:-0}"
# metadata: table only; download: table + FASTQs; align: table + existing FASTQs;
# all: download, update table, then align. No large downloads unless requested.
MODE="${MODE:-align}"
# Optional text file containing one selected SRR/ERR/DRR run accession per line.
# Use your bulk RNA-seq run list, excluding single-nucleus experiments.
RUN_ACCESSIONS="${RUN_ACCESSIONS:-}"

case "$MODE" in metadata|download|align|all) ;; *) echo "Invalid MODE: $MODE" >&2; exit 1;; esac
case "$SPECIES" in all|human|chimpanzee|bonobo|macaque) ;; *) echo "Invalid species: $SPECIES" >&2; exit 1;; esac
[[ "$MAX_RUNS" =~ ^[0-9]+$ ]] || { echo 'MAX_RUNS must be a nonnegative integer' >&2; exit 1; }
command -v python3 >/dev/null || { echo 'Missing python3' >&2; exit 1; }

ARGS=(--read-root "$READ_ROOT" --species "$SPECIES" --max-runs "$MAX_RUNS")
[[ -z "$RUN_ACCESSIONS" ]] || ARGS+=(--accessions "$RUN_ACCESSIONS")
if [[ "$MODE" == download || "$MODE" == all ]]; then ARGS+=(--download); fi
python3 "$SCRIPT_DIR/primate_sample_table.py" "${ARGS[@]}"
if [[ "$MODE" == metadata || "$MODE" == download ]]; then exit 0; fi

for cmd in STAR samtools gzip; do
    command -v "$cmd" >/dev/null || { echo "Missing $cmd" >&2; exit 1; }
done

get_reference() {
    case "$1" in
        human) echo "$REF_ROOT/H9T2Thap1_GCA_054883195.1_H9_T2T.hap1_genomic.fna";;
        chimpanzee) echo "$REF_ROOT/mPanTro3v20pri_GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna";;
        bonobo) echo "$REF_ROOT/mPanPan1v20pri_GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna";;
        macaque) echo "$REF_ROOT/Mmul10_GCF_003339765.1_Mmul_10_genomic.fna";;
    esac
}

if [[ "$SPECIES" == all ]]; then SPECIES_LIST=(human chimpanzee bonobo macaque); else SPECIES_LIST=("$SPECIES"); fi
for species in "${SPECIES_LIST[@]}"; do
    REF=$(get_reference "$species")
    READ_DIR="$READ_ROOT/$species"
    INDEX_DIR="$OUT_ROOT/indexes/$species"
    BAM_DIR="$OUT_ROOT/bam/$species"
    # Consume only complete pairs selected by the metadata helper.
    PLAN="$READ_ROOT/metadata/alignment_plan_${SPECIES}.tsv"
    RUNS=()
    while IFS=$'\t' read -r sp run sample experiment; do
        [[ "$sp" == "$species" ]] && RUNS+=("$run|$sample|$experiment")
    done < "$PLAN"
    if [[ ${#RUNS[@]} -eq 0 ]]; then echo "No complete selected pairs for $species; skipping."; continue; fi
    [[ -s "$REF" ]] || { echo "Reference missing: $REF" >&2; exit 1; }
    mkdir -p "$INDEX_DIR" "$BAM_DIR"
    if [[ ! -s "$INDEX_DIR/Genome" || ! -s "$INDEX_DIR/SA" ]]; then
        STAR --runMode genomeGenerate --runThreadN "$THREADS" \
            --genomeDir "$INDEX_DIR" --genomeFastaFiles "$REF"
    fi
    for item in "${RUNS[@]}"; do
        IFS='|' read -r RUN SAMPLE EXPERIMENT <<< "$item"
        R1="$READ_DIR/${RUN}_1.fastq.gz"
        R2="$READ_DIR/${RUN}_2.fastq.gz"
        # Metadata-only mode can run during an external download. Align only after
        # that downloader finishes: files still being written are not safe inputs.
        gzip -t "$R1" "$R2"
        RUN_DIR="$BAM_DIR/$RUN"
        mkdir -p "$RUN_DIR"
        PREFIX="$RUN_DIR/${RUN}_"
        BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
        COMPLETE="${PREFIX}alignment.complete"
        if [[ -s "$BAM" && -f "$COMPLETE" ]] && samtools quickcheck -q "$BAM"; then
            echo "SKIP: $species / $RUN"
        else
            # A BAM quickcheck alone cannot establish that STAR finished.
            rm -f -- "$COMPLETE"
            STAR --runThreadN "$THREADS" --genomeDir "$INDEX_DIR" \
                --readFilesIn "$R1" "$R2" --readFilesCommand gzip -cd \
                --twopassMode Basic --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM "$SORT_RAM" \
                --outSAMattrRGline "ID:$RUN" "SM:$SAMPLE" "LB:$EXPERIMENT" "PL:ILLUMINA" \
                --outFileNamePrefix "$PREFIX"
            samtools quickcheck -v "$BAM"
            touch "$COMPLETE"
        fi
        [[ -s "${BAM}.bai" ]] || samtools index -@ "$THREADS" "$BAM"
        samtools flagstat -@ "$THREADS" "$BAM" > "${PREFIX}flagstat.txt"
    done
done
echo "Finished. Sample table: $READ_ROOT/metadata/samples.tsv"
echo "Alignments: $OUT_ROOT"
