# FragIncRNA Bash Scripts

This directory contains helper scripts for downloading primate genomes, running
human lncRNA searches against those genomes, summarizing BAM alignment quality,
and finding shared lncRNA transcripts or local sequence regions.

Most scripts contain hard-coded default paths. Check and edit the variables near
the top of each script before running them on a different machine or dataset.

## Typical Workflow

1. Download primate reference genomes with `download_primates.sh`.
2. Align human lncRNA transcripts to each primate genome with
   `run_minimap2.sh`.
3. Put the generated sorted BAM paths in `bam_list.txt`.
4. Generate per-sample alignment QC with `lncrna_alignment_qc.sh`.
5. Export per-transcript alignment metrics and status matrices with
   `export_lncRNA_bam_alignment_info.sh`.
6. Find transcripts passing QC in all genomes with
   `find_shared_lncRNA_transcripts.sh`.
7. Optionally run the k-mer workflow with `run_FraglncRNA_k15.sh` and/or find
   shared local sequence regions with `find_shared_lncRNA_regions.sh`.

## Scripts

### `download_primates.sh`

Downloads selected primate genome FASTA files from NCBI FTP using `wget`.

Inputs:

- No input files.
- Network access to NCBI FTP.

Downloaded references:

- Bonobo:
  - `GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.fna.gz`
  - `GCF_000258655.2_panpan1.1_genomic.fna.gz`
  - `GCF_013052645.1_Mhudiblu_PPA_v0_genomic.fna.gz`
- Chimpanzee:
  - `GCF_002880755.1_Clint_PTRv2_genomic.fna.gz`
  - `GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna.gz`
- Human:
  - `GCA_054883195.1_H9_T2T.hap1_genomic.fna.gz`
  - `GCA_054883265.1_H9_T2T.hap2_genomic.fna.gz`
  - `GCA_018503265.2_NA19240_pat_hprc_f2_genomic.fna.gz`
  - `GCA_018503275.2_NA19240_mat_hprc_f2_genomic.fna.gz`
- Rhesus macaque:
  - `GCF_003339765.1_Mmul_10_genomic.fna.gz`
  - `GCF_049350105.2_T2T-MMU8v2.0_genomic.fna.gz`

Generated files:

- The `.fna.gz` genome files listed above, written to the current working
  directory where the script is executed.

Dependencies:

- `wget`

### `run_minimap2.sh`

Builds a temporary minimap2 index for each primate genome, aligns human lncRNA
transcripts to that genome, converts the SAM output to sorted BAM, indexes the
BAM, and removes temporary SAM/index files.

Default inputs:

- Reference genome directory: `./data/all_primates`
- Reference file extensions searched: `*.fasta.gz`, `*.fa.gz`, `*.fna.gz`,
  `*.fas.gz`
- Query transcript FASTA: `./gencode.v49.lncRNA_transcripts.fa`
- Minimap2 executable: `./minimap2-2.30_x64-linux/minimap2`

Generated files:

- Output directory: `./output`
- Log file: `./output/run_all_primates_YYYYMMDD_HHMMSS.log`
- For each reference genome:
  - `./output/human_lncRNA_vs_<prefix>.sorted.bam`
  - `./output/human_lncRNA_vs_<prefix>.sorted.bam.bai`

Temporary files removed by the script:

- `./output/<reference_basename>.mmi`
- `./output/human_lncRNA_vs_<prefix>.sam`

Notes:

- `<prefix>` is taken from the reference basename before the first underscore.
- Uses minimap2 with `-ax splice:hq`, `--secondary=no`, and k-mer size `15`.
- Default thread settings are `THREADS_MINIMAP=30` and `THREADS_SAMTOOLS=8`.

Dependencies:

- `minimap2`
- `samtools`
- `tee`, `find`, `wc`

### `run_FraglncRNA_k15.sh`

Creates a TOML configuration file and runs the compiled `lncrna_mers` binary to
search lncRNA transcript k-mers against reference genomes using an IBF or HIBF
index.

Default inputs:

- Reference genome directory: `/mnt/d/primates_bmc/genomes/`
- Query transcript FASTA:
  `/mnt/d/primates_bmc/lncrna/gencode.v49.lncRNA_transcripts.fa`
- Executable: `./build/lncrna_mers`

Generated files:

- Output directory: `/mnt/d/primates_bmc/output_kmers`
- Generated config:
  `/mnt/d/primates_bmc/output_kmers/config_k15_per_reference.toml`
- Main configured result file:
  `/mnt/d/primates_bmc/output_kmers/results.tsv`
- Main configured log file:
  `/mnt/d/primates_bmc/output_kmers/ibf_run.log`
- Per-reference result files in the output directory.
- Unique k-mer files in:
  `/mnt/d/primates_bmc/output_kmers/unique_mers`

Important default parameters:

- `INDEX_METHOD=ibf`
- `KMER_SIZE=15`
- `FRAGMENT_SIZE=80000`
- `HIT_THRESHOLD=15`
- `THREADS=1`
- `STORE_FRAGMENTS=false`
- `STORE_INDEX=false`
- `CLEANUP_INDEX=true`

Dependencies:

- Built and executable `./build/lncrna_mers`
- Reference genomes in the configured `REF_DIR`
- Query FASTA at the configured `QUERY_FILE`

### `lncrna_alignment_qc.sh`

Summarizes alignment quality for each sorted BAM in a BAM list. It reports
mapping rates, MAPQ distribution, mean coverage, mean identity, soft clipping,
SA-tag frequency, pass counts under basic and strict QC thresholds, and the
percentage of alignments on large contigs.

Usage:

```bash
./lncrna_alignment_qc.sh [BAM_LIST] [OUT_DIR]
```

Default inputs:

- BAM list: `bam_list.txt`
- Output directory: `bam_comparison/qc`

Input file format:

- `BAM_LIST` is a plain text file with one BAM path per line.
- Blank lines and lines beginning with `#` are ignored.
- BAM files should be sorted BAMs, usually the files generated by
  `run_minimap2.sh`.

Generated files:

- `bam_comparison/qc/lncrna_alignment_sample_qc.tsv`

Configurable environment variables:

- `MIN_COV` default `80`
- `MIN_ID` default `80`
- `MIN_MAPQ` default `10`
- `STRICT_MIN_ID` default `90`
- `STRICT_MIN_MAPQ` default `30`
- `LARGE_CONTIG_BP` default `10000000`

Dependencies:

- `samtools`
- `awk`

### `export_lncRNA_bam_alignment_info.sh`

Exports per-transcript alignment metrics from each BAM and builds matrices that
compare transcript status, query coverage, and identity across samples/genomes.

Usage:

```bash
./export_lncRNA_bam_alignment_info.sh [BAM_LIST] [OUT_DIR]
```

Default inputs:

- BAM list: `bam_list.txt`
- Output directory: `bam_comparison/qc`

Input file format:

- `BAM_LIST` is a plain text file with one BAM path per line.
- Blank lines and lines beginning with `#` are ignored.
- Transcript names are expected to use GENCODE-style pipe-delimited FASTA
  headers, for example `transcript_id|gene_id|...|transcript_name|gene_name|...`.

Generated files:

- `bam_comparison/qc/lncrna_transcript_alignment_qc.tsv`
- `bam_comparison/qc/lncrna_transcript_status_matrix.tsv`
- `bam_comparison/qc/lncrna_transcript_coverage_matrix.tsv`
- `bam_comparison/qc/lncrna_transcript_identity_matrix.tsv`

Optional generated file:

- If `EXPORT_SUPPLEMENTARY=1`:
  `bam_comparison/qc/lncrna_supplementary_alignment_segments.tsv`

Status values in `lncrna_transcript_status_matrix.tsv`:

- `PASS_STRICT`: mapped, coverage is at least `MIN_COV`, identity is at least
  `STRICT_MIN_ID`, and MAPQ is at least `STRICT_MIN_MAPQ`.
- `PASS_BASIC`: mapped, coverage is at least `MIN_COV`, identity is at least
  `MIN_ID`, and MAPQ is at least `MIN_MAPQ`.
- `LOW_QC`: mapped but did not pass basic QC.
- `UNMAPPED`: transcript is unmapped.
- `MISSING`: transcript/sample combination was not present in the exported data.

Configurable environment variables:

- `MIN_COV` default `80`
- `MIN_ID` default `80`
- `MIN_MAPQ` default `10`
- `STRICT_MIN_ID` default `90`
- `STRICT_MIN_MAPQ` default `30`
- `EXPORT_SUPPLEMENTARY` default `0`

Dependencies:

- `samtools`
- `awk`

### `find_shared_lncRNA_transcripts.sh`

Uses the transcript status matrix to identify lncRNA transcripts that pass QC in
all genomes/samples.

Usage:

```bash
./find_shared_lncRNA_transcripts.sh [STATUS_MATRIX] [OUT_DIR]
```

Default inputs:

- Status matrix:
  `bam_comparison/qc/lncrna_transcript_status_matrix.tsv`
- Output directory:
  `bam_comparison/qc/shared_lncRNA_transcripts`

Input file format:

- A tab-separated status matrix generated by
  `export_lncRNA_bam_alignment_info.sh`.
- First four columns are transcript metadata:
  `transcript_id`, `gene_id`, `transcript_name`, `gene_name`.
- Remaining columns are per-sample status values.

Generated files:

- `bam_comparison/qc/shared_lncRNA_transcripts/shared_lncRNA_transcripts_strict_all_genomes.tsv`
- `bam_comparison/qc/shared_lncRNA_transcripts/shared_lncRNA_transcripts_basic_all_genomes.tsv`
- `bam_comparison/qc/shared_lncRNA_transcripts/lncRNA_transcript_pass_counts.tsv`
- `bam_comparison/qc/shared_lncRNA_transcripts/shared_lncRNA_transcripts_summary.tsv`

Output meaning:

- Strict shared transcripts are transcripts with `PASS_STRICT` in every
  sample/genome column.
- Basic shared transcripts are transcripts with either `PASS_STRICT` or
  `PASS_BASIC` in every sample/genome column.
- Pass counts report the number of strict and basic passes per transcript.
- The summary reports total transcripts and the number/percentage shared under
  strict and basic rules.

Dependencies:

- `awk`

### `find_shared_lncRNA_regions.sh`

Finds local sequence regions shared among human lncRNA transcripts using
MMseqs2. This is an all-vs-all local similarity workflow, not a global multiple
sequence alignment.

Usage:

```bash
./find_shared_lncRNA_regions.sh [options]
```

Main options:

- `-i FASTA`: input lncRNA FASTA, default
  `gencode.v49.lncRNA_transcripts.fa`
- `-o OUTDIR`: output directory, default `lncRNA_shared_regions_mmseqs`
- `-t THREADS`: number of threads, default is detected from the system
- `-p FLOAT`: minimum local sequence identity, default `0.70`
- `-l INT`: minimum local alignment length in bp, default `80`
- `-e EVALUE`: maximum e-value, default `1e-10`
- `-m INT`: MMseqs2 `--max-seqs`, default `1000`
- `-c FLOAT`: transcript-level clustering coverage, default `0.70`
- `-P FLOAT`: transcript-level clustering minimum identity, default `0.80`
- `-g`: keep same-gene local block hits; by default same-gene hits are removed
- `-s FLOAT`: MMseqs2 sensitivity, default `7.5`
- `-L INT`: MMseqs2 maximum sequence length, default `1000000`
- `-h`: show help

Input file format:

- FASTA of lncRNA transcripts.
- Headers are expected to be GENCODE-style pipe-delimited identifiers with the
  transcript ID in field 1 and gene ID in field 2.

Generated files:

- Transcript-level cluster outputs:
  - `OUTDIR/transcript_clusters/transcript_cluster_cluster.tsv`
  - `OUTDIR/transcript_clusters/transcript_cluster_rep_seq.fasta`
  - `OUTDIR/transcript_clusters/transcript_cluster_all_seqs.fasta`
- Local-region outputs:
  - `OUTDIR/local_regions/all_vs_all.tsv`
  - `OUTDIR/local_regions/shared_blocks.filtered.tsv`
  - `OUTDIR/local_regions/shared_region_clusters.tsv`
  - `OUTDIR/local_regions/shared_region_members.bed`
  - `OUTDIR/local_regions/shared_region_sequences.fa`
  - `OUTDIR/local_regions/run_summary.txt`
- Methods text:
  - `OUTDIR/METHODS_workflow.txt`
- Temporary working files:
  - `OUTDIR/tmp/input.singleline.fa`
  - MMseqs2 temporary directories under `OUTDIR/tmp`

Output meaning:

- `all_vs_all.tsv` contains raw MMseqs2 local alignments.
- `shared_blocks.filtered.tsv` removes self hits and, by default, same-gene
  hits.
- `shared_region_clusters.tsv` summarizes connected shared-region families.
- `shared_region_members.bed` lists merged member intervals. Coordinates use
  0-based start and 1-based exclusive end.
- `shared_region_sequences.fa` contains extracted sequences for each member
  interval.
- `run_summary.txt` records parameters and output counts.
- `METHODS_workflow.txt` provides a methods-ready description of the workflow.

Dependencies:

- `mmseqs`
- `awk`
- `python3`

## Common Input Files

- `gencode.v49.lncRNA_transcripts.fa`: human lncRNA transcript FASTA used as the
  query sequence set by the minimap2, k-mer, and MMseqs2 region workflows.
- Primate genome FASTA files: gzipped reference genomes downloaded by
  `download_primates.sh` or placed in the directories configured in the scripts.
- `bam_list.txt`: one sorted BAM path per line, usually pointing to BAM files
  generated by `run_minimap2.sh`.

## Common Generated Directories

- `output/`: minimap2 BAMs, BAM indexes, and run logs.
- `bam_comparison/qc/`: sample-level QC, transcript-level QC, and comparison
  matrices.
- `bam_comparison/qc/shared_lncRNA_transcripts/`: shared-transcript tables and
  summaries.
- `lncRNA_shared_regions_mmseqs/`: MMseqs2 transcript clusters and local shared
  region outputs.
- `/mnt/d/primates_bmc/output_kmers/`: default output location for the k-mer
  workflow.
