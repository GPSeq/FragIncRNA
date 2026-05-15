#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Find shared local sequence regions among human lncRNA transcripts with MMseqs2.

Outputs:
  OUTDIR/transcript_clusters/
    transcript_cluster_cluster.tsv          transcript-level similarity clusters
    transcript_cluster_rep_seq.fasta        representative transcript sequences
    transcript_cluster_all_seqs.fasta       all transcript sequences in clusters

  OUTDIR/local_regions/
    all_vs_all.tsv                          raw MMseqs2 local alignments
    shared_blocks.filtered.tsv              filtered pairwise shared blocks
    shared_region_clusters.tsv              connected components of shared local regions
    shared_region_members.bed               BED-like member intervals per region cluster
    shared_region_sequences.fa              extracted sequence for each member interval
    run_summary.txt                         run settings and counts
  OUTDIR/METHODS_workflow.txt               methods-ready workflow description

Usage:
  scripts/find_shared_lncRNA_regions.sh [options]

Options:
  -i FASTA       Input lncRNA FASTA
                 default: gencode.v49.lncRNA_transcripts.fa
  -o OUTDIR      Output directory
                 default: lncRNA_shared_regions_mmseqs
  -t THREADS     Threads
                 default: all available logical CPUs where detectable
  -p FLOAT       Minimum sequence identity for local blocks, 0-1
                 default: 0.70
  -l INT         Minimum local alignment length in bp
                 default: 80
  -e EVALUE      Maximum e-value
                 default: 1e-10
  -m INT         MMseqs2 --max-seqs for local all-vs-all search
                 default: 1000
  -c FLOAT       Coverage for transcript-level clustering, 0-1
                 default: 0.70
  -P FLOAT       Minimum sequence identity for transcript-level clustering, 0-1
                 default: 0.80
  -g             Keep same-gene local block hits
                 default: remove same-gene hits
  -s FLOAT       MMseqs2 sensitivity
                 default: 7.5
  -L INT         MMseqs2 maximum sequence length
                 default: 1000000
  -h             Show this help

Examples:
  scripts/find_shared_lncRNA_regions.sh

  scripts/find_shared_lncRNA_regions.sh \
    -i gencode.v49.lncRNA_transcripts.fa \
    -o results_shared_regions \
    -t 12 \
    -p 0.85 \
    -l 100 \
    -m 5000

Notes:
  - This finds local shared regions, not one global multiple sequence alignment.
  - Low-complexity masking is enabled in MMseqs2. Repeat-derived lncRNA sequence
    can still create many hits. Increase -p/-l or lower -m if output explodes.
  - FASTA headers are expected to look like GENCODE transcript headers:
    transcript_id|gene_id|...
USAGE
}

fasta="gencode.v49.lncRNA_transcripts.fa"
run_command="$0 $*"
outdir="lncRNA_shared_regions_mmseqs"
threads=""
min_seq_id="0.70"
min_aln_len="80"
evalue="1e-10"
max_seqs="1000"
cluster_cov="0.70"
cluster_min_seq_id="0.80"
keep_same_gene="0"
sensitivity="7.5"
max_seq_len="1000000"

while getopts ":i:o:t:p:l:e:m:c:P:gs:L:h" opt; do
  case "${opt}" in
    i) fasta="${OPTARG}" ;;
    o) outdir="${OPTARG}" ;;
    t) threads="${OPTARG}" ;;
    p) min_seq_id="${OPTARG}" ;;
    l) min_aln_len="${OPTARG}" ;;
    e) evalue="${OPTARG}" ;;
    m) max_seqs="${OPTARG}" ;;
    c) cluster_cov="${OPTARG}" ;;
    P) cluster_min_seq_id="${OPTARG}" ;;
    g) keep_same_gene="1" ;;
    s) sensitivity="${OPTARG}" ;;
    L) max_seq_len="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) echo "Missing argument for -${OPTARG}" >&2; usage >&2; exit 2 ;;
    \?) echo "Unknown option: -${OPTARG}" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "${threads}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    threads="$(nproc)"
  elif command -v sysctl >/dev/null 2>&1; then
    threads="$(sysctl -n hw.logicalcpu)"
  else
    threads="4"
  fi
fi

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: required command not found on PATH: $1" >&2
    exit 1
  fi
}

require_cmd mmseqs
require_cmd awk
require_cmd python3

if [[ ! -s "${fasta}" ]]; then
  echo "ERROR: input FASTA does not exist or is empty: ${fasta}" >&2
  exit 1
fi

mkdir -p "${outdir}/transcript_clusters" "${outdir}/local_regions" "${outdir}/tmp"

cluster_prefix="${outdir}/transcript_clusters/transcript_cluster"
analysis_fasta="${outdir}/tmp/input.singleline.fa"
cluster_tmp="${outdir}/tmp/transcript_cluster_tmp"
search_tmp="${outdir}/tmp/all_vs_all_tmp"
raw_hits="${outdir}/local_regions/all_vs_all.tsv"
filtered_hits="${outdir}/local_regions/shared_blocks.filtered.tsv"
region_clusters="${outdir}/local_regions/shared_region_clusters.tsv"
region_members="${outdir}/local_regions/shared_region_members.bed"
region_fasta="${outdir}/local_regions/shared_region_sequences.fa"
summary="${outdir}/local_regions/run_summary.txt"
methods="${outdir}/METHODS_workflow.txt"

echo "[1/6] Normalizing FASTA records for MMseqs2"
awk '
  BEGIN { seq = "" }
  /^>/ {
    if (NR > 1) print seq
    print
    seq = ""
    next
  }
  NF > 0 {
    gsub(/[ \t\r]/, "")
    seq = seq $0
  }
  END {
    if (NR > 0) print seq
  }
' "${fasta}" > "${analysis_fasta}"

echo "[2/6] Transcript-level clustering with MMseqs2"
mmseqs easy-cluster \
  "${analysis_fasta}" \
  "${cluster_prefix}" \
  "${cluster_tmp}" \
  --dbtype 2 \
  --min-seq-id "${cluster_min_seq_id}" \
  -c "${cluster_cov}" \
  --cov-mode 0 \
  --cluster-mode 1 \
  --single-step-clustering 1 \
  --mask 1 \
  --shuffle 0 \
  --createdb-mode 1 \
  --max-seq-len "${max_seq_len}" \
  -s "${sensitivity}" \
  -v 2 \
  --threads "${threads}"

echo "[3/6] Local all-vs-all MMseqs2 search"
mmseqs easy-search \
  "${analysis_fasta}" \
  "${analysis_fasta}" \
  "${raw_hits}" \
  "${search_tmp}" \
  --search-type 3 \
  --min-seq-id "${min_seq_id}" \
  --min-aln-len "${min_aln_len}" \
  -e "${evalue}" \
  --max-seqs "${max_seqs}" \
  --mask 1 \
  --max-seq-len "${max_seq_len}" \
  -s "${sensitivity}" \
  --format-mode 4 \
  --format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov \
  -v 2 \
  --threads "${threads}"

echo "[4/6] Filtering self hits and same-gene hits"
awk -v keep_same_gene="${keep_same_gene}" '
  BEGIN { FS = OFS = "\t" }
  NR == 1 && $1 == "query" {
    print $0, "qgene", "tgene"
    next
  }
  {
    q = $1
    t = $2
    split(q, qa, "|")
    split(t, ta, "|")
    qgene = qa[2]
    tgene = ta[2]
    if (q == t) next
    if (!keep_same_gene && qgene != "" && tgene != "" && qgene == tgene) next
    print $0, qgene, tgene
  }
' "${raw_hits}" > "${filtered_hits}"

echo "[5/6] Clustering shared local intervals into connected region families"
python3 - "${filtered_hits}" "${analysis_fasta}" "${region_clusters}" "${region_members}" "${region_fasta}" <<'PY'
import sys
from collections import defaultdict

hits_path, fasta_path, clusters_path, members_path, regions_fa_path = sys.argv[1:6]

parent = []
rank = []
nodes = []
by_seq = defaultdict(list)

def make_node(seq_id, gene_id, start, end):
    if start > end:
        start, end = end, start
    idx = len(nodes)
    nodes.append((seq_id, gene_id, start, end))
    parent.append(idx)
    rank.append(0)
    by_seq[seq_id].append((start, end, idx))
    return idx

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra = find(a)
    rb = find(b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        parent[ra] = rb
    elif rank[ra] > rank[rb]:
        parent[rb] = ra
    else:
        parent[rb] = ra
        rank[ra] += 1

with open(hits_path) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    if not header or header[0] != "query":
        raise SystemExit("Expected a header line from MMseqs2 format-mode 4")
    col = {name: i for i, name in enumerate(header)}
    required = ["query", "target", "qstart", "qend", "tstart", "tend", "pident", "alnlen", "evalue", "bits", "qgene", "tgene"]
    missing = [name for name in required if name not in col]
    if missing:
        raise SystemExit("Missing required columns: " + ",".join(missing))
    for line in fh:
        if not line.strip():
            continue
        row = line.rstrip("\n").split("\t")
        q = row[col["query"]]
        t = row[col["target"]]
        qgene = row[col["qgene"]]
        tgene = row[col["tgene"]]
        qstart = int(row[col["qstart"]])
        qend = int(row[col["qend"]])
        tstart = int(row[col["tstart"]])
        tend = int(row[col["tend"]])
        qnode = make_node(q, qgene, qstart, qend)
        tnode = make_node(t, tgene, tstart, tend)
        union(qnode, tnode)

# Merge overlapping/near-identical local intervals on the same transcript.
# This joins multiple pairwise hits that point to the same underlying local region.
for seq_id, intervals in by_seq.items():
    intervals.sort(key=lambda x: (x[0], x[1]))
    active_end = -1
    active_node = None
    for start, end, idx in intervals:
        if active_node is not None and start <= active_end:
            union(active_node, idx)
            if end > active_end:
                active_end = end
                active_node = idx
        else:
            active_end = end
            active_node = idx

components = defaultdict(list)
for idx in range(len(nodes)):
    components[find(idx)].append(idx)

seqs = {}
current = None
parts = []
with open(fasta_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if current is not None:
                seqs[current] = "".join(parts)
            current = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip())
    if current is not None:
        seqs[current] = "".join(parts)

summaries = []
member_rows = []
component_items = sorted(components.items(), key=lambda kv: (-len({nodes[i][0] for i in kv[1]}), kv[0]))

for comp_num, (_root, idxs) in enumerate(component_items, start=1):
    comp_id = f"SR{comp_num:07d}"
    per_seq = defaultdict(list)
    genes = set()
    for idx in idxs:
        seq_id, gene_id, start, end = nodes[idx]
        per_seq[seq_id].append((start, end))
        if gene_id:
            genes.add(gene_id)
    merged_count = 0
    total_bp = 0
    rep = None
    for seq_id, intervals in sorted(per_seq.items()):
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            elif end > merged[-1][1]:
                merged[-1][1] = end
        gene_id = seq_id.split("|")[1] if "|" in seq_id else ""
        for start, end in merged:
            merged_count += 1
            span = end - start + 1
            total_bp += span
            if rep is None or span > rep[3]:
                rep = (seq_id, start, end, span)
            # BED uses 0-based start and 1-based exclusive end.
            member_rows.append((comp_id, seq_id, gene_id, start - 1, end, span))
    summaries.append((
        comp_id,
        len(idxs),
        merged_count,
        len(per_seq),
        len(genes),
        total_bp,
        rep[0] if rep else "",
        rep[1] if rep else "",
        rep[2] if rep else "",
        rep[3] if rep else "",
    ))

with open(clusters_path, "w") as out:
    out.write("region_cluster_id\tn_raw_intervals\tn_merged_member_intervals\tn_transcripts\tn_genes\ttotal_member_bp\trepresentative_transcript\trepresentative_start\trepresentative_end\trepresentative_len\n")
    for row in summaries:
        out.write("\t".join(map(str, row)) + "\n")

with open(members_path, "w") as out, open(regions_fa_path, "w") as fa:
    out.write("region_cluster_id\ttranscript_id\tgene_id\tstart0\tend1\tlength\n")
    for comp_id, seq_id, gene_id, start0, end1, span in member_rows:
        out.write(f"{comp_id}\t{seq_id}\t{gene_id}\t{start0}\t{end1}\t{span}\n")
        seq = seqs.get(seq_id)
        if seq is None:
            continue
        subseq = seq[start0:end1]
        fa.write(f">{comp_id}|{seq_id}|{gene_id}|{start0 + 1}-{end1}|len={len(subseq)}\n")
        for i in range(0, len(subseq), 60):
            fa.write(subseq[i:i + 60] + "\n")

print(f"region_clusters={len(summaries)}")
print(f"member_intervals={len(member_rows)}")
PY

echo "[6/6] Writing run summary and methods text"
{
  echo "input_fasta=${fasta}"
  echo "analysis_fasta=${analysis_fasta}"
  echo "outdir=${outdir}"
  echo "threads=${threads}"
  echo "mmseqs_version=$(mmseqs version)"
  echo "local_min_seq_id=${min_seq_id}"
  echo "local_min_aln_len=${min_aln_len}"
  echo "local_evalue=${evalue}"
  echo "local_max_seqs=${max_seqs}"
  echo "local_keep_same_gene=${keep_same_gene}"
  echo "cluster_min_seq_id=${cluster_min_seq_id}"
  echo "cluster_cov=${cluster_cov}"
  echo "sensitivity=${sensitivity}"
  echo "max_seq_len=${max_seq_len}"
  echo
  echo "input_transcripts=$(awk 'BEGIN{n=0} /^>/{n++} END{print n}' "${fasta}")"
  echo "raw_local_hits=$(awk 'END{n=NR-1; if(n<0)n=0; print n}' "${raw_hits}")"
  echo "filtered_local_hits=$(awk 'END{n=NR-1; if(n<0)n=0; print n}' "${filtered_hits}")"
  echo "shared_region_clusters=$(awk 'END{n=NR-1; if(n<0)n=0; print n}' "${region_clusters}")"
  echo "shared_region_member_intervals=$(awk 'END{n=NR-1; if(n<0)n=0; print n}' "${region_members}")"
} > "${summary}"

cat > "${methods}" <<METHODS
Methods workflow: identification of shared local sequence regions among human lncRNA transcripts

Analysis overview
This workflow was designed to identify local sequence regions shared among human long non-coding RNA (lncRNA) transcripts. The analysis does not perform one global multiple sequence alignment of all lncRNAs. Instead, it uses a local all-against-all nucleotide similarity search to identify pairwise shared blocks, filters trivial hits, and then clusters overlapping local intervals into shared-region families.

Input data
Input FASTA: ${fasta}
MMseqs2 working FASTA: ${analysis_fasta}
The workflow expects GENCODE-style transcript FASTA headers in which fields are separated by pipe characters, with the transcript identifier in field 1 and the gene identifier in field 2. For example:
transcript_id|gene_id|...

The input FASTA was normalized to a working FASTA containing one sequence line per record before MMseqs2 was run. This does not change sequence identifiers or sequence content; it avoids MMseqs2 automatic recomputation messages caused by multiline FASTA records.

Software
MMseqs2 version: $(mmseqs version)
Python version: $(python3 --version 2>&1)
awk was used for tabular filtering and summary counts.

Run command
${run_command}

Runtime parameters
Threads: ${threads}
MMseqs2 sensitivity (-s): ${sensitivity}
MMseqs2 maximum sequence length (--max-seq-len): ${max_seq_len}
Low-complexity masking: enabled with MMseqs2 --mask 1

Transcript-level clustering parameters
Minimum sequence identity (--min-seq-id): ${cluster_min_seq_id}
Coverage threshold (-c): ${cluster_cov}
Coverage mode (--cov-mode): 0, requiring coverage of both query and target sequences
Cluster mode (--cluster-mode): 1, connected-component clustering
Single-step clustering (--single-step-clustering): enabled
Nucleotide database mode (--dbtype): 2

Local all-against-all search parameters
Search type (--search-type): 3, nucleotide-vs-nucleotide search
Minimum local sequence identity (--min-seq-id): ${min_seq_id}
Minimum local alignment length (--min-aln-len): ${min_aln_len} bp
Maximum e-value (-e): ${evalue}
Maximum reported target sequences per query (--max-seqs): ${max_seqs}
Output format: tabular MMseqs2 alignment output with query, target, percent identity, alignment length, mismatch count, gap-open count, query start/end, target start/end, e-value, bit score, query length, target length, query coverage, and target coverage.

Step-by-step algorithm
1. Transcript-level clustering
   MMseqs2 easy-cluster was run on the normalized lncRNA FASTA in nucleotide mode. This step groups transcripts that are globally similar according to the transcript-level sequence identity and coverage thresholds above. Connected-component clustering was used so that transcripts connected by qualifying similarity edges are assigned to the same transcript cluster. These transcript clusters are useful for identifying highly similar transcript families or likely isoform-like/redundant sequence groups, but they are not the primary definition of local shared regions.

2. Local all-against-all nucleotide similarity search
   MMseqs2 easy-search was run with the normalized lncRNA FASTA as both query and target. This produces local pairwise alignments between lncRNA transcripts. Local alignment is used because lncRNAs can share short sequence blocks while being otherwise unrelated in full-length transcript structure.

3. Filtering of local hits
   The raw local alignment table was filtered to remove self-hits where the query and target transcript identifiers were identical. By default, hits between transcripts from the same gene were also removed using the second pipe-delimited FASTA header field as the gene identifier. This focuses the analysis on sequence sharing between different lncRNA genes. If the script is run with -g, same-gene hits are retained.

4. Conversion of pairwise hits to interval nodes
   Each retained pairwise local alignment contributes two transcript intervals: one interval on the query transcript and one interval on the target transcript. The two intervals from the same alignment are connected because they represent the same pairwise shared sequence block.

5. Clustering of local intervals into shared-region families
   A union-find connected-component algorithm was used to cluster local intervals. First, query and target intervals from each retained pairwise alignment were joined. Second, intervals that overlapped on the same transcript were also joined. This creates connected components representing shared-region families. Each family may contain intervals from multiple transcripts and genes.

6. Merging and summarizing member intervals
   Within each shared-region family, overlapping intervals on the same transcript were merged into non-overlapping member intervals. For each region family, the workflow reports the number of raw intervals, the number of merged member intervals, the number of transcripts, the number of genes, the total member base pairs, and the longest representative member interval.

7. Sequence extraction
   The merged member intervals were extracted from the normalized working FASTA and written to a FASTA file. The working FASTA has the same sequence identifiers and sequence content as the input FASTA. These extracted sequences can be used for downstream multiple sequence alignment, motif discovery, repeat annotation, or manual inspection of individual shared-region families.

Output files
${cluster_prefix}_cluster.tsv
Transcript-level MMseqs2 cluster table.

${cluster_prefix}_rep_seq.fasta
Representative transcript sequence for each transcript-level cluster.

${cluster_prefix}_all_seqs.fasta
All transcript sequences grouped by transcript-level cluster.

${raw_hits}
Raw MMseqs2 local all-against-all alignment table.

${filtered_hits}
Filtered local alignment table after removal of self-hits and, unless -g was used, same-gene hits.

${region_clusters}
Summary table of shared-region families.

${region_members}
BED-like table of merged member intervals for each shared-region family. Coordinates are 0-based start and 1-based exclusive end.

${region_fasta}
FASTA file containing the extracted member sequences for each shared-region family.

${summary}
Machine-readable run summary with key parameter values and output counts.

Interpretation notes
The resulting shared-region families represent local sequence similarity among lncRNA transcripts, not necessarily evolutionary conservation. Human lncRNAs often contain repeat-derived and low-complexity sequence. MMseqs2 low-complexity masking is enabled, but repeat-derived blocks can still contribute to shared-region clusters. For biological interpretation, shared-region families should be annotated against repeat databases and, where relevant, compared with independent evolutionary conservation tracks or orthology-aware alignments.

Suggested software citation
MMseqs2: Steinegger M. and Soding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology 35, 1026-1028 (2017). For clustering: Steinegger M. and Soding J. Clustering huge protein sequence sets in linear time. Nature Communications 9, 2542 (2018).

MMseqs2 messages addressed
The workflow normalizes the input FASTA to one sequence line per record and explicitly sets --shuffle 0 and --createdb-mode 1 during transcript-level clustering to avoid automatic recomputation warnings related to shuffled databases and multiline FASTA input. The workflow also enables --single-step-clustering 1 with connected-component clustering, following the MMseqs2 recommendation printed when connected-component mode is used.
METHODS

echo
echo "Done."
echo "Main outputs:"
echo "  ${cluster_prefix}_cluster.tsv"
echo "  ${filtered_hits}"
echo "  ${region_clusters}"
echo "  ${region_members}"
echo "  ${region_fasta}"
echo "  ${summary}"
echo "  ${methods}"
