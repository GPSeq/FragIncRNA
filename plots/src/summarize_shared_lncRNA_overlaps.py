#python ./summarize_shared_lncRNA_overlaps.py --input-dir lncrna/lncRNA_shared_regions_mmseqs --output lncrna/lncRNA_shared_regions_summary/shared_lncRNA_overlap_summary.tsv
import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median


@dataclass
class MemberAccumulator:
    transcript_ids: set[str] = field(default_factory=set)
    gene_ids: set[str] = field(default_factory=set)
    gene_symbols: Counter[str] = field(default_factory=Counter)
    lengths: list[int] = field(default_factory=list)
    total_bp: int = 0

    def add(self, transcript_header: str, gene_id: str, length: int) -> None:
        transcript_id, parsed_gene_id, _transcript_name, gene_symbol, _transcript_len = parse_gencode_header(
            transcript_header
        )
        self.transcript_ids.add(transcript_id)
        self.gene_ids.add(gene_id or parsed_gene_id)
        if gene_symbol:
            self.gene_symbols[gene_symbol] += 1
        self.lengths.append(length)
        self.total_bp += length


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize shared local lncRNA region clusters from "
            "lncRNA_shared_regions_mmseqs/local_regions outputs."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("lncrna/lncRNA_shared_regions_mmseqs"),
        help="MMseqs shared-region output directory. Default: lncrna/lncRNA_shared_regions_mmseqs",
    )
    parser.add_argument(
        "--clusters",
        type=Path,
        help="Optional explicit path to shared_region_clusters.tsv.",
    )
    parser.add_argument(
        "--members",
        type=Path,
        help="Optional explicit path to shared_region_members.bed.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("lncrna/lncRNA_shared_regions_summary/shared_lncRNA_overlap_summary.tsv"),
        help="Output enriched cluster summary TSV.",
    )
    parser.add_argument(
        "--top-genes",
        type=int,
        default=10,
        help="Number of most frequent gene symbols to include per shared-region cluster.",
    )
    return parser.parse_args()


def resolve_inputs(args: argparse.Namespace) -> tuple[Path, Path]:
    clusters = args.clusters or args.input_dir / "local_regions" / "shared_region_clusters.tsv"
    members = args.members or args.input_dir / "local_regions" / "shared_region_members.bed"
    if not clusters.is_file():
        raise FileNotFoundError(f"Cluster summary not found: {clusters}")
    if not members.is_file():
        raise FileNotFoundError(f"Member interval file not found: {members}")
    return clusters, members


def parse_gencode_header(header: str) -> tuple[str, str, str, str, int | None]:
    fields = header.split("|")
    transcript_id = fields[0] if len(fields) > 0 else header
    gene_id = fields[1] if len(fields) > 1 else ""
    transcript_name = fields[4] if len(fields) > 4 else ""
    gene_symbol = fields[5] if len(fields) > 5 else ""
    transcript_len = None
    if len(fields) > 6 and fields[6].isdigit():
        transcript_len = int(fields[6])
    return transcript_id, gene_id, transcript_name, gene_symbol, transcript_len


def read_cluster_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "region_cluster_id",
            "n_raw_intervals",
            "n_merged_member_intervals",
            "n_transcripts",
            "n_genes",
            "total_member_bp",
            "representative_transcript",
            "representative_start",
            "representative_end",
            "representative_len",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
        return list(reader)


def read_member_stats(path: Path) -> dict[str, MemberAccumulator]:
    stats: dict[str, MemberAccumulator] = defaultdict(MemberAccumulator)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"region_cluster_id", "transcript_id", "gene_id", "length"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")

        for row in reader:
            length = int(row["length"])
            stats[row["region_cluster_id"]].add(row["transcript_id"], row["gene_id"], length)
    return stats


def format_top_gene_symbols(counter: Counter[str], top_n: int) -> str:
    if not counter:
        return "NA"
    return ";".join(f"{gene}:{count}" for gene, count in counter.most_common(top_n))


def output_rows(
    cluster_rows: list[dict[str, str]],
    member_stats: dict[str, MemberAccumulator],
    top_genes: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for cluster in cluster_rows:
        cluster_id = cluster["region_cluster_id"]
        stats = member_stats.get(cluster_id, MemberAccumulator())
        lengths = stats.lengths
        rep_tx, rep_gene, rep_tx_name, rep_gene_symbol, rep_tx_len = parse_gencode_header(
            cluster["representative_transcript"]
        )

        mean_len = stats.total_bp / len(lengths) if lengths else 0
        row = {
            "region_cluster_id": cluster_id,
            "n_raw_intervals": cluster["n_raw_intervals"],
            "n_member_intervals": str(len(lengths)),
            "n_transcripts": str(len(stats.transcript_ids)),
            "n_genes": str(sum(1 for gene_id in stats.gene_ids if gene_id)),
            "n_gene_symbols": str(len(stats.gene_symbols)),
            "total_member_bp": str(stats.total_bp),
            "mean_member_len": f"{mean_len:.2f}",
            "median_member_len": f"{median(lengths):.2f}" if lengths else "0.00",
            "min_member_len": str(min(lengths)) if lengths else "0",
            "max_member_len": str(max(lengths)) if lengths else "0",
            "representative_transcript_id": rep_tx,
            "representative_gene_id": rep_gene,
            "representative_transcript_name": rep_tx_name,
            "representative_gene_symbol": rep_gene_symbol,
            "representative_transcript_length": str(rep_tx_len or ""),
            "representative_start": cluster["representative_start"],
            "representative_end": cluster["representative_end"],
            "representative_region_len": cluster["representative_len"],
            "top_gene_symbols": format_top_gene_symbols(stats.gene_symbols, top_genes),
            "source_n_member_intervals": cluster["n_merged_member_intervals"],
            "source_n_transcripts": cluster["n_transcripts"],
            "source_n_genes": cluster["n_genes"],
            "source_total_member_bp": cluster["total_member_bp"],
            "representative_transcript_header": cluster["representative_transcript"],
        }
        rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError("No rows to write.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    clusters_path, members_path = resolve_inputs(args)
    cluster_rows = read_cluster_rows(clusters_path)
    member_stats = read_member_stats(members_path)
    rows = output_rows(cluster_rows, member_stats, args.top_genes)
    write_tsv(args.output, rows)
    print(f"Wrote {len(rows)} shared-region cluster summaries to: {args.output}")


if __name__ == "__main__":
    main()
