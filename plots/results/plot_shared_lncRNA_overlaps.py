#python  ./plot_shared_lncRNA_overlaps.py --input-summary lncrna/lncRNA_shared_regions_summary/shared_lncRNA_overlap_summary.tsv --output-dir lncrna/lncRNA_shared_regions_summary --top-n 30

import argparse
import csv
import math
import os
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a publication-ready multi-panel plot from shared_lncRNA_overlap_summary.tsv."
    )
    parser.add_argument(
        "--input-summary",
        type=Path,
        default=Path("lncrna/lncRNA_shared_regions_summary/shared_lncRNA_overlap_summary.tsv"),
        help="Summary TSV produced by summarize_shared_lncRNA_overlaps.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lncrna/lncRNA_shared_regions_summary"),
        help="Directory where the plot files will be written.",
    )
    parser.add_argument(
        "--prefix",
        default="shared_lncRNA_overlap_figure",
        help="Output filename prefix. PNG and PDF files are written.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=30,
        help="Number of largest shared-region clusters to show in the top-cluster panel.",
    )
    return parser.parse_args()


def require_matplotlib(output_dir: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import FuncFormatter
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for plotting. Install it with: pip install matplotlib"
        ) from exc
    return plt, LogNorm, FuncFormatter


def read_summary(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Summary TSV not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "region_cluster_id",
            "n_member_intervals",
            "n_transcripts",
            "n_genes",
            "total_member_bp",
            "representative_region_len",
            "representative_gene_symbol",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No data rows found in: {path}")
    return rows


def as_int(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    if value in {"", "NA"}:
        return 0
    return int(float(value))


def human_number(value: float, _position: int | None = None) -> str:
    if value == 0:
        return "0"
    abs_value = abs(value)
    for threshold, suffix in ((1_000_000_000, "B"), (1_000_000, "M"), (1_000, "K")):
        if abs_value >= threshold:
            return f"{value / threshold:.1f}{suffix}"
    return f"{value:.0f}"


def cluster_size_bins(values: list[int]) -> tuple[list[str], list[int]]:
    bins = [
        ("2-5", 2, 5),
        ("6-20", 6, 20),
        ("21-100", 21, 100),
        ("101-1,000", 101, 1_000),
        (">1,000", 1_001, math.inf),
    ]
    labels = []
    counts = []
    for label, lower, upper in bins:
        labels.append(label)
        counts.append(sum(1 for value in values if lower <= value <= upper))
    return labels, counts


def short_cluster_label(row: dict[str, str]) -> str:
    gene_symbol = row.get("representative_gene_symbol", "")
    if gene_symbol and gene_symbol != "NA":
        return f"{row['region_cluster_id']} ({gene_symbol})"
    return row["region_cluster_id"]


def save_figure(fig, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")


def make_plot(rows: list[dict[str, str]], output_dir: Path, prefix: str, top_n: int) -> None:
    plt, LogNorm, FuncFormatter = require_matplotlib(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 150,
        }
    )

    n_transcripts = [max(as_int(row, "n_transcripts"), 1) for row in rows]
    n_genes = [max(as_int(row, "n_genes"), 1) for row in rows]
    total_bp = [max(as_int(row, "total_member_bp"), 1) for row in rows]
    rep_lens = [max(as_int(row, "representative_region_len"), 1) for row in rows]

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5))
    fig.suptitle("Shared local sequence regions among lncRNA transcripts", fontsize=13, y=0.995)

    ax = axes[0][0]
    bin_labels, bin_counts = cluster_size_bins(n_transcripts)
    ax.bar(bin_labels, bin_counts, color="#4477AA", edgecolor="#263C59", linewidth=0.5)
    ax.set_ylabel("Shared-region clusters")
    ax.set_xlabel("Transcripts per cluster")
    ax.set_title("Cluster-size distribution")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    for index, value in enumerate(bin_counts):
        ax.text(index, value, f"{value:,}", ha="center", va="bottom", fontsize=8)
    ax.text(-0.13, 1.08, "a", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    ax = axes[0][1]
    ranked = sorted(n_transcripts, reverse=True)
    ranks = list(range(1, len(ranked) + 1))
    ax.plot(ranks, ranked, color="#228833", linewidth=1.8)
    ax.scatter(ranks[:10], ranked[:10], color="#228833", s=12, zorder=3)
    ax.set_yscale("log")
    ax.set_xlabel("Shared-region cluster rank")
    ax.set_ylabel("Transcripts per cluster")
    ax.set_title("Ranked shared-region membership")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.text(-0.13, 1.08, "b", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    ax = axes[1][0]
    max_bp = max(total_bp)
    sizes = [25 + 260 * (math.log10(value + 1) / math.log10(max_bp + 1)) for value in total_bp]
    scatter = ax.scatter(
        n_transcripts,
        n_genes,
        c=total_bp,
        s=sizes,
        cmap="viridis",
        norm=LogNorm(vmin=max(min(total_bp), 1), vmax=max_bp),
        alpha=0.72,
        edgecolors="#222222",
        linewidths=0.25,
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Transcripts per cluster")
    ax.set_ylabel("Genes per cluster")
    ax.set_title("Transcript and gene sharing per region")
    ax.xaxis.set_major_formatter(FuncFormatter(human_number))
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Total member bp")
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.text(-0.13, 1.08, "c", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    ax = axes[1][1]
    top_rows = sorted(rows, key=lambda row: as_int(row, "n_transcripts"), reverse=True)[:top_n]
    labels = [short_cluster_label(row) for row in top_rows]
    values = [as_int(row, "n_transcripts") for row in top_rows]
    y_positions = list(range(len(top_rows)))
    ax.barh(y_positions, values, color="#CC6677", edgecolor="#66333B", linewidth=0.45)
    ax.set_yticks(y_positions, labels=labels)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Transcripts per cluster")
    ax.set_title(f"Top {len(top_rows)} largest shared-region clusters")
    ax.xaxis.set_major_formatter(FuncFormatter(human_number))
    for y_position, value in zip(y_positions, values):
        ax.text(value * 1.05, y_position, f"{value:,}", va="center", fontsize=7)
    ax.text(-0.13, 1.08, "d", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    fig.tight_layout(rect=(0, 0, 1, 0.975))
    save_figure(fig, output_dir / prefix)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rows = read_summary(args.input_summary)
    make_plot(rows, args.output_dir, args.prefix, args.top_n)
    print(f"Wrote plot files to: {args.output_dir}")


if __name__ == "__main__":
    main()
