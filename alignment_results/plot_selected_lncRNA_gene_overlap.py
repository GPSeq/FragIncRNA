
import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT = (
    "/mnt/d/primates_bmc/bam_comparison/qc/shared_lncRNA_transcripts/"
    "lncRNA_transcript_pass_counts.tsv"
)

GENE_ALIASES = {
    "MALAT": "MALAT1",
    "NEAT": "NEAT1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create one multi-panel figure describing alignment overlap for "
            "selected lncRNA gene names."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        default=DEFAULT_INPUT,
        help="Input lncRNA_transcript_pass_counts.tsv file.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for plots and summary tables.",
    )
    parser.add_argument(
        "--genes",
        default="HOTAIR,MALAT,XIST,NEAT,MIAT",
        help=(
            "Comma-separated gene names to plot. MALAT is treated as MALAT1 "
            "and NEAT as NEAT1 by default."
        ),
    )
    return parser.parse_args()


def savefig(fig: plt.Figure, path_prefix: Path) -> None:
    fig.savefig(path_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def short_sample_name(name: str) -> str:
    return name.replace("human_lncRNA_vs_", "")


def status_columns(df: pd.DataFrame) -> list[str]:
    required = {
        "transcript_id",
        "gene_id",
        "transcript_name",
        "gene_name",
        "strict_pass_count",
        "basic_pass_count",
        "genome_count",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SystemExit(f"Missing required column(s): {', '.join(missing)}")
    start = df.columns.get_loc("gene_name") + 1
    end = df.columns.get_loc("strict_pass_count")
    return list(df.columns[start:end])


def canonical_gene_list(raw: str) -> list[str]:
    genes = []
    for item in raw.split(","):
        gene = item.strip()
        if not gene:
            continue
        gene = GENE_ALIASES.get(gene.upper(), gene)
        genes.append(gene)
    seen = set()
    ordered = []
    for gene in genes:
        if gene not in seen:
            ordered.append(gene)
            seen.add(gene)
    return ordered


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.6)
    ax.set_axisbelow(True)


def matrix_for_status(
    df: pd.DataFrame,
    genes: list[str],
    status_cols: list[str],
    mode: str,
) -> pd.DataFrame:
    rows = []
    for gene in genes:
        sub = df[df["gene_name"] == gene]
        row = []
        for col in status_cols:
            if mode == "strict":
                row.append(int(sub[col].eq("PASS_STRICT").sum()))
            elif mode == "basic":
                row.append(int(sub[col].isin(["PASS_STRICT", "PASS_BASIC"]).sum()))
            else:
                raise ValueError(mode)
        rows.append(row)
    return pd.DataFrame(rows, index=genes, columns=status_cols)


def annotate_heatmap(ax: plt.Axes, data: pd.DataFrame) -> None:
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            value = int(data.iat[i, j])
            ax.text(
                j,
                i,
                str(value),
                ha="center",
                va="center",
                fontsize=7,
                color="white" if value > data.values.max() * 0.55 else "#252525",
            )


def plot_count_heatmap(
    ax: plt.Axes,
    data: pd.DataFrame,
    title: str,
    cmap: str,
    panel_label: str,
) -> None:
    im = ax.imshow(data.to_numpy(), aspect="auto", cmap=cmap)
    annotate_heatmap(ax, data)
    ax.set_title(title)
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index)
    ax.set_xticks(range(len(data.columns)))
    ax.set_xticklabels([short_sample_name(c) for c in data.columns], rotation=60, ha="right")
    ax.tick_params(axis="both", length=0)
    add_panel_label(ax, panel_label)
    return im


def build_gene_summary(
    df: pd.DataFrame,
    genes: list[str],
    status_cols: list[str],
) -> pd.DataFrame:
    genome_count = int(df["genome_count"].max())
    rows = []
    for gene in genes:
        sub = df[df["gene_name"] == gene].copy()
        strict_all = sub["strict_pass_count"].eq(genome_count)
        basic_all = sub["basic_pass_count"].eq(genome_count)
        strict_genomes = [
            col for col in status_cols if sub[col].eq("PASS_STRICT").any()
        ]
        basic_genomes = [
            col for col in status_cols if sub[col].isin(["PASS_STRICT", "PASS_BASIC"]).any()
        ]
        rows.append(
            {
                "gene_name": gene,
                "associated_gene_ids": ",".join(sorted(sub["gene_id"].dropna().unique())),
                "total_transcripts": int(sub["transcript_id"].nunique()),
                "strict_shared_transcripts_all_genomes": int(strict_all.sum()),
                "basic_shared_transcripts_all_genomes": int(basic_all.sum()),
                "strict_detected_genomes": len(strict_genomes),
                "basic_detected_genomes": len(basic_genomes),
                "strict_shared_gene_all_genomes": len(strict_genomes) == len(status_cols),
                "basic_shared_gene_all_genomes": len(basic_genomes) == len(status_cols),
            }
        )
    return pd.DataFrame(rows)


def plot_transcript_bars(ax: plt.Axes, summary: pd.DataFrame) -> None:
    x = np.arange(len(summary))
    width = 0.25
    bars1 = ax.bar(
        x - width,
        summary["total_transcripts"],
        width,
        label="All transcripts",
        color="#636363",
    )
    bars2 = ax.bar(
        x,
        summary["strict_shared_transcripts_all_genomes"],
        width,
        label="Strict shared",
        color="#2b8cbe",
    )
    bars3 = ax.bar(
        x + width,
        summary["basic_shared_transcripts_all_genomes"],
        width,
        label="Basic shared",
        color="#7bccc4",
    )
    for bars in (bars1, bars2, bars3):
        ax.bar_label(bars, padding=2, fontsize=7, rotation=90)
    ax.set_title("Transcript counts per selected gene", pad=30)
    ax.set_ylabel("Transcripts")
    ax.set_xticks(x)
    ax.set_xticklabels(summary["gene_name"])
    ax.legend(
        frameon=False,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.13),
        ncol=3,
        handlelength=1.1,
        columnspacing=1.0,
    )
    add_panel_label(ax, "c")
    style_axis(ax)


def plot_gene_presence(ax: plt.Axes, summary: pd.DataFrame) -> None:
    data = summary[
        [
            "strict_detected_genomes",
            "basic_detected_genomes",
            "strict_shared_gene_all_genomes",
            "basic_shared_gene_all_genomes",
        ]
    ].copy()
    data["strict_shared_gene_all_genomes"] = data[
        "strict_shared_gene_all_genomes"
    ].astype(int)
    data["basic_shared_gene_all_genomes"] = data[
        "basic_shared_gene_all_genomes"
    ].astype(int)

    x = np.arange(len(summary))
    width = 0.36
    strict_bars = ax.bar(
        x - width / 2,
        data["strict_detected_genomes"],
        width,
        color="#2b8cbe",
        label="Strict detected genomes",
    )
    basic_bars = ax.bar(
        x + width / 2,
        data["basic_detected_genomes"],
        width,
        color="#7bccc4",
        label="Basic detected genomes",
    )
    ax.axhline(11, color="#252525", linewidth=0.8, linestyle="--")
    ax.bar_label(strict_bars, padding=2, fontsize=7)
    ax.bar_label(basic_bars, padding=2, fontsize=7)
    ax.set_ylim(0, 12)
    ax.set_title("Genome overlap for selected gene names", pad=30)
    ax.set_ylabel("Genomes with >=1 passing transcript")
    ax.set_xticks(x)
    ax.set_xticklabels(summary["gene_name"])
    ax.legend(
        frameon=False,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.13),
        ncol=2,
        handlelength=1.1,
        columnspacing=1.0,
    )
    add_panel_label(ax, "d")
    style_axis(ax)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = canonical_gene_list(args.genes)
    df = pd.read_csv(args.input, sep="\t", dtype=str)
    status_cols = status_columns(df)

    for col in ["strict_pass_count", "basic_pass_count", "genome_count"]:
        df[col] = pd.to_numeric(df[col], errors="raise")

    missing = [gene for gene in genes if gene not in set(df["gene_name"])]
    if missing:
        print("WARNING: missing gene_name value(s): " + ", ".join(missing))

    selected = df[df["gene_name"].isin(genes)].copy()
    selected.to_csv(outdir / "selected_lncRNA_gene_transcripts.tsv", sep="\t", index=False)

    summary = build_gene_summary(df, genes, status_cols)
    summary.to_csv(outdir / "selected_lncRNA_gene_summary.tsv", sep="\t", index=False)

    strict_counts = matrix_for_status(df, genes, status_cols, "strict")
    basic_counts = matrix_for_status(df, genes, status_cols, "basic")
    strict_counts.to_csv(outdir / "selected_lncRNA_strict_counts_by_genome.tsv", sep="\t")
    basic_counts.to_csv(outdir / "selected_lncRNA_basic_counts_by_genome.tsv", sep="\t")

    fig = plt.figure(figsize=(15, 10.8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], hspace=0.62, wspace=0.25)

    ax_a = fig.add_subplot(gs[0, 0])
    im_a = plot_count_heatmap(
        ax_a,
        strict_counts,
        "Strict QC transcripts per genome",
        "Blues",
        "a",
    )
    fig.colorbar(im_a, ax=ax_a, fraction=0.035, pad=0.02, label="Transcript count")

    ax_b = fig.add_subplot(gs[0, 1])
    im_b = plot_count_heatmap(
        ax_b,
        basic_counts,
        "Basic-or-strict QC transcripts per genome",
        "GnBu",
        "b",
    )
    fig.colorbar(im_b, ax=ax_b, fraction=0.035, pad=0.02, label="Transcript count")

    ax_c = fig.add_subplot(gs[1, 0])
    plot_transcript_bars(ax_c, summary)

    ax_d = fig.add_subplot(gs[1, 1])
    plot_gene_presence(ax_d, summary)

    savefig(fig, outdir / "selected_lncRNA_gene_overlap")

    print(f"Wrote selected gene plot and tables to: {outdir}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
