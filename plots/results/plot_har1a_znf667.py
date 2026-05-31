#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


DEFAULT_KMER = "kmers_comparison/qc/kmer_transcript_pass_counts_by_k.tsv"
DEFAULT_ALIGNMENT = "FragIncRNA/alignment_results/lncRNA_transcript_pass_counts.tsv"
DEFAULT_OUTDIR = "HOR1_ZNF667_results"
GENES = ["HAR1A", "ZNF667-AS1"]

STATUS_ORDER = ["UNMAPPED", "LOW_QC", "PASS_BASIC", "PASS_STRICT"]
STATUS_COLORS = {
    "UNMAPPED": "#d7301f",
    "LOW_QC": "#fdae61",
    "PASS_BASIC": "#7bccc4",
    "PASS_STRICT": "#2b8cbe",
}
GENE_COLORS = {
    "HAR1A": "#2b8cbe",
    "ZNF667-AS1": "#756bb1",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a four-panel paper figure comparing HAR1A and ZNF667-AS1 "
            "alignment and k-mer evidence across primate genome assemblies."
        )
    )
    parser.add_argument(
        "--kmer",
        default=DEFAULT_KMER,
        help="Input kmer_transcript_pass_counts_by_k.tsv file.",
    )
    parser.add_argument(
        "--alignment",
        default=DEFAULT_ALIGNMENT,
        help="Input lncRNA_transcript_pass_counts.tsv file.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=DEFAULT_OUTDIR,
        help="Output directory for figures and summary tables.",
    )
    return parser.parse_args()


def short_genome_name(name: str) -> str:
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
        raise SystemExit(f"Missing alignment column(s): {', '.join(missing)}")
    start = df.columns.get_loc("gene_name") + 1
    end = df.columns.get_loc("strict_pass_count")
    cols = list(df.columns[start:end])
    if not cols:
        raise SystemExit("No genome status columns found in alignment table.")
    return cols


def load_alignment(path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(path, sep="\t")
    status_cols = status_columns(df)
    df = df[df["gene_name"].isin(GENES)].copy()
    if df.empty:
        raise SystemExit(f"No selected genes found in {path}")
    df["gene_name"] = pd.Categorical(df["gene_name"], categories=GENES, ordered=True)
    df = df.sort_values(["gene_name", "transcript_name"]).reset_index(drop=True)
    df["gene_name"] = df["gene_name"].astype(str)
    return df, status_cols


def load_kmer(path: Path) -> pd.DataFrame:
    usecols = [
        "k_size",
        "transcript_name",
        "gene_name",
        "gene_id",
        "transcript_id",
        "transcript_length",
        "total_kmers",
        "genome_count",
        "strict_pass_count",
        "basic_pass_count",
        "min_matched_kmer_ratio",
        "mean_matched_kmer_ratio",
        "max_matched_kmer_ratio",
    ]
    chunks = []
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=250_000):
        sub = chunk[chunk["gene_name"].isin(GENES)]
        if not sub.empty:
            chunks.append(sub.copy())
    if not chunks:
        raise SystemExit(f"No selected genes found in {path}")
    df = pd.concat(chunks, ignore_index=True)
    df["gene_name"] = pd.Categorical(df["gene_name"], categories=GENES, ordered=True)
    df = df.sort_values(["gene_name", "transcript_name", "k_size"]).reset_index(drop=True)
    df["gene_name"] = df["gene_name"].astype(str)
    return df


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left",
    )


def style_axis(ax: plt.Axes, xgrid: bool = False, ygrid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if ygrid:
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.7)
    if xgrid:
        ax.grid(axis="x", color="#e5e5e5", linewidth=0.7)
    ax.set_axisbelow(True)


def plot_alignment_heatmap(
    ax: plt.Axes,
    alignment_df: pd.DataFrame,
    status_cols: list[str],
) -> None:
    rows = []
    labels = []
    annotations = []
    for gene in GENES:
        sub = alignment_df[alignment_df["gene_name"] == gene]
        total = int(sub["transcript_id"].nunique())
        values = []
        notes = []
        for col in status_cols:
            strict = int(sub[col].eq("PASS_STRICT").sum())
            values.append(strict / total if total else 0)
            notes.append(f"{strict}/{total}")
        rows.append(values)
        labels.append(gene)
        annotations.append(notes)

    matrix = np.array(rows)
    cmap = LinearSegmentedColormap.from_list(
        "strict_fraction",
        ["#fff7bc", "#7fcdbb", "#2b8cbe"],
    )
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)
    ax.set_title("Strict alignment support per genome assembly", pad=10)
    ax.set_xticks(range(len(status_cols)))
    ax.set_xticklabels([short_genome_name(c) for c in status_cols], rotation=45, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontweight="bold")
    ax.tick_params(axis="both", length=0)
    ax.set_xlabel("Genome assembly")
    ax.set_ylabel("Gene")
    for i, row in enumerate(annotations):
        for j, text in enumerate(row):
            value = matrix[i, j]
            color = "white" if value > 0.65 else "#252525"
            ax.text(j, i, text, ha="center", va="center", fontsize=8, color=color)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Fraction strict")
    add_panel_label(ax, "A")


def transcript_method_summary(
    alignment_df: pd.DataFrame,
    kmer_df: pd.DataFrame,
) -> pd.DataFrame:
    kmer_summary = (
        kmer_df.groupby(["gene_name", "transcript_id", "transcript_name"], observed=True)
        .agg(
            kmer_strict_pass_median=("strict_pass_count", "median"),
            kmer_strict_pass_min=("strict_pass_count", "min"),
            kmer_mean_ratio_mean=("mean_matched_kmer_ratio", "mean"),
            transcript_length=("transcript_length", "first"),
        )
        .reset_index()
    )
    align_summary = alignment_df[
        [
            "gene_name",
            "transcript_id",
            "transcript_name",
            "strict_pass_count",
            "basic_pass_count",
            "genome_count",
        ]
    ].rename(
        columns={
            "strict_pass_count": "alignment_strict_pass_count",
            "basic_pass_count": "alignment_basic_pass_count",
            "genome_count": "alignment_genome_count",
        }
    )
    merged = align_summary.merge(
        kmer_summary,
        on=["gene_name", "transcript_id", "transcript_name"],
        how="left",
    )
    merged["gene_name"] = pd.Categorical(merged["gene_name"], categories=GENES, ordered=True)
    return merged.sort_values(["gene_name", "transcript_name"]).reset_index(drop=True)


def plot_alignment_vs_kmer_bars(ax: plt.Axes, summary: pd.DataFrame) -> None:
    rows = []
    for gene in GENES:
        sub = summary[summary["gene_name"] == gene]
        total = int(sub["transcript_id"].nunique())
        rows.extend(
            [
                {
                    "gene_name": gene,
                    "method": "Alignment strict",
                    "full": int(
                        sub["alignment_strict_pass_count"]
                        .eq(sub["alignment_genome_count"].max())
                        .sum()
                    ),
                    "partial": int(
                        sub["alignment_strict_pass_count"]
                        .lt(sub["alignment_genome_count"].max())
                        .sum()
                    ),
                    "total": total,
                },
                {
                    "gene_name": gene,
                    "method": "K-mer strict",
                    "full": int(
                        sub["kmer_strict_pass_min"]
                        .eq(sub["alignment_genome_count"].max())
                        .sum()
                    ),
                    "partial": int(
                        sub["kmer_strict_pass_min"]
                        .lt(sub["alignment_genome_count"].max())
                        .sum()
                    ),
                    "total": total,
                },
            ]
        )
    plot_df = pd.DataFrame(rows)
    ylabels = [f"{row.gene_name}\n{row.method}" for row in plot_df.itertuples()]
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["full"], color="#2b8cbe", label="Strict in all genomes")
    ax.barh(
        y,
        plot_df["partial"],
        left=plot_df["full"],
        color="#fdae61",
        label="Not strict in all genomes",
    )
    for i, row in plot_df.reset_index(drop=True).iterrows():
        ax.text(
            row["total"] + 0.35,
            i,
            f'{int(row["full"])}/{int(row["total"])}',
            va="center",
            ha="left",
            fontsize=9,
        )
    ax.set_title("Transcripts with full strict support", pad=10)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels)
    ax.set_xlim(0, int(plot_df["total"].max()) + 4)
    ax.set_xlabel("Transcript count")
    ax.invert_yaxis()
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    style_axis(ax, xgrid=True, ygrid=False)
    add_panel_label(ax, "B")


def plot_kmer_robustness(ax: plt.Axes, kmer_df: pd.DataFrame) -> None:
    for gene in GENES:
        sub = kmer_df[kmer_df["gene_name"] == gene]
        for _, tx in sub.groupby("transcript_name", observed=True):
            ax.plot(
                tx["k_size"],
                tx["mean_matched_kmer_ratio"],
                color=GENE_COLORS[gene],
                alpha=0.18,
                linewidth=0.9,
            )
        grouped = (
            sub.groupby("k_size", observed=True)["mean_matched_kmer_ratio"]
            .agg(["mean", "min", "max"])
            .reset_index()
        )
        ax.fill_between(
            grouped["k_size"].to_numpy(),
            grouped["min"].to_numpy(),
            grouped["max"].to_numpy(),
            color=GENE_COLORS[gene],
            alpha=0.14,
            linewidth=0,
        )
        ax.plot(
            grouped["k_size"],
            grouped["mean"],
            color=GENE_COLORS[gene],
            linewidth=2.3,
            marker="o",
            markersize=4,
            label=gene,
        )
    ax.set_title("K-mer matched ratio across k sizes", pad=10)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Mean matched k-mer ratio")
    ax.set_ylim(0.985, 1.002)
    ax.set_xticks(sorted(kmer_df["k_size"].unique()))
    ax.legend(frameon=False, fontsize=9, loc="lower left")
    style_axis(ax)
    add_panel_label(ax, "C")


def status_count_table(
    alignment_df: pd.DataFrame,
    status_cols: list[str],
) -> pd.DataFrame:
    rows = []
    for gene in GENES:
        sub = alignment_df[alignment_df["gene_name"] == gene]
        for col in status_cols:
            counts = sub[col].value_counts()
            for status in STATUS_ORDER:
                rows.append(
                    {
                        "gene_name": gene,
                        "genome": short_genome_name(col),
                        "status": status,
                        "transcript_count": int(counts.get(status, 0)),
                    }
                )
    return pd.DataFrame(rows)


def plot_status_composition(ax: plt.Axes, status_counts: pd.DataFrame) -> None:
    zsub = status_counts[status_counts["gene_name"] == "ZNF667-AS1"]
    genomes = zsub["genome"].drop_duplicates().tolist()
    y = np.arange(len(genomes))
    left = np.zeros(len(genomes))
    for status in STATUS_ORDER:
        values = []
        for genome in genomes:
            value = zsub[
                (zsub["genome"] == genome)
                & (zsub["status"] == status)
            ]["transcript_count"].iloc[0]
            values.append(value)
        ax.barh(
            y,
            values,
            left=left,
            color=STATUS_COLORS[status],
            label=status,
            edgecolor="white",
            linewidth=0.5,
        )
        left += np.array(values)

    ax.set_title("ZNF667-AS1 alignment-status composition", pad=34)
    ax.set_xlabel("Transcript count")
    ax.set_yticks(y)
    ax.set_yticklabels(genomes)
    ax.invert_yaxis()
    ax.set_ylim(len(genomes) - 0.5, -1.35)
    ax.legend(
        frameon=False,
        fontsize=10,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        handlelength=1.3,
        columnspacing=1.2,
        ncol=4,
    )
    style_axis(ax, xgrid=True, ygrid=False)
    add_panel_label(ax, "D")


def save_outputs(
    outdir: Path,
    alignment_df: pd.DataFrame,
    kmer_df: pd.DataFrame,
    method_summary: pd.DataFrame,
    status_counts: pd.DataFrame,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    alignment_df.to_csv(outdir / "HAR1A_ZNF667_alignment_selected.tsv", sep="\t", index=False)
    kmer_df.to_csv(outdir / "HAR1A_ZNF667_kmer_selected.tsv", sep="\t", index=False)
    method_summary.to_csv(
        outdir / "HAR1A_ZNF667_transcript_method_summary.tsv",
        sep="\t",
        index=False,
    )
    status_counts.to_csv(
        outdir / "HAR1A_ZNF667_alignment_status_counts.tsv",
        sep="\t",
        index=False,
    )


def build_figure(
    alignment_df: pd.DataFrame,
    status_cols: list[str],
    kmer_df: pd.DataFrame,
    method_summary: pd.DataFrame,
    status_counts: pd.DataFrame,
) -> plt.Figure:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(21, 18), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.05, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_alignment_heatmap(ax_a, alignment_df, status_cols)
    plot_alignment_vs_kmer_bars(ax_b, method_summary)
    plot_kmer_robustness(ax_c, kmer_df)
    plot_status_composition(ax_d, status_counts)

    fig.suptitle(
        "HAR1A and ZNF667-AS1 alignment and k-mer comparison",
        fontsize=20,
        fontweight="bold",
        y=1.02,
    )
    return fig


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    alignment_df, status_cols = load_alignment(Path(args.alignment))
    kmer_df = load_kmer(Path(args.kmer))
    method_summary = transcript_method_summary(alignment_df, kmer_df)
    counts = status_count_table(alignment_df, status_cols)

    save_outputs(outdir, alignment_df, kmer_df, method_summary, counts)
    fig = build_figure(alignment_df, status_cols, kmer_df, method_summary, counts)
    prefix = outdir / "HAR1A_ZNF667_alignment_kmer_four_panel"
    fig.savefig(prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {prefix.with_suffix('.png')}")
    print(f"Wrote {prefix.with_suffix('.pdf')}")
    print(f"Wrote summary tables to {outdir}")


if __name__ == "__main__":
    main()
