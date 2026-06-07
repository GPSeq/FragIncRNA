#!/usr/bin/env python3
"""Compare high-confidence alignment candidates with k-mer support."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, PercentFormatter


BLUE = "#2F6C9E"
TEAL = "#1B9E77"
GOLD = "#D89C28"
PURPLE = "#756BB1"
INK = "#252525"
LIGHT_GRID = "#E8E8E8"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare high-confidence alignment candidate orthologs with "
            "all-genome k-mer support across k sizes."
        )
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        default=Path(
            "alignment_results/high_confident_candidate_orthologs/"
            "high_confident_candidate_ortholog_pass_counts.tsv"
        ),
        help="High-confidence alignment pass-count table.",
    )
    parser.add_argument(
        "--kmer",
        type=Path,
        default=Path(
            "kmer_results/ibf_new/qc/kmer_transcript_pass_counts_by_k.tsv"
        ),
        help="K-mer transcript pass-count table.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("alignment_results/high_confident_candidate_orthologs/plots"),
        help="Output directory.",
    )
    return parser.parse_args()


def configure_matplotlib() -> None:
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
            "savefig.dpi": 450,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def require_columns(df: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = sorted(columns.difference(df.columns))
    if missing:
        raise SystemExit(f"Missing column(s) in {path}: {', '.join(missing)}")


def load_alignment(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    require_columns(
        df,
        {
            "transcript_id",
            "gene_id",
            "high_confident_pass_count",
            "strict_pass_count",
            "basic_pass_count",
            "genome_count",
        },
        path,
    )
    if df["transcript_id"].duplicated().any():
        raise SystemExit(f"Duplicate transcript IDs in {path}")
    return df


def load_kmer(path: Path) -> pd.DataFrame:
    columns = [
        "k_size",
        "transcript_id",
        "gene_id",
        "strict_pass_count",
        "basic_pass_count",
        "genome_count",
    ]
    df = pd.read_csv(path, sep="\t", usecols=columns)
    require_columns(df, set(columns), path)
    return df


def count_unique_genes(df: pd.DataFrame, mask: pd.Series) -> int:
    return int(df.loc[mask, "gene_id"].nunique())


def build_alignment_summary(alignment: pd.DataFrame) -> pd.DataFrame:
    genome_count = alignment["genome_count"]
    rows = []
    for label, column in [
        ("Basic detectability", "basic_pass_count"),
        ("Strict alignment", "strict_pass_count"),
        ("High-confidence candidates", "high_confident_pass_count"),
    ]:
        mask = alignment[column].eq(genome_count)
        rows.append(
            {
                "alignment_tier": label,
                "transcript_count": int(mask.sum()),
                "unique_gene_count": count_unique_genes(alignment, mask),
            }
        )
    return pd.DataFrame(rows)


def build_overlap_summary(
    alignment: pd.DataFrame,
    kmer: pd.DataFrame,
) -> pd.DataFrame:
    high = alignment[
        ["transcript_id", "gene_id", "high_confident_pass_count", "genome_count"]
    ].copy()
    high["alignment_high_all"] = high["high_confident_pass_count"].eq(
        high["genome_count"]
    )
    high = high.drop(columns=["high_confident_pass_count", "genome_count"])

    rows = []
    for k_size, group in kmer.groupby("k_size", sort=True):
        work = group.merge(
            high,
            on="transcript_id",
            how="inner",
            suffixes=("_kmer", "_alignment"),
            validate="one_to_one",
        )
        strict = work["strict_pass_count"].eq(work["genome_count"])
        basic = work["basic_pass_count"].eq(work["genome_count"])
        high_all = work["alignment_high_all"]
        both_strict = high_all & strict
        both_basic = high_all & basic
        rows.append(
            {
                "k_size": int(k_size),
                "alignment_high_confident_all": int(high_all.sum()),
                "kmer_strict_all": int(strict.sum()),
                "kmer_basic_all": int(basic.sum()),
                "high_and_kmer_strict": int(both_strict.sum()),
                "high_and_kmer_basic": int(both_basic.sum()),
                "high_only_vs_strict": int((high_all & ~strict).sum()),
                "kmer_strict_only": int((strict & ~high_all).sum()),
                "neither_high_nor_kmer_strict": int((~high_all & ~strict).sum()),
                "high_and_kmer_strict_unique_genes": int(
                    work.loc[both_strict, "gene_id_alignment"].nunique()
                ),
                "high_supported_by_kmer_strict_fraction": (
                    float(both_strict.sum() / high_all.sum())
                    if high_all.any()
                    else np.nan
                ),
                "high_supported_by_kmer_basic_fraction": (
                    float(both_basic.sum() / high_all.sum())
                    if high_all.any()
                    else np.nan
                ),
                "kmer_strict_supported_by_high_fraction": (
                    float(both_strict.sum() / strict.sum())
                    if strict.any()
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def style_axis(ax: plt.Axes, axis: str = "y") -> None:
    ax.grid(axis=axis, color=LIGHT_GRID, linewidth=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, color="#444444")


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=INK,
    )


def human_number(value: float, _position: int | None = None) -> str:
    if abs(value) >= 1_000:
        return f"{value / 1_000:.0f}K"
    return f"{value:.0f}"


def plot_alignment_tiers(ax: plt.Axes, summary: pd.DataFrame) -> None:
    x = np.arange(len(summary))
    bars = ax.bar(
        x,
        summary["transcript_count"],
        color=[TEAL, BLUE, PURPLE],
        width=0.68,
    )
    ax.set_xticks(x, summary["alignment_tier"])
    ax.set_ylabel("Transcripts shared across all genomes")
    ax.set_title("Alignment evidence tiers")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.margins(y=0.16)
    for bar, row in zip(bars, summary.itertuples(index=False)):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{row.transcript_count:,} transcripts\n{row.unique_gene_count:,} genes",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    style_axis(ax)
    add_panel_label(ax, "a")


def plot_counts_by_k(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    x = overlap["k_size"].to_numpy()
    ax.plot(
        x,
        overlap["kmer_basic_all"],
        color=TEAL,
        marker="o",
        linewidth=2,
        label="K-mer basic, all genomes",
    )
    ax.plot(
        x,
        overlap["kmer_strict_all"],
        color=GOLD,
        marker="s",
        linewidth=2,
        label="K-mer strict, all genomes",
    )
    ax.plot(
        x,
        overlap["high_and_kmer_strict"],
        color=PURPLE,
        marker="D",
        linewidth=2,
        label="High-confidence alignment ∩ k-mer strict",
    )
    ax.axhline(
        overlap["alignment_high_confident_all"].iloc[0],
        color=BLUE,
        linestyle="--",
        linewidth=1.8,
        label="High-confidence alignment, all genomes",
    )
    ax.set_xticks(x)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Transcripts")
    ax.set_title("Alignment and k-mer support")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.margins(y=0.1)
    ax.legend(loc="best", fontsize=8)
    style_axis(ax)
    add_panel_label(ax, "b")


def plot_overlap_composition(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    x = np.arange(len(overlap))
    both = overlap["high_and_kmer_strict"].to_numpy()
    high_only = overlap["high_only_vs_strict"].to_numpy()
    kmer_only = overlap["kmer_strict_only"].to_numpy()
    ax.bar(x, both, color=PURPLE, label="Both")
    ax.bar(x, high_only, bottom=both, color=BLUE, label="Alignment only")
    ax.bar(
        x,
        kmer_only,
        bottom=both + high_only,
        color=GOLD,
        label="K-mer only",
    )
    ax.set_xticks(x, overlap["k_size"].astype(str))
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Union of supported transcripts")
    ax.set_title("High-confidence alignment versus strict k-mers")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.legend(loc="best", fontsize=8)
    style_axis(ax)
    add_panel_label(ax, "c")


def plot_support_fraction(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    x = overlap["k_size"].to_numpy()
    ax.plot(
        x,
        overlap["high_supported_by_kmer_basic_fraction"],
        color=TEAL,
        marker="o",
        linewidth=2,
        label="Alignment candidates with basic k-mer support",
    )
    ax.plot(
        x,
        overlap["high_supported_by_kmer_strict_fraction"],
        color=PURPLE,
        marker="s",
        linewidth=2,
        label="Alignment candidates with strict k-mer support",
    )
    ax.plot(
        x,
        overlap["kmer_strict_supported_by_high_fraction"],
        color=GOLD,
        marker="D",
        linewidth=2,
        label="Strict k-mer set with high-confidence alignment",
    )
    ax.set_xticks(x)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Fraction")
    ax.set_title("Cross-method confirmation")
    ax.yaxis.set_major_formatter(PercentFormatter(xmax=1.0))
    ax.set_ylim(0, 1.03)
    ax.legend(loc="best", fontsize=8)
    style_axis(ax)
    add_panel_label(ax, "d")


def make_figure(
    alignment_summary: pd.DataFrame,
    overlap_summary: pd.DataFrame,
) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(14.2, 9.0))
    plot_alignment_tiers(axes[0, 0], alignment_summary)
    plot_counts_by_k(axes[0, 1], overlap_summary)
    plot_overlap_composition(axes[1, 0], overlap_summary)
    plot_support_fraction(axes[1, 1], overlap_summary)
    fig.suptitle(
        "High-confidence alignment candidates and k-mer support",
        fontsize=14,
        y=0.995,
    )
    fig.subplots_adjust(wspace=0.3, hspace=0.38, top=0.92, bottom=0.1)
    return fig


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    args.outdir.mkdir(parents=True, exist_ok=True)

    alignment = load_alignment(args.alignment)
    kmer = load_kmer(args.kmer)
    alignment_summary = build_alignment_summary(alignment)
    overlap_summary = build_overlap_summary(alignment, kmer)

    alignment_summary.to_csv(
        args.outdir / "alignment_evidence_tier_summary.tsv",
        sep="\t",
        index=False,
    )
    overlap_summary.to_csv(
        args.outdir / "high_confident_alignment_kmer_overlap_by_k.tsv",
        sep="\t",
        index=False,
    )

    prefix = args.outdir / "high_confident_alignment_vs_kmers"
    fig = make_figure(alignment_summary, overlap_summary)
    fig.savefig(prefix.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {prefix.with_suffix('.png')}")
    print(f"Wrote {prefix.with_suffix('.pdf')}")
    print(f"Wrote summary tables to {args.outdir}")


if __name__ == "__main__":
    main()
