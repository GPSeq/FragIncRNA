#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import PercentFormatter


DEFAULT_RESULTS = Path("kmer_results/multi_exact_mapper")
DEFAULT_OUTDIR = DEFAULT_RESULTS / "plots"
DEFAULT_GENE = "HAR1A"

GENOME_LABELS = {
    "GCF_002880755.1_Clint_PTRv2_genomic": "Clintptrv2",
    "GCA_054883195.1_H9_T2T.hap1_genomic": "H9T2Thap1",
    "GCA_054883265.1_H9_T2T.hap2_genomic": "H9T2Thap2",
    "GCF_013052645.1_Mhudiblu_PPA_v0_genomic": "MhudibluPPAv0",
    "GCF_003339765.1_Mmul_10_genomic": "Mmul10",
    "GCF_049350105.2_T2T-MMU8v2.0_genomic": "T2TMMU8v20",
    "GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic": "mPanPan1v20pri",
    "GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic": "mPanTro3v20pri",
    "GCA_018503275.2_NA19240_mat_hprc_f2_genomic": "matHPRCF2",
    "GCF_000258655.2_panpan1.1_genomic": "panpan11",
    "GCA_018503265.2_NA19240_pat_hprc_f2_genomic": "patHPRCf2",
}

GENOME_ORDER = [
    "Clintptrv2",
    "H9T2Thap1",
    "H9T2Thap2",
    "MhudibluPPAv0",
    "Mmul10",
    "T2TMMU8v20",
    "mPanPan1v20pri",
    "mPanTro3v20pri",
    "matHPRCF2",
    "panpan11",
    "patHPRCf2",
]

COLORS = [
    "#0072B2",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#E69F00",
    "#56B4E9",
    "#7F3C8D",
    "#11A579",
    "#3969AC",
    "#F2B701",
    "#E73F74",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot HAR1A exact k-mer recovery across genomes and k-mer sizes."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS,
        help="Directory containing results_XX summary folders.",
    )
    parser.add_argument(
        "--gene",
        default=DEFAULT_GENE,
        help="Safe gene prefix used in exact-mapper filenames.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory for figures and combined summary TSV.",
    )
    parser.add_argument(
        "--prefix",
        default="HAR1A_multi_exact_mapper",
        help="Output filename prefix.",
    )
    return parser.parse_args()


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.labelsize": 9.5,
            "axes.titlesize": 10.5,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "figure.dpi": 160,
            "savefig.dpi": 450,
        }
    )


def genome_from_path(path: Path, gene: str) -> str:
    pattern = rf"^{re.escape(gene)}_(.+)_k\d+_kmer_locations\.summary\.tsv$"
    match = re.match(pattern, path.name)
    if not match:
        raise ValueError(f"Cannot parse genome name from: {path}")
    genome = match.group(1)
    return GENOME_LABELS.get(genome, genome)


def load_summaries(results_dir: Path, gene: str) -> pd.DataFrame:
    paths = sorted(results_dir.glob("results_*/*_kmer_locations.summary.tsv"))
    if not paths:
        raise FileNotFoundError(f"No summary TSV files found under {results_dir}")

    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep="\t")
        frame["source_file"] = str(path)
        frame["genome"] = genome_from_path(path, gene)
        frames.append(frame)

    df = pd.concat(frames, ignore_index=True)
    required = {
        "transcript_name",
        "transcript_id",
        "k_size",
        "total_kmer_windows",
        "matched_kmer_windows",
        "matched_kmer_ratio",
        "unique_kmers",
        "unique_kmers_found",
        "unique_kmer_ratio",
        "reference_hit_sum",
        "genome",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing summary column(s): {', '.join(missing)}")

    numeric_cols = [
        "k_size",
        "total_kmer_windows",
        "matched_kmer_windows",
        "matched_kmer_ratio",
        "unique_kmers",
        "unique_kmers_found",
        "unique_kmer_ratio",
        "reference_hit_sum",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col])

    df["genome"] = pd.Categorical(df["genome"], categories=GENOME_ORDER, ordered=True)
    return df.sort_values(["genome", "transcript_name", "k_size"]).reset_index(drop=True)


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


def draw_heatmap(ax: plt.Axes, df: pd.DataFrame) -> None:
    heat = (
        df.groupby(["genome", "k_size"], observed=True)["matched_kmer_ratio"]
        .mean()
        .unstack("k_size")
        .reindex(GENOME_ORDER)
    )
    cmap = LinearSegmentedColormap.from_list(
        "exact_ratio", ["#f7f7f7", "#c7e9b4", "#41ab5d", "#005a32"]
    )
    im = ax.imshow(heat.to_numpy(), aspect="auto", vmin=0, vmax=1, cmap=cmap)
    ax.set_title("Mean exact k-mer recovery")
    ax.set_xlabel("k-mer size")
    ax.set_xticks(np.arange(len(heat.columns)))
    ax.set_xticklabels(heat.columns.astype(int))
    ax.set_yticks(np.arange(len(heat.index)))
    ax.set_yticklabels(heat.index)

    for y in range(heat.shape[0]):
        for x in range(heat.shape[1]):
            value = heat.iat[y, x]
            if pd.notna(value):
                color = "white" if value >= 0.58 else "#252525"
                ax.text(
                    x,
                    y,
                    f"{value:.0%}",
                    ha="center",
                    va="center",
                    fontsize=6.8,
                    color=color,
                )
    ax.tick_params(length=0)
    add_panel_label(ax, "A")
    return im


def draw_transcript_lines(ax: plt.Axes, df: pd.DataFrame, transcript: str, label: str) -> None:
    sub = df[df["transcript_name"].eq(transcript)].copy()
    for color, genome in zip(COLORS, GENOME_ORDER, strict=False):
        line = sub[sub["genome"].eq(genome)].sort_values("k_size")
        if line.empty:
            continue
        ax.plot(
            line["k_size"],
            line["matched_kmer_ratio"],
            marker="o",
            markersize=3.2,
            linewidth=1.2,
            color=color,
            label=genome,
        )

    ax.set_title(transcript)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Matched k-mer ratio")
    ax.set_ylim(-0.02, 1.02)
    ax.yaxis.set_major_formatter(PercentFormatter(1.0))
    ax.grid(axis="y", color="#e8e8e8", linewidth=0.8)
    ax.set_axisbelow(True)
    add_panel_label(ax, label)


def build_figure(df: pd.DataFrame) -> plt.Figure:
    transcripts = sorted(df["transcript_name"].dropna().unique())
    if len(transcripts) != 3:
        raise ValueError(
            f"Expected 3 HAR1A transcripts for 4-panel plot, found {len(transcripts)}"
        )

    fig = plt.figure(figsize=(13.5, 9.5), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, width_ratios=[1.35, 1], height_ratios=[1, 1, 1])
    ax_heat = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])
    ax_d = fig.add_subplot(grid[2, 1])

    im = draw_heatmap(ax_heat, df)
    draw_transcript_lines(ax_b, df, transcripts[0], "B")
    draw_transcript_lines(ax_c, df, transcripts[1], "C")
    draw_transcript_lines(ax_d, df, transcripts[2], "D")

    handles, labels = ax_b.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=4,
        fontsize=7.5,
        handlelength=1.6,
        columnspacing=1.1,
    )
    fig.colorbar(
        im,
        ax=ax_heat,
        fraction=0.035,
        pad=0.025,
        label="Mean matched k-mer ratio",
        format=PercentFormatter(1.0),
    )
    fig.suptitle("HAR1A exact k-mer recovery across primate genomes", fontsize=14)
    return fig


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    args.outdir.mkdir(parents=True, exist_ok=True)

    summary = load_summaries(args.results_dir, args.gene)
    table_path = args.outdir / f"{args.prefix}_summary.tsv"
    png_path = args.outdir / f"{args.prefix}_4panel.png"
    pdf_path = args.outdir / f"{args.prefix}_4panel.pdf"

    summary.to_csv(table_path, sep="\t", index=False)
    fig = build_figure(summary)
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {table_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
