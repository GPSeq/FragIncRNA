#!/usr/bin/env python3
"""
Example:
  python plots/results/plot_kmer_comparison_qc.py \
      --qc-dir ../kmers_comparison/qc \
      --outdir ../kmers_comparison/qc/paper_plots
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache")
os.environ.setdefault("ARROW_USER_SIMD_LEVEL", "NONE")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter, MaxNLocator


BLUE = "#2F6C9E"
TEAL = "#1B9E77"
GOLD = "#D89C28"
RED = "#C44E52"
INK = "#252525"
GRAY = "#6F6F6F"
LIGHT_GRID = "#E8E8E8"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-ready plots from summarize_kmer_comparison_qc.sh outputs."
        )
    )
    parser.add_argument(
        "--qc-dir",
        type=Path,
        default=Path("../kmers_comparison/qc"),
        help="Directory containing kmer_shared_summary_by_k.tsv and related QC files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        help="Output directory. Default: QC_DIR/paper_plots.",
    )
    parser.add_argument(
        "--prefix",
        default="kmer_qc",
        help="Filename prefix for generated figures and tables.",
    )
    return parser.parse_args()


def configure_matplotlib(outdir: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(outdir / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
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
            "axes.linewidth": 0.8,
            "figure.dpi": 160,
            "savefig.dpi": 450,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required QC file not found: {path}")
    return pd.read_csv(path, sep="\t", **kwargs)


def load_inputs(qc_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = read_tsv(qc_dir / "kmer_shared_summary_by_k.tsv")
    pass_counts = read_tsv(
        qc_dir / "kmer_transcript_pass_counts_by_k.tsv",
        usecols=[
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
        ],
    )
    status_counts = count_statuses_by_genome(qc_dir / "status_matrices")
    return summary, pass_counts, status_counts


def count_statuses_by_genome(status_dir: Path) -> pd.DataFrame:
    if not status_dir.is_dir():
        raise FileNotFoundError(f"Status matrix directory not found: {status_dir}")

    rows: list[dict[str, object]] = []
    for path in sorted(status_dir.glob("kmer_transcript_status_matrix_k*.tsv")):
        k_size = int(path.stem.rsplit("k", 1)[1])
        matrix = pd.read_csv(path, sep="\t")
        genome_cols = list(matrix.columns[4:])
        for genome in genome_cols:
            counts = matrix[genome].value_counts()
            rows.append(
                {
                    "k_size": k_size,
                    "genome": genome,
                    "PASS_STRICT_KMER": int(counts.get("PASS_STRICT_KMER", 0)),
                    "PASS_BASIC_KMER": int(counts.get("PASS_BASIC_KMER", 0)),
                    "LOW_KMER": int(counts.get("LOW_KMER", 0)),
                }
            )

    if not rows:
        raise ValueError(f"No status matrices found in: {status_dir}")
    return pd.DataFrame(rows).sort_values(["k_size", "genome"]).reset_index(drop=True)


def build_plot_summary(summary: pd.DataFrame, pass_counts: pd.DataFrame) -> pd.DataFrame:
    work = summary.copy()
    work["k_size"] = work["k_size"].astype(int)
    work = work.sort_values("k_size").reset_index(drop=True)

    grouped = pass_counts.groupby("k_size", sort=True)
    derived = grouped.apply(summarize_pass_rows, include_groups=False).reset_index()
    plot_summary = work.merge(derived, on="k_size", how="left")

    plot_summary["strict_retention_fraction"] = (
        plot_summary["strict_shared_transcripts_all_genomes"]
        / plot_summary["total_transcripts"]
    )
    plot_summary["basic_retention_fraction"] = (
        plot_summary["basic_shared_transcripts_all_genomes"]
        / plot_summary["total_transcripts"]
    )
    plot_summary["strict_not_retained"] = (
        plot_summary["total_transcripts"]
        - plot_summary["strict_shared_transcripts_all_genomes"]
    )
    plot_summary["basic_not_retained"] = (
        plot_summary["total_transcripts"]
        - plot_summary["basic_shared_transcripts_all_genomes"]
    )
    return plot_summary


def summarize_pass_rows(group: pd.DataFrame) -> pd.Series:
    genome_count = int(group["genome_count"].max())
    zero_kmer = group["total_kmers"].eq(0)
    perfect_all = group["min_matched_kmer_ratio"].eq(1.0)
    strict_all = group["strict_pass_count"].eq(genome_count)
    basic_all = group["basic_pass_count"].eq(genome_count)
    subperfect_retained = group["min_matched_kmer_ratio"].lt(1.0) & strict_all & ~zero_kmer

    return pd.Series(
        {
            "zero_kmer_transcripts": int(zero_kmer.sum()),
            "perfect_all_genomes": int(perfect_all.sum()),
            "subperfect_but_strict_shared": int(subperfect_retained.sum()),
            "strict_failed_transcripts": int((~strict_all).sum()),
            "basic_failed_transcripts": int((~basic_all).sum()),
            "lowest_nonzero_min_ratio": float(
                group.loc[~zero_kmer, "min_matched_kmer_ratio"].min()
            ),
        }
    )


def edge_case_rows(pass_counts: pd.DataFrame) -> pd.DataFrame:
    genome_count = pass_counts["genome_count"]
    mask = (
        pass_counts["total_kmers"].eq(0)
        | pass_counts["min_matched_kmer_ratio"].lt(1.0)
        | pass_counts["strict_pass_count"].lt(genome_count)
    )
    edge = pass_counts.loc[mask].copy()
    edge["edge_case"] = np.select(
        [
            edge["total_kmers"].eq(0),
            edge["strict_pass_count"].eq(edge["genome_count"])
            & edge["min_matched_kmer_ratio"].lt(1.0),
        ],
        ["no_valid_kmers", "subperfect_but_shared"],
        default="not_shared_all_genomes",
    )
    return edge.sort_values(
        ["k_size", "edge_case", "transcript_length", "transcript_name"]
    ).reset_index(drop=True)


def short_genome_label(name: str) -> str:
    clean = name.replace("_genomic", "")
    replacements = {
        "GCA_018503265.2_NA19240_pat_hprc_f2": "NA19240 pat",
        "GCA_018503275.2_NA19240_mat_hprc_f2": "NA19240 mat",
        "GCA_054883195.1_H9_T2T.hap1": "H9 hap1",
        "GCA_054883265.1_H9_T2T.hap2": "H9 hap2",
        "GCF_000258655.2_panpan1.1": "Bonobo panpan1.1",
        "GCF_002880755.1_Clint_PTRv2": "Chimp Clint",
        "GCF_003339765.1_Mmul_10": "Rhesus Mmul_10",
        "GCF_013052645.1_Mhudiblu_PPA_v0": "Orangutan PPA",
        "GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri": "Chimp mPanTro3",
        "GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri": "Bonobo mPanPan1",
        "GCF_049350105.2_T2T-MMU8v2.0": "Rhesus T2T-MMU8",
    }
    return replacements.get(clean, clean)


def human_number(value: float, _position: int | None = None) -> str:
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.0f}K"
    return f"{value:.0f}"


def percent_label(value: float, _position: int | None = None) -> str:
    return f"{value * 100:.4f}"


def style_axis(ax: plt.Axes, grid_axis: str = "y") -> None:
    ax.grid(axis=grid_axis, color=LIGHT_GRID, linewidth=0.7)
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


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_retention_panel(ax: plt.Axes, plot_summary: pd.DataFrame) -> None:
    x = plot_summary["k_size"].to_numpy()
    y = plot_summary["strict_retention_fraction"].to_numpy()
    ax.plot(x, y, color=BLUE, marker="o", markersize=4.4, linewidth=2.0)
    ax.scatter(x, y, color="white", edgecolor=BLUE, linewidth=1.2, s=28, zorder=3)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("All-genome retained transcripts (%)")
    ax.yaxis.set_major_formatter(FuncFormatter(percent_label))

    ymin = max(0.99996, y.min() - 0.000005)
    ax.set_ylim(ymin, 1.000002)
    ax.set_xticks(x)
    ax.set_title("Transcript retention by k-mer size")
    style_axis(ax)
    add_panel_label(ax, "a")

    retained = int(plot_summary["strict_shared_transcripts_all_genomes"].iloc[-1])
    total = int(plot_summary["total_transcripts"].iloc[-1])
    ax.text(
        0.02,
        0.06,
        f"k={x[-1]}: {retained:,} / {total:,} retained",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8,
        color=INK,
    )


def plot_edge_case_bars(ax: plt.Axes, plot_summary: pd.DataFrame) -> None:
    x = np.arange(len(plot_summary))
    zero = plot_summary["zero_kmer_transcripts"].to_numpy()
    subperfect = plot_summary["subperfect_but_strict_shared"].to_numpy()
    width = 0.72

    ax.bar(x, zero, width=width, color=RED, label="No valid k-mers")
    ax.bar(x, subperfect, width=width, bottom=zero, color=GOLD, label="Sub-perfect, retained")
    ax.set_xticks(x, labels=plot_summary["k_size"].astype(str))
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Transcripts")
    ax.set_title("Transcript edge-case classes")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(loc="upper left")
    style_axis(ax)
    add_panel_label(ax, "b")


def plot_status_heatmap(
    ax: plt.Axes,
    fig: plt.Figure,
    status_counts: pd.DataFrame,
) -> None:
    heat = status_counts.pivot(index="genome", columns="k_size", values="LOW_KMER")
    heat = heat.sort_index()
    labels = [short_genome_label(name) for name in heat.index]
    cmap = LinearSegmentedColormap.from_list("low_kmer", ["#F7FBFF", "#9ECAE1", BLUE])

    image = ax.imshow(heat.to_numpy(), aspect="auto", cmap=cmap, vmin=0, vmax=max(3, heat.to_numpy().max()))
    ax.set_xticks(np.arange(heat.shape[1]), labels=[str(k) for k in heat.columns])
    ax.set_yticks(np.arange(heat.shape[0]), labels=labels)
    ax.set_xlabel("k-mer size")
    ax.set_title("LOW_KMER counts by genome")
    ax.tick_params(axis="y", length=0)

    for y in range(heat.shape[0]):
        for x in range(heat.shape[1]):
            value = int(heat.iat[y, x])
            if value > 0:
                ax.text(x, y, str(value), ha="center", va="center", fontsize=7, color=INK)

    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("LOW_KMER transcripts")
    cbar.ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    add_panel_label(ax, "c")


def plot_edge_case_map(ax: plt.Axes, edge: pd.DataFrame, plot_summary: pd.DataFrame) -> None:
    keep = edge[
        edge["transcript_name"].isin(
            ["KIFC1-201", "CABIN1-220", "C4B-236", "KCNQ1DN-202", "ENST00000625083"]
        )
    ].copy()
    if keep.empty:
        ax.set_axis_off()
        return

    order = (
        keep[["transcript_name", "gene_name", "transcript_length"]]
        .drop_duplicates()
        .sort_values(["transcript_length", "transcript_name"])
    )
    y_lookup = {name: index for index, name in enumerate(order["transcript_name"])}
    y_labels = [
        f"{row.transcript_name}\n{int(row.transcript_length)} nt"
        for row in order.itertuples(index=False)
    ]

    style = {
        "no_valid_kmers": {
            "label": "No valid k-mers",
            "color": RED,
            "marker": "s",
            "size": 82,
        },
        "subperfect_but_shared": {
            "label": "Sub-perfect, retained",
            "color": GOLD,
            "marker": "o",
            "size": 72,
        },
    }
    for edge_case, spec in style.items():
        subset = keep[keep["edge_case"].eq(edge_case)]
        if subset.empty:
            continue
        ax.scatter(
            subset["k_size"],
            subset["transcript_name"].map(y_lookup),
            s=spec["size"],
            marker=spec["marker"],
            color=spec["color"],
            edgecolor="white",
            linewidth=0.7,
            label=spec["label"],
            zorder=3,
        )

    ax.set_xlabel("k-mer size")
    ax.set_yticks(range(len(y_labels)), labels=y_labels)
    ax.set_xticks(plot_summary["k_size"].to_numpy())
    ax.set_xlim(plot_summary["k_size"].min() - 0.7, plot_summary["k_size"].max() + 0.7)
    ax.set_ylim(-0.6, len(y_labels) - 0.4)
    ax.invert_yaxis()
    ax.set_title("Edge cases by transcript and k-mer size")
    ax.legend(loc="lower left")
    style_axis(ax, grid_axis="x")
    add_panel_label(ax, "d")


def make_main_figure(
    plot_summary: pd.DataFrame,
    edge: pd.DataFrame,
    status_counts: pd.DataFrame,
    outdir: Path,
    prefix: str,
) -> None:
    fig = plt.figure(figsize=(12.3, 8.2))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 1.15],
        height_ratios=[1.0, 1.1],
        wspace=0.33,
        hspace=0.38,
    )

    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    plot_retention_panel(ax_a, plot_summary)
    plot_edge_case_bars(ax_b, plot_summary)
    plot_status_heatmap(ax_c, fig, status_counts)
    plot_edge_case_map(ax_d, edge, plot_summary)

    fig.suptitle(
        "Matched k-mer transcript support across primate genome assemblies",
        fontsize=13,
        y=0.99,
    )
    save_figure(fig, outdir / f"{prefix}_paper_figure")


def make_supplement_figures(
    plot_summary: pd.DataFrame,
    status_counts: pd.DataFrame,
    outdir: Path,
    prefix: str,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    x = plot_summary["k_size"].to_numpy()
    ax.plot(
        x,
        plot_summary["strict_not_retained"],
        marker="o",
        color=RED,
        linewidth=2,
        label="Strict/basic not retained",
    )
    ax.plot(
        x,
        plot_summary["subperfect_but_strict_shared"],
        marker="s",
        color=GOLD,
        linewidth=2,
        label="Sub-perfect but retained",
    )
    ax.set_xticks(x)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Transcripts")
    ax.set_title("Rare non-perfect k-mer outcomes by k")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(loc="upper left")
    style_axis(ax)
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_rare_outcomes_by_k")

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    heat = status_counts.pivot(index="genome", columns="k_size", values="LOW_KMER")
    heat = heat.sort_index()
    plot_status_heatmap(ax, fig, status_counts)
    fig.tight_layout()
    save_figure(fig, outdir / f"{prefix}_low_kmer_heatmap")


def write_tables(
    plot_summary: pd.DataFrame,
    edge: pd.DataFrame,
    status_counts: pd.DataFrame,
    outdir: Path,
    prefix: str,
) -> None:
    plot_summary.to_csv(outdir / f"{prefix}_plot_summary.tsv", sep="\t", index=False)
    edge.to_csv(outdir / f"{prefix}_edge_cases.tsv", sep="\t", index=False)
    status_counts.to_csv(outdir / f"{prefix}_status_counts_by_genome.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    outdir = args.outdir or args.qc_dir / "paper_plots"
    outdir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib(outdir)

    summary, pass_counts, status_counts = load_inputs(args.qc_dir)
    plot_summary = build_plot_summary(summary, pass_counts)
    edge = edge_case_rows(pass_counts)

    write_tables(plot_summary, edge, status_counts, outdir, args.prefix)
    make_main_figure(plot_summary, edge, status_counts, outdir, args.prefix)
    make_supplement_figures(plot_summary, status_counts, outdir, args.prefix)

    print(f"Wrote k-mer QC plots and tables to: {outdir}")


if __name__ == "__main__":
    main()
