#!/usr/bin/env python3
"""
Compare bonobo alignment, k-mer, and gene-expression evidence.

Expression was quantified with featureCounts using gene_id, so expression
evidence is gene-level. Transcript plots therefore mean transcripts whose
parent gene is expressed, not direct transcript-isoform expression.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, MaxNLocator


ALIGNMENT_COLUMN = "human_lncRNA_vs_mPanPan1v20pri"
BONOBO_KMER_COLUMN = "GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic"

BLUE = "#2F6C9E"
TEAL = "#1B9E77"
GOLD = "#D89C28"
RED = "#C44E52"
PURPLE = "#756BB1"
INK = "#252525"
LIGHT_GRID = "#E8E8E8"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot bonobo alignment and k-mer support together with gene-level "
            "bonobo expression evidence."
        )
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        default=Path("alignment_results/lncRNA_transcript_pass_counts.tsv"),
        help="Alignment transcript status table.",
    )
    parser.add_argument(
        "--kmer",
        type=Path,
        default=Path(
            "kmer_results/ibf_new/qc/kmer_transcript_pass_counts_by_k.tsv"
        ),
        help="K-mer transcript summary table.",
    )
    parser.add_argument(
        "--status-dir",
        type=Path,
        help="Directory with per-k status matrices. Default: KMER parent/status_matrices.",
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("expression_data/bonobo/bonobo_lnc_fpkm_PASS_strict.txt"),
        help="Gene-by-sample FPKM matrix.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("bonobo_expression_support"),
        help="Output directory.",
    )
    parser.add_argument(
        "--k-size",
        type=int,
        help="k-mer size for detailed panels. Default: largest available k.",
    )
    parser.add_argument(
        "--expression-threshold",
        type=float,
        default=0.0,
        help="A sample is expressed when FPKM is greater than this value.",
    )
    parser.add_argument(
        "--min-expressed-samples",
        type=int,
        default=1,
        help="Minimum number of expressed samples required to call a gene expressed.",
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


def require_columns(df: pd.DataFrame, columns: list[str], source: Path) -> None:
    missing = sorted(set(columns).difference(df.columns))
    if missing:
        raise SystemExit(f"Missing column(s) in {source}: {', '.join(missing)}")


def read_alignment(path: Path) -> pd.DataFrame:
    columns = [
        "transcript_id",
        "transcript_name",
        "gene_id",
        "gene_name",
        ALIGNMENT_COLUMN,
    ]
    df = pd.read_csv(path, sep="\t", usecols=columns)
    require_columns(df, columns, path)
    if df["transcript_id"].duplicated().any():
        raise SystemExit(f"Duplicate transcript_id values in {path}")
    return df.rename(columns={ALIGNMENT_COLUMN: "alignment_status"})


def read_kmer_summary(path: Path) -> pd.DataFrame:
    columns = [
        "k_size",
        "transcript_id",
        "gene_id",
        "strict_pass_count",
        "genome_count",
    ]
    df = pd.read_csv(path, sep="\t", usecols=columns)
    require_columns(df, columns, path)
    df["k_size"] = df["k_size"].astype(int)
    df["all_genomes_kmer_strict"] = df["strict_pass_count"].eq(df["genome_count"])
    return df


def read_expression(
    path: Path,
    threshold: float,
    min_samples: int,
) -> pd.DataFrame:
    if min_samples < 1:
        raise SystemExit("--min-expressed-samples must be at least 1")

    matrix = pd.read_csv(path, sep="\t", index_col=0)
    if matrix.index.duplicated().any():
        raise SystemExit(f"Duplicate gene IDs in expression matrix: {path}")
    matrix = matrix.apply(pd.to_numeric, errors="raise")

    positive = matrix.gt(threshold)
    expressed_samples = positive.sum(axis=1)
    result = pd.DataFrame(
        {
            "gene_id": matrix.index.astype(str),
            "expression_available": True,
            "expressed_sample_count": expressed_samples.to_numpy(),
            "expressed_sample_fraction": (expressed_samples / matrix.shape[1]).to_numpy(),
            "mean_fpkm": matrix.mean(axis=1).to_numpy(),
            "median_fpkm": matrix.median(axis=1).to_numpy(),
            "max_fpkm": matrix.max(axis=1).to_numpy(),
        }
    )
    result["gene_expressed"] = result["expressed_sample_count"].ge(min_samples)
    return result


def status_matrix_paths(status_dir: Path) -> dict[int, Path]:
    paths: dict[int, Path] = {}
    for path in status_dir.glob("kmer_transcript_status_matrix_k*.tsv"):
        try:
            k_size = int(path.stem.rsplit("k", 1)[1])
        except ValueError:
            continue
        paths[k_size] = path
    if not paths:
        raise SystemExit(f"No k-mer status matrices found in {status_dir}")
    return dict(sorted(paths.items()))


def read_bonobo_kmer_status(path: Path) -> pd.DataFrame:
    columns = ["transcript_id", BONOBO_KMER_COLUMN]
    df = pd.read_csv(path, sep="\t", usecols=columns)
    require_columns(df, columns, path)
    return df.rename(columns={BONOBO_KMER_COLUMN: "bonobo_kmer_status"})


def build_transcript_table(
    alignment: pd.DataFrame,
    kmer_for_k: pd.DataFrame,
    bonobo_status: pd.DataFrame,
    expression: pd.DataFrame,
) -> pd.DataFrame:
    kmer_columns = [
        "transcript_id",
        "all_genomes_kmer_strict",
        "strict_pass_count",
        "genome_count",
    ]
    merged = alignment.merge(
        kmer_for_k[kmer_columns],
        on="transcript_id",
        how="left",
        validate="one_to_one",
    )
    merged = merged.merge(
        bonobo_status,
        on="transcript_id",
        how="left",
        validate="one_to_one",
    )
    merged = merged.merge(
        expression,
        on="gene_id",
        how="left",
        validate="many_to_one",
    )

    merged["alignment_strict"] = merged["alignment_status"].eq("PASS_STRICT")
    merged["bonobo_kmer_strict"] = merged["bonobo_kmer_status"].eq(
        "PASS_STRICT_KMER"
    )
    merged["all_genomes_kmer_strict"] = (
        merged["all_genomes_kmer_strict"].fillna(False).astype(bool)
    )
    merged["expression_available"] = (
        merged["expression_available"].fillna(False).astype(bool)
    )
    merged["gene_expressed"] = merged["gene_expressed"].fillna(False).astype(bool)
    for column in [
        "expressed_sample_count",
        "expressed_sample_fraction",
        "mean_fpkm",
        "median_fpkm",
        "max_fpkm",
    ]:
        merged[column] = merged[column].fillna(0)

    merged["aligned_and_gene_expressed"] = (
        merged["alignment_strict"] & merged["gene_expressed"]
    )
    merged["bonobo_three_way_support"] = (
        merged["aligned_and_gene_expressed"] & merged["bonobo_kmer_strict"]
    )
    merged["all_genomes_three_way_support"] = (
        merged["aligned_and_gene_expressed"]
        & merged["all_genomes_kmer_strict"]
    )
    return merged


def count_transcripts_and_genes(
    transcript_table: pd.DataFrame,
    mask: pd.Series,
) -> tuple[int, int]:
    selected = transcript_table.loc[mask]
    return len(selected), selected["gene_id"].nunique()


def build_k_summary(
    alignment: pd.DataFrame,
    kmer_summary: pd.DataFrame,
    expression: pd.DataFrame,
    matrix_paths: dict[int, Path],
) -> pd.DataFrame:
    base = alignment[["transcript_id", "gene_id", "alignment_status"]].copy()
    base["alignment_strict"] = base["alignment_status"].eq("PASS_STRICT")
    base = base.merge(
        expression[["gene_id", "gene_expressed"]],
        on="gene_id",
        how="left",
        validate="many_to_one",
    )
    base["gene_expressed"] = base["gene_expressed"].fillna(False).astype(bool)

    rows: list[dict[str, int]] = []
    for k_size, matrix_path in matrix_paths.items():
        kmer_for_k = kmer_summary[kmer_summary["k_size"].eq(k_size)][
            ["transcript_id", "all_genomes_kmer_strict"]
        ]
        bonobo = read_bonobo_kmer_status(matrix_path)
        work = base.merge(
            kmer_for_k,
            on="transcript_id",
            how="left",
            validate="one_to_one",
        ).merge(
            bonobo,
            on="transcript_id",
            how="left",
            validate="one_to_one",
        )
        work["all_genomes_kmer_strict"] = (
            work["all_genomes_kmer_strict"].fillna(False).astype(bool)
        )
        work["bonobo_kmer_strict"] = work["bonobo_kmer_status"].eq(
            "PASS_STRICT_KMER"
        )
        eligible = work["alignment_strict"] & work["gene_expressed"]
        bonobo_supported = eligible & work["bonobo_kmer_strict"]
        all_supported = eligible & work["all_genomes_kmer_strict"]
        rows.append(
            {
                "k_size": k_size,
                "eligible_transcripts": int(eligible.sum()),
                "eligible_genes": int(work.loc[eligible, "gene_id"].nunique()),
                "bonobo_supported_transcripts": int(bonobo_supported.sum()),
                "bonobo_supported_genes": int(
                    work.loc[bonobo_supported, "gene_id"].nunique()
                ),
                "all_genomes_supported_transcripts": int(all_supported.sum()),
                "all_genomes_supported_genes": int(
                    work.loc[all_supported, "gene_id"].nunique()
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("k_size").reset_index(drop=True)


def build_gene_table(transcript_table: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        transcript_table.groupby("gene_id", sort=False)
        .agg(
            gene_name=("gene_name", "first"),
            transcript_count=("transcript_id", "nunique"),
            alignment_strict_transcripts=("alignment_strict", "sum"),
            bonobo_kmer_strict_transcripts=("bonobo_kmer_strict", "sum"),
            all_genomes_kmer_strict_transcripts=("all_genomes_kmer_strict", "sum"),
            bonobo_three_way_transcripts=("bonobo_three_way_support", "sum"),
            all_genomes_three_way_transcripts=("all_genomes_three_way_support", "sum"),
            expression_available=("expression_available", "max"),
            gene_expressed=("gene_expressed", "max"),
            expressed_sample_count=("expressed_sample_count", "max"),
            expressed_sample_fraction=("expressed_sample_fraction", "max"),
            mean_fpkm=("mean_fpkm", "max"),
            median_fpkm=("median_fpkm", "max"),
            max_fpkm=("max_fpkm", "max"),
        )
        .reset_index()
    )
    grouped["alignment_strict"] = grouped["alignment_strict_transcripts"].gt(0)
    grouped["bonobo_kmer_strict"] = grouped["bonobo_three_way_transcripts"].gt(0)
    grouped["all_genomes_kmer_strict"] = grouped[
        "all_genomes_three_way_transcripts"
    ].gt(0)
    grouped["kmer_conservation_class"] = np.select(
        [
            grouped["all_genomes_kmer_strict"],
            grouped["bonobo_kmer_strict"],
        ],
        [
            "Strict across all genomes",
            "Strict in bonobo only",
        ],
        default="Not strict in bonobo",
    )
    return grouped


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
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.0f}K"
    return f"{value:.0f}"


def plot_evidence_funnel(
    ax: plt.Axes,
    transcript_table: pd.DataFrame,
    k_size: int,
) -> None:
    masks = [
        pd.Series(True, index=transcript_table.index),
        transcript_table["alignment_strict"],
        transcript_table["aligned_and_gene_expressed"],
        transcript_table["bonobo_three_way_support"],
        transcript_table["all_genomes_three_way_support"],
    ]
    labels = [
        "All annotated transcripts",
        "Strict alignment in bonobo",
        "Strict alignment + expressed gene",
        f"+ strict bonobo k-mer (k={k_size})",
        f"+ strict k-mer in all genomes (k={k_size})",
    ]
    counts = [count_transcripts_and_genes(transcript_table, mask) for mask in masks]
    transcript_counts = [count[0] for count in counts]
    colors = ["#A7A7A7", BLUE, TEAL, GOLD, PURPLE]
    y = np.arange(len(labels))

    ax.barh(y, transcript_counts, color=colors, height=0.68)
    ax.set_yticks(y, labels=labels)
    ax.invert_yaxis()
    ax.set_xlabel("Transcripts")
    ax.set_title("Nested evidence support")
    ax.xaxis.set_major_formatter(FuncFormatter(human_number))
    ax.margins(x=0.2)
    for index, (transcripts, genes) in enumerate(counts):
        ax.text(
            transcripts + max(transcript_counts) * 0.015,
            index,
            f"{transcripts:,} transcripts\n{genes:,} genes",
            va="center",
            fontsize=8,
            color=INK,
        )
    style_axis(ax, axis="x")
    add_panel_label(ax, "a")


def plot_k_size_support(ax: plt.Axes, summary: pd.DataFrame) -> None:
    x = summary["k_size"].to_numpy()
    ax.plot(
        x,
        summary["bonobo_supported_transcripts"],
        color=TEAL,
        marker="o",
        linewidth=2.2,
        label="Strict k-mer in bonobo",
    )
    ax.plot(
        x,
        summary["all_genomes_supported_transcripts"],
        color=PURPLE,
        marker="s",
        linewidth=2.2,
        label="Strict k-mer in all genomes",
    )
    ax.set_xticks(x)
    ax.set_xlabel("k-mer size")
    ax.set_ylabel("Supported transcripts")
    ax.set_title("Expression-supported transcripts across k sizes")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.margins(y=0.1)
    ax.legend(loc="best")
    style_axis(ax)
    add_panel_label(ax, "b")


def expression_groups(gene_table: pd.DataFrame) -> tuple[list[str], list[pd.DataFrame]]:
    order = [
        "Not strict in bonobo",
        "Strict in bonobo only",
        "Strict across all genomes",
    ]
    expressed = gene_table[
        gene_table["alignment_strict"] & gene_table["gene_expressed"]
    ]
    groups = [
        expressed[expressed["kmer_conservation_class"].eq(label)] for label in order
    ]
    return order, groups


def styled_boxplot(
    ax: plt.Axes,
    groups: list[pd.DataFrame],
    column: str,
    colors: list[str],
) -> None:
    values = [
        group[column].to_numpy() if not group.empty else np.array([np.nan])
        for group in groups
    ]
    artists = ax.boxplot(
        values,
        patch_artist=True,
        showfliers=False,
        widths=0.6,
        medianprops={"color": INK, "linewidth": 1.5},
        whiskerprops={"color": "#555555"},
        capprops={"color": "#555555"},
    )
    for patch, color in zip(artists["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)


def plot_expression_breadth(ax: plt.Axes, gene_table: pd.DataFrame) -> None:
    labels, groups = expression_groups(gene_table)
    colors = [RED, GOLD, PURPLE]
    styled_boxplot(ax, groups, "expressed_sample_fraction", colors)
    tick_labels = [
        f"{label}\n(n={len(group):,} genes)" for label, group in zip(labels, groups)
    ]
    ax.set_xticks(range(1, len(labels) + 1), labels=tick_labels)
    ax.set_ylabel("Fraction of 98 samples with expression")
    ax.set_title("Expression breadth by k-mer conservation class")
    ax.set_ylim(-0.03, 1.03)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:.0%}"))
    style_axis(ax)
    add_panel_label(ax, "c")


def plot_expression_abundance(ax: plt.Axes, gene_table: pd.DataFrame) -> None:
    labels, groups = expression_groups(gene_table)
    colors = [RED, GOLD, PURPLE]
    plot_groups = [group.assign(log_mean=np.log10(group["mean_fpkm"] + 0.01)) for group in groups]
    styled_boxplot(ax, plot_groups, "log_mean", colors)
    tick_labels = [
        f"{label}\n(n={len(group):,} genes)" for label, group in zip(labels, groups)
    ]
    ax.set_xticks(range(1, len(labels) + 1), labels=tick_labels)
    ax.set_ylabel("Mean FPKM (log10 scale)")
    ax.set_title("Expression abundance by k-mer conservation class")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=6))
    style_axis(ax)
    add_panel_label(ax, "d")


def make_figure(
    transcript_table: pd.DataFrame,
    gene_table: pd.DataFrame,
    k_summary: pd.DataFrame,
    k_size: int,
) -> plt.Figure:
    fig, axes = plt.subplots(2, 2, figsize=(14.8, 9.2))
    plot_evidence_funnel(axes[0, 0], transcript_table, k_size)
    plot_k_size_support(axes[0, 1], k_summary)
    plot_expression_breadth(axes[1, 0], gene_table)
    plot_expression_abundance(axes[1, 1], gene_table)
    fig.suptitle(
        "Bonobo lncRNA support from alignment, k-mers, and gene expression",
        fontsize=14,
        y=0.995,
    )
    fig.subplots_adjust(wspace=0.34, hspace=0.38, top=0.92, bottom=0.1)
    return fig


def write_readme(
    path: Path,
    args: argparse.Namespace,
    k_size: int,
    sample_count: int,
) -> None:
    prefilter_note = ""
    if "PASS_strict" in args.expression.name:
        prefilter_note = (
            "The expression matrix name indicates that it was prefiltered for "
            "strict alignment support, so expression is not independent of the "
            "alignment set.\n"
        )
    text = f"""Bonobo alignment/k-mer/expression comparison

Expression input: {args.expression}
Expression is gene-level because featureCounts used gene_id.
A gene is called expressed when FPKM > {args.expression_threshold:g} in at least
{args.min_expressed_samples} of {sample_count} samples.
{prefilter_note}

Alignment support is transcript-level PASS_STRICT for {ALIGNMENT_COLUMN}.
Bonobo k-mer support is transcript-level PASS_STRICT_KMER for
{BONOBO_KMER_COLUMN}.
All-genome k-mer support requires strict_pass_count == genome_count.

Detailed expression-distribution panels use k={k_size}.
The phrase "transcripts whose gene is expressed" must not be interpreted as
transcript-isoform-specific expression.
"""
    path.write_text(text, encoding="ascii")


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    args.outdir.mkdir(parents=True, exist_ok=True)

    status_dir = args.status_dir or args.kmer.parent / "status_matrices"
    matrix_paths = status_matrix_paths(status_dir)
    available_k = sorted(matrix_paths)
    k_size = args.k_size if args.k_size is not None else max(available_k)
    if k_size not in matrix_paths:
        choices = ", ".join(map(str, available_k))
        raise SystemExit(f"k={k_size} not found in {status_dir}; available: {choices}")

    alignment = read_alignment(args.alignment)
    kmer_summary = read_kmer_summary(args.kmer)
    expression = read_expression(
        args.expression,
        args.expression_threshold,
        args.min_expressed_samples,
    )
    kmer_for_k = kmer_summary[kmer_summary["k_size"].eq(k_size)].copy()
    bonobo_status = read_bonobo_kmer_status(matrix_paths[k_size])

    transcript_table = build_transcript_table(
        alignment,
        kmer_for_k,
        bonobo_status,
        expression,
    )
    gene_table = build_gene_table(transcript_table)
    k_summary = build_k_summary(
        alignment,
        kmer_summary,
        expression,
        matrix_paths,
    )

    prefix = args.outdir / "bonobo_alignment_kmer_expression"
    transcript_table.to_csv(
        args.outdir / f"bonobo_transcript_evidence_k{k_size}.tsv",
        sep="\t",
        index=False,
    )
    transcript_table[transcript_table["bonobo_three_way_support"]].to_csv(
        args.outdir / f"bonobo_three_way_supported_transcripts_k{k_size}.tsv",
        sep="\t",
        index=False,
    )
    gene_table.to_csv(
        args.outdir / f"bonobo_gene_evidence_k{k_size}.tsv",
        sep="\t",
        index=False,
    )
    k_summary.to_csv(
        args.outdir / "bonobo_expression_supported_counts_by_k.tsv",
        sep="\t",
        index=False,
    )

    fig = make_figure(transcript_table, gene_table, k_summary, k_size)
    fig.savefig(prefix.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    sample_count = pd.read_csv(args.expression, sep="\t", nrows=0).shape[1] - 1
    write_readme(args.outdir / "README.txt", args, k_size, sample_count)

    print(f"Wrote {prefix.with_suffix('.png')}")
    print(f"Wrote {prefix.with_suffix('.pdf')}")
    print(f"Wrote merged evidence tables to {args.outdir}")


if __name__ == "__main__":
    main()
