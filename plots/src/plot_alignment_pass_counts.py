
import argparse
import os
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STATUS_COLORS = {
    "PASS_STRICT": "#2b8cbe",
    "PASS_BASIC": "#7bccc4",
    "LOW_QC": "#fdae61",
    "UNMAPPED": "#d7301f",
    "MISSING": "#bdbdbd",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate alignment-result plots from "
            "lncRNA_transcript_pass_counts.tsv."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        default=(
            "/mnt/d/primates_bmc/bam_comparison/qc/shared_lncRNA_transcripts/"
            "lncRNA_transcript_pass_counts.tsv"
        ),
        help="Input lncRNA_transcript_pass_counts.tsv file.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory for figures and summary tables.",
    )
    parser.add_argument(
        "--top-intersections",
        type=int,
        default=25,
        help=(
            "Number of largest gene_name intersections to show in each UpSet "
            "plot. Use 0 to show all nonzero intersections."
        ),
    )
    parser.add_argument(
        "--strip-gene-version",
        action="store_true",
        help="Remove version suffixes from gene_name values, e.g. ENSG... .1.",
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
    cols = list(df.columns[start:end])
    if not cols:
        raise SystemExit("No genome/sample status columns found.")
    return cols


def count_unique_genes(series: pd.Series) -> int:
    return int(series.dropna().astype(str).replace("", np.nan).dropna().nunique())


def join_unique(values: pd.Series) -> str:
    unique = sorted(
        {
            str(value)
            for value in values.dropna()
            if str(value) != ""
        }
    )
    return ",".join(unique)


def build_gene_summary(df: pd.DataFrame, status_cols: Iterable[str]) -> pd.DataFrame:
    genome_count = int(df["genome_count"].max())
    strict_shared = df["strict_pass_count"].eq(genome_count)
    basic_shared = df["basic_pass_count"].eq(genome_count)
    unmapped_all = df[list(status_cols)].eq("UNMAPPED").all(axis=1)

    work = df[
        [
            "gene_id",
            "gene_name",
            "transcript_id",
            "strict_pass_count",
            "basic_pass_count",
        ]
    ].copy()
    work["strict_shared_transcript"] = strict_shared
    work["basic_shared_transcript"] = basic_shared
    work["unmapped_all_genomes_transcript"] = unmapped_all

    summary = (
        work.groupby("gene_name", dropna=False)
        .agg(
            associated_gene_ids=("gene_id", join_unique),
            n_transcripts=("transcript_id", "nunique"),
            strict_shared_transcripts=("strict_shared_transcript", "sum"),
            basic_shared_transcripts=("basic_shared_transcript", "sum"),
            unmapped_all_genomes_transcripts=("unmapped_all_genomes_transcript", "sum"),
            max_strict_pass_count=("strict_pass_count", "max"),
            max_basic_pass_count=("basic_pass_count", "max"),
        )
        .reset_index()
    )
    summary["strict_shared_gene"] = summary["strict_shared_transcripts"].gt(0)
    summary["basic_shared_gene"] = summary["basic_shared_transcripts"].gt(0)
    return summary


def write_summary_tables(
    df: pd.DataFrame,
    gene_summary: pd.DataFrame,
    status_cols: list[str],
    outdir: Path,
) -> None:
    genome_count = int(df["genome_count"].max())
    strict_shared = df["strict_pass_count"].eq(genome_count)
    basic_shared = df["basic_pass_count"].eq(genome_count)
    unmapped_all = df[status_cols].eq("UNMAPPED").all(axis=1)

    overview_rows = [
        ("total_transcripts", int(df["transcript_id"].nunique())),
        ("total_unique_gene_names", count_unique_genes(df["gene_name"])),
        ("strict_shared_transcripts_all_genomes", int(strict_shared.sum())),
        ("basic_shared_transcripts_all_genomes", int(basic_shared.sum())),
        (
            "strict_shared_unique_gene_names",
            int(gene_summary["strict_shared_gene"].sum()),
        ),
        (
            "basic_shared_unique_gene_names",
            int(gene_summary["basic_shared_gene"].sum()),
        ),
        ("transcripts_unmapped_in_all_genomes", int(unmapped_all.sum())),
        (
            "unique_gene_names_with_transcripts_unmapped_in_all_genomes",
            count_unique_genes(df.loc[unmapped_all, "gene_name"]),
        ),
        ("genome_count", genome_count),
    ]
    pd.DataFrame(overview_rows, columns=["metric", "value"]).to_csv(
        outdir / "alignment_plot_summary.tsv", sep="\t", index=False
    )

    gene_summary.to_csv(outdir / "unique_gene_summary.tsv", sep="\t", index=False)
    df.loc[unmapped_all].to_csv(
        outdir / "transcripts_unmapped_in_all_genomes.tsv",
        sep="\t",
        index=False,
    )

    status_counts = []
    for col in status_cols:
        counts = df[col].value_counts()
        for status in STATUS_COLORS:
            status_counts.append(
                {
                    "genome": col,
                    "status": status,
                    "transcript_count": int(counts.get(status, 0)),
                }
            )
    pd.DataFrame(status_counts).to_csv(
        outdir / "per_genome_status_counts.tsv", sep="\t", index=False
    )


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.6)
    ax.set_axisbelow(True)


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


def plot_transcript_result_sets(
    ax: plt.Axes,
    df: pd.DataFrame,
    status_cols: list[str],
) -> None:
    genome_count = int(df["genome_count"].max())
    strict_shared = df["strict_pass_count"].eq(genome_count)
    basic_shared = df["basic_pass_count"].eq(genome_count)
    unmapped_all = df[status_cols].eq("UNMAPPED").all(axis=1)

    shared_counts = [
        int(strict_shared.sum()),
        int(basic_shared.sum()),
        int(unmapped_all.sum()),
    ]
    bars = ax.bar(
        ["Strict\nshared", "Basic\nshared", "Unmapped\nall"],
        shared_counts,
        color=["#2b8cbe", "#7bccc4", "#d7301f"],
    )
    ax.set_title("Transcript-level result sets")
    ax.set_ylabel("Transcripts")
    ax.bar_label(bars, padding=3, fontsize=8)
    add_panel_label(ax, "a")
    style_axis(ax)


def plot_unique_gene_counts(
    ax: plt.Axes,
    df: pd.DataFrame,
    gene_summary: pd.DataFrame,
    status_cols: list[str],
) -> None:
    unmapped_all = df[status_cols].eq("UNMAPPED").all(axis=1)

    gene_counts = [
        count_unique_genes(df["gene_name"]),
        int(gene_summary["strict_shared_gene"].sum()),
        int(gene_summary["basic_shared_gene"].sum()),
        count_unique_genes(df.loc[unmapped_all, "gene_name"]),
    ]
    bars = ax.bar(
        ["All", "Strict\nshared", "Basic\nshared", "Unmapped\nall"],
        gene_counts,
        color=["#636363", "#2b8cbe", "#7bccc4", "#d7301f"],
    )
    ax.set_title("Unique gene names associated with transcripts")
    ax.set_ylabel("Unique gene names")
    ax.bar_label(bars, padding=3, fontsize=8)
    add_panel_label(ax, "b")
    style_axis(ax)


def gene_membership(
    df: pd.DataFrame,
    status_cols: list[str],
    mode: str,
) -> pd.DataFrame:
    if mode == "strict":
        transcript_pass = df[status_cols].eq("PASS_STRICT")
    elif mode == "basic":
        transcript_pass = df[status_cols].isin(["PASS_STRICT", "PASS_BASIC"])
    else:
        raise ValueError(f"Unknown mode: {mode}")

    gene_pass = transcript_pass.copy()
    gene_pass["gene_name"] = df["gene_name"].astype(str)
    return gene_pass.groupby("gene_name", dropna=False)[status_cols].any()


def intersection_table(
    membership: pd.DataFrame,
    status_cols: list[str],
) -> pd.DataFrame:
    rows = []
    for pattern, count in membership.groupby(status_cols).size().items():
        if not isinstance(pattern, tuple):
            pattern = (pattern,)
        if not any(pattern):
            continue
        genomes = [col for col, present in zip(status_cols, pattern) if present]
        rows.append(
            {
                "intersection": "&".join(genomes),
                "n_genomes": len(genomes),
                "gene_count": int(count),
                **{col: bool(present) for col, present in zip(status_cols, pattern)},
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=["intersection", "n_genomes", "gene_count", *status_cols]
        )
    table = pd.DataFrame(rows)
    return table.sort_values(
        ["gene_count", "n_genomes", "intersection"],
        ascending=[False, False, True],
    ).reset_index(drop=True)


def select_upset_rows(
    table: pd.DataFrame,
    status_cols: list[str],
    top_n: int,
) -> pd.DataFrame:
    """Keep the all-genomes shared intersection first, then largest others."""
    if table.empty:
        return table.copy()

    plot_df = table.copy()
    plot_df["all_genomes_shared"] = plot_df[status_cols].all(axis=1)

    if top_n <= 0:
        return plot_df.sort_values(
            ["all_genomes_shared", "gene_count", "n_genomes"],
            ascending=[False, False, False],
        ).reset_index(drop=True)

    all_shared = plot_df[plot_df["all_genomes_shared"]]
    remaining = plot_df[~plot_df["all_genomes_shared"]].sort_values(
        ["gene_count", "n_genomes", "intersection"],
        ascending=[False, False, True],
    )
    if not all_shared.empty:
        all_shared = all_shared.sort_values("gene_count", ascending=False).head(1)
        remaining = remaining.head(max(top_n - 1, 0))
        return pd.concat([all_shared, remaining], ignore_index=True)

    return remaining.head(top_n).reset_index(drop=True)


def plot_upset(
    table: pd.DataFrame,
    status_cols: list[str],
    outdir: Path,
    mode: str,
    top_n: int,
) -> None:
    table.to_csv(
        outdir / f"gene_upset_intersections_{mode}.tsv", sep="\t", index=False
    )
    if table.empty:
        return

    plot_df = select_upset_rows(table, status_cols, top_n)
    n_patterns = len(plot_df)
    n_sets = len(status_cols)
    labels = [short_sample_name(c) for c in status_cols]

    fig = plt.figure(figsize=(max(10, n_patterns * 0.45), 6.8))
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=[2.2, 1.45],
        hspace=0.05,
    )
    ax_bar = fig.add_subplot(gs[0])
    ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(n_patterns)
    counts = plot_df["gene_count"].to_numpy()
    color = "#2b8cbe" if mode == "strict" else "#7bccc4"
    colors = np.repeat(color, n_patterns).astype(object)
    colors[plot_df["all_genomes_shared"].to_numpy(dtype=bool)] = "#08519c"
    bars = ax_bar.bar(x, counts, color=colors)
    ax_bar.set_ylabel("Unique gene names")
    all_shared_count = int(
        plot_df.loc[plot_df["all_genomes_shared"], "gene_count"].iloc[0]
    ) if plot_df["all_genomes_shared"].any() else 0
    ax_bar.set_title(
        f"Gene-name intersections across genomes "
        #f"all genomes shared = {all_shared_count:,}"
    )
    ax_bar.bar_label(bars, padding=2, fontsize=7, rotation=90)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.grid(axis="y", color="#e5e5e5", linewidth=0.6)
    ax_bar.set_axisbelow(True)
    plt.setp(ax_bar.get_xticklabels(), visible=False)

    for y, col in enumerate(status_cols):
        ax_matrix.scatter(
            x,
            np.full(n_patterns, y),
            s=22,
            color="#d9d9d9",
            zorder=1,
        )
        present = plot_df[col].to_numpy(dtype=bool)
        ax_matrix.scatter(
            x[present],
            np.full(present.sum(), y),
            s=34,
            color="#252525",
            zorder=2,
        )
        for xi, is_present in zip(x, present):
            if is_present:
                members = np.where(plot_df.iloc[xi][status_cols].to_numpy(dtype=bool))[0]
                if len(members) > 1:
                    ax_matrix.plot(
                        [xi, xi],
                        [members.min(), members.max()],
                        color="#252525",
                        linewidth=1,
                        zorder=1.5,
                    )

    ax_matrix.set_yticks(range(n_sets))
    ax_matrix.set_yticklabels(labels)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xlabel("Genome-set intersections")
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels(
        ["All\ngenomes" if shared else "" for shared in plot_df["all_genomes_shared"]],
        fontsize=8,
    )
    ax_matrix.tick_params(axis="x", length=0)
    ax_matrix.spines["top"].set_visible(False)
    ax_matrix.spines["right"].set_visible(False)
    ax_matrix.spines["left"].set_visible(False)
    ax_matrix.tick_params(axis="y", length=0)

    savefig(fig, outdir / f"gene_upset_{mode}")


def draw_upset_panel(
    fig: plt.Figure,
    outer_spec,
    table: pd.DataFrame,
    status_cols: list[str],
    mode: str,
    top_n: int,
    panel_label: str,
) -> None:
    inner = outer_spec.subgridspec(
        2,
        1,
        height_ratios=[2.0, 1.55],
        hspace=0.04,
    )
    ax_bar = fig.add_subplot(inner[0])
    ax_matrix = fig.add_subplot(inner[1], sharex=ax_bar)

    if table.empty:
        ax_bar.text(0.5, 0.5, "No gene intersections", ha="center", va="center")
        ax_bar.set_axis_off()
        ax_matrix.set_axis_off()
        add_panel_label(ax_bar, panel_label)
        return

    plot_df = select_upset_rows(table, status_cols, top_n)
    n_patterns = len(plot_df)
    n_sets = len(status_cols)
    labels = [short_sample_name(c) for c in status_cols]
    x = np.arange(n_patterns)
    counts = plot_df["gene_count"].to_numpy()
    color = "#2b8cbe" if mode == "strict" else "#7bccc4"
    colors = np.repeat(color, n_patterns).astype(object)
    colors[plot_df["all_genomes_shared"].to_numpy(dtype=bool)] = "#08519c"

    bars = ax_bar.bar(x, counts, color=colors)
    ax_bar.set_ylabel("Unique gene names")
    all_shared_count = int(
        plot_df.loc[plot_df["all_genomes_shared"], "gene_count"].iloc[0]
    ) if plot_df["all_genomes_shared"].any() else 0
    ax_bar.set_title(
        f"{mode.capitalize()} gene-name intersections"
        #f"all genomes shared = {all_shared_count:,}"
    )
    ax_bar.bar_label(bars, padding=2, fontsize=6, rotation=90)
    style_axis(ax_bar)
    plt.setp(ax_bar.get_xticklabels(), visible=False)
    add_panel_label(ax_bar, panel_label)

    for y, col in enumerate(status_cols):
        ax_matrix.scatter(
            x,
            np.full(n_patterns, y),
            s=15,
            color="#d9d9d9",
            zorder=1,
        )
        present = plot_df[col].to_numpy(dtype=bool)
        ax_matrix.scatter(
            x[present],
            np.full(present.sum(), y),
            s=24,
            color="#252525",
            zorder=2,
        )
        for xi, is_present in zip(x, present):
            if not is_present:
                continue
            members = np.where(plot_df.iloc[xi][status_cols].to_numpy(dtype=bool))[0]
            if len(members) > 1:
                ax_matrix.plot(
                    [xi, xi],
                    [members.min(), members.max()],
                    color="#252525",
                    linewidth=0.9,
                    zorder=1.5,
                )

    ax_matrix.set_yticks(range(n_sets))
    ax_matrix.set_yticklabels(labels, fontsize=7)
    ax_matrix.invert_yaxis()
    ax_matrix.set_xlabel("Genome-set intersections")
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels(
        ["All\ngenomes" if shared else "" for shared in plot_df["all_genomes_shared"]],
        fontsize=7,
    )
    ax_matrix.tick_params(axis="x", length=0)
    ax_matrix.spines["top"].set_visible(False)
    ax_matrix.spines["right"].set_visible(False)
    ax_matrix.spines["left"].set_visible(False)
    ax_matrix.tick_params(axis="y", length=0)


def plot_combined_paper_figure(
    df: pd.DataFrame,
    gene_summary: pd.DataFrame,
    status_cols: list[str],
    strict_table: pd.DataFrame,
    basic_table: pd.DataFrame,
    outdir: Path,
    top_n: int,
) -> None:
    fig = plt.figure(figsize=(15, 11.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.55], hspace=0.35, wspace=0.28)

    ax_a = fig.add_subplot(gs[0, 0])
    plot_transcript_result_sets(ax_a, df, status_cols)

    ax_b = fig.add_subplot(gs[0, 1])
    plot_unique_gene_counts(ax_b, df, gene_summary, status_cols)

    draw_upset_panel(fig, gs[1, 0], strict_table, status_cols, "strict", top_n, "c")
    draw_upset_panel(fig, gs[1, 1], basic_table, status_cols, "basic", top_n, "d")

    savefig(fig, outdir / "alignment_results_overview")


def main() -> None:
    args = parse_args()
    in_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_path, sep="\t", dtype=str)
    status_cols = status_columns(df)

    for col in ["strict_pass_count", "basic_pass_count", "genome_count"]:
        df[col] = pd.to_numeric(df[col], errors="raise")

    if args.strip_gene_version:
        df["gene_name"] = df["gene_name"].str.replace(r"\.[0-9]+$", "", regex=True)

    gene_summary = build_gene_summary(df, status_cols)
    write_summary_tables(df, gene_summary, status_cols, outdir)

    intersection_tables = {}
    for mode in ["strict", "basic"]:
        membership = gene_membership(df, status_cols, mode)
        table = intersection_table(membership, status_cols)
        intersection_tables[mode] = table
        plot_upset(table, status_cols, outdir, mode, args.top_intersections)

    plot_combined_paper_figure(
        df,
        gene_summary,
        status_cols,
        intersection_tables["strict"],
        intersection_tables["basic"],
        outdir,
        args.top_intersections,
    )

    print(f"Wrote plots and summary tables to: {outdir}")
    print(f"Input transcripts: {df['transcript_id'].nunique()}")
    print(f"Unique gene names: {df['gene_name'].nunique()}")


if __name__ == "__main__":
    main()
