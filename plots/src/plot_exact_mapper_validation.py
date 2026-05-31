#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_EXACT = Path("../../kmer_results/exact_mapper/ZNF667_AS1_Mmul10_k35_kmer_locations.summary.tsv")
DEFAULT_OUTDIR = Path("ZNF667_results")
DEFAULT_GENE = "ZNF667-AS1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an exact-matching validation plot for ZNF667-AS1 against Mmul10."
    )
    parser.add_argument(
        "--exact-summary",
        type=Path,
        default=DEFAULT_EXACT,
        help="Transcript-level exact-matching summary TSV.",
    )
    parser.add_argument(
        "--gene",
        default=DEFAULT_GENE,
        help="Gene name to plot.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory for plots and TSV.",
    )
    return parser.parse_args()


def load_exact(path: Path, gene: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
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
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SystemExit(f"Missing exact summary column(s): {', '.join(missing)}")

    sub = df.copy()
    if gene:
        sub = sub[sub["transcript_name"].str.startswith(gene)].copy()
    if sub.empty:
        raise SystemExit(f"No rows for gene {gene!r} found in {path}")

    sub = sub.sort_values(
        ["matched_kmer_ratio", "transcript_name"],
        ascending=[False, True],
    ).reset_index(drop=True)
    return sub


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="left",
    )


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)


def annotate_bars(ax: plt.Axes, values: pd.Series, labels: list[str], x_offset: float) -> None:
    for idx, (value, label) in enumerate(zip(values, labels, strict=True)):
        ax.text(
            value + x_offset,
            idx,
            label,
            va="center",
            ha="left",
            fontsize=9,
        )


def build_figure(summary: pd.DataFrame, gene: str) -> plt.Figure:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(
        figsize=(13, max(8, 0.35 * len(summary) + 3)),
        constrained_layout=True,
    )

    y = np.arange(len(summary))
    ax.barh(
        y,
        summary["matched_kmer_ratio"],
        color="#d95f0e",
        edgecolor="white",
        linewidth=0.5,
        label="Exact matched k-mer ratio",
    )
    labels = [
        f"{int(matched)}/{int(total)}"
        for matched, total in zip(
            summary["matched_kmer_windows"],
            summary["total_kmer_windows"],
            strict=True,
        )
    ]
    annotate_bars(ax, summary["matched_kmer_ratio"], labels, 0.002)
    ax.set_yticks(y)
    ax.set_yticklabels(summary["transcript_name"])
    ax.invert_yaxis()
    ax.set_xlim(0, max(0.12, float(summary["matched_kmer_ratio"].max()) * 1.25))
    ax.set_xlabel("Matched k-mer ratio from exact string search")
    #ax.set_title(f"{gene} exact-match validation against Mmul10", pad=10)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    style_axis(ax)
    add_panel_label(ax, "A")

    fig.suptitle(
        f"{gene} exact-matching validation against Mmul10",
        fontsize=18,
        fontweight="bold",
        y=1.02,
    )
    return fig


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    summary = load_exact(args.exact_summary, args.gene)

    outstem = args.outdir / f"{args.gene.replace('-', '_')}_exact_validation"
    summary_path = outstem.with_suffix(".tsv")
    summary.to_csv(summary_path, sep="\t", index=False)

    fig = build_figure(summary, args.gene)
    fig.savefig(outstem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outstem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {outstem.with_suffix('.png')}")
    print(f"Wrote {outstem.with_suffix('.pdf')}")
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
