#!/usr/bin/env python3
"""
Create transcript- and gene-level synteny summary outputs.

The input transcript-level synteny table can be large, so this script streams it
and aggregates one row per lncRNA gene per assembly.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_INPUT_DIR = Path("alignment_vs_kmers/synteny_support")
ASSEMBLY_ORDER = [
    "Clintptrv2",
    "mPanTro3v20pri",
    "panpan11",
    "mPanPan1v20pri",
    "MhudibluPPAv0",
    "Mmul10",
    "T2TMMU8v20",
]

ASSEMBLY_LABELS = {
    "Clintptrv2": "Chimp Clint",
    "mPanTro3v20pri": "Chimp mPanTro3",
    "panpan11": "Bonobo panpan1.1",
    "mPanPan1v20pri": "Bonobo mPanPan1",
    "MhudibluPPAv0": "Bonobo Mhudiblu",
    "Mmul10": "Rhesus Mmul10",
    "T2TMMU8v20": "Rhesus T2T-MMU8",
}

SPECIES_COLORS = {
    "chimpanzee": "#4477AA",
    "bonobo": "#228833",
    "macaque": "#AA3377",
}


@dataclass
class GeneStats:
    assembly: str
    species_group: str
    gene_id: str
    gene_name: str
    transcript_rows: int = 0
    transcripts_with_expected_neighbors: int = 0
    transcripts_with_target_neighbors: int = 0
    transcripts_with_any_match: int = 0
    transcripts_with_two_or_more_matches: int = 0
    max_matched_neighbor_count: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create gene-level synteny summaries and a two-panel summary figure."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help="Directory containing high_conf_lncRNA_synteny_support.tsv.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        help="Output directory. Default: INPUT_DIR/gene_level_summary.",
    )
    parser.add_argument(
        "--prefix",
        default="synteny_support",
        help="Output filename prefix.",
    )
    return parser.parse_args()


def truthy(value: str) -> bool:
    return value.strip().lower() == "true"


def percent(numerator: int | float, denominator: int | float) -> float:
    if denominator == 0:
        return 0.0
    return 100.0 * numerator / denominator


def assembly_sort_key(assembly: str) -> tuple[int, str]:
    try:
        return ASSEMBLY_ORDER.index(assembly), assembly
    except ValueError:
        return len(ASSEMBLY_ORDER), assembly


def load_transcript_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df["any_match_percent"] = (
        100.0 * df["rows_with_any_synteny_match"] / df["high_conf_transcripts"]
    )
    df["two_or_more_percent"] = (
        100.0
        * df["rows_with_two_or_more_synteny_matches"]
        / df["high_conf_transcripts"]
    )
    df["target_neighbor_percent"] = (
        100.0 * df["rows_with_target_pc_neighbors"] / df["high_conf_transcripts"]
    )
    df["expected_neighbor_percent"] = (
        100.0
        * df["rows_with_expected_one2one_neighbors"]
        / df["high_conf_transcripts"]
    )
    df["assembly_order"] = df["assembly"].map(lambda value: assembly_sort_key(value)[0])
    return df.sort_values(["assembly_order", "assembly"]).reset_index(drop=True)


def build_gene_summary(path: Path) -> pd.DataFrame:
    genes: dict[tuple[str, str, str], GeneStats] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["assembly"], row["gene_id"], row["gene_name"])
            stats = genes.get(key)
            if stats is None:
                stats = GeneStats(
                    assembly=row["assembly"],
                    species_group=row["species_group"],
                    gene_id=row["gene_id"],
                    gene_name=row["gene_name"],
                )
                genes[key] = stats
            stats.transcript_rows += 1
            if int(row["expected_neighbor_count"]):
                stats.transcripts_with_expected_neighbors += 1
            if row["target_pc_neighbors"]:
                stats.transcripts_with_target_neighbors += 1
            if truthy(row["synteny_supported_any"]):
                stats.transcripts_with_any_match += 1
            if truthy(row["synteny_supported_two_or_more"]):
                stats.transcripts_with_two_or_more_matches += 1
            stats.max_matched_neighbor_count = max(
                stats.max_matched_neighbor_count,
                int(row["matched_neighbor_count"]),
            )

    rows = []
    for stats in genes.values():
        rows.append(
            {
                "assembly": stats.assembly,
                "species_group": stats.species_group,
                "gene_id": stats.gene_id,
                "gene_name": stats.gene_name,
                "transcript_rows": stats.transcript_rows,
                "transcripts_with_expected_neighbors": (
                    stats.transcripts_with_expected_neighbors
                ),
                "transcripts_with_target_neighbors": (
                    stats.transcripts_with_target_neighbors
                ),
                "transcripts_with_any_match": stats.transcripts_with_any_match,
                "transcripts_with_two_or_more_matches": (
                    stats.transcripts_with_two_or_more_matches
                ),
                "any_match_transcript_fraction": (
                    stats.transcripts_with_any_match / stats.transcript_rows
                ),
                "two_or_more_transcript_fraction": (
                    stats.transcripts_with_two_or_more_matches / stats.transcript_rows
                ),
                "max_matched_neighbor_count": stats.max_matched_neighbor_count,
                "gene_supported_any_transcript": (
                    stats.transcripts_with_any_match > 0
                ),
                "gene_supported_two_or_more_any_transcript": (
                    stats.transcripts_with_two_or_more_matches > 0
                ),
                "gene_supported_any_majority_transcripts": (
                    stats.transcripts_with_any_match / stats.transcript_rows >= 0.5
                ),
                "gene_supported_two_or_more_majority_transcripts": (
                    stats.transcripts_with_two_or_more_matches
                    / stats.transcript_rows
                    >= 0.5
                ),
            }
        )

    df = pd.DataFrame(rows)
    df["assembly_order"] = df["assembly"].map(lambda value: assembly_sort_key(value)[0])
    return df.sort_values(["assembly_order", "gene_name", "gene_id"]).drop(
        columns=["assembly_order"]
    )


def summarize_gene_by_assembly(gene_df: pd.DataFrame) -> pd.DataFrame:
    grouped = gene_df.groupby(["assembly", "species_group"], sort=False)
    summary = grouped.agg(
        high_conf_genes=("gene_id", "nunique"),
        gene_rows=("gene_id", "size"),
        genes_with_any_synteny_match=("gene_supported_any_transcript", "sum"),
        genes_with_two_or_more_synteny_match=(
            "gene_supported_two_or_more_any_transcript",
            "sum",
        ),
        genes_majority_transcripts_any_match=(
            "gene_supported_any_majority_transcripts",
            "sum",
        ),
        genes_majority_transcripts_two_or_more=(
            "gene_supported_two_or_more_majority_transcripts",
            "sum",
        ),
        median_any_match_transcript_fraction=(
            "any_match_transcript_fraction",
            "median",
        ),
        median_two_or_more_transcript_fraction=(
            "two_or_more_transcript_fraction",
            "median",
        ),
    ).reset_index()
    summary["genes_with_any_synteny_match_percent"] = (
        100.0
        * summary["genes_with_any_synteny_match"]
        / summary["high_conf_genes"]
    )
    summary["genes_with_two_or_more_synteny_match_percent"] = (
        100.0
        * summary["genes_with_two_or_more_synteny_match"]
        / summary["high_conf_genes"]
    )
    summary["genes_majority_transcripts_any_match_percent"] = (
        100.0
        * summary["genes_majority_transcripts_any_match"]
        / summary["high_conf_genes"]
    )
    summary["genes_majority_transcripts_two_or_more_percent"] = (
        100.0
        * summary["genes_majority_transcripts_two_or_more"]
        / summary["high_conf_genes"]
    )
    summary["assembly_order"] = summary["assembly"].map(
        lambda value: assembly_sort_key(value)[0]
    )
    return summary.sort_values(["assembly_order", "assembly"]).drop(
        columns=["assembly_order"]
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.1,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        ha="left",
        va="top",
    )


def style_percent_axis(ax: plt.Axes, xlabel: str) -> None:
    ax.set_xlim(0, 105)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_summary(
    transcript_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
    outbase: Path,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(13.5, 6.2),
        gridspec_kw={"wspace": 0.45},
    )

    y = np.arange(len(transcript_summary))
    labels = [
        ASSEMBLY_LABELS.get(assembly, assembly)
        for assembly in transcript_summary["assembly"]
    ]
    colors = [
        SPECIES_COLORS.get(species, "#777777")
        for species in transcript_summary["species_group"]
    ]

    ax = axes[0]
    ax.barh(
        y - 0.18,
        transcript_summary["any_match_percent"],
        height=0.34,
        color=colors,
        alpha=0.9,
        label=">=1 matched neighbor",
    )
    ax.barh(
        y + 0.18,
        transcript_summary["two_or_more_percent"],
        height=0.34,
        color=colors,
        alpha=0.45,
        label=">=2 matched neighbors",
    )
    for index, row in transcript_summary.iterrows():
        ax.text(
            row["any_match_percent"] + 1.0,
            index - 0.18,
            f"{row['any_match_percent']:.1f}",
            va="center",
            fontsize=7,
        )
        ax.text(
            row["two_or_more_percent"] + 1.0,
            index + 0.18,
            f"{row['two_or_more_percent']:.1f}",
            va="center",
            fontsize=7,
        )
    ax.set_yticks(y, labels=labels)
    ax.invert_yaxis()
    ax.set_title("Transcript-level synteny support")
    style_percent_axis(ax, "Supported transcripts (%)")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        frameon=False,
        fontsize=8,
        ncol=2,
    )
    add_panel_label(ax, "a")

    ax = axes[1]
    y = np.arange(len(gene_summary))
    labels = [
        ASSEMBLY_LABELS.get(assembly, assembly)
        for assembly in gene_summary["assembly"]
    ]
    colors = [
        SPECIES_COLORS.get(species, "#777777")
        for species in gene_summary["species_group"]
    ]
    ax.barh(
        y - 0.18,
        gene_summary["genes_with_any_synteny_match_percent"],
        height=0.34,
        color=colors,
        alpha=0.9,
        label=">=1 transcript supported",
    )
    ax.barh(
        y + 0.18,
        gene_summary["genes_majority_transcripts_two_or_more_percent"],
        height=0.34,
        color=colors,
        alpha=0.45,
        label="majority transcripts, >=2 neighbors",
    )
    for index, row in gene_summary.iterrows():
        ax.text(
            row["genes_with_any_synteny_match_percent"] + 1.0,
            index - 0.18,
            f"{row['genes_with_any_synteny_match_percent']:.1f}",
            va="center",
            fontsize=7,
        )
        ax.text(
            row["genes_majority_transcripts_two_or_more_percent"] + 1.0,
            index + 0.18,
            f"{row['genes_majority_transcripts_two_or_more_percent']:.1f}",
            va="center",
            fontsize=7,
        )
    ax.set_yticks(y, labels=labels)
    ax.invert_yaxis()
    ax.set_title("Gene-level synteny support")
    style_percent_axis(ax, "Supported genes (%)")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        frameon=False,
        fontsize=8,
        ncol=2,
    )
    add_panel_label(ax, "b")

    fig.suptitle(
        "Synteny support for high-confidence alignment/k-mer lncRNA candidates",
        y=0.995,
    )
    fig.savefig(outbase.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    outdir = args.outdir or input_dir / "gene_level_summary"
    outdir.mkdir(parents=True, exist_ok=True)

    transcript_summary_path = input_dir / "high_conf_lncRNA_synteny_summary_by_assembly.tsv"
    transcript_table_path = input_dir / "high_conf_lncRNA_synteny_support.tsv"
    transcript_summary = load_transcript_summary(transcript_summary_path)
    gene_df = build_gene_summary(transcript_table_path)
    gene_summary = summarize_gene_by_assembly(gene_df)

    gene_table_path = outdir / f"{args.prefix}_gene_level.tsv"
    gene_summary_path = outdir / f"{args.prefix}_gene_level_by_assembly.tsv"
    transcript_summary_out = outdir / f"{args.prefix}_transcript_level_by_assembly.tsv"
    figure_base = outdir / f"{args.prefix}_summary_ab"

    gene_df.to_csv(gene_table_path, sep="\t", index=False)
    gene_summary.to_csv(gene_summary_path, sep="\t", index=False)
    transcript_summary.drop(columns=["assembly_order"], errors="ignore").to_csv(
        transcript_summary_out,
        sep="\t",
        index=False,
    )
    plot_summary(transcript_summary, gene_summary, figure_base)

    print(f"Wrote {gene_table_path}")
    print(f"Wrote {gene_summary_path}")
    print(f"Wrote {transcript_summary_out}")
    print(f"Wrote {figure_base.with_suffix('.png')}")
    print(f"Wrote {figure_base.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
