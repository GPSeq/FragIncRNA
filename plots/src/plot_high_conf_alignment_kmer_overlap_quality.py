#!/usr/bin/env python3
"""Plot overlap between high-confidence alignments and k-mer support."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


BLUE = "#2F6C9E"
TEAL = "#1B9E77"
GOLD = "#D89C28"
RED = "#C44E52"
PURPLE = "#756BB1"
INK = "#252525"
GRID = "#E8E8E8"


ASSEMBLY_META = {
    "H9T2Thap1": ("Human", "T2T/H9 hap1", "high"),
    "H9T2Thap2": ("Human", "T2T/H9 hap2", "high"),
    "matHPRCF2": ("Human", "HPRC maternal", "high"),
    "patHPRCf2": ("Human", "HPRC paternal", "high"),
    "mPanTro3v20pri": ("Chimpanzee", "T2T/primary", "high"),
    "Clintptrv2": ("Chimpanzee", "Clint PTRv2", "reference"),
    "mPanPan1v20pri": ("Bonobo", "T2T/primary", "high"),
    "MhudibluPPAv0": ("Bonobo", "Mhudiblu PPA", "reference"),
    "panpan11": ("Bonobo", "panpan1.1", "reference"),
    "T2TMMU8v20": ("Rhesus macaque", "T2T-MMU8", "high"),
    "Mmul10": ("Rhesus macaque", "Mmul_10", "reference"),
}

ASSEMBLY_TO_KMER_COLUMN = {
    "patHPRCf2": "GCA_018503265.2_NA19240_pat_hprc_f2_genomic",
    "matHPRCF2": "GCA_018503275.2_NA19240_mat_hprc_f2_genomic",
    "H9T2Thap1": "GCA_054883195.1_H9_T2T.hap1_genomic",
    "H9T2Thap2": "GCA_054883265.1_H9_T2T.hap2_genomic",
    "panpan11": "GCF_000258655.2_panpan1.1_genomic",
    "Clintptrv2": "GCF_002880755.1_Clint_PTRv2_genomic",
    "Mmul10": "GCF_003339765.1_Mmul_10_genomic",
    "MhudibluPPAv0": "GCF_013052645.1_Mhudiblu_PPA_v0_genomic",
    "mPanTro3v20pri": "GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic",
    "mPanPan1v20pri": "GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic",
    "T2TMMU8v20": "GCF_049350105.2_T2T-MMU8v2.0_genomic",
}

SPECIES_ORDER = {
    "Human": 0,
    "Chimpanzee": 1,
    "Bonobo": 2,
    "Rhesus macaque": 3,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create overlap and assembly-quality plots for high-confidence "
            "alignment candidates and all-genome k-mer support."
        )
    )
    parser.add_argument(
        "--alignment-pass-counts",
        type=Path,
        default=Path(
            "alignment_results/high_confident_candidate_orthologs/"
            "high_confident_candidate_ortholog_pass_counts.tsv"
        ),
        help="High-confidence candidate ortholog pass-count table.",
    )
    parser.add_argument(
        "--alignment-by-genome",
        type=Path,
        default=Path(
            "alignment_results/high_confident_candidate_orthologs/"
            "high_confident_candidate_ortholog_summary_by_genome.tsv"
        ),
        help="Per-genome high-confidence candidate ortholog summary table.",
    )
    parser.add_argument(
        "--alignment-status-matrix",
        type=Path,
        default=Path(
            "alignment_results/high_confident_candidate_orthologs/"
            "high_confident_candidate_ortholog_status_matrix.tsv"
        ),
        help="Per-transcript high-confidence candidate ortholog status matrix.",
    )
    parser.add_argument(
        "--kmer-pass-counts",
        type=Path,
        default=Path("kmer_results/ibf_new/qc/kmer_transcript_pass_counts_by_k.tsv"),
        help="K-mer transcript pass-count table.",
    )
    parser.add_argument(
        "--kmer-status-matrix",
        type=Path,
        help=(
            "Per-transcript k-mer status matrix. Default: "
            "kmer_results/ibf_new/qc/status_matrices/kmer_transcript_status_matrix_k<K>.tsv"
        ),
    )
    parser.add_argument(
        "--genome-n-content",
        type=Path,
        default=Path("../genomes/genome_N_content.tsv"),
        help="Genome N-content table from check_genome_N_content.sh.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("alignment_vs_kmers/high_conf"),
        help="Output directory for plots and summary tables.",
    )
    parser.add_argument(
        "--k-size",
        type=int,
        default=24,
        help="k-mer size to use for the Venn-style overlap panel.",
    )
    parser.add_argument(
        "--kmer-mode",
        choices=("strict", "basic"),
        default="strict",
        help="All-genome k-mer support definition for the overlap panel.",
    )
    return parser.parse_args()


def configure_matplotlib(outdir: Path) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(outdir / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
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
            "figure.dpi": 160,
            "savefig.dpi": 450,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_tsv(path: Path, **kwargs) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required file not found: {path}")
    return pd.read_csv(path, sep="\t", **kwargs)


def require_columns(df: pd.DataFrame, columns: set[str], path: Path) -> None:
    missing = sorted(columns.difference(df.columns))
    if missing:
        raise ValueError(f"Missing column(s) in {path}: {', '.join(missing)}")


def assembly_name(genome_column: str) -> str:
    return genome_column.replace("human_lncRNA_vs_", "")


def n_content_by_assembly(path: Path) -> dict[str, float]:
    if not path.is_file():
        return {}
    table = read_tsv(path)
    require_columns(table, {"file", "n_percent"}, path)
    lookup = {}
    for row in table.itertuples(index=False):
        assembly = str(row.file).split("_G", 1)[0]
        lookup[assembly] = float(row.n_percent)
    return lookup


def load_alignment_sets(path: Path) -> tuple[pd.DataFrame, set[str]]:
    columns = {
        "transcript_id",
        "gene_id",
        "gene_name",
        "high_confident_pass_count",
        "genome_count",
    }
    df = read_tsv(path)
    require_columns(df, columns, path)
    high_mask = df["high_confident_pass_count"].eq(df["genome_count"])
    return df, set(df.loc[high_mask, "transcript_id"].astype(str))


def load_kmer_sets(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "k_size",
        "transcript_id",
        "gene_id",
        "gene_name",
        "strict_pass_count",
        "basic_pass_count",
        "genome_count",
    ]
    df = read_tsv(path, usecols=columns)
    require_columns(df, set(columns), path)
    rows = []
    for k_size, group in df.groupby("k_size", sort=True):
        genome_count = group["genome_count"]
        strict = group["strict_pass_count"].eq(genome_count)
        basic = group["basic_pass_count"].eq(genome_count)
        rows.append(
            {
                "k_size": int(k_size),
                "kmer_strict_all": int(strict.sum()),
                "kmer_basic_all": int(basic.sum()),
                "kmer_strict_unique_genes": int(group.loc[strict, "gene_id"].nunique()),
                "kmer_basic_unique_genes": int(group.loc[basic, "gene_id"].nunique()),
            }
        )
    return df, pd.DataFrame(rows)


def kmer_supported_transcripts(kmer: pd.DataFrame, k_size: int, mode: str) -> set[str]:
    group = kmer[kmer["k_size"].eq(k_size)]
    if group.empty:
        raise ValueError(f"No k-mer rows found for k={k_size}")
    pass_column = "strict_pass_count" if mode == "strict" else "basic_pass_count"
    mask = group[pass_column].eq(group["genome_count"])
    return set(group.loc[mask, "transcript_id"].astype(str))


def build_overlap_by_k(
    alignment: pd.DataFrame,
    alignment_high_set: set[str],
    kmer: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    total = int(alignment["transcript_id"].nunique())
    high = set(alignment_high_set)
    for k_size, group in kmer.groupby("k_size", sort=True):
        strict = set(
            group.loc[
                group["strict_pass_count"].eq(group["genome_count"]), "transcript_id"
            ].astype(str)
        )
        basic = set(
            group.loc[
                group["basic_pass_count"].eq(group["genome_count"]), "transcript_id"
            ].astype(str)
        )
        for mode, kmer_set in (("strict", strict), ("basic", basic)):
            overlap = high & kmer_set
            union = high | kmer_set
            rows.append(
                {
                    "k_size": int(k_size),
                    "kmer_mode": mode,
                    "total_transcripts": total,
                    "alignment_high_confident_all": len(high),
                    "kmer_all_genomes": len(kmer_set),
                    "overlap": len(overlap),
                    "alignment_only": len(high - kmer_set),
                    "kmer_only": len(kmer_set - high),
                    "neither": total - len(union),
                    "overlap_fraction_of_alignment": len(overlap) / len(high),
                    "overlap_fraction_of_kmer": len(overlap) / len(kmer_set)
                    if kmer_set
                    else np.nan,
                    "jaccard_index": len(overlap) / len(union) if union else np.nan,
                }
            )
    return pd.DataFrame(rows)


def build_assembly_summary(path: Path, n_lookup: dict[str, float]) -> pd.DataFrame:
    df = read_tsv(path)
    require_columns(
        df,
        {
            "genome",
            "total_transcripts",
            "high_confident_transcripts",
            "high_confident_percent",
            "high_confident_unique_genes",
        },
        path,
    )
    df["assembly"] = df["genome"].map(assembly_name)
    df["species"] = df["assembly"].map(
        lambda x: ASSEMBLY_META.get(x, ("Other", x, "other"))[0]
    )
    df["display_name"] = df["assembly"].map(
        lambda x: ASSEMBLY_META.get(x, ("Other", x, "other"))[1]
    )
    df["quality_class"] = df["assembly"].map(
        lambda x: ASSEMBLY_META.get(x, ("Other", x, "other"))[2]
    )
    df["n_percent"] = df["assembly"].map(n_lookup).fillna(0.0)
    df["species_order"] = df["species"].map(SPECIES_ORDER).fillna(99).astype(int)
    df["quality_order"] = df["quality_class"].map({"high": 0, "reference": 1}).fillna(2)
    return df.sort_values(
        ["species_order", "quality_order", "high_confident_percent"],
        ascending=[True, True, False],
    ).reset_index(drop=True)


def representative_species_summary(assembly: pd.DataFrame) -> pd.DataFrame:
    return (
        assembly.sort_values(
            ["species_order", "quality_order", "high_confident_percent"],
            ascending=[True, True, False],
        )
        .groupby("species", sort=False, as_index=False)
        .head(1)
        .sort_values("species_order")
        .reset_index(drop=True)
    )


def build_per_assembly_overlap(
    alignment_status_path: Path,
    kmer_status_path: Path,
    n_lookup: dict[str, float],
    mode: str,
) -> pd.DataFrame:
    alignment = read_tsv(alignment_status_path)
    kmer = read_tsv(kmer_status_path)
    require_columns(
        alignment,
        {"transcript_id", "gene_id", "gene_name"},
        alignment_status_path,
    )
    require_columns(kmer, {"transcript_id"}, kmer_status_path)

    alignment_cols = [col for col in alignment.columns if col.startswith("human_lncRNA_vs_")]
    if not alignment_cols:
        raise ValueError(f"No alignment status columns found in {alignment_status_path}")

    merged = alignment.merge(
        kmer,
        on="transcript_id",
        how="inner",
        suffixes=("_alignment", "_kmer"),
        validate="one_to_one",
    )
    total = int(merged["transcript_id"].nunique())
    rows = []
    for align_col in alignment_cols:
        assembly = assembly_name(align_col)
        if assembly not in ASSEMBLY_TO_KMER_COLUMN:
            continue
        kmer_col = ASSEMBLY_TO_KMER_COLUMN[assembly]
        if kmer_col not in merged.columns:
            raise ValueError(f"Missing k-mer status column for {assembly}: {kmer_col}")
        species, display_name, quality_class = ASSEMBLY_META.get(
            assembly, ("Other", assembly, "reference")
        )
        alignment_high = merged[align_col].eq("HIGH_CONFIDENT")
        if mode == "strict":
            kmer_pass = merged[kmer_col].eq("PASS_STRICT_KMER")
        else:
            kmer_pass = merged[kmer_col].isin(["PASS_STRICT_KMER", "PASS_BASIC_KMER"])
        overlap = alignment_high & kmer_pass
        rows.append(
            {
                "genome": align_col,
                "assembly": assembly,
                "species": species,
                "display_name": display_name,
                "quality_class": quality_class,
                "total_transcripts": total,
                "alignment_high_confident_transcripts": int(alignment_high.sum()),
                "kmer_pass_transcripts": int(kmer_pass.sum()),
                "overlap_transcripts": int(overlap.sum()),
                "overlap_percent_total": 100 * float(overlap.mean()),
                "overlap_unique_genes": int(merged.loc[overlap, "gene_id_alignment"].nunique()),
                "n_percent": n_lookup.get(assembly, 0.0),
            }
        )
    summary = pd.DataFrame(rows)
    summary["species_order"] = summary["species"].map(SPECIES_ORDER).fillna(99).astype(int)
    summary["quality_order"] = summary["quality_class"].map({"high": 0, "reference": 1}).fillna(2)
    return summary.sort_values(
        ["species_order", "quality_order", "overlap_transcripts"],
        ascending=[True, True, False],
    ).reset_index(drop=True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.12,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
        color=INK,
    )


def human_number(value: float, _position: int | None = None) -> str:
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.0f}K"
    return f"{value:.0f}"


def draw_assembly_panel(ax: plt.Axes, assembly: pd.DataFrame, k_size: int) -> None:
    colors = {
        "Human": "#3B6EA8",
        "Chimpanzee": "#4FA35B",
        "Bonobo": "#E88C2A",
        "Rhesus macaque": "#A66AA0",
    }
    x = np.arange(len(assembly))
    bars = []
    for xpos, row in zip(x, assembly.itertuples(index=False)):
        quality_is_high = row.quality_class == "high"
        bar = ax.bar(
            xpos,
            row.overlap_transcripts,
            color=colors.get(row.species, "#999999"),
            alpha=0.95 if quality_is_high else 0.48,
            edgecolor="#222222" if quality_is_high else "#555555",
            linewidth=1.5 if quality_is_high else 1.0,
            hatch="" if quality_is_high else "///",
            zorder=2,
        )[0]
        bars.append(bar)
    ax.set_xticks(x, assembly["display_name"], rotation=45, ha="right")
    ax.set_ylabel("Alignment-k-mer overlap transcripts")
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    ax.set_ylim(0, assembly["overlap_transcripts"].max() * 1.16)
    
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_axisbelow(True)

    label_offset = assembly["overlap_transcripts"].max() * 0.012
    for bar, row in zip(bars, assembly.itertuples(index=False)):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + label_offset,
            f"{int(row.overlap_transcripts):,}",
            ha="center",
            va="bottom",
            fontsize=7,
            rotation=90,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.4},
            zorder=5,
        )

    ax_n = ax.twinx()
    ax_n.plot(
        x,
        assembly["n_percent"],
        color=RED,
        marker="D",
        linewidth=1.1,
        markersize=3.5,
        alpha=0.78,
        label="N content",
        zorder=1,
    )
    ax_n.set_ylabel("Assembly N content (%)", color=RED)
    ax_n.tick_params(axis="y", colors=RED)
    ax_n.set_ylim(0, max(35, assembly["n_percent"].max() * 1.6))
    ax_n.spines["top"].set_visible(False)

    species_handles = [
        Patch(facecolor=color, edgecolor="none", label=species)
        for species, color in colors.items()
        if species in set(assembly["species"])
    ]
    quality_handles = [
        Patch(facecolor="white", edgecolor="#222222", linewidth=1.5, label="High-quality/T2T/HPRC"),
        Patch(facecolor="white", edgecolor="#555555", linewidth=1.0, hatch="///", label="Non-T2T/reference"),
    ]
    ax.legend(
        handles=species_handles + quality_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.34),
        fontsize=7,
        frameon=False,
        ncol=3,
    )
    #add_panel_label(ax)


def overlap_transcript_ids_per_assembly(
    alignment_status_path: Path,
    kmer_status_path: Path,
    mode: str,
) -> dict[str, pd.DataFrame]:
    alignment = read_tsv(alignment_status_path)
    kmer = read_tsv(kmer_status_path)
    require_columns(
        alignment,
        {"transcript_id", "gene_id", "gene_name"},
        alignment_status_path,
    )
    require_columns(kmer, {"transcript_id"}, kmer_status_path)

    alignment_cols = [col for col in alignment.columns if col.startswith("human_lncRNA_vs_")]
    if not alignment_cols:
        raise ValueError(f"No alignment status columns found in {alignment_status_path}")

    merged = alignment.merge(
        kmer,
        on="transcript_id",
        how="inner",
        suffixes=("_alignment", "_kmer"),
        validate="one_to_one",
    )

    outputs: dict[str, pd.DataFrame] = {}
    for align_col in alignment_cols:
        assembly = assembly_name(align_col)
        if assembly not in ASSEMBLY_TO_KMER_COLUMN:
            continue
        kmer_col = ASSEMBLY_TO_KMER_COLUMN[assembly]
        if kmer_col not in merged.columns:
            raise ValueError(f"Missing k-mer status column for {assembly}: {kmer_col}")
        alignment_high = merged[align_col].eq("HIGH_CONFIDENT")
        if mode == "strict":
            kmer_pass = merged[kmer_col].eq("PASS_STRICT_KMER")
        else:
            kmer_pass = merged[kmer_col].isin(["PASS_STRICT_KMER", "PASS_BASIC_KMER"])
        overlap = alignment_high & kmer_pass
        outputs[assembly] = (
            merged.loc[overlap, ["transcript_id", "gene_id_alignment", "gene_name_alignment"]]
            .rename(
                columns={
                    "gene_id_alignment": "gene_id",
                    "gene_name_alignment": "gene_name",
                }
            )
            .sort_values(["gene_name", "transcript_id"])
            .reset_index(drop=True)
        )
    return outputs


def write_outputs(
    outdir: Path,
    overlap_by_k: pd.DataFrame,
    assembly_overlap: pd.DataFrame,
    overlap_ids: dict[str, pd.DataFrame],
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    tables_dir = outdir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    overlap_by_k.to_csv(outdir / "high_conf_alignment_kmer_overlap_by_k.tsv", sep="\t", index=False)
    assembly_overlap.drop(columns=["species_order", "quality_order"]).to_csv(
        outdir / "high_conf_alignment_kmer_assembly_overlap_summary.tsv",
        sep="\t",
        index=False,
    )
    for assembly, table in overlap_ids.items():
        table.to_csv(tables_dir / f"{assembly}.tsv", sep="\t", index=False)


def make_figure(
    assembly_overlap: pd.DataFrame,
    k_size: int,
) -> plt.Figure:
    fig, axes = plt.subplots(
        1,
        1,
        figsize=(16.4, 6.1),
    )
    draw_assembly_panel(axes, assembly_overlap, k_size)
    fig.suptitle(f"Assembly-wise high-confidence alignment and k-mer overlap (k={k_size})", fontsize=13, y=0.98)
    return fig


def main() -> None:
    args = parse_args()
    configure_matplotlib(args.outdir)

    alignment, alignment_high_set = load_alignment_sets(args.alignment_pass_counts)
    kmer, _kmer_counts = load_kmer_sets(args.kmer_pass_counts)
    overlap_by_k = build_overlap_by_k(alignment, alignment_high_set, kmer)

    selected = overlap_by_k[
        overlap_by_k["k_size"].eq(args.k_size)
        & overlap_by_k["kmer_mode"].eq(args.kmer_mode)
    ]
    if selected.empty:
        available = overlap_by_k[["k_size", "kmer_mode"]].drop_duplicates()
        raise SystemExit(
            f"No overlap row for k={args.k_size}, mode={args.kmer_mode}. "
            f"Available combinations:\n{available.to_string(index=False)}"
        )
    n_lookup = n_content_by_assembly(args.genome_n_content)
    kmer_status_matrix = args.kmer_status_matrix or Path(
        f"kmer_results/ibf_new/qc/status_matrices/kmer_transcript_status_matrix_k{args.k_size}.tsv"
    )
    assembly_overlap = build_per_assembly_overlap(
        args.alignment_status_matrix,
        kmer_status_matrix,
        n_lookup,
        args.kmer_mode,
    )
    overlap_ids = overlap_transcript_ids_per_assembly(
        args.alignment_status_matrix,
        kmer_status_matrix,
        args.kmer_mode,
    )

    write_outputs(
        args.outdir,
        overlap_by_k,
        assembly_overlap,
        overlap_ids,
    )
    fig = make_figure(
        assembly_overlap,
        args.k_size,
    )
    prefix = args.outdir / "high_conf_alignment_kmer_overlap_quality"
    fig.savefig(prefix.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {prefix.with_suffix('.png')}")
    print(f"Wrote {prefix.with_suffix('.pdf')}")
    print(f"Wrote summary tables to {args.outdir}")


if __name__ == "__main__":
    main()
