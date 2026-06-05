import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_INPUT = Path("../../negative_set_results/negative_set.tsv")
DEFAULT_OUTDIR = Path("../../negative_set_results")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot negative-set generation and IBF search results."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="negative_set.tsv produced by the negative_set executable.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Directory where figures and the summary TSV will be written.",
    )
    return parser.parse_args()


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 15,
            "axes.labelsize": 12,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_input(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {
        "id",
        "chromosome",
        "start",
        "end",
        "length",
        "total_kmers",
        "matched_kmer_positions",
        "match_count",
        "match_fraction",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SystemExit(f"Missing required column(s): {', '.join(missing)}")

    numeric_cols = [
        "start",
        "end",
        "length",
        "total_kmers",
        "matched_kmer_positions",
        "match_count",
        "match_fraction",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if df[numeric_cols].isna().any().any():
        raise SystemExit("Input TSV contains non-numeric values in expected numeric columns.")

    return df


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    summary = pd.DataFrame(
        [
            {
                "metric": "rows",
                "value": int(len(df)),
            },
            {
                "metric": "mean_length",
                "value": float(df["length"].mean()),
            },
            {
                "metric": "mean_total_kmers",
                "value": float(df["total_kmers"].mean()),
            },
            {
                "metric": "mean_matched_positions",
                "value": float(df["matched_kmer_positions"].mean()),
            },
            {
                "metric": "mean_match_count",
                "value": float(df["match_count"].mean()),
            },
            {
                "metric": "mean_match_fraction",
                "value": float(df["match_fraction"].mean()),
            },
            {
                "metric": "max_match_fraction",
                "value": float(df["match_fraction"].max()),
            },
        ]
    )
    return summary


def plot_results(df: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 6.5), constrained_layout=True)

    scatter = ax.scatter(
        df["length"],
        df["match_fraction"],
        s=24,
        alpha=0.65,
        c=df["matched_kmer_positions"],
        cmap="magma",
        edgecolors="none",
    )
    ax.set_xlabel("sampled sequence length")
    ax.set_ylabel("fraction of k-mer positions with at least one IBF hit")
    ax.set_title("Negative sequences: length versus IBF hit fraction")
    ax.grid(True, color="#e6e6e6", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.colorbar(scatter, ax=ax, label="matched k-mer positions")

    #fig.suptitle("Negative set IBF search results", fontsize=16, y=1.01)

    fig.savefig(outdir / "negative_set_results.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "negative_set_results.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    df = load_input(args.input)
    summary = build_summary(df)
    summary.to_csv(args.output_dir / "negative_set_plot_summary.tsv", sep="\t", index=False)
    plot_results(df, args.output_dir)

    print(f"Wrote {args.output_dir / 'negative_set_results.png'}")
    print(f"Wrote {args.output_dir / 'negative_set_results.pdf'}")
    print(f"Wrote {args.output_dir / 'negative_set_plot_summary.tsv'}")


if __name__ == "__main__":
    main()
