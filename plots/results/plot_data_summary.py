# python plot_data_summary.py --input-genomes genomes --input-lncrna lncrna --output-dir data_summary_plots
# or if we have the summary file then: 
# python ./plot_data_summary.py --input-summary data_summary_plots/data_summary_counts.tsv --output-dir data_summary_plots
import argparse
import gzip
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


LNCRNA_CANDIDATE_NAMES = (
    "gencode.v49.lncRNA_transcripts.fa",
    "gencode.v49.lncRNA_transcripts.fa.gz",
    "gencode.v49.lncRNA_transcript.fa",
    "gencode.v49.lncRNA_transcript.fa.gz",
)

FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fas",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fas.gz",
)


@dataclass(frozen=True)
class FastaStats:
    path: Path
    label: str
    records: int
    total_bp: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create publication-ready summary plots and a TSV for genome FASTA "
            "references and a GENCODE lncRNA transcript FASTA."
        )
    )
    parser.add_argument(
        "--input-genomes",
        type=Path,
        help="Directory containing genome FASTA files, optionally gzip-compressed.",
    )
    parser.add_argument(
        "--input-lncrna",
        type=Path,
        help=(
            "Path to the lncRNA FASTA file, or a directory containing "
            "gencode.v49.lncRNA_transcripts.fa."
        ),
    )
    parser.add_argument(
        "--input-summary",
        type=Path,
        help=(
            "Existing data_summary_counts.tsv/data_summary_count.tsv file to reuse. "
            "If this file exists and is not empty, FASTA inputs are not required."
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory where plots and the TSV summary will be written.",
    )
    return parser.parse_args()


def open_maybe_gzip(path: Path) -> TextIO:
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def is_fasta_path(path: Path) -> bool:
    if path.name.startswith("._"):
        return False
    name = path.name.lower()
    return any(name.endswith(suffix) for suffix in FASTA_SUFFIXES)


def discover_genome_fastas(input_genomes: Path) -> list[Path]:
    if not input_genomes.exists():
        raise FileNotFoundError(f"Genome input does not exist: {input_genomes}")
    if input_genomes.is_file():
        if not is_fasta_path(input_genomes):
            raise ValueError(f"Genome input file does not look like FASTA: {input_genomes}")
        return [input_genomes]
    fastas = sorted(path for path in input_genomes.rglob("*") if path.is_file() and is_fasta_path(path))
    if not fastas:
        raise ValueError(f"No FASTA files found in genome directory: {input_genomes}")
    return fastas


def resolve_lncrna_path(input_lncrna: Path) -> Path:
    if not input_lncrna.exists():
        raise FileNotFoundError(f"lncRNA input does not exist: {input_lncrna}")
    if input_lncrna.is_file():
        if input_lncrna.name.startswith("._"):
            raise ValueError(f"Refusing AppleDouble metadata file as lncRNA FASTA: {input_lncrna}")
        return input_lncrna

    for candidate in LNCRNA_CANDIDATE_NAMES:
        path = input_lncrna / candidate
        if path.exists() and path.is_file():
            return path

    fastas = sorted(path for path in input_lncrna.rglob("*") if path.is_file() and is_fasta_path(path))
    if not fastas:
        raise ValueError(f"No lncRNA FASTA file found in directory: {input_lncrna}")
    return fastas[0]


def display_label(path: Path) -> str:
    name = path.name
    for suffix in (".fna.gz", ".fasta.gz", ".fas.gz", ".fa.gz", ".fna", ".fasta", ".fas", ".fa"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_genomic"):
        name = name[: -len("_genomic")]
    for marker in ("_GCF_", "_GCA_"):
        if marker in name:
            name = name.split(marker, 1)[0]
            break
    return name


def count_fasta(path: Path) -> FastaStats:
    records = 0
    total_bp = 0
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if not line:
                continue
            if line.startswith(">"):
                records += 1
            else:
                total_bp += len(line.strip())
    return FastaStats(path=path, label=display_label(path), records=records, total_bp=total_bp)


def write_summary_tsv(output_path: Path, genome_stats: list[FastaStats], lncrna_stats: FastaStats) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("category\tname\tcount_type\tcount\ttotal_bp\tsource_file\n")
        handle.write(
            "summary\ttotal_genome_references\tgenomes\t"
            f"{len(genome_stats)}\tNA\t{common_parent(genome_stats)}\n"
        )
        handle.write(
            "summary\ttotal_lncRNA_transcripts\tlncRNA_transcripts\t"
            f"{lncrna_stats.records}\t{lncrna_stats.total_bp}\t{lncrna_stats.path}\n"
        )
        for stat in genome_stats:
            handle.write(
                "genome\t"
                f"{stat.label}\tchromosomes_or_records\t{stat.records}\t{stat.total_bp}\t{stat.path}\n"
            )


def write_summary_outputs(output_dir: Path, genome_stats: list[FastaStats], lncrna_stats: FastaStats) -> None:
    write_summary_tsv(output_dir / "data_summary_counts.tsv", genome_stats, lncrna_stats)
    write_summary_tsv(output_dir / "data_summary_count.tsv", genome_stats, lncrna_stats)


def read_summary_tsv(input_path: Path) -> tuple[list[FastaStats], FastaStats]:
    import csv

    genome_stats: list[FastaStats] = []
    lncrna_stats: FastaStats | None = None

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"category", "name", "count_type", "count", "total_bp", "source_file"}
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            raise ValueError(
                f"Summary TSV is missing required columns: {', '.join(sorted(required_columns))}"
            )

        for row in reader:
            category = row["category"]
            count_type = row["count_type"]
            if category == "summary" and count_type == "lncRNA_transcripts":
                lncrna_stats = FastaStats(
                    path=Path(row["source_file"]),
                    label=row["name"],
                    records=parse_int_field(row["count"], row["name"], "count"),
                    total_bp=parse_optional_int(row["total_bp"]),
                )
            elif category == "genome":
                genome_stats.append(
                    FastaStats(
                        path=Path(row["source_file"]),
                        label=row["name"],
                        records=parse_int_field(row["count"], row["name"], "count"),
                        total_bp=parse_optional_int(row["total_bp"]),
                    )
                )

    if lncrna_stats is None:
        raise ValueError(f"No lncRNA summary row found in: {input_path}")
    if not genome_stats:
        raise ValueError(f"No genome rows found in: {input_path}")
    return genome_stats, lncrna_stats


def parse_int_field(value: str, row_name: str, column_name: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"Invalid integer in {column_name!r} for row {row_name!r}: {value!r}") from exc


def parse_optional_int(value: str) -> int:
    if value in {"", "NA", "nan", "NaN"}:
        return 0
    return int(value)


def common_parent(stats: Iterable[FastaStats]) -> str:
    paths = [str(stat.path.parent) for stat in stats]
    if not paths:
        return "NA"
    return os.path.commonpath(paths)


def require_matplotlib():
    try:
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required for plotting. Install it with: pip install matplotlib"
        ) from exc
    return plt, FuncFormatter


def human_number(value: float, _position: int | None = None) -> str:
    if value == 0:
        return "0"
    abs_value = abs(value)
    for threshold, suffix in ((1_000_000_000, "B"), (1_000_000, "M"), (1_000, "K")):
        if abs_value >= threshold:
            return f"{value / threshold:.1f}{suffix}"
    return f"{value:.0f}"


def save_figure(fig, output_base: Path) -> None:
    fig.savefig(output_base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")


def make_summary_plot(output_dir: Path, genome_stats: list[FastaStats], lncrna_stats: FastaStats) -> None:
    plt, FuncFormatter = require_matplotlib()

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
            "figure.dpi": 150,
        }
    )

    sorted_by_name = sorted(genome_stats, key=lambda stat: stat.label.lower())
    figure_height = max(4.6, 0.26 * len(sorted_by_name) + 2.2)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(11.5, figure_height),
        gridspec_kw={"width_ratios": [0.8, 1.6], "wspace": 0.42},
        constrained_layout=False,
    )

    ax = axes[0]
    overview_labels = ["Genome\nreferences", "lncRNA\ntranscripts"]
    overview_counts = [len(genome_stats), lncrna_stats.records]
    colors = ["#4477AA", "#228833"]
    ax.bar(overview_labels, overview_counts, color=colors, width=0.62)
    ax.set_ylabel("Count")
    ax.set_title("Dataset overview", pad=10)
    ax.yaxis.set_major_formatter(FuncFormatter(human_number))
    if max(overview_counts) / max(1, min(count for count in overview_counts if count > 0)) > 50:
        ax.set_yscale("log")
        ax.set_ylabel("Count (log scale)")
    for index, value in enumerate(overview_counts):
        ax.text(index, value, human_number(value), ha="center", va="bottom", fontsize=8)
    ax.text(-0.18, 1.05, "a", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    ax = axes[1]
    labels = [stat.label for stat in sorted_by_name]
    values = [stat.records for stat in sorted_by_name]
    y_positions = range(len(sorted_by_name))
    ax.barh(y_positions, values, color="#66CCEE", edgecolor="#334E5C", linewidth=0.45)
    ax.set_yticks(list(y_positions), labels=labels)
    ax.invert_yaxis()
    ax.set_xlabel("Chromosomes/FASTA records")
    ax.set_title("Reference composition", pad=10)
    ax.xaxis.set_major_formatter(FuncFormatter(human_number))
    max_value = max(values) if values else 0
    ax.set_xlim(0, max_value * 1.12 if max_value else 1)
    for y_position, value in zip(y_positions, values):
        ax.text(value + max(max_value * 0.01, 0.5), y_position, f"{value:,}", va="center", fontsize=7)
    ax.text(-0.14, 1.05, "b", transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")

    fig.suptitle("Genome reference and lncRNA transcript summary", fontsize=12, y=0.995)
    save_figure(fig, output_dir / "data_summary_ab")
    plt.close(fig)


def make_basepair_plot(output_dir: Path, genome_stats: list[FastaStats]) -> None:
    plt, FuncFormatter = require_matplotlib()

    sorted_by_bp = sorted(genome_stats, key=lambda stat: stat.total_bp, reverse=True)
    sorted_by_bp = [stat for stat in sorted_by_bp if stat.total_bp > 0]
    if not sorted_by_bp:
        print("Skipping base-pair plot because the summary has no genome total_bp values.")
        return
    figure_height = max(4.2, 0.27 * len(sorted_by_bp) + 1.8)
    fig, ax = plt.subplots(figsize=(8.6, figure_height))
    labels = [stat.label for stat in sorted_by_bp]
    values = [stat.total_bp for stat in sorted_by_bp]
    y_positions = range(len(sorted_by_bp))

    ax.barh(y_positions, values, color="#CC6677", edgecolor="#5E2B33", linewidth=0.45)
    ax.set_yticks(list(y_positions), labels=labels)
    ax.invert_yaxis()
    ax.set_xlabel("Total base pairs")
    ax.set_title("", pad=10)
    ax.xaxis.set_major_formatter(FuncFormatter(human_number))
    max_value = max(values) if values else 0
    ax.set_xlim(0, max_value * 1.14 if max_value else 1)
    for y_position, value in zip(y_positions, values):
        ax.text(
            value + max(max_value * 0.01, 1),
            y_position,
            human_number(value),
            va="center",
            fontsize=7,
        )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save_figure(fig, output_dir / "reference_base_pairs")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(args.output_dir / ".matplotlib"))
    Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    using_existing_summary = (
        args.input_summary is not None
        and args.input_summary.exists()
        and args.input_summary.stat().st_size > 0
    )
    if using_existing_summary:
        genome_stats, lncrna_stats = read_summary_tsv(args.input_summary)
        print(f"Reused counts from: {args.input_summary}")
    else:
        if args.input_summary is not None and args.input_summary.exists() and args.input_summary.stat().st_size == 0:
            print(f"Summary TSV is empty, counting FASTA inputs instead: {args.input_summary}")
        if args.input_genomes is None or args.input_lncrna is None:
            raise SystemExit(
                "Provide --input-summary with a non-empty TSV, or provide both "
                "--input-genomes and --input-lncrna to count from FASTA files."
            )
        genome_paths = discover_genome_fastas(args.input_genomes)
        lncrna_path = resolve_lncrna_path(args.input_lncrna)
        genome_stats = [count_fasta(path) for path in genome_paths]
        lncrna_stats = count_fasta(lncrna_path)

    write_summary_outputs(args.output_dir, genome_stats, lncrna_stats)
    make_summary_plot(args.output_dir, genome_stats, lncrna_stats)
    make_basepair_plot(args.output_dir, genome_stats)

    print(f"Genome references counted: {len(genome_stats)}")
    print(f"lncRNA transcripts counted: {lncrna_stats.records}")
    print(f"Outputs written to: {args.output_dir}")


if __name__ == "__main__":
    main()
