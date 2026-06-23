#!/usr/bin/env python3
"""
Summarize synteny support for high-confidence aligned lncRNA transcripts.

The report uses human protein-coding neighbors as anchors, maps those anchors
through BioMart one-to-one orthology, then checks whether the expected target
gene symbols are near the BAM-mapped lncRNA locus in each target assembly.
"""


import argparse
import bisect
import csv
import gzip
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam


DEFAULT_HIGH_CONF_DIR = Path("alignment_vs_kmers/high_conf/tables")
DEFAULT_BAM_DIR = Path("../output_minimap2")
DEFAULT_HUMAN_LNCRNA_GTF = Path("../lncrna/gencode.v49.long_noncoding_RNAs.gtf")
DEFAULT_HUMAN_FEATURES = Path("../homology/features.txt")
DEFAULT_HOMOLOGUES = Path("../homology/Homologues.txt")
DEFAULT_TARGET_ANNOTATIONS_DIR = Path("../homology/assemblies_gtf")
DEFAULT_OUTDIR = Path("alignment_vs_kmers/synteny_support")

TARGET_ANNOTATIONS = {
    "Clintptrv2": (
        "chimpanzee",
        "GCF_002880755.1_Clint_PTRv2_genomic.gtf.gz",
    ),
    "MhudibluPPAv0": (
        "bonobo",
        "GCF_013052645.1_Mhudiblu_PPA_v0_genomic.gtf.gz",
    ),
    "Mmul10": (
        "macaque",
        "GCF_003339765.1_Mmul_10_genomic.gtf.gz",
    ),
    "T2TMMU8v20": (
        "macaque",
        "GCF_049350105.2_T2T-MMU8v2.0_genomic.gtf.gz",
    ),
    "mPanPan1v20pri": (
        "bonobo",
        "GCF_029289425.2_NHGRI_mPanPan1-v2.1_pri_genomic.gff.gz",
    ),
    "mPanTro3v20pri": (
        "chimpanzee",
        "GCF_028858775.2_NHGRI_mPanTro3-v2.1_pri_genomic.gtf.gz",
    ),
    "panpan11": (
        "bonobo",
        "GCF_000258655.2_panpan1.1_genomic.gtf.gz",
    ),
}

SPECIES_COLUMNS = {
    "chimpanzee": (
        "Chimpanzee gene stable ID",
        "Chimpanzee gene name",
        "Chimpanzee homology type",
    ),
    "bonobo": (
        "Bonobo gene stable ID",
        "Bonobo gene name",
        "Bonobo homology type",
    ),
    "macaque": (
        "Macaque gene stable ID",
        "Macaque gene name",
        "Macaque homology type",
    ),
}


@dataclass(frozen=True)
class Gene:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str


@dataclass(frozen=True)
class LncRna:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str


@dataclass(frozen=True)
class Neighbor:
    gene: Gene
    relation: str
    distance: int


@dataclass(frozen=True)
class AlignmentHit:
    chrom: str
    start: int
    end: int
    strand: str
    mapq: int
    cigar: str
    aligned_bp: int
    alignment_score: int


class GeneIndex:
    def __init__(self, genes: list[Gene]) -> None:
        self.by_chrom_start: dict[str, list[Gene]] = defaultdict(list)
        self.by_chrom_end: dict[str, list[Gene]] = defaultdict(list)
        for gene in genes:
            self.by_chrom_start[gene.chrom].append(gene)
            self.by_chrom_end[gene.chrom].append(gene)

        self.starts: dict[str, list[int]] = {}
        self.ends: dict[str, list[int]] = {}
        for chrom, chrom_genes in self.by_chrom_start.items():
            chrom_genes.sort(key=lambda gene: (gene.start, gene.end, gene.gene_name))
            self.starts[chrom] = [gene.start for gene in chrom_genes]
        for chrom, chrom_genes in self.by_chrom_end.items():
            chrom_genes.sort(key=lambda gene: (gene.end, gene.start, gene.gene_name))
            self.ends[chrom] = [gene.end for gene in chrom_genes]

    def nearest(self, chrom: str, start: int, end: int, count: int) -> list[Neighbor]:
        if chrom not in self.by_chrom_start:
            return []

        starts = self.starts[chrom]
        ends = self.ends[chrom]
        by_start = self.by_chrom_start[chrom]
        by_end = self.by_chrom_end[chrom]
        neighbors: list[Neighbor] = []
        seen: set[str] = set()

        overlap_limit = bisect.bisect_right(starts, end)
        overlaps = []
        for index in range(overlap_limit - 1, max(-1, overlap_limit - 501), -1):
            gene = by_start[index]
            if gene.end >= start:
                overlaps.append(gene)
            if len(overlaps) >= count:
                break
        overlaps.sort(
            key=lambda gene: (
                0 if gene.start <= start and gene.end >= end else 1,
                min(abs(gene.start - start), abs(gene.end - end)),
                gene.gene_name,
            )
        )
        for gene in overlaps[:count]:
            neighbors.append(Neighbor(gene, "overlap", 0))
            seen.add(gene.gene_id)

        upstream_limit = bisect.bisect_left(ends, start)
        for gene in reversed(by_end[:upstream_limit]):
            if gene.gene_id in seen:
                continue
            neighbors.append(Neighbor(gene, "upstream", start - gene.end))
            seen.add(gene.gene_id)
            if sum(item.relation == "upstream" for item in neighbors) >= count:
                break

        downstream_start = bisect.bisect_right(starts, end)
        for gene in by_start[downstream_start:]:
            if gene.gene_id in seen:
                continue
            neighbors.append(Neighbor(gene, "downstream", gene.start - end))
            seen.add(gene.gene_id)
            if sum(item.relation == "downstream" for item in neighbors) >= count:
                break

        return sorted(
            neighbors,
            key=lambda item: (
                {"upstream": 0, "overlap": 1, "downstream": 2}[item.relation],
                item.distance,
                item.gene.gene_name,
            ),
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a synteny-support table for high-confidence lncRNA "
            "alignments using nearby protein-coding genes."
        )
    )
    parser.add_argument("--high-conf-dir", type=Path, default=DEFAULT_HIGH_CONF_DIR)
    parser.add_argument("--bam-dir", type=Path, default=DEFAULT_BAM_DIR)
    parser.add_argument("--human-lncrna-gtf", type=Path, default=DEFAULT_HUMAN_LNCRNA_GTF)
    parser.add_argument("--human-features", type=Path, default=DEFAULT_HUMAN_FEATURES)
    parser.add_argument("--homologues", type=Path, default=DEFAULT_HOMOLOGUES)
    parser.add_argument(
        "--target-annotations-dir",
        type=Path,
        default=DEFAULT_TARGET_ANNOTATIONS_DIR,
    )
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--neighbor-count",
        type=int,
        default=2,
        help="Number of upstream and downstream protein-coding genes to inspect.",
    )
    parser.add_argument(
        "--genes",
        help="Optional comma-separated lncRNA gene names for a focused report.",
    )
    parser.add_argument(
        "--assemblies",
        help="Optional comma-separated assembly aliases to process.",
    )
    return parser.parse_args()


def smart_open(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_gtf_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for match in re.finditer(r'([A-Za-z0-9_]+)\s+"([^"]*)"', text):
        attrs[match.group(1)] = match.group(2)
    return attrs


def parse_gff_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def parse_attributes(text: str) -> dict[str, str]:
    if "=" in text and not re.search(r'\w+\s+"', text):
        return parse_gff_attributes(text)
    return parse_gtf_attributes(text)


def strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def normalize_name(name: str) -> str:
    return name.strip().upper()


def normalize_human_chrom(chrom: str) -> str:
    if chrom.startswith("chr"):
        return chrom
    if chrom == "MT":
        return "chrM"
    if chrom in {"X", "Y"} or chrom.isdigit():
        return f"chr{chrom}"
    return chrom


def load_human_protein_genes(path: Path) -> list[Gene]:
    genes: dict[str, Gene] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = row["Gene stable ID"]
            if gene_id in genes:
                continue
            genes[gene_id] = Gene(
                chrom=normalize_human_chrom(row["Chromosome/scaffold name"]),
                start=int(row["Gene start (bp)"]),
                end=int(row["Gene end (bp)"]),
                strand="+" if row["Strand"] == "1" else "-",
                gene_id=gene_id,
                gene_name=row["Gene name"],
            )
    return list(genes.values())


def load_homologues(path: Path) -> dict[str, dict[str, tuple[str, str]]]:
    homologues: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            human_gene_id = row["Gene stable ID"]
            for species, (id_col, name_col, type_col) in SPECIES_COLUMNS.items():
                if row[type_col] != "ortholog_one2one":
                    continue
                target_gene_id = row[id_col].strip()
                target_gene_name = row[name_col].strip()
                if target_gene_id and target_gene_name:
                    homologues[human_gene_id][species] = (
                        target_gene_id,
                        target_gene_name,
                    )
    return homologues


def load_human_lncRNAs(path: Path) -> dict[str, LncRna]:
    transcripts: dict[str, LncRna] = {}
    with smart_open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue
            transcripts[transcript_id] = LncRna(
                chrom=fields[0],
                start=int(fields[3]),
                end=int(fields[4]),
                strand=fields[6],
                gene_id=attrs.get("gene_id", ""),
                gene_name=attrs.get("gene_name", ""),
            )
    return transcripts


def target_gene_name(attrs: dict[str, str]) -> str:
    return (
        attrs.get("gene")
        or attrs.get("gene_name")
        or attrs.get("Name")
        or attrs.get("gene_id")
        or attrs.get("ID")
        or ""
    )


def target_gene_id(attrs: dict[str, str], gene_name: str) -> str:
    return attrs.get("gene_id") or attrs.get("ID") or gene_name


def load_target_genes(path: Path) -> list[Gene]:
    genes: list[Gene] = []
    with smart_open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = parse_attributes(fields[8])
            biotype = (
                attrs.get("gene_biotype")
                or attrs.get("gene_type")
                or attrs.get("gbkey")
                or ""
            )
            if biotype != "protein_coding":
                continue
            gene_name = target_gene_name(attrs)
            if not gene_name:
                continue
            genes.append(
                Gene(
                    chrom=fields[0],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    gene_id=target_gene_id(attrs, gene_name),
                    gene_name=gene_name,
                )
            )
    return genes


def load_high_conf(path: Path, genes: set[str] | None) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if genes is not None and row["gene_name"] not in genes:
                continue
            rows[row["transcript_id"]] = row
    return rows


def alignment_score(read: pysam.AlignedSegment) -> int:
    try:
        return int(read.get_tag("AS"))
    except KeyError:
        return 0


def aligned_reference_bp(read: pysam.AlignedSegment) -> int:
    if read.reference_start is None or read.reference_end is None:
        return 0
    return int(read.reference_end - read.reference_start)


def better_hit(new: AlignmentHit, old: AlignmentHit | None) -> bool:
    if old is None:
        return True
    return (
        new.mapq,
        new.alignment_score,
        new.aligned_bp,
    ) > (
        old.mapq,
        old.alignment_score,
        old.aligned_bp,
    )


def extract_best_hits(bam_path: Path, transcript_ids: set[str]) -> dict[str, AlignmentHit]:
    hits: dict[str, AlignmentHit] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            transcript_id = read.query_name.split("|", 1)[0]
            if transcript_id not in transcript_ids:
                continue
            if read.reference_start is None or read.reference_end is None:
                continue
            hit = AlignmentHit(
                chrom=bam.get_reference_name(read.reference_id),
                start=int(read.reference_start) + 1,
                end=int(read.reference_end),
                strand="-" if read.is_reverse else "+",
                mapq=int(read.mapping_quality),
                cigar=read.cigarstring or "",
                aligned_bp=aligned_reference_bp(read),
                alignment_score=alignment_score(read),
            )
            if better_hit(hit, hits.get(transcript_id)):
                hits[transcript_id] = hit
    return hits


def format_neighbors(neighbors: list[Neighbor]) -> str:
    return ";".join(
        f"{item.gene.gene_name}|{item.gene.gene_id}|{item.relation}|{item.distance}"
        for item in neighbors
    )


def expected_orthologs(
    neighbors: list[Neighbor],
    species: str,
    homologues: dict[str, dict[str, tuple[str, str]]],
) -> list[tuple[str, str, str, str]]:
    expected = []
    seen: set[str] = set()
    for item in neighbors:
        ortholog = homologues.get(item.gene.gene_id, {}).get(species)
        if ortholog is None:
            continue
        target_gene_id, target_gene_name = ortholog
        key = normalize_name(target_gene_name)
        if key in seen:
            continue
        seen.add(key)
        expected.append(
            (
                item.gene.gene_name,
                item.gene.gene_id,
                target_gene_name,
                target_gene_id,
            )
        )
    return expected


def format_expected(expected: list[tuple[str, str, str, str]]) -> str:
    return ";".join(
        f"{human_name}|{human_id}|{target_name}|{target_id}"
        for human_name, human_id, target_name, target_id in expected
    )


def matched_expected(
    expected: list[tuple[str, str, str, str]],
    target_neighbors: list[Neighbor],
) -> list[str]:
    observed = {normalize_name(item.gene.gene_name) for item in target_neighbors}
    return [
        target_name
        for _human_name, _human_id, target_name, _target_id in expected
        if normalize_name(target_name) in observed
    ]


def assembly_from_table(path: Path) -> str:
    return path.stem


def selected_assemblies(args: argparse.Namespace) -> list[str]:
    assemblies = list(TARGET_ANNOTATIONS)
    if args.assemblies:
        requested = [item.strip() for item in args.assemblies.split(",") if item.strip()]
        unknown = sorted(set(requested).difference(TARGET_ANNOTATIONS))
        if unknown:
            raise SystemExit(f"Unknown or unsupported assemblie(s): {', '.join(unknown)}")
        assemblies = requested
    return assemblies


def write_report(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    genes = None
    if args.genes:
        genes = {item.strip() for item in args.genes.split(",") if item.strip()}

    print("Loading human protein-coding genes...")
    human_gene_index = GeneIndex(load_human_protein_genes(args.human_features))
    print("Loading human-to-primate orthologues...")
    homologues = load_homologues(args.homologues)
    print("Loading human lncRNA transcript coordinates...")
    lncRNAs = load_human_lncRNAs(args.human_lncrna_gtf)

    report_path = args.outdir / "high_conf_lncRNA_synteny_support.tsv"
    summary_path = args.outdir / "high_conf_lncRNA_synteny_summary_by_assembly.tsv"

    fieldnames = [
        "assembly",
        "species_group",
        "gene_name",
        "gene_id",
        "transcript_id",
        "human_locus",
        "target_locus",
        "target_mapq",
        "target_strand",
        "human_pc_neighbors",
        "expected_target_one2one_orthologs",
        "target_pc_neighbors",
        "matched_expected_target_genes",
        "expected_neighbor_count",
        "matched_neighbor_count",
        "synteny_supported_any",
        "synteny_supported_two_or_more",
        "status",
    ]

    summaries: list[dict[str, object]] = []
    with report_path.open("w", newline="") as report_handle:
        writer = csv.DictWriter(report_handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()

        for assembly in selected_assemblies(args):
            species, annotation_name = TARGET_ANNOTATIONS[assembly]
            table_path = args.high_conf_dir / f"{assembly}.tsv"
            bam_path = args.bam_dir / f"human_lncRNA_vs_{assembly}.sorted.bam"
            annotation_path = args.target_annotations_dir / annotation_name

            if not table_path.is_file():
                print(f"Skipping {assembly}: missing high-confidence table {table_path}")
                continue
            if not bam_path.is_file():
                print(f"Skipping {assembly}: missing BAM {bam_path}")
                continue
            if not annotation_path.is_file():
                print(f"Skipping {assembly}: missing annotation {annotation_path}")
                continue

            print(f"Processing {assembly}...")
            high_conf = load_high_conf(table_path, genes)
            target_gene_index = GeneIndex(load_target_genes(annotation_path))
            hits = extract_best_hits(bam_path, set(high_conf))

            stats = {
                "assembly": assembly,
                "species_group": species,
                "high_conf_transcripts": len(high_conf),
                "primary_alignments_found": len(hits),
                "rows_with_human_lncRNA_coordinates": 0,
                "rows_with_target_pc_neighbors": 0,
                "rows_with_expected_one2one_neighbors": 0,
                "rows_with_any_synteny_match": 0,
                "rows_with_two_or_more_synteny_matches": 0,
            }

            for transcript_id, high_conf_row in sorted(
                high_conf.items(),
                key=lambda item: (item[1]["gene_name"], item[0]),
            ):
                lncRNA = lncRNAs.get(transcript_id)
                hit = hits.get(transcript_id)
                status = "ok"

                human_neighbors: list[Neighbor] = []
                target_neighbors: list[Neighbor] = []
                expected: list[tuple[str, str, str, str]] = []
                matched: list[str] = []
                human_locus = ""
                target_locus = ""
                target_mapq = ""
                target_strand = ""

                if lncRNA is None:
                    status = "missing_human_lncRNA_coordinates"
                else:
                    stats["rows_with_human_lncRNA_coordinates"] += 1
                    human_locus = (
                        f"{lncRNA.chrom}:{lncRNA.start}-{lncRNA.end}:{lncRNA.strand}"
                    )
                    human_neighbors = human_gene_index.nearest(
                        lncRNA.chrom,
                        lncRNA.start,
                        lncRNA.end,
                        args.neighbor_count,
                    )
                    expected = expected_orthologs(human_neighbors, species, homologues)

                if hit is None:
                    status = (
                        "missing_primary_alignment"
                        if status == "ok"
                        else f"{status};missing_primary_alignment"
                    )
                else:
                    target_locus = f"{hit.chrom}:{hit.start}-{hit.end}:{hit.strand}"
                    target_mapq = str(hit.mapq)
                    target_strand = hit.strand
                    target_neighbors = target_gene_index.nearest(
                        hit.chrom,
                        hit.start,
                        hit.end,
                        args.neighbor_count,
                    )
                    if target_neighbors:
                        stats["rows_with_target_pc_neighbors"] += 1

                if expected:
                    stats["rows_with_expected_one2one_neighbors"] += 1
                if expected and target_neighbors:
                    matched = matched_expected(expected, target_neighbors)
                    if matched:
                        stats["rows_with_any_synteny_match"] += 1
                    if len(matched) >= 2:
                        stats["rows_with_two_or_more_synteny_matches"] += 1

                writer.writerow(
                    {
                        "assembly": assembly,
                        "species_group": species,
                        "gene_name": high_conf_row["gene_name"],
                        "gene_id": high_conf_row["gene_id"],
                        "transcript_id": transcript_id,
                        "human_locus": human_locus,
                        "target_locus": target_locus,
                        "target_mapq": target_mapq,
                        "target_strand": target_strand,
                        "human_pc_neighbors": format_neighbors(human_neighbors),
                        "expected_target_one2one_orthologs": format_expected(expected),
                        "target_pc_neighbors": format_neighbors(target_neighbors),
                        "matched_expected_target_genes": ";".join(matched),
                        "expected_neighbor_count": len(expected),
                        "matched_neighbor_count": len(matched),
                        "synteny_supported_any": bool(matched),
                        "synteny_supported_two_or_more": len(matched) >= 2,
                        "status": status,
                    }
                )

            summaries.append(stats)

    with summary_path.open("w", newline="") as summary_handle:
        writer = csv.DictWriter(
            summary_handle,
            delimiter="\t",
            fieldnames=[
                "assembly",
                "species_group",
                "high_conf_transcripts",
                "primary_alignments_found",
                "rows_with_human_lncRNA_coordinates",
                "rows_with_target_pc_neighbors",
                "rows_with_expected_one2one_neighbors",
                "rows_with_any_synteny_match",
                "rows_with_two_or_more_synteny_matches",
            ],
        )
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Wrote {report_path}")
    print(f"Wrote {summary_path}")


def main() -> None:
    write_report(parse_args())


if __name__ == "__main__":
    main()
