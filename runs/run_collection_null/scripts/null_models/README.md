# Human null-model dataset preparation

```bash
#Used
python3 scripts/null_models/prepare_null_models.py \
  --input-dir . \
  --output-dir generated_3000 \
  --cache-dir /tmp/fraglncrna_null_models_cache \
  --n-transcripts 3000 \
  --shuffles-per-transcript 10 \
  --intergenic-candidates 4096 \
  --seed 42 \
  --threads 15 \
  --k-values 15,20,32 \
  --allow-partial
```
Prepare a single human query collection for later FragIncRNA comparisons across
primate assemblies. This command only prepares inputs. It never downloads data,
changes the originals, builds an IBF, or starts genome searches.

Python 3.11+ and NumPy 1.24+ are required. Tests require pytest 7+.
Use an existing environment or explicitly install the small requirements file:

```bash
python3 -m pip install -r scripts/null_models/requirements.txt
python3 -m pytest -q -o cache_dir=/tmp/null_models_pytest_cache tests/test_null_models.py
```

After preparation, run the separate standard-library auditor if desired:

```bash
python3 tests/check_null_models_output.py \
  --dataset generated_smoke_25 \
  --annotation-gtf gencode.v49.annotation.gtf.gz --flank 10000
```

It imports none of the preparer and checks sliding-window counts, shuffle
invariants and parental replicate counts, transcript-BED mask reconstruction,
all-gene/flank exclusion, intergenic nonoverlap, and hard length/GC constraints.

The supplied files in this workspace are directly in `null_models/`, so use
`--input-dir .` here. In a checkout with an `input/` subdirectory, point there.
The supporting Python files stay beside the entry point. The same
`scripts/null_models/` and `tests/test_null_models.py` layout can be copied into
the FragIncRNA reporting checkout without modifying its C++ code or grid runner.

## Inputs and validation

| Input override | Default basename in input directory | Purpose |
| --- | --- | --- |
| `--lncrna-fasta` | `gencode.v49.lncRNA_transcripts.fa` | Original complete lncRNA transcripts |
| `--lncrna-gtf` | `gencode.v49.long_noncoding_RNAs.gtf` | lncRNA exon structures |
| `--annotation-gtf` | `gencode.v49.annotation.gtf` | Comprehensive all-gene exclusion and coding exons |
| `--coding-fasta` | `gencode.v49.pc_transcripts.fa` | Complete protein-coding transcripts, including UTRs |
| `--genome-fasta` | `GRCh38.primary_assembly.genome.fa` | Human source reference only |
| `--repeats-bed` | `GRCh38.repeats.bed` | Human RepeatMasker coordinates |

The optional seventh file, `rmsk.txt`, is unnecessary because the converted BED is
available. It is never read. No raw-UCSC conversion is performed implicitly.
Each default basename is searched exactly with no suffix, `.gz`, `.bgz`, or
`.bgzf`. Multiple matches fail with candidates and require an explicit override.
Compression is detected by gzip magic bytes, including hidden/misleading suffixes.
Normal preparation fully reads inputs and checks gzip CRC/truncation.

Versioned transcript and gene IDs are authoritative; symbols are display fields.
Duplicate FASTA IDs, malformed DNA, header-length discrepancies, missing eligible
lncRNA annotation IDs, incompatible annotated releases/assemblies, FASTA/GTF gene-ID
differences, and conflicting comprehensive/lncRNA exon structures fail clearly.
GTF quoted attributes and GFF3 `key=value` attributes are supported. GFF3 must
provide authoritative transcript/gene IDs and biotypes; multi-parent exons must
be expanded before use. Missing assembly/version comments cannot independently
certify a release; ID agreement and complete exon-spliced validation still apply.

Gene and transcript biotypes are distinct. The default lncRNA population requires
`gene_type=lncRNA`, with accepted transcript types
`lncRNA,retained_intron,processed_transcript` (set explicitly with
`--lncrna-transcript-types`). GENCODE v49's lncRNA export also contains TEC entries;
those are excluded and listed by ID. Noncanonical accepted lncRNA transcript
types carry quality flags. Coding positives require both gene and transcript
type `protein_coding`; retained-intron, nonsense-mediated-decay and other coding
gene isoforms are excluded. Biotype counts and every excluded source transcript
are reported before selection, so these filters cannot shrink an already
selected population silently.

The coding FASTA has ALL-region scope, while `gencode.v49.annotation.gtf` has
chromosome-only scope; these are different official exports of the same release.
See the [GENCODE v49 file descriptions](https://www.gencodegenes.org/human/release_49.html).
By default `--missing-coding-annotation exclude` lists every coding FASTA ID
without a supplied coordinate model and excludes it before matching. No ALT
coordinate or chromosome identity is guessed for unresolved records. Use
`--missing-coding-annotation error` for strict complete-source coverage. This
policy does not admit any unannotated coding control: every eligible positive
still requires matching versioned IDs, coding biotypes, and exact splicing.

GTF coordinates are converted from 1-based inclusive to 0-based half-open.
Only primary autosomes 1–22 are eligible by default. `chr1`/`1` aliases are
explicitly mapped to reference names and recorded; mixed ambiguous aliases fail.
Mitochondrial, ALT, and unlocalized contigs are excluded and counted. To opt into
sex chromosomes, explicitly pass, for example,
`--chromosomes 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y`.
The user must establish sample/assembly compatibility before doing this.

Every eligible lncRNA and coding sequence is compared with the reference spliced
through its annotated exons, including reverse complementation on minus strands.
Exon lengths must match exactly. The default substitution tolerance is zero;
`--max-splice-mismatch FRACTION` permits an explicit recorded substitution fraction
but never repairs coordinates or exon lengths. Sequences without a valid 32-mer
are excluded before sampling and counted. Short sequences below 200 bases are
flagged. Eligibility always requires a valid 32-mer, independent of `--k-values`.

## Commands

From this workspace, inspect paths and input headers without writing anything:

```bash
python3 scripts/null_models/prepare_null_models.py \
  --input-dir . --output-dir generated \
  --n-transcripts 5000 --shuffles-per-transcript 10 \
  --seed 42 --threads 15 --k-values 15,20,32 --dry-run
```

The dry run explicitly reports its limited header-probe scope and approximate
resources; it does not claim to have validated all records or matching feasibility.

Small real-input smoke run:

```bash
python3 scripts/null_models/prepare_null_models.py \
  --input-dir . --output-dir generated_smoke_25 \
  --cache-dir /tmp/fraglncrna_null_models_cache \
  --n-transcripts 25 --shuffles-per-transcript 2 \
  --seed 42 --threads 4 --k-values 15,20,32
```

Full preparation command, **for the user to launch after reviewing smoke QC**:

```bash
python3 scripts/null_models/prepare_null_models.py \
  --input-dir /mnt/d/primates_bmc/review/null_models \
  --output-dir /mnt/d/primates_bmc/review/null_models/generated \
  --cache-dir /tmp/fraglncrna_null_models_cache \
  --n-transcripts 5000 --shuffles-per-transcript 10 \
  --seed 42 --threads 15 --k-values 15,20,32
```

This is preparation only. Full-scale matching is not guaranteed by a 25-parent
smoke run. If matching fails, review the sibling `*.matching_failure.json` file.
There is no automatic relaxation or replacement of hard-to-match selected parents.
Explicitly increase `--intergenic-candidates` (default 512), change a documented
tolerance, or use `--allow-partial` in a new run. Partial mode retains the fixed
baseline population, records every unmatched parent, and names each control
FASTA with its actual record count, including zero.

`--resume` verifies input/code/config provenance and all output SHA-256 hashes for
a completed run. Worker count and output location are not scientific fingerprint
fields. A compatible interrupted staging directory is restarted deterministically;
partially written FASTAs are never appended to. `--force` replaces only a completed
directory owned by this tool and keeps the previous output until its replacement
is complete. Neither flag overwrites arbitrary directories. Normal errors remove
unpublished staging data; explicit caches and matching-failure diagnostics remain.
A hard kill can leave `.OUTPUT.preparing` and lock files. Confirm the recorded
PID is dead before removing stale locks, then use `--resume` with unchanged inputs
and settings. Source changes require a new preparation cache/output.

## Frozen selection and matching

1. Each lncRNA gene contributes at most one seeded representative by default.
   `--max-per-gene N` changes that cap; zero disables it. This defines the
   gene-limited sampling population before stratification.
2. Proportional stratified sampling uses Hamilton allocation and deterministic
   seeded tie-breaking. Fixed cut points are length
   `[200,500,1000,2000,4000,8000,16000]`, GC fraction `[.35,.45,.55,.65]`, repeat
   fraction `[.01,.25,.5,.75]`, and exon count `[1,2,4,8]`. Cut points belong to
   the lower bin. `qc/strata.tsv` records source, gene-limited and sample counts.
3. Coding assignment processes parents with the fewest feasible candidates first.
   It minimizes the sum of squared standardized differences in log length, GC,
   repeat fraction and log(1+exon count), plus 0.1 for a chromosome mismatch.
   Standard deviations come from the eligible coding pool. Exact ties use sorted
   versioned transcript IDs. Default hard limits are relative length error ≤0.15
   and absolute GC fraction error ≤0.05 (five percentage points).
   `--repeat-tolerance` adds a hard repeat-fraction limit. A coding transcript is
   used once and a coding gene once, unless `--coding-isoforms` is explicitly set.
   Assignment is greedy, not a global minimum-cost solution.
4. Intergenic availability subtracts **all comprehensive annotated gene bodies**
   plus a 10 kb flank (`--flank`), and every non-ACGT run. Sampling is uniform over
   available starts of exact parent-length windows across the entire eligible
   human universe, not just neighborhoods of lncRNA genes. There are at most 512
   draws per parent. Feasible candidates meet the same GC and optional repeat
   limits. Their cost is `(GC_delta/.1)^2 + (repeat_delta/.25)^2`, plus 0.1 for a
   chromosome mismatch. Cost/chromosome/start break ties deterministically.
   Longer parents go first; each selected interval is removed to enforce
   nonoverlap. Available starts, draws and acceptable candidates are recorded per
   parent. Bounded failure is not proof that no feasible locus exists.

Repeat fractions use **union coverage**, including overlapping RepeatMasker
annotations. Within-class union fractions are reported separately and can sum
above the total because classes can overlap. The supplied BED layout is
`chrom,start,end,name,class,family`; BED3 is accepted as unclassified. Standard
numeric BED5 scores are not treated as repeat classes.

Repeat fraction, chromosome and exon count are soft matching features by default.
No repeat-class or exon-structure equality is guaranteed. Intergenic controls are
contiguous loci; no artificial block-matched variant is implemented.

## Shuffling and masking

Shuffling constructs an Euler trail through the dinucleotide multigraph, using
Wilson's loop-erased random-walk last-exit tree and randomized outgoing edge
orders. This is the Euler/tree approach described in the
[uShuffle paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/), implemented
locally for A/C/G/T runs without adding a shuffle-library dependency. Each run's
length, mononucleotides and all 16 dinucleotide counts are checked exactly on
every attempt. Ambiguous bases remain at their exact original positions, run
endpoints are retained, and lowercase status is retained at each position.

SHA-256-derived seed namespaces include the global seed, complete parent query
ID and replicate. Worker scheduling and Python's randomized `hash()` are never
used. Up to 32 attempts per replicate reject unchanged/duplicate sequences.
Unavoidable duplicates are kept with explicit flags after that cap. Rejection
conditions the randomization distribution; this is not an assertion of uniform
independent draws over distinct sequences. Ten shuffles remain ten paired
replicates of one transcript, not ten independent transcripts.

Genome repeats are projected through transcript exon order, converting each
minus-strand overlap to offsets from the corresponding exon end. Selected
sequences are N-masked only after splicing has been validated. Baseline sequence
characters/case and source headers are preserved (FASTA wrapping is normalized
to 80 columns). The masked query's length and all unmasked bases remain unchanged.

For each requested k, metadata records `max(length-k+1,0)`, the actual ACGT-only
window denominator, and windows lost due to query masking. One masked base can
invalidate up to k windows. `--k-values` changes only QC, not the FASTAs.

## Artifacts, provenance and resource use

The default full output has the five requested files under `fasta/`:
`real_lncrna_5000.fa.gz`, `dinuc_shuffled_50000.fa.gz`,
`intergenic_matched_5000.fa.gz`, `coding_matched_5000.fa.gz`, and
`real_lncrna_5000_repeat_masked.fa.gz`.

`metadata/` contains all requested selected/shuffled/paired/mask tables,
matching diagnostics, exclusion counts, valid-window counts, and
`transcript_repeat_intervals.bed.gz`. This BED uses **masked FASTA query IDs and
transcript coordinates**, with union intervals. Repeat-class coverage is in
`transcript_repeat_classes.tsv.gz`. `query_metadata.tsv.gz` is a unified table
covering every FASTA record, including sequence SHA-256. It records fraction
units (multiply by 100 for percentages), versioned source IDs, parental links,
replicate, exon/interval coordinates, signed deltas, costs and seed namespaces.
Shuffled sequences have no genomic repeat annotation: their repeat fraction and
source coordinates are blank; use the parent table for the source's repeat burden.
`ineligible_transcripts.tsv.gz` lists excluded source IDs, reasons and original/
valid windows for k=15,20,32 (plus any additional requested k).

`qc/` contains `summary.json`, `report.md`, `strata.tsv`, and the per-pair
`length_gc_repeat_matching.tsv`. Plots are optional and are not generated.
`provenance/` contains resolved TOML, input/output SHA-256 tables, code hashes,
scientific configuration, actual counts, Python version and worker count.
Input hashes cover original file bytes, including compression. Output receipt
hashes cover every final artifact except the receipt itself.

All gzip output has zero timestamps and no embedded output filename. FASTAs and
metadata are byte-identical for the same inputs/seed/config, independent of
threads, process restart and output directory. Provenance records runtime paths
and worker count and therefore is not expected to be byte-identical across them.
Every output is reparsed and its gzip CRC, record IDs, sequence hashes and window
counts are verified before atomic publication.

The human reference streams into a disk-backed uppercase byte file with indexed
chromosome offsets. It is never loaded as a whole-genome RAM string. Transcript
sequences live in SQLite; repeat coordinates use compact arrays and prefix maxima
for overlap lookup. Genome and transcript caches have SHA-256 integrity checks,
input/code/config fingerprints and exclusive locks. Allow roughly 4–5 GB of local
cache disk and several GB of RAM for comprehensive annotations and features.
Prefer `--cache-dir /tmp/...` or a native Linux disk on WSL. Without that option,
the temporary cache is removed after preparation. No list of all possible genomic
candidate windows is created. `--threads` parallelizes shuffling using processes;
parsing, stratification and assignment are ordered and serial.
