# FragIncRNA

The selectable production workflow is now available as schema 3. Use
`--mode ibf --scoring-method per_bin_verified` for IBF candidates followed by
exact verification within indexed candidate intervals. The main recovery and
status then use candidate-local verified support. It is a lower bound on
whole-reference exact recovery; it is not the same biological quantity.
See the [production audit](docs/PER_BIN_PRODUCTION_AUDIT.md),
[current method](docs/PER_BIN_IBF_METHOD.md),
[measured benchmark](docs/PER_BIN_PRODUCTION_BENCHMARK.md), and
[acceptance status](docs/PER_BIN_PRODUCTION_STATUS.md).

Ready-to-run fixture examples, after the measured build below:

```bash
python3 scripts/measure_command.py --output /tmp/frag-configure.json -- cmake -S . -B build-production -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLNCRNA_MERS_ENABLE_HIBF=OFF
python3 scripts/measure_command.py --output /tmp/frag-build.json -- cmake --build build-production -j 4
python3 scripts/measure_command.py --output /tmp/frag-tests.json -- ctest --test-dir build-production --output-on-failure
python3 scripts/measure_command.py --output /tmp/frag-production.json -- build-production/lncrna_mers_v2 --config configs/per_bin_production.toml
python3 scripts/measure_command.py --output /tmp/frag-validate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_production.toml --mode validate --output /tmp/frag-validate.tsv
python3 scripts/measure_command.py --output /tmp/frag-calibrate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_production.toml --mode calibrate --output /tmp/frag-calibrate.tsv
```

The fixture uses k=4/4-base bins and a one-hit threshold solely for executable
checks. For research, set the supplied v49 query/reference paths, k=20,
fragment_size=80000 and per_bin_min_hits=2; the 5%/top-64 rule is heuristic.
Run one reference per invocation. The [bounded runner](scripts/run_per_bin_production_benchmark.py)
processes 100 then 5,000 queries and never launches the full manuscript grid.
The data-preparation script stages and validates an extraction before exposing it.

After preparing the supplied v49 subsets and complete first contig, these use
the same real-data settings as the benchmark (outputs are new example files):

```bash
python3 scripts/measure_command.py --output /tmp/v49-production.json -- build-production/lncrna_mers_v2 --config configs/per_bin_v49_bounded.toml --load-index results/per_bin_production_validation/complete_chr1/k20_nominal.ibf --reference-cache results/per_bin_production_validation/complete_chr1/reference.bases --output /tmp/v49-production.tsv
python3 scripts/measure_command.py --output /tmp/v49-validate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_v49_bounded.toml --mode validate --load-index results/per_bin_production_validation/complete_chr1/k20_nominal.ibf --reference-cache results/per_bin_production_validation/complete_chr1/reference.bases --output /tmp/v49-validate.tsv
python3 scripts/measure_command.py --output /tmp/v49-calibrate.json -- build-production/lncrna_mers_v2 --config configs/per_bin_v49_bounded.toml --mode calibrate --load-index results/per_bin_production_validation/complete_chr1/k20_nominal.ibf --reference-cache results/per_bin_production_validation/complete_chr1/reference.bases --output /tmp/v49-calibrate.tsv
```

Omit `--load-index` on first use to build instead. The indexed reference here
is one complete macaque contig, not all 22 contigs. Validation/calibration scan
the full supplied reference independently; production does not. **Do not use
the default candidate filter as a drop-in global-recovery replacement:** on
the 5,000-query sample it missed 36.34% of exact-positive windows and 3,005 of
4,482 queries with any exact support. Biological accuracy is not established.

Choose `global_ibf` for historical any-bin recovery, `per_bin_ibf` for
approximate support in selected bins, `per_bin_verified` for exact membership
in selected intervals, or `global_exact` for independent full-reference
membership. Validation reports all four, even if exact remains primary.
Old unspecified exact/global invocations retain schema 2; per-bin IBF
production and explicit scoring use schema 3 with method/scope/provenance.
The merger checks source, scoring, schema, thresholds and candidate settings.

`--ibf-auto-tune off|recommend|apply` diagnoses settings and writes
`<output>.ibf_tuning.json`. Explicit off/recommend retains settings; apply
can budget a global-IBF objective using `--genome-fpr`. Per-bin methods retain
their design settings and recommend measured candidate-recall validation.
Unspecified legacy global-IBF invocations preserve their previous automatic
whole-reference budget; this compatibility rule is explicitly logged.
`--ibf-max-index-memory-gb` preflights raw index bytes; it is not a total
process RSS cap. Dry-run uses `--ibf-dry-run`. Loaded archives are checked and
never resized/rebuilt silently. Fragment alternatives are estimates, never
automatic scientific-setting changes. Per-bin Bloom FPR, approximate union
FPR, transcript classification errors and biological FDR are distinct.

`--threads 15` uses OpenMP for bounded query and interval work; actual workers
are capped by available resources and recorded. Default is 1. Explicit threads
wins over OMP_NUM_THREADS. Index building/global exact scans remain serial.
Without OpenMP the binary reports serial fallback. A normalized original-base
cache enables random access for plain/ordinary gzip FASTA; use
`--reference-cache FILE` for reuse. Original-reference inspection and cache
checksum cost remain in total time. The candidate-only byte counter excludes
this separately reported provenance/preparation I/O.

Persistent UTC text logs and matching JSONL default to `<output>.log` and
`<output>.log.jsonl`. Levels are trace/debug/info/warning/error/off;
`--no-log-to-console` disables console copies. Progress uses
`--progress-interval-seconds` and bounded `--progress-every-queries` batches.
Debug adds per-query scoring events; `--log-per-query`/`--log-per-bin` explicitly
enable expensive detail. Lines flush immediately. Signals/worker exceptions
produce terminal error events where graceful termination is possible.
Schema-3 main TSV/metrics finalize by rename; failed sidecars/partial files
must not be treated as successful runs. SIGKILL cannot produce a footer.

## Output reference: files, columns and interpretation

This reference describes the **current writers**, not just historical example
files. Start with the main TSV, read `scoring_method`/`primary_match_source`,
then use the per-bin tables for candidate detail and the metrics/logs for run
provenance. A verified k-mer is not a verified ortholog or a reconstructed
transcript. `fdr_control_per_bin` is a historical option name, not a statistical
FDR guarantee. FPR is FP/(FP+TN); the observed false-discovery proportion is
FP/(TP+FP). They have different denominators.

Contents: [file inventory](#output-file-inventory),
[main TSV](#main-result-tsv), [candidate tables](#per-bin-output-tables),
[optional reports](#optional-coordinate-repeat-and-calibration-reports),
[JSON/logs](#json-metadata-and-logs),
[Python utilities](#python-utility-output-files),
[legacy/negative-set outputs](#legacy-and-negative-set-outputs).

### Output file inventory

Let `O` be the exact `--output` path, including its extension. For example,
`--output results/run.tsv` produces `results/run.tsv.metrics.json`, NOT
`results/run.metrics.json`. Paths explicitly supplied on the CLI override
defaults. Relative paths are relative to the **working directory**, not the
TOML file's directory. V2 analyzes one reference file per invocation; a file
may contain multiple contigs.

| File | Created when | Contents / row unit |
|---|---|---|
| `O` | Every completed non-dry run | Main TSV; one row per input query, including zero-hit/invalid queries |
| `O.metrics.json` | Every completed non-dry run; override `--metrics` | Run settings, hashes, aggregate counts, timings and process HWM |
| `O.log` and `O.log.jsonl` | Logging enabled; override `--log-file` | Human-readable log and one JSON object per event; `--log-level off` disables files |
| `O.ibf_tuning.json` | IBF modes, including `--ibf-dry-run`; override `--ibf-tuning-report` | Requested/effective design, formulas, recommendations and alternatives; absent in exact mode |
| `P = O.per_bin.tsv` | Per-bin screening enabled; override `--per-bin-output P` | One paired candidate/truth summary per query |
| `P.candidates.tsv` | Same as `P` | One row per nominated candidate locus; no rows for a query with no candidates |
| `O.bins.tsv` | `--mode calibrate`; override `--calibration-bins` | One row per real reference bin, with occupancy and empirical membership errors |
| Path from `--fragments` | Option supplied | One row per real reference bin with original coordinates; available in exact mode too |
| Path from `--positions` | Option supplied in exact/validate/calibrate | One row per exact query-window/reference occurrence, NOT a unique-window count; `.gz` enables gzip |
| `O.repeat_classes.tsv` | `--repeat annotate`; override `--repeat-classes` | Query-side repeat-class window counts and selected-score support |
| Path from `--reference-hit-classes` | Option supplied with repeat annotation and exact-capable mode | Reference-side repeat attribution of globally exact occurrences |
| Path from `--discover-loci` | Option supplied in exact-capable mode | Loci clustered from ALL global exact occurrences; different from IBF-nominated candidates |
| `O.expected_loci.tsv` | `--expected-loci-bed` supplied; override `--expected-loci-output` | One row per supplied expected interval whose query occurs in the FASTA |
| No separate exon output file | `--exon-blocks` supplied | Exon-aware results are columns in `O` and aggregates in metrics |
| Path from `--save-index` plus `.metadata.json` | Index saved | Checked binary flat-IBF archive and readable metadata; no TSV columns in the binary |
| Path from `--reference-cache` | Candidate verification enabled | Reusable binary original-base cache, one byte/base; not a FASTA or an IBF. No sidecar: contig offsets are reconstructed from the original reference |
| Temporary cache in system temp directory | Verification enabled without explicit cache path | Same original-base access, normally removed at exit |
| `*.partial.<unique>` | Staged writes / interrupted run | Incomplete files, NEVER valid final results or indexes |

`--ibf-dry-run` inspects the reference and writes tuning/log information, but
does not allocate/query the index or emit a main result/metrics report.
There is no automatic compression of main/candidate TSVs merely because their
names end in `.gz`; only the positions writer and the Python merger implement
their documented gzip output behavior.

Use a **new output prefix/directory for each run**. Existing paths can be
overwritten; main TSV and metrics are staged separately, not as a transaction
covering all sidecars. A failed rerun can leave older final files beside newer
partial sidecars/logs. Check exit status, matching `run_id` in logs/metrics,
and terminal `run_completed`; neither file existence nor a candidate file
alone proves success. SIGKILL can leave staging files without a footer.

### Shared conventions and score definitions

- TSVs have a tab-separated header. `NA` means unavailable/undefined, not zero.
  JSON uses `null`. An available count of zero is written as `0`.
- Query IDs are the first whitespace-delimited FASTA header token. Embedded
  `|` characters are preserved. `reference_id` is the reference filename stem
  with FASTA/gzip suffix removed; `contig` is an actual FASTA record ID.
- `query_index` and `bin_index` are zero-based; candidate `rank` is one-based.
  Reference/query interval coordinates are **zero-based, half-open** `[start,end)`.
  The exception is positions `query_window_index`, a **one-based original
  window-start index**, including gaps left by invalid windows.
- Fractions are 0–1, not percentages. Query duplicates retain positional
  multiplicity; repeated reference occurrences do not multiply recovery.
  Multibin sums and occurrence counts can exceed the number of query windows.
- Let `N = valid_windows`, `G =` windows positive in any IBF bin, `C =` windows
  positive in any selected bin, `V =` windows exactly supported in the union
  of candidate loci, and `E =` independently globally exact-positive windows.
  The four recoveries are `G/N`, `C/N`, `V/N`, `E/N`. `V <= E` is expected
  under identical preprocessing, but `V` need not be bounded by `C`: a merged
  locus includes intervening sequence between nominated bins.
- Global error fields compare **G against E**, not the selected method.
  Candidate error fields compare **V against E**. Candidate-IBF error fields
  compare **C against E**. Error counts require independent truth.
- Main/metrics use schema 3 for explicit scoring or per-bin IBF production.
  Unspecified legacy exact/global invocations retain schema 2. Most optional
  side tables still label their own schema `2`; their number is independent
  of the main schema. Current candidate tables add `source_bin_ids` when the
  main run uses schema 3. Always inspect actual headers when joining old runs.

### Main result TSV

The current schema-3 header contains the schema-2 prefix plus method-specific
columns. Optional-report columns remain present and become `NA` when disabled.
The following tables cover every main column; groups follow their meaning
rather than requiring a particular column-number lookup.

| Identity / configuration column | Meaning |
|---|---|
| `schema_version` | Main result schema, 2 or 3 |
| `query_index` | Zero-based input-record order |
| `query_id` | Stable FASTA token, including embedded pipes |
| `reference_id` | Reference file ID, not the name of a selected contig |
| `k` | K-mer length |
| `strand_mode` | `forward`, `both` or `canonical` |
| `repeat_policy` | `include`, `exclude` or `annotate` |
| `ambiguity_policy` | `skip` or `error`; bases are never deleted to bridge N-runs |
| `mode` | `exact`, `ibf`, `validate` or `calibrate` |
| `fragment_size` | Reference bin start-range length in bases |
| `hash_functions` | Configured index hash count; exact mode does not build an index |
| `configured_per_bin_fpr` | Effective design setting, not empirical FPR; legacy prefix retains this configured value even in exact mode |
| `basic_threshold`, `strict_threshold` | Actual inclusive recovery thresholds used to classify this query |
| `threshold_source` | `cli` (including ordinary TOML settings) or `table` |
| `training_control_set`, `evaluation_control_set` | Supplied control-set labels, otherwise `NA`; labels alone do not prove biological independence |

| Window / primary-score column | Meaning |
|---|---|
| `total_windows` | `max(query_length-k+1,0)`, before ambiguity/repeat exclusions |
| `valid_windows` | `N`: valid A/C/G/T windows remaining under query masking rules |
| `invalid_windows` | `total_windows - valid_windows`; includes excluded repeat windows |
| `matched_windows` | Number supported by the SELECTED method, each original window at most once |
| `unmatched_windows` | `N - matched_windows`, not the number of reference mismatches |
| `distinct_query_kmers` | Distinct canonical identities among valid windows, even when `strand_mode=forward` |
| `distinct_matched_kmers` | Distinct canonical identities among selected-method matched windows |
| `total_bin_hits` | Sum of IBF-positive bin counts over valid windows; can exceed `N`; zero in exact-only runs because no IBF ran |
| `recovery_fraction` | `matched_windows/N`, or `NA` for `N=0` |
| `status` | Classification from ONLY `recovery_fraction`, using the rules below |
| `ibf_recovery` | `G/N`; `NA` in exact-only mode |
| `exact_recovery` | `E/N`; `NA` in IBF-only production |
| `false_positive_windows` | Windows with global IBF positive and global exact negative; NOT verified-local FP |
| `false_negative_windows` | Windows with global exact positive and global IBF negative; NOT candidate misses |
| `recovery_difference` | `(G-E)/N`, signed; global error fields are `NA` outside validate/calibrate |

Status priority is `NO_VALID_KMERS` for `N=0`, `NO_MATCH` for zero support,
then `PASS_STRICT` when recovery >= strict threshold, `PASS_BASIC` when
recovery >= basic threshold, otherwise `LOW_SUPPORT`. Passing strict also
satisfies basic, but only the stricter status string is emitted. Zero support
remains `NO_MATCH` even if a configured threshold is zero.

Precision warning: the current writer prints thresholds and many fractions
to six decimals, while classification uses their full internal values. A
fitted threshold just above 0.75 can print as `0.750000` yet correctly reject
an exact 0.75 score. Do not reconstruct boundary calls from rounded TSV cells.
Metrics preserve more threshold precision; downstream comparison precision is
an open technical finding in the review linked below. Default 0.75/0.85 runs
are not the demonstrated near-boundary fitted-threshold case.

| Schema-3 method / validation column | Meaning |
|---|---|
| `scoring_method` | Selected `global_ibf`, `per_bin_ibf`, `per_bin_verified` or `global_exact` |
| `primary_match_source` | Respectively `ibf`, `ibf_candidates`, `candidate_exact`, `exact` |
| `global_ibf_recovery` | Explicitly named alias of `ibf_recovery` (`G/N`) |
| `per_bin_candidate_count` | Number of candidate LOCI, not the number of bins |
| `per_bin_candidate_window_support` | `C`: distinct original windows with a positive probe in at least one selected bin |
| `per_bin_ibf_recovery` | `C/N`, approximate and unverified |
| `per_bin_verified_local_recovery` | `V/N`, union over all verified candidate intervals, not best-locus recovery |
| `per_bin_verified_locus_count` | Candidate loci with at least one exact-supported window; not necessarily passing basic/strict |
| `exact_global_recovery` | Explicitly named alias of `exact_recovery` (`E/N`) |
| `verification_scope` | Scope of SELECTED score: `none`, `candidate_intervals`, `whole_reference`; validation may compute additional scopes |
| `result_certainty` | `approximate` or `exact_membership`; exact membership does not imply biological homology |
| `index_fpr_per_bin` | Effective per-bin design value; `NA` without an IBF |
| `index_fpr_genome_target` | Requested global membership target when applicable; not observed FPR, not biological FDR |
| `exact_global_status`, `global_ibf_status`, `per_bin_ibf_status`, `per_bin_verified_status` | Each method independently classified with the same thresholds; `NA` if not evaluated |
| `candidate_window_recall` | `V/E`; `NA` without verification/truth or when `E=0` |
| `candidate_false_positive_windows` | `V`-positive / `E`-negative windows, counted by overlap, not by count subtraction |
| `candidate_false_negative_windows` | `E`-positive windows unsupported in candidates; equals `E-V` when verification is correct |
| `candidate_recovery_difference` | `(V-E)/N`, signed, requiring independent truth |
| `bins_rejected_min_hits` | Positive-hit bins failing the minimum hit count; bins with zero observed hits are not enumerated |
| `bins_rejected_fraction` | Positive-hit bins passing min-hits but failing minimum fraction; rejection categories are sequential, not overlapping |
| `per_bin_ibf_false_positive_windows` | `C`-positive / `E`-negative windows |
| `per_bin_ibf_false_negative_windows` | `E`-positive / `C`-negative windows |

| Query-repeat column | Meaning |
|---|---|
| `repeat_fraction` | Repeat-overlapping A/C/G/T query windows divided by all base-valid A/C/G/T query windows, BEFORE repeat exclusion; not a fraction of bases |
| `valid_nonrepeat_windows` | Valid query windows not overlapping softmask/BED repeats |
| `matched_nonrepeat_windows` | Such windows supported by the selected method |
| `nonrepeat_recovery` | `matched_nonrepeat_windows/valid_nonrepeat_windows` |
| `repeat_matched_windows` | Selected-method supported query windows overlapping repeats; zero if excluded on query side |

The following fixed/sliding localization columns require `--local-window L`
and global exact occurrences. Their denominator is `N`, not candidate support.
They are **not** the per-bin verified candidate summary.

| Localization column | Meaning |
|---|---|
| `global_recovery` | `E/N`, populated only when local-window reporting is enabled; use `exact_global_recovery` for the unambiguous truth field |
| `best_region_recovery`, `best_region_matched_windows` | Recovery/count of original query windows in the best fixed-width region; windows assigned by reference k-mer START |
| `number_of_supported_regions` | Fixed regions with any exact support, NOT a count of homologous loci |
| `candidate_contig`, `candidate_start`, `candidate_end` | Best FIXED region; `NA` if no supported region |
| `local_concentration` | Best fixed-region matched windows divided by `E`; `NA` when `E=0` |
| `off_region_fraction` | `1-local_concentration` |
| `multimapped_query_windows` | Original windows with at least two exact occurrence rows; this can include distinct strands at one genomic start |
| `sliding_best_region_recovery`, `sliding_best_region_matched_windows` | Best exact support in an L-base sliding start interval |
| `sliding_candidate_contig`, `sliding_candidate_start`, `sliding_candidate_end` | Best sliding interval, starting at an observed hit start |
| `boundary_gain_windows` | Sliding best count minus fixed best count; shows loss caused by arbitrary fixed-bin boundaries |

| Additional optional main column | Meaning / enabling option |
|---|---|
| `exon_eligible_windows` | Valid windows wholly contained within a supplied query exon block (`--exon-blocks`) |
| `exon_matched_windows` | Eligible windows with an exact occurrence at the block's mapped reference position and strand |
| `exon_aware_recovery` | `exon_matched_windows/exon_eligible_windows`; not global recovery |
| `number_of_discovered_loci` | Exact-hit clusters from `--discover-loci`; not IBF candidates |
| `best_locus_matched_windows`, `best_locus_recovery` | Best discovered exact locus support and count/N; no discovered loci gives count `0` and recovery `NA` |
| `best_locus_contig`, `best_locus_start`, `best_locus_end` | Best discovered exact locus interval, otherwise `NA` |
| `expected_locus_count` | Number of intervals supplied for this query by `--expected-loci-bed` |
| `expected_locus_matched_windows`, `expected_locus_recovery` | Distinct windows with a FULL k-mer occurrence contained in any expected interval, and that union count/N |
| `off_expected_locus_windows` | `E - expected_locus_matched_windows`; windows may match both inside and outside expected loci |
| `expected_locus_hit_rows` | Exact occurrence rows wholly contained in expected loci; can exceed unique supported windows |

Queries absent from supplied exon/expected-locus annotations receive `NA` for
those annotation-specific main fields, not invented zeros. Exact positions
and derived reports are subject to `--max-position-rows` PER BATCH, even if
the positions file itself was not requested. Hitting the cap aborts the run;
it does not silently truncate occurrence reports.

### Per-bin output tables

`P` has one row per query. Common identity fields are `schema_version` (sidecar
schema 2), `query_index`, `query_id`, `reference_id`, `k`, `valid_windows`.
Every other column is defined here:

| Paired column | Meaning |
|---|---|
| `exact_global_matched_windows`, `exact_global_recovery` | `E` and `E/N`; `NA` without independent truth |
| `ibf_union_matched_windows`, `ibf_union_recovery` | `G` and `G/N`, not the selected candidate score |
| `max_bin_hits`, `max_bin_support` | Highest single-bin `H_b` and `H_b/N` BEFORE filtering/top-K |
| `candidate_bins_before_top_k` | Bins passing BOTH min-hits and min-fraction |
| `candidate_bin_count` | Selected bins after top-K |
| `candidate_locus_count` | Selected bins clustered by contig/gap/span rules |
| `top_k_dropped_bins` | Eligible bins discarded by top-K; not all rejected bins |
| `best_locus_exact_matched_windows`, `best_locus_exact_recovery` | MAXIMUM single-candidate support and count/N; NOT the union `V` |
| `best_locus_exact_concentration` | Best single-candidate support divided by `E`; `NA` without global truth |
| `candidate_union_exact_matched_windows` | `V`: unique supported windows in the candidate union |
| `global_exact_windows_outside_candidates` | `E-V`; unavailable in production without global truth |
| `verified_candidate_loci`, `rejected_candidate_loci` | Candidate loci with nonzero versus zero exact support; `NA` if unverified |
| `union_status` | Basic/strict classification of GLOBAL any-bin `G/N` |
| `per_bin_status` | Discovery state, not basic/strict score: `NO_VALID_KMERS`, `NO_CANDIDATE`, `CANDIDATE_UNVERIFIED`, `VERIFIED_CANDIDATE`, `CANDIDATE_REJECTED` |
| `exact_status` | Basic/strict classification of `E/N`, otherwise `NA` |

`P.candidates.tsv` contains the nominated intervals, including technically
empty ones. A `VERIFIED` row means **at least one exact-supported window**,
not that the locus reaches either manuscript threshold or is biologically true.

| Candidate column | Meaning |
|---|---|
| `schema_version`, `query_index`, `query_id`, `reference_id`, `k` | Shared identities; sidecar schema remains 2 |
| `rank` | One-based deterministic rank: descending maximum bin hits, then summed bin hits, then contig/start/end |
| `contig`, `start`, `end` | Candidate interval for reference WINDOW STARTS; includes gaps between merged bins |
| `sequence_end` | `min(end+k-1,contig_length)`; fetched sequence extends beyond the start interval to complete boundary k-mers |
| `candidate_bin_count` | Selected bins belonging to THIS locus |
| `max_bin_hits`, `sum_bin_hits` | Maximum/sum of contributing `H_b`; sum can double-count query windows and is NOT recovery |
| `max_bin_support` | `max_bin_hits/N` |
| `exact_local_matched_windows`, `exact_local_recovery` | Unique original query windows with an exact start in this locus, and count/N |
| `exact_concentration` | This locus's exact-supported windows divided by `E`; `NA` without truth |
| `verification_status` | `VERIFIED` (support >0), `REJECTED` (support=0), or `UNVERIFIED` |
| `source_bin_ids` | Schema-3-run addition: comma-separated global zero-based bin IDs that formed the locus, in selected-bin rank order |

Bins rank by descending hits then reference contig order/start/bin ID. Top-K
precedes locus clustering. Candidate loci never cross contigs. Their gaps can
bridge introns, but merging cannot rescue bins removed by earlier thresholds.

### Optional coordinate, repeat and calibration reports

Repeated identity columns retain the definitions above. Most optional TSVs
have sidecar `schema_version=2`; positions have no schema column.

| Report | All columns, in emitted order; meaning |
|---|---|
| Fragment coordinates | `schema_version`, `bin_index`, `contig`, `start`, `end`, `sequence_end`, `inserted_reference_windows`, `distinct_reference_kmers`. `[start,end)` assigns k-mer STARTS; sequence extends to `sequence_end`. Insertions include repeated valid reference occurrences, while distinct counts are per-bin indexed identities. Both counts are `NA` in exact mode. |
| Exact positions | `query_index`, `query_id`, `reference_id`, `contig`, `start`, `end`, `strand`, `query_window_index`. `[start,end)` is the full k-mer; `strand` is query-relative `+`/`-`, or `.` for reverse-complement palindromes. Window index is one-based original query start. Multiple reference occurrences produce multiple rows. |
| Query repeat classes | `schema_version`, `query_index`, `query_id`, `repeat_class`, `valid_windows`, `matched_windows`, `recovery_fraction`. Counts are WITHIN this query repeat class; matched uses the selected primary method, and recovery is class matched/class valid. Nonrepeat class `NONE` is omitted; `SOFTMASK` and `MULTIPLE_CLASSES` can appear. |
| Reference hit classes | `schema_version`, `query_index`, `query_id`, `reference_id`, `reference_repeat_class`, `matched_query_windows`, `exact_occurrences`, `exact_matched_query_windows`, `fraction_of_exact_matched_windows`. First count deduplicates original windows within a class, occurrences do not; denominator is global `E`. Class `NONE` can appear. A window may match multiple reference classes, so class fractions need not sum to 1. |
| Discovered exact loci | `schema_version`, `query_index`, `query_id`, `reference_id`, `rank`, `contig`, `start`, `end`, `matched_windows`, `exact_occurrences`, `valid_query_windows`, `locus_recovery`, `locus_concentration`. Locus recovery=count/N; concentration=count/E. Occurrences are deduplicated `(contig,start,query-window)` events (strand is not a separate event). Loci rank by supported-window count, then coordinates; span includes complete k-mers. |
| Expected loci | `schema_version`, `query_index`, `query_id`, `reference_id`, `contig`, `start`, `end`, `label`, `valid_query_windows`, `matched_windows`, `locus_recovery`, `exact_occurrences`. Counts require the full k-mer inside the supplied interval; recovery=count/N. An annotated query/interval with no hits still gets a zero-count row; label defaults to `.`. |

Calibration `--mode calibrate` produces the following per-bin columns:

| Calibration column | Meaning |
|---|---|
| `schema_version`, `bin_index`, `contig`, `start`, `end` | Sidecar schema 2 and zero-based start-assignment coordinates |
| `inserted_reference_windows`, `distinct_reference_kmers` | Repeated insertions versus distinct indexed keys in this bin |
| `bin_bits` | Bits in ONE bin, not total index bytes |
| `occupied_bits`, `bit_occupancy_fraction` | Set bits in this bin and `occupied_bits/bin_bits` |
| `insertions_per_bit` | `inserted_reference_windows/bin_bits`; not an FPR |
| `tested_absent_query_windows` | Original valid query windows with no exact match in THIS bin, summed over queries |
| `false_positive_bin_pairs` | Such absent window/bin pairs for which the IBF said positive |
| `false_negative_bin_pairs` | Exact-present window/bin pairs for which the IBF said negative |
| `empirical_per_bin_fpr` | False-positive pairs/tested absent query windows for THIS bin; `NA` when denominator=0 |

Exon-block INPUT format is a TSV with header `query_id`, `query_start`,
`query_end`, `contig`, `ref_start`, `ref_end`, `strand`. Query/reference spans
must have equal positive lengths and query blocks must not overlap. Expected
locus INPUT is BED4/5: `contig start end query_id [label]` (tabs). Repeat INPUT
is BED3/4: `id start end [class]`. These are supplied annotations, NOT generated
homology evidence. Interval validation preserves original coordinates/Ns.

### JSON metadata and logs

`O.metrics.json` is an object, not a per-query table. Below, groups list exact
top-level key names. All `*_seconds` fields are wall seconds unless explicitly
CPU; all `*_bytes` fields are bytes. Process RSS is a lifetime high-water mark,
NOT current RSS or a per-stage incremental allocation. Stage timings overlap
in some cases and should not be blindly summed. Some disabled-stage counters
remain zero; interpret them with mode/options and use `null` where emitted.

| Metrics keys | Meaning |
|---|---|
| `schema_version`, `software`, `software_version`, `compiler_version`, `git_sha` | Result schema and build provenance; a Git SHA alone does not identify a dirty working tree |
| `semantic_fasta_hash_encoding`, `reference_sequence_sha256`, `query_sequence_sha256`, `reference_sequence_fnv1a64`, `query_sequence_fnv1a64` | Semantic FASTA hashes include IDs, original characters/case and record boundaries but not wrapping/header descriptions; FNV is supplementary, not a cryptographic integrity guarantee |
| `config_path`, `config_sha256`, `reference`, `queries` | Input/config paths and raw config hash; absent config gives empty strings |
| `reference_repeat_bed_sha256`, `query_repeat_bed_sha256`, `exon_blocks_sha256`, `expected_loci_bed_sha256`, `threshold_table_sha256` | Raw auxiliary-input SHA-256; empty string when not supplied |
| `mode`, `strand_mode`, `repeat_policy`, `mask_target`, `k`, `fragment_size`, `hash_functions` | Actual analysis settings; `mask_target` chooses query/reference/both for exclusion |
| `reference_repeat_bed`, `query_repeat_bed`, `repeat_classes`, `reference_hit_classes`, `exon_blocks`, `discovered_loci`, `expected_loci_bed`, `expected_loci_output`, `positions`, `fragments`, `calibration_bins`, `per_bin_output` | Auxiliary input/output paths; empty strings mean disabled |
| `max_repeat_intervals`, `max_exon_blocks`, `max_locus_rows`, `max_expected_loci`, `max_query_bases` | Input/output caps; exceeding a cap aborts rather than silently truncating |
| `locus_gap`, `max_locus_span`, `local_window`, `max_position_rows_per_batch` | Exact-locus/localization settings and occurrence cap; max span also constrains candidate clustering |
| `fdr_control_per_bin`, `per_bin_exact_verify`, `per_bin_min_hits`, `per_bin_min_fraction`, `per_bin_top_k`, `per_bin_max_locus_gap` | Effective candidate settings; the first flag is NOT statistical FDR control |
| `threshold_source`, `threshold_table`, `training_control_set`, `evaluation_control_set`, `basic_threshold`, `strict_threshold` | Threshold values and supplied provenance labels |
| `index_format_version`, `index_provenance`, `index_path`, `index_file_sha256`, `disk_index_bytes` | Archive version, `built`/`loaded`/disabled provenance, path/hash/size; no archive written means zero disk bytes, not a zero-memory index |
| `bin_count`, `bin_bits`, `configured_per_bin_fpr`, `requested_per_bin_fpr_cap`, `target_genome_wide_fpr` | Real bins, bits per bin, effective/requested design and applicable global target |
| `nominal_any_bin_fpr_independent_bins`, `observed_occupancy_union_estimate` | Theoretical independent-bin design estimate and occupancy-based union estimate; neither is biological FDR nor measured transcript FPR |
| `physical_occupied_bits_all_modes`, `distinct_bin_kmers_sum`, `occupied_bin_bits_sum` | Physical set-bit count (includes padded storage), sum of distinct per-bin keys (not distinct genome-wide), and real-bin set-bit sum in calibrate mode |
| `physical_occupancy_seconds`, `bit_occupancy_measurement_seconds`, `peak_rss_after_bit_occupancy_bytes` | Occupancy scan cost and process HWM afterwards; calibration-specific fields can be `null` |
| `max_calibration_pairs_per_batch`, `max_ibf_bin_pairs_per_batch` | Caps on exact key/bin associations and stored positive query-window/bin pairs; not total-process memory limits |
| `scoring_method`, `primary_match_source`, `status_counts` | Selected score/source; mapping from main status strings to query counts |
| `threads`, `requested_threads`, `parallel_backend` | Actual observed team size, requested workers and compiled `openmp`/`serial` backend; backend can be OpenMP with one worker |
| `run_id`, `log_file`, `log_level`, `warning_count`, `ibf_tuning_report` | Log join ID, settings, warning count and diagnostics path; empty tuning path in exact mode |
| `candidate_reference_bases_fetched`, `candidate_reference_bytes_read`, `candidate_interval_reads` | Logical candidate-cache bytes/bases and read tasks; first two are equal for one-byte/base cache and INCLUDE repeated k-1 task padding. Not physical disk I/O, and EXCLUDE provenance/cache preprocessing reads |
| `candidate_repeated_reference_keys_skipped`, `candidate_key_cache_insertions`, `candidate_key_cache_limit_per_worker` | Matched repeated reference keys skipped, bounded cache insertions over tasks, and per-task/worker limit (8192) |
| `candidate_verification_seconds`, `candidate_selection_seconds`, `per_bin_comparison_seconds` | Interval verification, ranking/selected-bin support, and independent per-bin truth comparison times |
| `whole_reference_exact_scans`, `whole_reference_exact_bases_visited` | Independent full truth scans and bases actually visited. Verified production should have zero; reference inspection/cache checking is a separate pass |
| `reference_cache_preprocessing_bases`, `reference_cache_reused` | Bases processed during cache/provenance preparation, and whether an explicit existing cache was checksum-validated/reused; count is not necessarily newly written bytes |
| `query_batches`, `total_reference_bases`, `reference_ambiguous_bases`, `reference_softmasked_bases`, `valid_reference_windows`, `query_count`, `valid_query_windows` | Batch/input QC totals; reference masking counts are bases, query window totals are windows |
| `per_bin_queries_with_candidates`, `per_bin_candidate_bins`, `per_bin_candidate_loci`, `per_bin_verified_candidate_loci`, `per_bin_rejected_candidate_loci`, `per_bin_top_k_dropped_bins` | Candidate workload totals; loci verified/rejected mean any exact support versus none, not basic/strict calls |
| `per_bin_missed_global_exact_queries`, `per_bin_global_exact_windows_outside_candidates` | Exact-supported queries with no verified candidate support, and missed exact-positive windows; `null` without candidate verification and global truth |
| `discovered_locus_rows`, `expected_locus_matched_windows`, `expected_locus_recovery` | Optional global-exact locus/expected-locus totals; expected recovery pools only annotated queries' valid windows |
| `exon_eligible_windows`, `exon_matched_windows`, `exon_aware_recovery` | Optional eligible/mapped exon support and pooled matched/eligible recovery |
| `tested_absent_windows`, `exact_positive_windows` | Aggregate `N-E` and `E`; `null` without independent truth |
| `false_positive_windows`, `false_negative_windows`, `measured_false_negative_rate`, `measured_genome_wide_fpr` | GLOBAL any-bin errors relative to truth, and respectively FN/E and FP/(N-E); NOT selected verified-local error. `null` without an evaluated comparison or denominator |
| `tested_absent_bin_pairs`, `false_positive_bin_pairs`, `false_negative_bin_pairs`, `measured_per_bin_fpr` | Pooled per-bin comparison totals; FPR is sum FP pairs/sum absent pairs, NOT the unweighted mean of bin FPRs |
| `reference_inspection_seconds`, `peak_rss_after_reference_inspection_bytes`, `index_construction_seconds`, `index_load_seconds`, `index_save_seconds`, `peak_rss_after_index_bytes` | Reference/preparation and index phase wall costs/HWM; cold construction and loaded-index query costs are different experiments |
| `occupancy_guard_seconds`, `peak_rss_after_occupancy_guard_bytes`, `fragment_output_seconds`, `peak_rss_after_fragment_output_bytes` | Additional guard/report phase costs and lifetime HWM checkpoints |
| `query_seconds`, `exact_scan_seconds`, `result_output_seconds` | Query-window/IBF processing, independent global exact scans, and scoring/report emission times; candidate timing is separate |
| `reference_valid_window_throughput_per_second`, `query_valid_window_throughput_per_second` | Valid windows divided by construction/query phase time; zero when that phase has zero duration, not an end-to-end throughput |
| `total_seconds`, `peak_resident_bytes`, `user_cpu_seconds`, `system_cpu_seconds` | In-process wall/HWM/CPU observations taken before final metrics close, rename and final output checksum/logging; external command wall can be larger and CPU sum can exceed wall under parallelism |

The v2 metrics are NOT the same as `scripts/measure_command.py` reports. Its
chosen output JSON contains `command` (argv array), `exit_code`, `wall_seconds`,
`peak_child_resident_bytes`, `reaped_child_peak_resident_bytes`,
`linux_sampled_descendant_peak_bytes`, `rss_sampling_interval_seconds`,
`rss_measurement_scope`. The reported peak is maximum single-process HWM,
not summed process-tree RSS. Linux sampling supplements reaped-child rusage;
it can miss short-lived/orphaned processes. Old forced-kill measurements can
omit an unreaped child: the failed whole-assembly pilot's 13.7-MB wrapper
value is explicitly NOT its program peak. JSON `exit_code=-9` denotes a
signal-killed measured command; nested wrapper shell exit codes may differ.

Tuning JSON keys and units:

| Tuning keys | Meaning |
|---|---|
| `scoring_method`, `auto_tune`, `legacy_budget_rule`, `changed` | Objective/policy, whether historical automatic budgeting applied, whether effective FPR changed |
| `bin_count`, `contig_count`, `fragment_size`, `bin_bits`, `estimated_index_bytes`, `max_index_memory_gb` | Real layout and raw index estimate; GB budget uses 1e9 bytes, NOT total process RSS |
| `requested_per_bin_fpr`, `effective_per_bin_fpr`, `genome_target`, `p_union_independent`, `union_bound` | Requested/effective designs and theoretical one-key global estimates. `genome_target` is the configured/default objective, not proof it was applied to a per-bin method |
| `maximum_lookup_orientations`, `query_window_union_bound` | Both-strand windows can probe two orientations; window-level budget includes that factor |
| `unsafe_global_settings`, `formula`, `assumptions`, `load_estimator`, `recommendation` | Risk flag, exact formula/assumptions, conservative load model and actionable advice |
| `hash_sensitivity` | Array with `hashes`, `same_p_bin_bytes`, `same_p_bin_error`, `budgeted_p_bin`, `budgeted_bytes`, `budgeted_error`; unrepresentable UNCHOSEN alternatives use `null` plus explanatory error |
| `fragment_alternatives` | Array with `fragment_size`, `bins`, `estimated_bytes`, `error`, `applied`; overflowing alternatives can include `scale` instead of `bins`. Alternatives are estimates, never silently applied |

Index metadata JSON keys are `index_format_version`,
`semantic_fasta_hash_encoding`, `software_version`, `compiler_version`,
`seqan_version`, `size_t_bytes`, `byte_order`, `reference_sequence_sha256`,
`reference_repeat_bed_sha256`, `strand_mode`, `ambiguity_policy`, `repeat_policy`,
`mask_target`, `k`, `fragment_size`, `fragment_overlap`, `bin_count`, `bin_bits`,
`hash_functions`, `configured_per_bin_fpr`, `valid_reference_windows`,
`distinct_bin_kmers_sum`, `disk_index_bytes`, `index_file_sha256`,
`requested_per_bin_fpr`, `ibf_auto_tune`, `tuning_report`,
`tuning_report_sha256`, and `contigs` (array of `id`, `length`). These describe
the saved archive, its ABI/settings/reference identity and the tuning report
used for construction. `fragment_overlap=k-1`; the metadata is NOT a portable
FASTA index. Legacy `.ibf` files are a different format and cannot be loaded
as checked v2 archives.

JSONL log fields are `timestamp` (UTC ISO-8601 with milliseconds), `run_id`,
`level`, `phase`, `event`, `elapsed_seconds`, `phase_seconds`, `message`.
Events include `stage_started`, `stage_progress`, `stage_completed`, `warning`,
`stage_failed`, `run_completed`, `run_failed`. Counts embedded in `message`
remain message text, not independent JSON numeric fields. A text line has the
same information in human-readable form. Messages are thread-safe and flushed;
`trace`/`debug` add query detail, and explicit per-bin logging can be expensive.

### Python utility output files

These are helper-generated artifacts, not additional v2 main columns.

| Utility | Output files and schemas |
|---|---|
| `sample_fasta.py` | Requested FASTA (original headers/order, sampled with `seed`) plus default `<sample>.metrics.json` or `--metrics`: `source`, `output`, `total_records`, `sampled_records`, `seed`, `wall_seconds`, `peak_resident_bytes` |
| `prepare_production_benchmark.py` | Validated `reference_first_contig.fa` and `input_inventory.json`; staging uses `.partial`. Inventory keys: `created`, `host`, `git_sha`, `query_release`, `selected_assembly`, `reference_scope`, `reference_base_cap`, `genomes`, `queries`, `benchmark_reference`. Per-input QC: `path`, `file_bytes`, `raw_sha256`, `type`, `record_count`, `total_bases`, `contigs` (`id`, `length`, `ambiguous`, `lowercase`), `contig_details_complete`, `ambiguous_bases`, `lowercase_bases`, `extracted_bases`. Benchmark reference has `path`, `raw_sha256`, `original_contig_id_preserved`. Default extraction is ONE complete contig, not the entire assembly. |
| `run_per_bin_production_benchmark.py` | Per-case v2 outputs, `<case>.command.json`, `<case>.config.json` (`command`, `timeout_seconds`, `actual_query_records`, `executable_sha256`, `source_sha256`); named cached `.ibf`/`.bases`; `serial_parallel_checksums.json` or `optimized_serial_parallel_checksums.json` maps worker counts to main/sidecar SHA-256. Whole stage adds `benchmark_scope.json`: `reference_scope`, `reference`, `raw_sha256`, `query_count_cap`, `k`. Sweep adds `development.fa`, `heldout.fa`, `threshold_split.json`, `threshold_selection.json` as described below. Use a new experiment directory; completed-name reuse is not a fresh provenance validation. |
| `analyze_per_bin_production.py` | `method_comparison.json`, `method_comparison.tsv`, `timings.tsv`, `checksums.json`; includes only successful command records with metrics. Checksums map relative filenames to raw SHA-256, excluding reference copies and the checksums file itself. JSON contains additional class/length detail beyond the TSV. |
| `merge_kmer_v2.py` | Requested merged TSV (same main columns; query-major/reference-ID-sorted rows; optional deterministic gzip) and `--summary` JSON. Summary: `schema_version`, `query_count`, `reference_count`, `query_sequence_sha256`, `primary_match_source`, `ibf_target_genome_wide_fpr`, `ibf_requested_per_bin_fpr_cap`, `inputs`, `references`, `wall_seconds`, `peak_resident_bytes`. Each reference entry: `query_count`, `valid_windows`, `matched_windows`, `false_positive_windows`, `false_negative_windows`, `statuses`, `window_weighted_recovery`. Recovery uses pooled selected matched/valid counts. Current merger initializes error aggregates to zero even when all source errors are `NA`; those zeros MUST NOT be interpreted as measured zero FP/FN. |
| `fit_kmer_thresholds.py` | Requested threshold TSV: `reference_id`, `k`, `basic_threshold`, `strict_threshold`, `training_control_set`; plus summary keys `schema_version`, `training`, `evaluation`, `training_label`, `evaluation_label`, `target_basic_null_pass_rate`, `target_strict_null_pass_rate`, `groups`, `wall_seconds`, `peak_resident_bytes`. Each group has `reference_id`, `k`, `training_count`, `evaluation_count`, both thresholds and `training_basic_pass_count`, `training_strict_pass_count`, `evaluation_basic_pass_count`, `evaluation_strict_pass_count`. This tool currently accepts schema-2 controls only, NOT schema-3 production TSVs. Empirical null pass rate is not FDR control. |

For the production analyzer, all `method_comparison.tsv` columns are:

| Columns | Meaning |
|---|---|
| `run`, `k`, `queries`, `method` | Case filename, k, total input query count and evaluated prediction method |
| `valid_windows`, `exact_positive_windows`, `matched_windows` | Pooled `N`, truth `E`, and method support |
| `false_positive_windows`, `false_negative_windows`, `window_recall`, `window_fpr` | THAT method versus independent truth; FP/(N-E), recall=(E-FN)/E; `NA` if unavailable |
| `classification_agreement` | Fraction of classifiable queries with the same five-class label as exact truth; not statistical confidence |
| `mean_absolute_recovery_error`, `median_absolute_recovery_error`, `p95_absolute_recovery_error` | Per-query absolute recovery error summary; p95 uses the implementation's sorted integer-rank selection, not interpolation |
| `any_match_tp`, `any_match_fp`, `any_match_fn`, `any_match_tn`, `any_match_sensitivity`, `any_match_specificity`, `any_match_precision`, `any_match_f1` | Endpoint is AT LEAST ONE matched window, not a homologous transcript |
| `basic_tp`, `basic_fp`, `basic_fn`, `basic_tn`, `basic_sensitivity`, `basic_specificity`, `basic_precision`, `basic_f1` | Endpoint is recovery >= per-row basic threshold |
| `strict_tp`, `strict_fp`, `strict_fn`, `strict_tn`, `strict_sensitivity`, `strict_specificity`, `strict_precision`, `strict_f1` | Endpoint is recovery >= per-row strict threshold |

For every confusion endpoint: sensitivity=TP/(TP+FN), specificity=TN/(TN+FP),
precision=TP/(TP+FP), F1=2TP/(2TP+FP+FN). Zero denominators give `NA`/`null`,
not 100%. Queries with zero valid windows are excluded from these fractions.

The JSON report is keyed by result filename. Per-run keys are `path`,
`reference_scope`, `k`, `biological_truth` (false), `query_count`,
`zero_valid_queries`, `valid_windows`, `exact_positive_windows`,
`candidate_window_recall`, `candidate_false_negative_windows`,
`candidate_false_positive_windows`, `global_ibf_false_positive_windows`,
`top_k_dropped_bins`, `methods`. Per-method keys include the TSV statistics
above (JSON calls the agreement `classification_agreement`), `query_count`
(classifiable queries), `weighted_recovery`, `any_match_confusion`,
`threshold_confusions` (`basic`, `strict`), `class_confusion` (counts keyed
`truth -> prediction`) and `length_strata_by_valid_windows`. Length groups
`lt500`, `500to2000`, `ge2000` contain `queries` and
`mean_absolute_recovery_error`. Full exact is the truth/reference row, not an
independent prediction-validation success. Currently the analyzer expects
the complete four-method validation column set; it is not a general TSV reader.

`timings.tsv` columns: `run`, `queries`, `k`, `method`, `threads`, `wall_seconds`,
`peak_rss_bytes`, `index_build_seconds`, `index_load_seconds`, `index_bytes`,
`query_seconds`, `candidate_selection_seconds`, `verification_seconds`,
`exact_seconds`, `candidate_bases`, `global_exact_scans`,
`requested_per_bin_fpr`, `effective_per_bin_fpr`. Wall is external command
time; RSS here comes from the binary's own `peak_resident_bytes`; candidate
bases are logical cache reads. Index save time remains in metrics but is not
a separate timings.tsv column.

`threshold_split.json` records `seed`, `grouping`,
`biological_family_independence_established`, and `development`/`heldout`
objects containing `query_ids`, `gene_ids`. Gene-disjoint is not automatically
family-independent. `threshold_selection.json` records `objective`,
`development_results`, `selected`; each result contains `configuration`
(`[label,min_hits,min_fraction,top_k,gap]`), `window_recall`, `candidates`.
Only development scores select parameters; held-out outputs remain separate.

Historical utilities are **not needed for new production runs**:

- `correct_per_bin_union_status.py` intentionally rewrites old paired TSVs,
  keeps requested backups and writes `union_status_corrections.json` with
  `schema_version`, `records`; each record has `path`, `backup`, `rows`,
  `changed_union_status_rows`, `raw_sha256`, `corrected_sha256`. It changes
  only historical miswritten `union_status` cells. Do not run on new outputs.
- `analyze_per_bin_ibf.py` consumes the fixed historical A0/A1/B0/B1/B2/T case
  names. It writes `summary.tsv`, `threshold_grid.tsv` (subset of the same
  columns), `paired_100.tsv`, `occupancy_summary.tsv`, `aggregate.json`.
  Its explicit historical column dictionary is below; do not substitute these
  legacy endpoint labels for the current four-method validation report.

<details>
<summary>Historical analyzer column dictionary</summary>

`summary.tsv` / `threshold_grid.tsv`: `run`, `k`, `queries`, `valid_windows`,
`exact_present_windows`, `exact_absent_windows` are case and input/truth counts;
`window_union_fp`, `window_union_fpr`, `window_union_fn` are global union error;
`candidate_queries`, `candidate_bins`, `candidate_loci`, `verified_loci`,
`rejected_loci`, `rejected_locus_fraction`, `missed_exact_queries`,
`exact_windows_outside_candidates`, `top_k_dropped_bins` describe candidate
burden/loss. `technical_endpoint`, `TP`, `FP`, `FN`, `TN`, `technical_precision`,
`technical_recall`, `technical_specificity` compare any globally exact support
against any verified candidate. `wall_seconds`, `peak_rss_bytes`,
`queries_per_second`, `valid_windows_per_second` are external resource/rate
fields; `index_load_seconds_internal`, `index_build_seconds_internal`,
`index_save_seconds_internal`, `query_seconds_internal`,
`exact_scan_seconds_internal` are internal phase costs. `index_bytes`,
`bin_count`, `bin_bits`, `effective_per_bin_fpr`, `target_whole_reference_fpr`,
`index_sha256`, `reference_semantic_sha256`, `query_semantic_sha256` record
layout and provenance.

`paired_100.tsv`: `query_id`, `k`, `exact_global_recovery`,
`exact_global_matched_windows`, `exact_status`; for `A0`/`A1`, suffixes
`_ibf_union_recovery`, `_union_status`; for `B0`/`B1`/`B2`, suffixes
`_max_bin_support`, `_candidate_bin_count`, `_candidate_locus_count`,
`_best_locus_exact_recovery`, `_exact_windows_missed`, `_per_bin_status`.
A0=nominal historical global, A1=budgeted global, B0=unverified nominal bins,
B1=verified nominal bins, B2=verified budgeted bins. Values have the same
support/count definitions as the paired table, not biological accuracy.

`occupancy_summary.tsv`: `index` (high/low design), `bins`, `occupied_bits`,
`total_bin_bits`, `aggregate_occupancy_fraction`,
`median_bin_occupancy_fraction`, `max_bin_occupancy_fraction`,
`occupancy_scan_seconds_internal`, `run_wall_seconds`, `peak_rss_bytes`,
`index_bytes`, `index_sha256`. The aggregate fraction pools bits, not bin FPRs.

`aggregate.json`: `schema_version`, `runs`, `paired_queries`, `source`,
`method_ids`, `run_manifest`. Each manifest entry has `id`, `command_json`,
`result_tsv`, `metrics_json`, `per_bin_tsv`, `index_sha256`,
`reference_semantic_sha256`, `query_semantic_sha256`.

</details>

### Legacy and negative-set outputs

These belong to `lncrna_mers` / `negative_set`, NOT the v2 main TSV. Legacy
FASTA parsing/fragment layout and raw archives differ from v2; never mix their
columns, thresholds or binary indexes as if interchangeable.

| Legacy file | Contents |
|---|---|
| `output_dir/output_file` | When `single_results_writer=true`: one query row, repeated reference-specific column groups |
| `output_dir/results_<ref>.tsv` | When `single_results_writer=false`: the same columns for one reference |
| `output_dir/unique_mers/<ref>.tsv` | `query_index`, `ibf_unique_kmer`: one distinct query sequence word whose probe hits EXACTLY ONE bin; empty cell row when none. Not a genome-unique or exact-verified sequence. |
| `output_dir/<ref>_kmers.tsv` | `query_index`, `kmer_indices`, `kmer_counts`: slash-separated parallel lists of one-based query-window starts and nonzero positive-bin counts; empty lists if no hits |
| `output_dir/<ref>_fragments.fasta` | With `store_fragments`: emitted sequence fragments, IDs `<ref>_frag<zero-based index>`; no coordinate TSV |
| `output_dir/<ref>.ibf` or `.hibf` | Raw legacy archive when storage enabled; not checked v2 format. Multi-index IBF writes `<ref>_multi_<index>.ibf` and current workflow cleans these temporary archives after use. |
| `output_dir/log_file` | Legacy text log (local-time timestamps), not v2 JSONL |

Legacy main columns are `query_index` (zero-based) and, for every reference:

| Legacy column | Actual formula |
|---|---|
| `<ref>_count` | SUM of positive-bin counts over query windows; not the count of matched windows |
| `<ref>_ibf_unique_kmer` | Number of DISTINCT query sequence words having a probe positive in exactly one bin |
| `<ref>_pass` | 1 when `<ref>_count >= hit_threshold`, otherwise 0; not a recovery threshold |
| `<ref>_pct` | `<ref>_count / max(query_length-k+1,0)` (zero for no windows), stored as a fraction with four decimals, NOT multiplied by 100; can exceed 1 and must not be interpreted as v2 recovery |

`query_index` is the only legacy query key; join back to FASTA order rather
than guessing IDs. Reference enumeration is not currently sorted/uniqueness
checked: avoid filenames that collapse to the same reference ID. Flat-IBF
ambiguity-safe canonical matching does not make these legacy scores equivalent
to v2 selected-window scoring. Optional HIBF remains outside the validated
production path and is not being integrated by this documentation change.

`negative_set` writes configurable `output_fasta` (default
`negative_set.fasta`), `output_tsv` (default `negative_set.tsv`) and `log_file`
under `output_dir`. Generated FASTA IDs are `negative_<zero-based index>`.
The shuffled samples originate from selected input chromosomes; they are not
guaranteed globally exact-negative, and do not automatically establish a
biological null distribution.

| Negative-set TSV column | Meaning |
|---|---|
| `id`, `chromosome` | Generated sample ID and original source chromosome ID |
| `start`, `end`, `length` | Original sampled zero-based half-open span and shuffled sequence length; coordinates identify the UNshuffled source, not an alignment of the shuffled output |
| `total_kmers` | `max(length-k+1,0)` |
| `matched_kmer_positions` | Shuffled query windows with any positive IBF bin |
| `match_count` | Sum of positive-bin counts, potentially above `total_kmers` |
| `match_fraction` | `matched_kmer_positions/total_kmers`, or zero for no windows; bounded 0–1, still approximate IBF membership |

For the source-based technical review and prioritized next steps, see
[technical code review](docs/TECHNICAL_REVIEW_2026-10-18.md). No technical
recommendation here means that biological FDR or manuscript readiness has
been established.

The CTest `output_documentation_contract` runs a tiny synthetic calibration
with every optional v2 writer enabled and checks that all emitted TSV names,
metrics/tuning/index-metadata keys and JSONL fields are documented here. Run
it alone with `python3 scripts/measure_command.py --output /tmp/output-docs.json
-- python3 tests/test_output_documentation.py build-production/lncrna_mers_v2`
(as one command). This checks documentation coverage and basic fixture values,
not biological validity or a proof that every implementation is bug-free.

The older schema-2 documentation below is retained for compatibility/history;
its full-reference candidate verifier and IBF-only restrictions are superseded.
[![CMake tests](https://github.com/GPSeq/FragIncRNA/actions/workflows/cmake-single-platform.yml/badge.svg?branch=reporting)](https://github.com/GPSeq/FragIncRNA/actions/workflows/cmake-single-platform.yml)
`FragIncRNA` builds per-reference fragment indexes from FASTA references, queries lncRNA sequences against them, and writes match summaries plus unique matching k-mers. The index backend can be either a flat `ibf` or a hierarchical `hibf`.

## Corrected k-mer analysis (schema 2 preview)

`lncrna_mers_v2` is the corrected analysis path. It reports per-window recovery in [0,1], default canonical strand matching, independent exact validation, and wall-time/peak-memory measurements. Its default `validate` mode now classifies from exact verified windows and reports Bloom-filter recovery separately. The TOML driven `lncrna_mers [config.toml]` retains its original output columns, but its flat-IBF path now shares v2's position-preserving, ambiguity-safe, canonical k-mer encoder and uses `k-1` fragment overlap. `negative_set` shares those flat-IBF rules too. Historical flat-IBF result values can therefore change when rerun. The legacy `*_pct` still sums bin hits and can exceed 1; do not treat it as a recovery fraction. The v2 path supports flat IBF for one reference per invocation. See [methods](docs/KMER_METHODS.md), [audit](docs/KMER_AUDIT.md), [resource measurements](docs/KMER_BENCHMARKS.md), [real chr1 run](docs/KMER_REAL_DATA_FIRST_RUN.md), and [implementation status](docs/KMER_IMPLEMENTATION_STATUS.md) before using it for biological conclusions.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLNCRNA_MERS_ENABLE_HIBF=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j 4
ctest --test-dir build --output-on-failure
./build/lncrna_mers_v2 --reference tests/fixtures/kmer_v2/reference.fa --queries tests/fixtures/kmer_v2/queries.fa --output /tmp/fraglncrna-v2.tsv --k 4 --fragment-size 4 --mode validate
./build/lncrna_mers_v2 --config configs/example.toml
```

The CLI example writes `/tmp/fraglncrna-v2.tsv` and its metrics JSON. The fully commented [TOML example](configs/example.toml) covers every v2 option; `--config` can be followed by CLI overrides, and the config file's SHA-256 is recorded in the run JSON. Legacy `lncrna_mers` continues to use its existing `configs/config_*.toml` layout. Run `./build/lncrna_mers_v2 --help` for options. Plain and `.gz` FASTA are accepted. The output is one row per query and reference, with stable FASTA IDs and explicit zero-hit rows. `--mode exact` skips IBF construction; use it for the supplied chr1 setting. IBF-only mode now sizes the per-bin filter for a default whole-reference FPR target of 0.01 (`--genome-fpr`) and checks observed bit occupancy before classifying. `--mode calibrate` measures genome-wide per-window and per-bin false positives against exact truth and writes a `.bins.tsv` report with distinct k-mer counts and actual bit occupancy. `--save-index reference.ibf` stores the flat IBF with checked compatibility metadata; `--load-index reference.ibf` reuses it. The metrics JSON contains semantic FASTA SHA-256 hashes, auxiliary input hashes, archive size, wall time, and peak RSS. `--fragments fragments.tsv` writes explicit reference bin and overlap coordinates. `--local-window 8000 --positions /tmp/verified-hits.tsv.gz` adds exact positions, fixed-region support, and sliding-window sensitivity with a bounded occurrence cap. `--discover-loci candidates.tsv` groups exact hits into ranked locus candidates; `--expected-loci-bed expected.bed` writes a separate expected-locus report and per-query recovery columns. `--exon-blocks blocks.tsv` reports eligible-window recovery at supplied splice-alignment coordinates. `--reference-repeat-bed`, `--query-repeat-bed`, and `--mask-target` control repeat annotation and exclusion; `--repeat annotate` writes query class counts and `--reference-hit-classes hits.tsv` writes exact reference-side class attribution. A threshold table keyed by reference/k can be fitted on separate training and evaluation controls with `scripts/fit_kmer_thresholds.py`; `scripts/merge_kmer_v2.py` streams verified per-reference runs into a query-major comparison. Exact validation rescans the reference for each query batch, so benchmark a subset before planning a full-genome run.

For the first real-data check, use the [measured 5,000-transcript runbook](docs/KMER_REAL_DATA_FIRST_RUN.md). `scripts/sample_fasta.py` makes a seeded subset and writes its own time and peak-memory JSON.

### Optional per-bin coherence screen

`--fdr-control-per-bin` adds candidate fragment scores and exact candidate-region
verification without changing the existing main TSV columns. The requested
flag name does **not** imply statistical FDR control: it measures how many
original query windows hit each IBF bin, groups nearby candidate bins, and
checks selected regions against the original reference. The paired per-query
table defaults to `<output>.per_bin.tsv`; candidate regions go to
`<output>.per_bin.tsv.candidates.tsv`. Use `--per-bin-min-hits`,
`--per-bin-min-fraction`, `--per-bin-top-k`, and
`--per-bin-max-locus-gap` to set the heuristic screen. The same Boolean and
numeric options are available in `[v2]` TOML; see
[the complete example](configs/example.toml) and [method definition](docs/PER_BIN_IBF_METHOD.md).
For an unverified historical-index diagnostic, add `--mode ibf
--no-per-bin-exact-verify`; the separate candidate status must not be read as
global recovery. The [comparison report](docs/PER_BIN_IBF_BENCHMARK.md) records
measured chr1 outcomes and resource costs.

## Requirements

- CMake 3.20+
- A C++20 compiler

SeqAn3 is fetched automatically by CMake during configure. This project currently pins SeqAn3 `3.4.1`.
HIBF support is optional and currently requires a newer compiler toolchain than the one used for the older flat IBF-only setup. In practice, use GCC 12+ and CMake 3.25+ for HIBF builds.

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

If sdsl-lite is causing any errors from seqan3 build, use:

```bash
git clone https://github.com/xxsds/sdsl-lite.git
cd sdsl-lite
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/.local
cmake --build build -j
cmake --install build
```

This produces the executable `build/lncrna_mers`.

## HIBF Toolchain Notes

If CMake prints:

```text
HIBF support disabled: GCC 11.4.0 is too old. HIBF requires GCC >= 12.
```

then your current compiler is too old for the HIBF library. The flat `ibf` backend can still be built, but `hibf` needs GCC 12 or newer.

On Ubuntu or WSL Ubuntu, install GCC 12:

```bash
sudo apt update
sudo apt install gcc-12 g++-12
```

Check that the compiler is available:

```bash
which g++-12
g++-12 --version
```

If CMake reports:

```text
The CMAKE_CXX_COMPILER:
  g++-12
is not a full path and was not found in the PATH.
```

then either `g++-12` is not installed, or it is not visible on `PATH`. In that case, use the full compiler paths:

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-12 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-12
cmake --build build -j
```

You can also set the environment variables instead:

```bash
export CC=/usr/bin/gcc-12
export CXX=/usr/bin/g++-12
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j
```

If your CMake version is very new and configure fails while processing the fetched SeqAn3 project, add:

```bash
-DCMAKE_POLICY_VERSION_MINIMUM=3.5
```

Example:

```bash
rm -rf build
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j
```

After building with a supported compiler, enable HIBF in [config.toml](config.toml):

```toml
[general]
index_method = "hibf"
```

## Run

The program now reads its settings from a TOML file instead of command-line flags.

1. Edit [config.toml](config.toml).
2. Run:

```bash
./build/lncrna_mers
```

By default it loads `./config.toml`. You can also pass a different TOML file path:

```bash
./build/lncrna_mers /path/to/config.toml
```

Example `config.toml`:

```toml
[general]
index_method = "ibf"

ref_dir = "/path/to/references"
query_file = "/path/to/queries.fa"

output_dir = "./ibf_results"
output_file = "results.tsv"
log_file = "ibf_run.log"

fragment_size = 8000
kmer_size = 15
hit_threshold = 13
threads = 8

store_fragments = false
store_index = false
cleanup_index = false
single_results_writer = false

[ibf]
hash_functions = 3
fpr = 0.01
multi_indexing = false

[hibf]
hash_functions = 2
maximum_fpr = 0.05
relaxed_fpr = 0.30
threads = 8
sketch_bits = 12
tmax = 0
empty_bin_fraction = 0.0
alpha = 1.2
max_rearrangement_ratio = 0.5
disable_estimate_union = false
disable_rearrangement = false
```

## Config Keys

The TOML file is split into three sections:

- `[general]`
- `[ibf]`
- `[hibf]`

General keys:

- `general.index_method`: `"ibf"` or `"hibf"`
- `general.ref_dir`: directory containing `.fa`, `.fna`, `.fasta`, and gzipped variants
- `general.query_file`: FASTA file with query sequences
- `general.fragment_size`: fragment length used for reference splitting; must be at least `4`
- `general.kmer_size`: k-mer length used for indexing and querying; valid range `1..32`
- `general.hit_threshold`: minimum total k-mer hits required for a passing match
- `general.threads`: number of references to process in parallel; use `0` to auto-detect hardware concurrency
- `general.output_dir`: directory for logs, results, stored indexes, and optional fragment FASTA files
- `general.output_file`: combined TSV file name when `single_results_writer = true`
- `general.log_file`: log file name written under `output_dir`
- `general.store_fragments`: when `true`, writes `<reference>_fragments.fasta`
- `general.store_index`: when `true`, stores one serialized index file per reference
- `general.cleanup_index`: when `true`, removes the stored index file after processing
- `general.single_results_writer`: `true` writes one combined TSV; `false` writes one `results_<reference>.tsv` per reference

IBF-specific keys:

- `ibf.hash_functions`: number of IBF hash functions; valid range `1..32`
- `ibf.fpr`: target false positive rate per bin; must be in `(0, 0.5]`
- `ibf.multi_indexing`: when `true`, builds one flat IBF per reference fragment, splits each fragment into near-equal sub-bins, stores the `.ibf` files in `output_dir` while searching, and removes them after that reference is processed

HIBF-specific keys:

- `hibf.hash_functions`: number of HIBF hash functions; valid range `1..5`
- `hibf.maximum_fpr`: maximum false positive rate; must be in `(0, 1)`
- `hibf.relaxed_fpr`: relaxed FPR for merged bins; must be in `(0, 1)` and `>= maximum_fpr`
- `hibf.threads`: internal HIBF build threads; `0` falls back to `general.threads`
- `hibf.sketch_bits`: HyperLogLog sketch bits; valid range `5..32`
- `hibf.tmax`: maximum technical bins for each IBF level; `0` lets HIBF choose
- `hibf.empty_bin_fraction`: reserved empty-bin fraction; valid range `[0.0, 1.0)`
- `hibf.alpha`: HIBF alpha parameter; must be positive
- `hibf.max_rearrangement_ratio`: valid range `[0.0, 1.0]`
- `hibf.disable_estimate_union`: disables union estimation during layout
- `hibf.disable_rearrangement`: disables rearrangement during layout

## Output

Combined mode writes a TSV with these columns per reference:

- `*_count`: total IBF hits across bins for the query/reference pair
- `*_ibf_unique_kmer`: number of distinct query k-mers that
  1. have IBF count exactly `1` for that reference
  2. are counted once per distinct k-mer, even if they appear multiple times in the query
- `*_pass`: `1` if `count >= hit_threshold`, else `0`
- `*_pct`: `count / number_of_query_kmers`

Unique kmers are also written under `unique_mers/<reference>.tsv`, with
`query_index` in the first column using the 0-based order from the query FASTA file.
All matching query k-mers are additionally written to `<output_dir>/<reference>_kmers.tsv`
with three columns: `query_index`, slash-separated 1-based `kmer_indices`, and the
matching slash-separated `kmer_counts` for those indices.
The main results TSV also writes `query_index` instead of the full FASTA query name.

## Tests

Configure and build as usual, then run:

```bash
ctest --test-dir build --output-on-failure
```

To run only the unit test binary directly:

```bash
./build/lncrna_mers_tests
```
