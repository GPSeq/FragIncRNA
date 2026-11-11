# FragIncRNA output reference

This document lists every file, TSV column, JSON key and log field written by
`lncrna_mers_v2`. The CTest `output_documentation_contract` checks that every
emitted name is documented here or in the README.

This reference describes the **current writers**, not just historical example
files. Start with the main TSV, read `scoring_method`/`primary_match_source`,
then use the per-bin tables for candidate detail and the metrics/logs for run
provenance. A verified k-mer is not a verified ortholog or a reconstructed
transcript. `fdr_control_per_bin` is a historical option name, not a statistical
FDR guarantee. FPR is FP/(FP+TN); the observed false-discovery proportion is
FP/(TP+FP). They have different denominators.

Contents: file inventory,
main TSV, candidate tables,
optional reports,
JSON/logs,
Python utilities,
legacy outputs.

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
Metrics preserve more threshold precision.

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

### Legacy `lncrna_mers` outputs

These belong to the legacy `lncrna_mers`, NOT the v2 main TSV. Legacy
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

### Sharded index, production and window-report fields

Explicit layout adds TSV provenance fields (existing schema-3 meanings are kept):
`ibf_layout`, `shard_count`, `max_bins_per_ibf`, `global_bin_count`, `index_manifest`,
`manifest_sha256`, `per_bin_fpr_requested`, `per_bin_fpr_effective`,
`candidate_exact_verification`, `global_error_budget_policy`.
The metrics include those fields, plus `shards_built`, `shards_reused`,
`shard_query_loads`, `shard_query_bytes_read`, `shard_query_load_seconds`,
`max_resident_shards`, and `query_scheduler`.
Shard count is a physical-layout property; biological interpretation still depends
on `scoring_method`, `primary_match_source` and `verification_scope`.

The Linux worker counters `process_io_read_bytes` and `process_io_write_bytes`
measure storage I/O, while `process_io_read_characters` and
`process_io_write_characters` include cached read/write syscall traffic.
Shard build insertion uses OpenMP over disjoint 64-local-bin storage words,
as permitted by the SeqAn3 flat-IBF thread-safety contract. Each shard log
records its actual insertion thread limit; encoding, streams and distinct-key
sets belong to individual workers. Single-layout insertion remains serial.
After disk-mode release the log records current worker RSS and zero resident
shards. Total worker peak RSS also includes query/reference and archive overhead.

`global_error_budget_policy` is `union_bound` when the conservative global
budget rule was applied, `per_bin` when sizing remains per-bin, or
`requested_only` when an explicit genome target was recorded without applying
it under the selected tuner/method. Explicit recommend/off preserves FPR;
apply changes global_ibf. Supplying genome_fpr without an explicit tuner keeps
the historical explicit-target budgeting rule. Target columns are requests,
not achieved error-rate guarantees. `shard_query_bytes_read` sums serialized
archive sizes per load; checksum passes cause additional read traffic captured
by the worker read-character and physical-I/O counters.

Production cache and window-report metrics are additive. `window_output`,
`window_scope`, `window_rows`, and `kmer_fasta_output` describe optional exports.
`canonical_kmer_encoding` specifies the canonical ID representation.
`ibf_shard_cache_budget_bytes` is the effective cache bound;
`shard_cache_hits`, `shard_cache_misses`, `shard_cache_evictions`, and
`shard_cache_peak_bytes` measure deserialized-object reuse. `shard_checksum_count`,
`shard_checksum_bytes`, and `shard_checksum_seconds` measure whole-archive
verification once per shard per execution. `shard_first_pass_seconds` and
`shard_second_pass_seconds` include their shard visits; `shard_global_selection_seconds`
is global ranking time. `shard_query_load_seconds` measures query-time
metadata parsing/deserialization without whole-payload hashing. These timings overlap;
do not sum load time and pass time as disjoint stages.

`reference_page_cache_hits`, `reference_page_cache_misses`, and
`reference_page_cache_budget_bytes` describe the bounded shared reference cache.
`candidate_logical_requested_bases` sums requests before cross-query merging;
`candidate_unique_start_bases` sums the merged union of eligible genomic start
coordinates per batch; `candidate_unique_intervals` counts those merged intervals.
The legacy `candidate_reference_bytes_read` and `candidate_reference_bases_fetched`
both retain their old meaning: decoded task bases including k−1 lookahead.
`candidate_reference_userspace_read_bytes` measures actual bytes read from the
normalized reference file by the verifier (including page overfetch).
`process_io_read_bytes` is the OS-reported physical IO for the entire process,
including input/index verification; it is not verifier-only IO or OS cache hits.

`ibf_local_cache_directory` records an explicitly requested optional disk copy.
`shard_local_cache_copied_bytes` and `shard_local_cache_checksum_count` record its
copy work and integrity checks separately from canonical-archive checksums.

`logical_ibf_positive_pairs` is a checked 64-bit logical count, independent of
materialized calibration limits. `unique_query_lookup_keys` sums per-query
unique lookup keys (it does not deduplicate across transcripts).
`query_encoding_seconds` measures encoding and grouping separately for sharded
runs; it is null for single IBFs where the existing encoder also probes.
`candidate_interval_grouping_seconds` measures owner/mask setup and interval
normalization, a subset of `candidate_verification_seconds`.
`candidate_reference_read_seconds` measures userspace seek/read wall time;
`reference_page_cache_lock_wait_seconds` sums worker wait time for the page
cache lock. These are aggregate worker times, not isolated kernel IO wait.

The atomic `run_manifest.tsv` records reference paths/file hashes, k, settings
hashes, index/output directories, status, runtime, query count and validation.
The entire assembly output directory is staged, files are flushed, a
`COMPLETED.sha256.tsv` ledger is sealed, then the directory is atomically
published. Failed/interrupted staging directories are retained for diagnosis.
A retry verifies completed ledgers and skips their genomes. Completed results
are never replaced. The fingerprint includes effective arguments, the config,
input bytes and executable hash: changed settings or binaries require a new
output directory. A partial genome restarts from its query beginning while
reusing valid index checkpoints. There is no query-row append checkpoint.
`continue_on_reference_error` defaults false; enabling it processes later
assemblies after a failure and still returns nonzero. Shared directory-run
inputs must be regular files; descriptors detect changes between assemblies.

Window reports are disabled by default. Use `window_output`/`--window-output`
with TSV or .tsv.gz; `window_scope` selects all, verified, candidate or failed
(valid windows lacking an exact candidate hit). `window_query_ids` filters
using one exact query ID per line, and unknown IDs fail. `max_window_rows`
defaults 10 million across the run; exceeding it fails with staged reports
instead of publishing a silently truncated report. A large all-window export
can be much larger than transcript summaries. Directory runs apply the cap
independently to each assembly.

Window schema 1 includes `k_index` (the original zero-based sliding ordinal),
`query_start`, `query_end` (half-open), `valid`, `canonical_kmer_id`,
`ibf_any_bin_positive`, `ibf_candidate_positive`,
`exact_verified_in_candidates`, and `matching_candidate_regions_count`.
Duplicate positions remain separate; invalid/masked positions retain their
indices and use NA for sequence/ID and unavailable predicates. Canonical IDs
are dna2bit-v1:k:hex: A=0,C=1,G=2,T=3, first base most significant, minimum of
forward and reverse-complement encodings, hex without required leading zeroes.
The k is explicit, so this is an exact encoding, not a probabilistic hash.
`window_include_sequence` adds `kmer_sequence` with the actual query substring
(including original case) for valid windows. Rows also carry scoring, strand,
repeat, masking and verification-scope provenance. `kmer_fasta_output` streams
one valid original occurrence per record, optionally gzip, with percent-escaped
IDs carrying reference/query/index/window coordinates; duplicates retain IDs.

`candidate_prefilter` enables a small rejection filter before the verifier's
exact owner map (default true; `--no-candidate-prefilter` disables it for
measurement). `verification_prefilter_rejections` counts definite absences;
`verifier_full_key_lookups` counts subsequent exact map probes. Every owner key
sets its filter bits, so this shortcut cannot remove a true key. False filter
admissions only cost an extra exact lookup and do not alter reported IBF FPR,
candidate recall, verification scope, or exact window masks. Reference pages
are 16 KiB; an earlier 64 KiB trial overfetched on sparse 80-kb loci.
