# Analysis catalog — the source of truth

One row per analysis: what it answers, what to import, required columns/roles, the
statistic, an example file to try it on, and what it outputs.

**How these analyses were checked.** During development every statistic in this catalog was
matched numerically against a MATLAB builtin on fixed synthetic input — `ttest` / `ttest2` and
`corr` / `partialcorr` for the legacy tier, `fitlm` / `anova1` / `anovan` for the GLM presets,
`fitlme` for the mixed model, and independent references for the circular statistics — and each
analysis was run end to end on data with a known planted cluster, paired with negative controls
(a null design, and a deliberately mismatched permutation scheme, must *not* produce a cluster).
The measured error rates quoted throughout this page come from those runs. The verification
harness itself is internal and is not distributed, so what you can reproduce from a clone is the
cluster-recovery half: run any file from
[`example_data/`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/example_data/README.md) and
confirm the planted cluster near E129 comes back.

## Input sources
Two ways to feed any analysis:
- **CSV files** — wide (one row per subject, channel columns `E1..E256`/`Cz`) for the
  legacy t-tests/correlation and the circular tier; a single combined table for the GLM presets;
  long (one row per trial) for `mixedmodel`. Excel (`.xlsx`+sheet) or `.csv`.
- **Spectral folder** (multi-folder) — one EEG_processor folder per cohort/condition; each
  folder is a level of the factor. Covers `pairedT/onesampleT/unpairedT/anova1/rmanova`
  only (see `scripts/README_scripts.md`). Two cohort folders with **different subject counts**
  are fine on the `unpairedT` route (that combination used to throw before reaching a statistic).

Channel columns must use real montage labels (`E1..E256`/`Cz` for EGI, `Fp1..` for
Compumedics); set Recording system accordingly. Anything not a channel label is metadata.

**Source space (2447 cortical voxels).** Recording system `Source 2447 (cortical voxels)`
(`source2447`) analyses reconstructed cortical band power instead of scalp electrodes. The
inverse solution is computed **upstream** (GeoSource); the toolbox ingests a subjects×2447
matrix whose columns are **`src0001..src2447` in canonical graph-node order** (enforced by
`snpm_assert_source` — column *i*, neighbour-graph row *i* and coordinate *i* must all be the
same voxel) and runs the same permutation-t / TFCE / cluster-extent engines (the scalp t-test's
`snpm_cluster_analysis`, reused unchanged) on the source adjacency graph
(`NeighborMatrix_Sources_2447_Full.mat`). Values must be non-negative **band-
power magnitudes**, not signed current density — the `absolute` group t-test **rejects negative
input** (`snpm:source:signedMagnitude`). See the how-to
[Run a source-level (2447-voxel) analysis](../how-to/run-a-source-level-analysis.md) for the
export contract and the **placeholder-coordinate** precondition.

## Catalog

Every path in the **Example data** column is relative to
[`example_data/`](https://github.com/TancyKao/hdEEG-SnPM/blob/main/example_data/README.md) at the
repo root, which is in the repository — there is nothing to generate.

| Analysis (key) | Research question | Import | Required cols / roles | Statistic (ref) | Example data |
|---|---|---|---|---|---|
| **Paired t** `pairedT` | Same people, 2 conditions differ? | 2 wide files (A, B), row-matched | — | paired `ttest` | `paired_condA.csv` + `paired_condB.csv` |
| **One-sample** `onesampleT` | Within-subject A−B ≠ 0? | 2 wide files (A, B) | — | `ttest`(A−B) | `onesample_condA.csv` + `onesample_condB.csv`; or `onesample_change.csv` alone (single file vs 0) |
| **Unpaired t** `unpairedT` | 2 independent groups differ? | 2 wide files (A, B), **group sizes may differ** | — | `ttest2` | `unpaired_groupA.csv` + `unpaired_groupB.csv` |
| **Correlation P** `correlationP` | Linear EEG↔measure assoc? | 2 wide files, each w/ `Subject` | Subject col | Pearson `corr` | `corr_eeg.csv` + `corr_behavior.csv` |
| **Correlation S** `correlationS` | Monotonic assoc (robust)? | as above | Subject col | Spearman `corr` | as above |
| **Circular phase, 2 groups** `circ_phase_group` | Do 2 independent groups couple at a different phase? | 2 wide angle files in **radians** (**disjoint** Subject ids) **+ 2 wide event-count files** (GUI: *Event Count File*) | `circ_units`, `circ_convention`, `count1_file`, `count2_file` | covariate-adjusted Hotelling T² on (cos,sin), F(2, N−q−3) | `circ_phase_groupA.csv`, `circ_phase_groupB.csv` + `circ_phase_countsA.csv`, `circ_phase_countsB.csv` |
| **Circular phase, 2 groups (U²)** `circ_phase_group_u2` | As above, omnibus over the whole distribution | as above (**no whole-degree input**) | as above | Watson's U² (no covariate slot) | as above; also `circ_conc_groupA.csv`, `circ_conc_groupB.csv` + `circ_conc_countsA.csv`, `circ_conc_countsB.csv` |
| **Circular–linear** `circ_corrAngLinear` | Does a behavioural measure vary with phase? | 1 wide angle file + 1 measure file, **matched** Subject ids | `circ_units`, `circ_convention`, `measure_col` | circular–linear correlation (Mardia 1976), model F | `circ_anglinear_angles.csv` + `circ_anglinear_measure.csv` |
| **One-way ANOVA** `anova1` | 3+ groups differ? | 1 file | `group_col` (≥2 lvls) | F vs `anova1` | `glm_anova1.csv` |
| **ANCOVA** `ancova` | Groups differ adjusting covariates? | 1 file | `group_col` + `covariate_cols` | t (2 groups) / F (≥3, Freedman–Lane) | `glm_ancova.csv` (2 groups → t); `glm_ancova3.csv` (3 groups → F + pairwise post-hoc) |
| **Regression** `regression` | Scales with a predictor? | 1 file | `predictor_col` (+ `covariate_cols`) | t on slope vs `fitlm` | `glm_regression.csv` |
| **RM-ANOVA** `rmanova` | 3+ within conditions differ? | 1 file, long | `subject_col` + `condition_col` | F (within perm) | `glm_rmanova.csv` |
| **Two-way mixed ANOVA** `mixed2way` | Condition effect differs between groups? | 1 file, long | `group_col` + `subject_col` + `condition_col` | F on interaction | `glm_mixed2way.csv` (2×2); `glm_mixed2way_2x3.csv` (2×3) |
| **Mixed model** `mixedmodel` | Trial-level effect (many obs/subject)? | 1 long file | `Subject` + DV + fixed + effect | `fitlme` t/F (Stephan 2021) | `lmm_long.csv` |
| **Source-level t** (any wide t-test with Recording system `source2447`) | Which cortical voxels differ? | wide file(s), cols `src0001..src2447` (magnitudes) | — (same as the scalp t-test) | permutation t + TFCE + cluster-**extent**, on the source graph | none ships — source input comes from your own GeoSource export (see the [how-to](../how-to/run-a-source-level-analysis.md)) |

**Circular tier (`core_snpm_circ`).** Routed from `core_snpm_analysis` by an early guard, like the
GLM and LMM tiers — the three keys above are **not** `compstring` cases and must not be added to
the legacy stat-function switches. Hard constraints, all enforced with named errors:

| Constraint | Identifier |
|---|---|
| `circ_units` (`rad`/`deg`) required, no default | `core_snpm:circUnitsRequired` |
| `circ_convention` required (`literature_uppeak0` / `yasa_uppeak0` / `luna_zerocross0` / `turtlewave_pre_v4` / `custom`) | — |
| Pooled frontal grand mean in the inverted window ⇒ reject | `core_snpm:circPhaseConventionInverted` |
| Whole-degree angles on the U² path ⇒ reject (warn on Hotelling) | `core_snpm:circResolutionTooCoarse` / `…Coarse` |
| Precision files required for both group keys | `core_snpm:circCountsRequired` |
| `datatype` must be `absolute` | `core_snpm:circDatatypeNotSupported` |
| `tail` must be `both` (all three statistics are non-negative, one-tailed) | `core_snpm:circTailNotSupported` |

Permutation scheme is **Freedman–Lane on the stacked (cos,sin) response**, `free` (whole-subject
relabelling) — between-subject only; there is no within-subject circular scheme yet. Cluster
statistic is **mass**. There is **no global/omnibus test** (arithmetic channel averaging is
meaningless for angles); a **descriptive**, explicitly no-inference circular panel is emitted
instead. This tier was checked the same way as the rest: the three statistics against independent
references, and end-to-end planted-phase-shift recovery with null-design and
permutation-scheme negative controls, a precision-confound control, and the
convention/inversion guards. Background:
[About circular statistics for phase](../explanation/circular-statistics-for-phase.md); task
guide: [Run a circular (phase) analysis](../how-to/run-a-circular-phase-analysis.md).

**Precision covariate.** The nuisance covariate for the two group keys is **log(Rayleigh Z)**
(Z = nR², per subject per channel), fitted **per channel**. It subsumes event count, which is
retained only as a fallback where Z is unavailable; the two must never both enter one model. Z is
a covariate when the outcome is *phase* and an outcome when the analysis is about *coupling
strength* — never both in one model. There is **no Rayleigh mask**: weakly-coupled cells are kept
and adjusted, not deleted. *Current build:* the engine reads these two files as integer event
counts and rejects non-integer input (`core_snpm:circCountsInvalid`); the log(Z) path lands with
the queued engine change.

**Source path coverage.** The source2447 system is wired into `core_snpm_analysis` (the wide
t-tests — `pairedT`/`onesampleT`/`unpairedT`) and `core_snpm_lmm` (`mixedmodel`). The GLM
presets (`anova1`/`ancova`/`regression`/`rmanova`/`mixed2way`, via `core_snpm_glm`) are **not
yet source-wired** — run those on scalp montages. TFCE is the primary correction; the cluster
statistic on the source **t-test** path is cluster **extent** (size), the same as the scalp
t-tests (`snpm_cluster_analysis`). Cluster **mass** applies only to the source **LMM** path
(`core_snpm_lmm`, mean-Wald).

## Permutation p-values and sampling

- **Minimum-bias p everywhere.** All tiers (legacy t-tests/correlation, GLM, LMM) report the
  permutation p as `(b + 1) / (N + 1)`, where `b` = permutations at least as extreme as the
  observed statistic and `N` = permutation count. The p is never exactly 0; the smallest
  reportable value is `~1/(N+1)` (the legacy tiers force the observed labeling into the null,
  so their floor is `~2/(N+1)`). A result cannot clear α unless `N` is large enough — e.g.
  `N = 100` floors p near 0.01–0.02, so use `permutations ≥ 5000–10000` for final inference.
- **Exact vs Monte-Carlo (legacy `unpairedT` and correlation).** When the requested
  `permutations` meets or exceeds the exchangeability-group size — `nchoosek(nSubj, nGrp)` for
  the unpaired two-group relabeling, `nSubj!` for correlation — the engine runs an **exact**
  test over the full set of labelings (denominator = the true group size); otherwise it samples
  with replacement (Monte-Carlo). This also removes a former hang when more permutations than
  distinct labelings were requested. *Small-n caveat:* an exact test has a coarse p-grid — e.g.
  unpaired `nGrp = 3, nSubj = 6` → only 20 permutations → minimum p ≈ 0.048 — so significance
  can be unreachable with very small samples regardless of the requested count.
- **`unpairedT` accepts unequal group sizes.** Group A and group B may contain different numbers
  of subjects. (Before the 2026-08 fix this never worked at all: `core_snpm_analysis` combined an
  n1×1 and an n2×1 logical unconditionally, so any between-subjects t-test with differently-sized
  groups threw a dimension error before it reached a statistic.) This matters most on the
  **spectral-folder** route, where two cohort folders routinely differ in subject count — that
  configuration now runs end to end.
- **Missing data: a channel is tested only if it is complete.** The correlation paths now follow
  the same rule as the GLM presets. A channel is evaluable only if **every analysed subject has a
  value there**; the mask is computed **once**, before the permutation loop, and applied to the
  observed map and to every permuted map. Excluded channels appear in the `excludedChannels`
  Excel sheet and the HTML banner (see Outputs below). This is deliberately conservative — **one
  missing cell costs the whole channel** — and the reason is that the previous per-channel
  pairwise deletion recomputed the complete-pair set *inside* the permutation loop, so the
  observed statistic and its null were computed at different sample sizes. Measured family-wise
  error at 40% disjoint missingness was **0.288** against a nominal 0.05; it is now **0.050**.
  Complete-data results are bit-identical to before, so no existing complete-data analysis
  changes.
- **`per_channel_n` changed meaning.** It is now computed **after** the mask: a retained channel
  reports the full analysed n, an excluded channel reports 0. That is the n actually used in the
  statistic, but it is not the old pairwise count — if you are comparing output from before and
  after the 2026-08 fix, this is why the numbers differ.
- **The t-tests exclude channels that are degenerate under relabelling** (2026-08-03). This is a
  narrower rule than the complete-column rule above: `pairedT`/`onesampleT` drop a channel with
  **fewer than 2 usable pairs** (df = 0 returns ±Inf under every sign flip, and an infinite value
  also breaks the TFCE threshold ladder with `MATLAB:pmaxsize`); `unpairedT` drops a channel with
  **fewer than `max(n1,n2) + 2` usable subjects pooled**, the threshold below which some
  relabelling leaves one group holding a single subject — MATLAB's `var` of a scalar is 0, not
  NaN, so `ttest2` returns a *finite* t from a pooled SD that ignores that group (measured max
  |t| = 7.4 on pure noise with n1 = 1, n2 = 10) and that value sets the max-statistic threshold
  for the whole montage. Every other channel is still tested on the pairs it has. Measured cost of
  not doing this: a planted 9-channel cluster was TFCE-detected 93.4% of the time with no
  degenerate channel present and **16.4%** when four channels were available in only 3 of 20
  subjects; masking those four restored it to 94.1%. Excluded channels are reported exactly like
  the correlation exclusions. No-op on complete data.
- **The whole-head omnibus averages over a common channel set** (2026-08-03). `global_stat_test`
  averages only the channels **finite in every analysed row of both arms**, instead of averaging
  each unit over whatever channels it happened to have. Scalp power has a strong spatial gradient
  (roughly 4× front-to-back in NREM delta), so averaging two units over different channel sets
  puts a per-unit offset into the contrast that no sample size removes. Measured at nominal 0.05
  with no true effect: an unpaired test where 4 of 10 subjects lost an occipital band ran at
  0.112, and a paired test with *different* regions missing in the two conditions ran at
  **1.000**; under the rule all configurations return to 0.047–0.053. The reported `n_channels` /
  `channels_kept` say what the average was taken over. No-op for the correlation path (its
  complete-column mask already makes the two sets identical) and no-op on complete data. The
  property that was checked is the one the old code violated, rather than a golden number: given
  the retained set, the global statistic is bit-identical whether or not the dropped cells were
  ever dropped. When no channel is common to both arms the statistic is
  **NaN with a `global_stat_test:noCommonChannel` warning**, not a number — the channel-wise maps
  are unaffected.
- **Repeated-measures guard.** The between-subject GLM presets (`anova1`, `ancova`,
  `regression`, which use `free` permutation) **error** (`core_snpm_glm:repeatedMeasures`) if a
  subject appears in more than one row, and direct you to `rmanova` or `mixedmodel`. Feed them
  one row per subject.

## Outputs (every analysis, via `core_snpm_analysis`)
- `<base>_<timestamp>.mat` — `results_struct` with `T` (incl. `real_T`), `p`, `Clusters`,
  `uncorrsigch`, `correctTFCEsigch` (TFCE-corrected), `SnPMsigch` (cluster-corrected), `chanlocs`.
- `<base>_<timestamp>.xlsx` — significant-channel table (`func_genSnpmTable`).
- **Excluded channels.** A channel that cannot be evaluated in *every* permutation is not tested
  at all: the GLM presets and the circular tier require a channel to be complete (no missing cell)
  because the Freedman–Lane nuisance fit works on whole columns, so one missing cell would make
  the column all-NaN inside each permuted fit while leaving the observed map finite. Excluded
  channels are listed in an **`excludedChannels`** sheet in the `.xlsx` and in a banner in the
  HTML report, and in `results_struct.excluded_channels`. The behaviour this replaced was
  **anti-conservative** — measured family-wise error **0.066** against a nominal 0.05, now
  **0.054**. Read the excluded list, and *where* the exclusions fall, before reading the topoplot.
- **Correlation runs also emit** each channel's effective N as `results_struct.per_channel_n` and
  a guarded **`effectiveN`** sheet in the `.xlsx`; the partial-correlation global p uses
  `df = n − 2 − k` (`k` = number of covariates). Since the 2026-08 fix `per_channel_n` is the n
  **after** the completeness mask — full n for a retained channel, 0 for an excluded one — not the
  old pairwise matched-pair count.
- `<base>_<timestamp>_report.html` — one self-contained report with a **TFCE/Cluster toggle**, stat-aware (t-map for the t-tests, **F-map + post-hoc** table for the omnibus ANOVA presets), the mean topographies and both significance maps. (Replaced the former separate `_TFCE_report.html` / `_Cluster_report.html`.)
- Topoplot PNGs.

**Source-space runs (`source2447`) instead of a scalp topoplot:** a 2-D scalp disc is
meaningless for cortical voxels, so `core_snpm_analysis`/`core_snpm_lmm` **bypass every EEGLAB
topoplot** (zero scalp PNGs) and emit a significant-voxel list via `write_source_voxel_table`:
`<base>_<timestamp>_sigvoxels.csv` (one row per voxel — `idx, label, X, Y, Z, stat, p, pTFCE,
sig_uncorr, sig_TFCE, in_sig_cluster, cluster_id, cluster_p`) plus the same table as a
`sigVoxels` sheet in the `.xlsx`; `results_struct.is_source` is set true. **Coordinates are
PLACEHOLDER** (graph-Laplacian layout in `source2447_coords.mat`) until the real GeoSource MNI
export is dropped in — the statistics, significance flags and cluster membership are exact; only
the `X/Y/Z` columns are placeholder, and the source path was checked on stats, flags and cluster
membership only, never on the coordinates.

**Reading the cluster output:** the GLM/LMM presets score a cluster by **mean Wald**
(`sum(Wald)/n_channels`, height-weighted) while the legacy t-tests score by **extent** (size),
so cluster statistics are **not comparable across the two tiers**; and a significant cluster is
a **regional** claim, not a per-electrode one. See
[Interpreting cluster results](../explanation/interpreting-cluster-results.md).

Specialized reports: GLM faceted report `export_glm_report.m` (anova/t templates in
`templates/`); spectral dashboard `export_report.m`; spectral sweep `SWEEP_grid_*.csv`.

## Choosing an analysis

Design trade-offs (categorical vs continuous predictors, unequal-n / dropout, and the
LMM-vs-GLM cost decision) are discussed in the explanation
[Choosing an analysis: design, missing data, and cost](../explanation/choosing-an-analysis.md).
That page also carries the sleep-question → preset cheat-sheet and the method-grounding references.

## Known gaps
- **No within-subject circular scheme.** The circular tier is between-subject only: the
  permutation relabels whole subjects. Two angle files containing the *same* subjects is a
  within-subject design and is not yet supported — use the signed linearised measure with
  `pairedT`/`rmanova` instead.
- **The TurtleWave event bridge has no end-to-end check.** The `db_to_group_table.py` →
  `anova1` route is not covered by an end-to-end recovery check: the only such check ever built
  used channel indices that did not line up with the analysed 178-channel subset, so it never
  demonstrated anything. The `anova1` engine underneath it is checked (against MATLAB's
  `anova1`), but the *bridge* — that the exported table's channel columns land where you think
  they do on the montage — is not. Verify the column ordering yourself against
  `hdeeg_scalpchannels` before trusting an event-bridge topography.
- **Source coordinates are PLACEHOLDER**: `source2447_coords.mat` is a deterministic graph-
  Laplacian layout, **not real MNI space**. The permutation statistics, significance flags and
  cluster membership are exact and were checked directly, but the `X/Y/Z` columns in
  `_sigvoxels.csv` are placeholder until the real GeoSource MNI export (2447 rows: `label,X,Y,Z`
  in graph-node order) replaces the asset via `build_source2447_coords`. Report **regional**, not
  per-voxel, results — sLORETA leakage / low spatial resolution and volumetric adjacency that
  bridges sulci and hemispheres make single-voxel claims unsafe.
