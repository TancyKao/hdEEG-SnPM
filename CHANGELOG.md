# Changelog

All notable changes to **hdEEG-SnPM** are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] — 2026-08-08

First public release. Everything below is the content of 1.0.0; the entries are grouped by what
they do rather than by the development commit they arrived in, and the dates mark when the work
landed during private development.

### The toolbox

Permutation-based (TFCE and cluster-extent) statistics on high-density EEG topographies, with
channel-adjacency structure, following Maris & Oostenveld (2007) and Nichols & Holmes (2001).
It ships the core statistics, the GLM and LMM tiers, the App Designer GUI, the reporting stack,
the montage assets, the test suites and the documentation. Vendored EEGLAB, regenerable
`test_data/` outputs, reference PDFs and raw EEG are deliberately not tracked.

### Added

- **Circular (phase) tier** (2026-08-03). New `core_snpm_circ.m`, routed by an early guard in
  the same way as the GLM and LMM tiers, with eleven `snpm_circ_*` helpers. Three comparisons:
  `circ_phase_group` (covariate-adjusted Hotelling T² on the cos/sin embedding),
  `circ_phase_group_u2` (Watson's U²) and `circ_corrAngLinear` (circular–linear correlation).
  The angle matrix is never permuted — Freedman–Lane permutes the response through one index
  applied to the stacked `[cos, sin]` matrix, so cos and sin cannot desynchronise — and the
  estimability gate is computed once from the observed labels and held fixed, because
  recomputing it per permutation would let the analysed channel set vary across the null. There
  is deliberately no global omnibus for this tier: arithmetic channel averaging is meaningless
  for angles.
- **CircStat vendored** (2026-08-03). A partial, unmodified copy of CircStat (Berens 2009) —
  eleven functions plus the upstream licence and readme, with `README_VENDORED.md` recording
  source, version and vendoring date — under `dependencies/circstat/`, in its own subfolder to
  avoid `Contents.m` / `license.txt` / `readme.txt` collisions with the repo root. This closes
  the long-standing gap where the circular analyses were gated behind a toolbox that was never
  present, so their checks silently reported SKIP. `circ_wwtest.m` is included for byte-identity
  with upstream but is **not** used; the Wheeler–Watson path it implements is retired.
- **Source-level (2447-voxel) ingest path** (2026-08-03) for the t-test and LMM tiers, taking a
  subjects × voxels cortical table as input. Ingest only: no inverse solution is computed here,
  and the GLM presets reject source input outright. The voxel coordinates in
  `source2447_coords.mat` are **placeholder** (a graph embedding, not anatomy); they affect only
  the X/Y/Z columns of the emitted voxel table, and every statistic, cluster membership and
  significance flag is exact and independent of them.
- **Diátaxis documentation** (2026-07-16) under `docs/` — tutorials, how-to guides, reference
  and explanation, with an index README and `docs/archive/` for frozen project history. Adds a
  getting-started tutorial, a GUI how-to, a comparison-key reference, and the
  choosing-an-analysis and interpreting-cluster-results explanations. Documents the
  `(b+1)/(N+1)` p-value convention, exact-versus-Monte-Carlo sampling, the repeated-measures
  guard and the `effectiveN` output. Documents the two-program pipeline with the upstream
  detector [TurtleWave-hdEEG](https://github.com/TancyKao/TurtleWave-hdEEG) — the
  `db_to_group_table.py` hand-off, which analysis each event measure belongs in, and the
  pre-v4.0 180° preferred-phase inversion — while stating that TurtleWave is a supported source
  and not a dependency.
- **Test suites that had never executed before 1.0.0**: `test_perm_invariant`,
  `test_lmm_invariant`, `test_global_common_channels`, `test_circ_snpm`, `test_circ_stats` and
  `test_cluster_enhancement_identity`, plus `test_source_snpm` for the source path.
- `test_all('skip', {...})`, which reports the requested skips as SKIP and ends in
  `VERIFY: INCOMPLETE` — never `ALL GREEN`.

### Changed

- **Missing data: the complete-column rule** (2026-08-03). Correlation now tests a channel only
  if it is finite across every analysed subject. Pairwise deletion inside a permutation loop is
  an inference bug: the permutation re-pairs *x_i* with *y_j*, so the complete-pair set is a
  function of the permutation and the observed statistic is scored against a null built at a
  different *n*. Measured family-wise error at a nominal 0.05 before the fix: **0.288** (TFCE,
  40% disjoint missingness) and **0.168** with zeros-imputed covariates.
  The legacy t-test maps are deliberately **not** complete-column masked: what is invariant
  there is the analysed *subject set* per channel, because missingness travels with the row —
  not the effective *n*, which for `unpairedT` is not invariant at all (the n1/n2 split moves
  the standard error by up to 34%). That costs sensitivity, not validity; measured family-wise
  error stayed 0.045–0.055 across 20–40% missingness, including group-confounded missingness.
- **Degeneracy guard under relabelling** (2026-08-03) for the t-test maps: at least 2 usable
  pairs (paired/one-sample), at least `max(n1,n2)+2` (unpaired). MATLAB's `var` of a scalar is
  0, not NaN, so a group reduced to one subject returns a finite inflated t that sets the
  max-null for the whole montage — measured planted-cluster power **0.038 unguarded versus
  0.868 guarded** on EGI-178, family-wise error nominal in both. A single usable pair returned
  `-Inf` and crashed `ClusterEnhancement` with `MATLAB:pmaxsize`.
- **The global omnibus now averages every analysed unit over one common channel set**
  (2026-08-03). `global_stat_test` previously averaged each subject over its own channel set;
  with a different region missing in each condition — two nights, different bad channels — the
  paired test rejected on 100% of null data. Measured type-I error **1.0000 → 0.0512** (paired)
  and **0.0975 → 0.0555** (unpaired), bitwise identical on complete data and bit-identical for
  correlation, where the upstream mask makes the rule idempotent. Applied per effect in
  `core_snpm_glm`'s `mixed2way`; `core_snpm_lmm` was already correct.
- **Permutation p-values unified to the minimum-bias `(b+1)/(N+1)` form** across all tiers —
  legacy, GLM and LMM (2026-07-16).
- Legacy paired/one-sample sign-flip null no longer freezes high-index subjects (it was
  enumerating the first *K* integer bit patterns): exact when 2^n ≤ requested, otherwise a
  random fair coin with the observed labelling included (2026-07-16).
- Legacy unpaired and correlation permutation: exact enumeration when the requested count is at
  least the group size, otherwise Monte Carlo with replacement (2026-07-16).
- GLM between-subject presets (`anova1`, `ancova`, `regression`) now error on repeated-measures
  data, and correlation runs surface per-channel effective N (2026-07-16).
- The `checkcode` gate escalates the unused-value class it used to discard — the class that
  would have caught the partial-Spearman defect below at authoring time (2026-08-03).

### Removed

- `circ_wheeler_watson_Test` and `circ_WatsonsU2Test`, together with the
  `snpm_single_threshold_with_TFCE_circ` and `snpm_cluster_analysis_circ` engines. Both keys now
  raise a named error pointing at their replacements in the circular tier (2026-08-03).

### Fixed

- Left-tailed paired and unpaired **global** p-values were computed with the wrong tail and
  returned the complement; the partial-correlation global p now uses df = n − 2 − k
  (2026-07-16).
- An infinite hang in the legacy unpaired/correlation permutation path, and a phantom reference
  to a non-existent `snpm_enumerate_combinations` (2026-07-16).
- **Partial Spearman** ranked the data, computed a `stat_comparison`, then never passed it to
  the engines, so it silently reported Spearman of rank residuals. It now matches
  `partialcorr(..., 'Type', 'Spearman')` to 3.7e-13 (2026-08-03).
- `per_channel_n` was read but never assigned, which threw on every legacy analysis, and
  `snpm_exclusion_spatial_profile` was called but did not exist (2026-08-03).
- The channel mask is unified on `isfinite`, so a `-Inf` from `log10(0)` is reported rather than
  silently dropped (2026-08-03).

### Verification

`matlab -batch "test_all"` returns `VERIFY: ALL GREEN` — 16 checks, 0 skipped, 1 known xfail
(the event-bridge `synthetic_events` `planted_truth` channel-index mismatch).

[1.0.0]: https://github.com/TancyKao/hdEEG-SnPM/releases/tag/v1.0.0
