# Changelog

All notable changes to **hdEEG-SnPM** are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Behavioural single-column Data 2 for the correlations** (`correlationP`/`correlationS`):
  Data 2 may now be a `Subject` column plus exactly ONE behavioural measure column (e.g. a
  reaction time), broadcast across Data 1's channels — no more copying the value into a fake
  subjects×channels table. The legacy two-montage-file format is unchanged and bit-identical.
  Format is classified from the column names (`dependencies/snpm_detect_behavior_table.m`,
  shared by engine and GUI) into three recorded matching modes
  (`results_struct.data2_matched`): `label` (verified montage match), `behavior`
  (single-column broadcast), `positional` (unrecognised column names of Data 1's width, matched
  by position with a warning + report caveat; a differing width errors
  `snpm:correlation:widthMismatch` naming the columns). A non-numeric or constant measure is a
  named error, and a subject with a missing value is dropped whole with a warning
  (`snpm:behaviorNaNRows`) and a count in
  `results_struct.subject_matching.behavior_nan_removed`. Only the FIRST `subject`-containing
  column is the ID, so a `SubjectiveSleepiness` column is a valid measure. New example file
  `example_data/corr_behavior_single.csv`. **Disclosure:** prefer `datatype='absolute'` for
  EEG-vs-behaviour runs — under `normalize` the z-scored EEG rows sum to zero, so the broadcast
  r-map satisfies Σ r·sd = 0 and a global association cannot appear uniformly signed; for any
  correlation under `normalize` the whole-head global is analytically zero and is flagged
  (`results_struct.global_not_interpretable`, `snpm:globalNotInterpretable`).
  Behaviour mode + `normalize` also warns (`snpm:behavior:normalizeBroadcast`) and surfaces the
  caveat as `results_struct.normalize_broadcast_caveat`, a `CAVEAT:` results-text line and an
  amber **Data Type caveat** banner in the HTML report. Only the FIRST `subject`-containing
  column is stripped as the ID; any other non-numeric column that survives errors
  `snpm:nonNumericDataColumn` naming it.

- **Recording hook for permutation indices** (`dependencies/snpm_perm_tap.m`) and an export
  script that uses it (`scripts/export_parity_fixtures.m`). The hook is off by default. When
  it is switched on it copies out the permutation indices each analysis has already drawn
  (sign-flip patterns, group relabellings, row orders and the Freedman–Lane residual index),
  plus the matrices handed to the legacy engines and the sorted cluster-size null. It draws no
  random number and changes no result, whether on or off, and nothing in the toolbox reads the
  log back. The export script seeds the generator, runs `core_snpm_analysis` on the files in
  `example_data/` (and on a few inputs it builds from them) for the t-tests, the correlations
  and the five General Linear Model presets, and saves one `.mat` file per case holding the
  inputs, the indices and the results, with a `README.md` that lists every field. The files go
  to `test_fixtures/parity/`, which is not tracked. The purpose is to check a separate
  implementation of the same analyses against MATLAB using identical permutations rather than
  only an identical seed. For the General Linear Model presets the export script replays the
  recorded indices and confirms that they reproduce the pipeline's own permutation null before
  it saves a file. The circular and mixed-model draw sites are also recorded by the hook, but
  the export script does not cover those tiers yet, so no reference files exist for them and
  their recorded indices have not been checked by replay.

- **The GUI checks a `mixed2way` table for balance before Run** (`SnPMAnalysisGui.m`). Once the
  file and the Subject and Condition roles are set, the readiness list applies the same rule as
  the engine: Run stays disabled and the list names the subjects without exactly one row per
  condition, the subjects in more than one group, and any empty Condition cells. A file the GUI
  cannot read does not block Run, so the engine's own error still reports it.

- **`example_data/`**: ready-made synthetic input, tracked in the repository: one CSV (or file
  pair) per analysis, each with a 17-channel neighbour cluster planted near E129 on the EGI 256
  montage, plus a `README.md` mapping every file to its analysis, its column roles and its GUI
  settings. This replaces the generate-it-first step: the data is simply there.

### Fixed

- **Mixed model: an infinite value in a model column is dropped like a missing one**
  (`mixedmodel`; `dependencies/snpm_lmm_drop_missing_meta.m`). `fitlme` refuses a fit whose
  response contains `Inf`, and the per-channel catch then left every channel NaN; a trial with
  ±`Inf` in a numeric model column (dependent variable, a formula column, the subject column, the
  group column under `group_label`) is now dropped up front with the other missing trials, with the
  same `snpm:lmmMetaNaNRows` warning and `droppedRows` reporting.
- **Mixed model: an effect that names nothing in the model is refused** (`mixedmodel`;
  `core_snpm_lmm.m`, new `dependencies/snpm_lmm_check_effect.m`). An `lmm_effect` that named no
  coefficient (continuous) or term (factor) of the fixed part produced an all-NaN map without an
  error. It is now refused from the metadata, before the channel setup and any fit, with
  `core_snpm_lmm:unknownEffect`, e.g. "Effect 'POWERX' is not a coefficient of the fixed part
  'POWER + cov' (coefficients: (Intercept), cov, POWER)" or "Effect 'POWER:time' is not a term of
  the fixed part 'POWER*group' (terms: group, POWER, group:POWER)". A factor term written in
  another variable order (`POWER:group`) is now resolved to the name `fitlme` uses
  (`group:POWER`) instead of mapping NaN.
- **Mixed model: the `group_label` effect check no longer accepts any `<group>_…` name**
  (`core_snpm_lmm:effectNotPermuted`). A token now names the group only when it is the group column
  or `<group>_<level>` for a level present in the data (`group_B`, also inside an interaction such
  as `group_B:POWER`); an unrelated column such as `group_size` was accepted before and is now
  refused with the same message. A `<group>_<x>` for a level not in the data (`group_C`) is
  refused as an unknown effect, whose message lists the real coefficients. The names are taken
  from one ordinary least-squares fit of the fixed part to the metadata prepared exactly as the
  per-channel fits prepare it (new `dependencies/snpm_lmm_categorize.m`, shared with
  `snpm_lmm_fit`), so they carry `fitlme`'s reference levels: a text column the fit does not
  convert keeps its first-appearing level as the reference, and a logical group column is named
  `group_true`.
- **Mixed model: a within-subject shuffle on an effect constant within subject is refused**
  (`mixedmodel`; `core_snpm_lmm.m`, new `dependencies/snpm_lmm_effect_within_share.m`). With
  `lmm_perm = 'within_subject'` and an effect whose variables are all constant within every
  subject (the group, `group_B`, a subject-level covariate), the shuffle leaves each subject's mean
  response in place and gives no null: on simulated data with no group effect and POWER predicting
  the response, the false-positive rate at a nominal 0.05 was 0.35. Such a run is now refused
  before any fit with `core_snpm_lmm:effectConstantWithinSubject`. Effects involving POWER or a
  within-subject variable (`group:time`) run as before.
- **Mixed model, `group_label` with a metadata response** (`core_snpm_lmm.m`). A metadata
  response constant within every subject is refused (`core_snpm_lmm:dvConstantWithinSubject`,
  new message) instead of producing a map of perfect fits. A channel with one value on every trial
  is now excluded and listed ("same value in every subject"); before, `fitlme` refused its fit
  (rank-deficient fixed design) at every labelling and the channel came back NaN without being
  listed. With a random POWER slope by subject, channels constant within each subject are
  excluded ("predictor constant within each subject: its random slope by subject cannot be
  estimated"). Without such a slope those channels are kept, as before: under subject
  relabelling they are a valid subject-level covariate.
- **Mixed model: the zero-variance refit trigger is unit-free** (`dependencies/snpm_lmm_fitlme.m`).
  Each random-effect standard deviation is now multiplied by the root mean square of the design
  column it multiplies before it is compared with 0.001 times the residual SD, so the decision no
  longer depends on the units of a random-slope predictor (a slope on x and on 1000·x used to
  decide differently). An intercept has scale 1: intercept-only models are bit-identical (1600
  fits over the `lmm_boundary_stall_rate` designs compared with the previous code, 0 differ), and
  `scripts/lmm_boundary_stall_rate.m` prints the same table. The refit is still kept only on a
  log-likelihood gain above 1e-6 (absolute, deliberately: a log-likelihood difference does not
  depend on units).
- **Circular–linear exclusion reason named groups** (`circ_corrAngLinear`; `core_snpm_circ.m`).
  The analysis is paired by subject, but its excluded channels carried the group analyses' reason
  "not estimable: fewer than 8 subjects per group, incomplete column under a per-channel
  covariate, or log(count) collinear with group". It now reads "not estimable: fewer than 8
  subjects with both an angle and a measure value". The two group analyses are unchanged.
- **Legacy permutation p counted the observed labelling twice** (`pairedT`, `unpairedT`,
  `onesampleT`, `correlationP`, `correlationS`, including correlation with covariates;
  `snpm_single_threshold_with_TFCE.m`, `snpm_cluster_analysis.m`). The observed labelling is
  already one of the N labellings these engines evaluate (row 1 of a Monte-Carlo draw, or a
  member of an exact enumeration) and is counted in b, yet the p added 1 to numerator and
  denominator again. The max-statistic, TFCE and cluster p are now b/N, the exact permutation p
  under full enumeration and the Phipson & Smyth estimator under Monte-Carlo. Every legacy
  corrected p falls by less than 1/(N + 1): negligible at 10000 permutations, decisive at small
  n with exact enumeration (paired, 5 subjects, one-sided: smallest p 1/32 = 0.031, was
  2/33 = 0.061, so significance was unreachable). The GLM, circular and mixed-model tiers draw
  labellings without the observed one and keep (b + 1)/(N + 1). Also corrects the documented
  3-versus-3 minimum p (1/20 one-sided, 2/20 two-sided, not "about 0.048").
- **Mixed model: settings the permutation test cannot handle are refused up front**
  (`mixedmodel`; `core_snpm_lmm.m`). Affects `mixedmodel` runs with any of three settings, which
  used to run and return a map that could not be read as a test. (1) `tail` `left` or `right`:
  for `left` the null was built from signed maxima, so nothing could ever be significant; an F
  has no side; and clusters always used the two-sided model p. Now an error,
  `core_snpm_lmm:oneSidedTail`; `mixedmodel` takes `tail = 'both'` only. (2) The `group_label`
  permutation with an effect that does not involve the group column (the effect's variables,
  split on `:` or `*`, must include `lmm_group` or a coefficient of it such as `group_B`): relabelling subjects then leaves every permuted
  map equal to the observed one, so every corrected p was 1. Now an error,
  `core_snpm_lmm:effectNotPermuted`. (3) The `within_subject` permutation with a dependent
  variable that is not a metadata column (for example `POWER`): it used to fail only inside the
  permutation loop. Now an error before anything is read, `core_snpm_lmm:dvNotMetadata`. Runs
  with valid settings are unchanged.

- **Mixed model: channels constant within each subject, and fits with no residual variation**
  (`mixedmodel`; `core_snpm_lmm.m`, `dependencies/snpm_lmm_fit.m`, new
  `dependencies/snpm_lmm_within_ratio.m` and `dependencies/snpm_lmm_degenerate.m`). Affects
  `mixedmodel` runs with a channel whose value is constant, or almost constant, within each
  subject (including a flat channel), and runs whose dependent variable is constant within each
  subject. Before: when the channel was the response (e.g. `POWER ~ group + (1|Subject)`),
  `fitlme` did not fail on such a channel but returned a finite F from rounding noise (20.6,
  p = 2e-5, on a flat channel in one example), the same in every permutation, so it could raise
  the maximum-statistic null for the whole map; above a subject/residual SD ratio of about 1e4
  its fit is also off the optimum (35% off in F at 1e5). When the channel was the predictor under
  `within_subject`, a channel constant within subject has a between-subject slope only, which
  the within-subject shuffle never moves: its t sat unchanged in every permuted map, and in a
  test example it raised a planted effect's corrected p from 0.016 to 1.0. Now:
  - before any fit, once and the same for every permutation, a channel whose within-subject
    share of variation (within-subject sum of squares over total, grouped by subject only) is at
    most 1e-6 is excluded with a reason, through the usual excluded-channel list, Excel sheet and
    HTML banner: when the channel is the response, "channel value constant within each subject:
    no within-subject variation for the mixed model to estimate; analyse subject means with
    anova1/unpairedT"; when it is the `within_subject` predictor, "predictor constant within each
    subject: its effect is between-subject only, which the within-subject shuffle cannot test (use
    group_label or regression on subject means)";
  - a `within_subject` run whose dependent variable has a within-subject share of at most 1e-6 is
    refused (`core_snpm_lmm:dvConstantWithinSubject`): the shuffle has nothing to permute;
  - as a backstop at every fit, a fit with no usable residual variation (the response fitted
    exactly by the fixed effects, or by the fixed and random effects when they leave some residual
    degrees of freedom, to 1e-20 of its uncentred
    sum of squares, or, when the only random effect is a subject intercept, a subject/residual SD
    ratio of at least 1e3 estimated by the method of moments from those sums) has no statistic
    (stat NaN, Wald 0, p NaN); an observed channel caught this way is listed as "perfect fit,
    statistic undefined";
  - if the exclusions leave no channel, the run stops with `core_snpm_lmm:noEvaluableChannels`.
  Channels with ordinary within-subject variation are unchanged bit for bit.

- **Mixed model: a fit that stopped at zero subject variance** (`mixedmodel`;
  `dependencies/snpm_lmm_fit.m`, new `dependencies/snpm_lmm_fitlme.m`). Affects every
  `mixedmodel` run in which any channel fit, observed or permuted, came back with a random-effect
  standard deviation (SD) of 0, usually the subject SD. `fitlme` estimates each random-effect SD
  relative to the residual SD, and its likelihood is symmetric about 0, so 0 is always a flat
  point; the default optimiser sometimes stopped there when a positive SD fitted the data clearly
  better. The statistic was then that of a model with no subject effect. On the 3-group example in
  `test_lmm_snpm` (Tier 5) the group F was 20.75 where the best fit gives 8.59, with a
  log-likelihood about 6 higher. For a between-group effect this makes the statistic too large:
  in the observed map that is anti-conservative at that channel, in a permuted map it inflates
  the null and makes the test conservative. For a within-subject effect the statistic is still
  wrong but can move either way. The stall is common when the subject SD is small next to the
  trial-to-trial SD: on simulated data with no effect and a subject SD 0.4 times the residual SD,
  142 of 200 fits stopped at 0 although the best fit was elsewhere
  (`scripts/lmm_boundary_stall_rate.m` reproduces the counts). Now, when any random-effect SD comes back at 0 (at most 0.001
  times the residual SD), the model is refitted once from a random start drawn from a private
  random stream, so the permutation draws are unchanged, and the refit is kept only if its
  likelihood is higher; a refit that fails is discarded. Fits away from 0 are returned exactly as
  before, so a run with no fit at 0 gives identical results. Affected runs change at the affected
  channels and permutations.

- **Mixed model: trials with missing model variables** (`mixedmodel`; `core_snpm_lmm.m`, new
  `dependencies/snpm_lmm_drop_missing_meta.m`). Affects `mixedmodel` runs whose table has a blank
  or `NaN` in the dependent variable (when it is a metadata column; a missing channel value is
  still handled by the complete-column rule), a column named in the fixed or random formula, the
  subject column, or (for the `group_label` permutation) the group column. `fitlme` used to drop those
  trials silently at every channel, so the reported trial and subject counts included rows that
  were never fitted. Worse, under `group_label` each subject's group was taken from its first
  trial, so a group label missing on a later trial, or a group that changed within a subject,
  meant even the unpermuted relabelling did not reproduce the data, and the null was built on a
  different design from the observed statistic (an invalid test). Under `within_subject` a
  missing dependent variable was shuffled onto other trials, so the set of fitted trials
  changed from one permutation to the next. Now those trials are dropped once, before
  anything else, with a warning (`snpm:lmmMetaNaNRows`) giving the count per column and any
  subject that lost all its trials; the counts are recorded in `results_struct.lmm.dropped_rows`
  and shown in the results text, the HTML model card and an Excel `droppedRows` sheet, and the
  reported trial and subject counts are those analysed. Columns not in the model never drop a
  row. The design is then checked: fewer than 3 subjects (`core_snpm_lmm:tooFewSubjects`), a group
  that varies within a subject under `group_label` (`core_snpm_lmm:groupVariesWithinSubject`, also
  refused when nothing is missing), or fewer than 2 groups or a group with fewer than 2 subjects
  under `group_label` (`core_snpm_lmm:groupLevelTooFewSubjects`) is an error, as is a subject or
  (under `group_label`) group column not listed in `lmm_meta_cols`
  (`core_snpm_lmm:noSubjectColumn`, `core_snpm_lmm:noGroupColumn`); a subject left with one trial under
  `within_subject` is a warning (`snpm:lmmSingleTrialSubject`). Dropping assumes the missingness
  is unrelated to the effect tested (for example, awakenings with no report). Runs with complete
  model variables are unchanged; a `group_label` run whose group varied within a subject now
  stops with an error instead of producing an invalid map.

- **Circular–linear correlation with a missing measure value** (`circ_corrAngLinear`;
  `core_snpm_circ.m`). Affects runs where the measure file has a blank or non-finite value for
  any subject who also has angles. Such a subject used to stay in the analysis: its missing
  value travelled with each permutation, so every permuted map dropped a *different* subject,
  and the reported per-channel n counted the subject as if it had been analysed. The test was
  not anti-conservative, but it was conservative at small n and could stop with a
  `core_snpm:circGateDrift` error at random. Now a subject without a finite measure is dropped
  whole before anything is computed, with a warning (`core_snpm:circMeasureNaNRows`) naming
  them, and is listed in `results_struct.circ.dropped_subjects`, the results text and a report
  notice; fewer than 8 remaining subjects is an error (`core_snpm:circTooFewSubjects`). Missing
  angles are still deleted per channel, now on a fixed set of subjects. The reported n is
  correct, and results of affected runs can differ from before (small-n runs were
  conservative). Runs with a complete measure are unchanged.

- **Correlation with covariates now uses a Freedman–Lane permutation** (`correlationP` /
  `correlationS` with `use_covariates`; `core_snpm_analysis.m`, `snpm_single_threshold_with_TFCE.m`,
  `snpm_cluster_analysis.m`, new `dependencies/snpm_fl_corr.m`). The old scheme residualized
  both measures on `[1, covariates]` once and then permuted the EEG residuals against the fixed
  behavioural residuals. A permuted residual is no longer orthogonal to the covariates, so the
  null was too narrow whenever the covariates explain much of the data: simulated family-wise
  error at nominal 0.05 was 0.856 (max statistic) / 0.503 (TFCE) / 0.235 (cluster) with 10
  subjects and 4 covariates, and 0.140 / 0.118 / 0.080 with 20 subjects and 3 covariates. Each
  permuted map is now `corr(R_Z * x_res(order,:), y_res)`: the permuted residuals are projected
  off the covariates again before correlating (Freedman & Lane 1983; Winkler et al. 2014). The
  observed map is unchanged (it still equals `partialcorr`), the permutation p-values change.
  The test is approximate, not exact, and is reported that way in the results text and the
  HTML report ("Freedman–Lane, approximate"). Simulated family-wise error at nominal 0.05 is
  about 0.04–0.06 for the per-channel max statistic and at or below nominal for TFCE and
  cluster inference. With fewer than 5 error degrees of freedom (n − 2 − k < 5) the max
  statistic may be mildly liberal (0.056 at nominal 0.05, 3 error df) while TFCE and cluster
  inference are conservative, and the run warns (`snpm:covariateConservative`). With no error
  degrees of freedom left the run stops (`snpm:covariateNoErrorDF`). The covariate count used
  in every degrees-of-freedom statement is now the rank of the covariate design
  (`k = rank([1, covariates]) − 1`), so a duplicated covariate no longer costs a degree of
  freedom (`snpm:covariateRankDeficient` warns). New `results_struct.covariates` fields:
  `k_eff`, `error_df`, `method`, `spearman_caveat`, `conservative_note`. Runs without
  covariates are bit-identical.
- **Uncorrected per-channel p of a partial correlation used the wrong degrees of freedom**.
  `p.real` came from `corr()` on the residuals, i.e. df = n − 2 as if nothing had been fitted;
  it now uses df = n − 2 − k and equals the `partialcorr` p (the global test already did this).
  Without covariates `p.real` is unchanged.
- **Partial Spearman caveat documented.** `correlationS` with covariates ranks the data and
  the covariates and then partials the covariates out linearly, which removes a covariate only
  as far as its relation to the ranks is linear. The null it tests is that the rank residuals
  have no linear association once the ranked covariates are removed linearly, not that the two
  measures are unrelated given the covariates. A covariate whose effect is linear on the raw
  scale is not fully removed (simulated max-statistic family-wise error 0.10–0.13 with 3 such
  covariates, under any permutation scheme). The results text and the report now say so and
  point to `correlationP` or the `regression` preset when covariate control matters.

- **Event bridge: one detection run per band, and density over analysed time**
  (`db_to_group_table.py`). The bridge pooled every detection run of a band, so a per-stage
  run and a later joint run of the same night counted most events twice. It now reads only the
  runs whose stages lie inside the new `--stages` option and stops with an error naming the
  database and runs when two runs cover the same stage or two detection methods share a band.
  Density is now events per minute of TurtleWave's artefact-free analysed time for the stages
  the run searched, the same rule as TurtleWave's own density code, instead of raw hypnogram
  minutes. Databases older than TurtleWave 4.3 have no analysed time; they are refused unless
  `--allow-hypnogram-denominator` is given, and a table that mixes denominator types or stage
  sets is refused. New `--bands` option; a new `recordings.csv` records, per subject and band,
  the denominator used, its minutes, the stages and the interpolated channels. The database is
  opened read-only. **Who is affected:** anyone who ran the bridge on a database with more than
  one run of a band (counts and densities were inflated), and every density table, whose
  values change with the denominator. **What to re-run:** the bridge, then any analysis of
  density. Amplitude and duration change only where runs were pooled.
- **Per-channel Excel sheets are headed by channel labels** (`dependencies/func_genSnpmTable.m`).
  The `pVal_or_pCritVal` sheet always came out with headers `Var1..VarN`, because a length check
  compared a one-cell array with the channel count. `statVal`, `pVal_or_pCritVal`, `pValTFCE` and
  `effectiveN` now carry the channel labels; a width mismatch warns
  (`func_genSnpmTable:labelMismatch`) and keeps generic headers. No value changes.

- **`ClusterEnhancement` returned the wrong number of values for a map with no channel above
  the first threshold** (`dependencies/ClusterEnhancement.m`). Threshold-free cluster
  enhancement (TFCE) adds up a contribution per channel over a ladder of thresholds that
  starts at the step size, 0.1 by default. When every channel's |statistic| was below that
  first step the ladder had a single row, MATLAB summed along the row instead of down the
  column, and the function returned one number (or, when the map had missing channels, a
  vector that stopped at the last missing channel) instead of one value per channel. It now
  always returns one value per channel: 0 where the channel is finite and NaN where it is
  missing. **Who is affected:** only maps in which every channel's |statistic| is below 0.1
  (or below the step, for the Watson U² path, which uses 0.005). Any map with at least one
  channel at or above the step is unchanged; 204 of the 206 stored reference maps are
  bit-identical before and after, and the two that changed are the all-zero map and the
  all-zero map with missing channels. **What it looked like before**, as far as the calling
  code shows, depends on the analysis. In the t-tests and the Pearson correlation, a permuted
  map of this kind without missing channels came back as a single zero that MATLAB copied
  across the row, which is the same as the correct answer, so the run completed with the same
  numbers; if such a map also had a missing channel that was not the last channel, the short
  vector did not fit the row and the run stopped with a size-mismatch error. In the Spearman
  correlation an older length patch appended one NaN to the single zero, the two values did
  not fit the row, and the run stopped with a size-mismatch error whenever any permuted map
  was entirely below the step, with or without missing channels. In the General Linear Model
  presets and the circular tier, a permuted map of this kind contributed a maximum of 0 to
  the null both before and after, so the run completed with the same numbers; an observed map
  entirely below the step stopped the run with an index error. The mixed-model path already
  padded the short result back to full length. **What this means for results:** in every run
  that completed, the permutation null, the corrected p-values and the clusters are the same
  before and after the fix. The defect could only stop a run, in the three situations above.
  Runs that stopped this way now complete, and some of them may contain significant channels.

- **Pairwise post-hoc comparisons are now each a separate two-level test** (`core_snpm_glm.m`,
  `export_glm_report.m`, new `dependencies/snpm_glm_pair.m`). After an effect with three or
  more levels, the comparison "level a against level b" is now tested on the rows of those two
  levels only, with the preset's design rebuilt on those rows and its own permutations. Before,
  each pair was tested as a contrast inside the model for all levels, and the effects of the
  other levels stayed in the residuals that were permuted. **What the old scheme did**, measured
  by simulation with a and b equal and a nominal rate of 0.05: when a third level was raised
  across the whole head, the cluster test for "a against b" rejected 0.17 to 0.59 of the time
  and the threshold-free cluster enhancement (TFCE) test about 0.02 of the time; when the third
  level was only less variable than the other two, both tests rejected 0.19 to 0.33 of the
  time. With the new scheme all of these rates are between 0.036 and 0.060. **Who is affected:**
  anyone who has reported a post-hoc map from `anova1`, `ancova` or `rmanova` with three or
  more levels, from the group, condition or simple-effect post-hocs of `mixed2way`, or from
  `export_glm_report.m`. Analyses with two levels have no post-hoc maps and are not affected.
  **What to re-run:** every such analysis whose post-hoc maps you rely on. The omnibus F maps,
  their corrected p-values and their clusters are unchanged; the post-hoc maps will differ from
  earlier runs. **What it costs:** the error term of a pair no longer pools the other levels. In
  the simulated designs this lowered TFCE power by about 0.05 for between-subject comparisons
  and about 0.15 for within-subject comparisons. **Other consequences:** in an `ancova`
  post-hoc the covariate slopes are estimated from the two groups compared, not from all
  groups. In a within-subject post-hoc a subject who has only one of the two levels is left out
  of that comparison, with a warning that names the subject. A pair that cannot be estimated is
  skipped with a warning, so a report can show fewer post-hoc maps than there are pairs.
  Post-hoc maps remain corrected across channels and are not corrected across pairs.

- **`mixed2way` now requires a balanced table and uses one channel set for the whole report**
  (`core_snpm_glm.m`). Every subject must have exactly one row per condition and belong to one
  group. If not, the run stops before any statistic is computed with the error
  `core_snpm_glm:mixed2wayUnbalanced`, which names the subjects and what each is missing or
  has twice, and points to `mixedmodel` for unbalanced data. Every effect in the report (the
  group, condition and interaction effects, the simple effects, every post-hoc and the three
  whole-head tests) now uses the same channels: those with a value in every row of the table.
  Before, the subject means for the group effect were taken over whatever conditions a subject
  had. **What the old behaviour did**, measured on data with no group effect: with a condition
  effect present and 7 of 20 subjects missing a condition, the group effect rejected 0.61 of
  the time at a nominal 0.05. **Who is affected:** `mixed2way` runs on a table in which some
  subject was missing a condition or had a condition twice, and `mixed2way` runs on a table
  with missing values at some channels. **What to re-run:** for the first kind, the old run
  cannot be repeated as it was; remove or complete the subjects the error names and run again,
  or use `mixedmodel`. For the second kind, re-run: the group effect and the simple effects are
  now tested only on channels that are complete in every row of the table, which can be fewer
  than before, so their maps, the group whole-head test and the excluded-channel lists can
  differ. A balanced table with no missing values gives the same group and interaction F maps
  as before (their post-hoc maps change under the item above, and the condition effect changes
  under the item below).

- **The `mixed2way` condition main effect is now tested in the full model**
  (`dependencies/snpm_glm_design.m`). The condition effect is tested with the
  group-by-condition columns in the model as nuisance terms, with group coded to sum to zero.
  Before, it was tested in the reduced model `[1, Ds, Dc]` (intercept, subject and condition
  columns only), which left any real interaction in the residuals that were permuted and could
  hide the condition effect completely: in simulation with a crossover interaction a planted
  condition effect was detected 0.00 of the time under the old model and 0.93 of the time under
  the new one. The condition F now equals the F of a standard mixed analysis of variance; on a
  check table it reproduces MATLAB's `ranova`, F(2,40) = 30.9004. Its error degrees of freedom
  are (subjects − groups) × (conditions − 1). With unequal group sizes it is the condition
  effect averaged over groups with equal weight per group, which is the default of standard
  software. **Who is affected:** every `mixed2way` run. **What to re-run:** any `mixed2way`
  analysis whose condition main effect, or condition post-hoc maps, you rely on. The condition
  statistic, its degrees of freedom, its corrected p-values and its whole-head test will differ
  from earlier runs. The interaction effect and the simple effects are not changed by this
  item.

- **Channels with no variation no longer stop the run, and are reported instead of silently
  tested** (`core_snpm_analysis.m`, `core_snpm_glm.m`, `dependencies/snpm_single_threshold_with_TFCE.m`,
  `dependencies/snpm_cluster_analysis.m`, `dependencies/snpm_glm_stat.m`,
  `dependencies/snpm_glm_fl_context.m`, new `dependencies/snpm_zero_ss.m`,
  `dependencies/snpm_ttest_novar.m`, `dependencies/snpm_excluded_add.m`). This fixes the crash
  found and deferred in August 2026. One rule now applies to every tier. A channel with nothing
  to permute (every paired difference zero, the same value in every subject, a correlation
  column that is constant, or no variation left after the General Linear Model's nuisance model)
  is excluded once, before the permutations, and listed with a reason. At any labelling whose
  error sum of squares is zero (at most 10⁻²⁰ of the sum of squared values, a relative standard
  deviation of 10⁻¹⁰) the statistic is defined as 0, for the observed labelling and every
  permuted one alike, which keeps the permutation test exact. A channel whose observed data
  separate perfectly is shown as NaN and listed as "perfect separation, t undefined" (t-tests)
  or "perfect fit, statistic undefined" (General Linear Model presets), and is left untested.
  In the General Linear Model presets this is decided per test, so a post-hoc pair can exclude a
  channel its omnibus tested. **What happened before:** when every paired difference was the
  same non-zero value, or two groups were each constant at different values, `ttest`/`ttest2`
  returned an infinite or 10¹⁶-scale t and threshold-free cluster enhancement (TFCE) stopped
  with an out-of-memory or array-size error naming no channel. This also happened when only one
  permuted labelling was degenerate, for example a channel whose paired differences were all +1
  or −1; on sparse count data with 6 subjects per group, up to 36% of simulated data sets hit
  it. A channel constant at a non-zero value was kept and tested with t = 0 or r = 0 without
  being mentioned. Family-wise error under the rule, on null count-like data with about 78%
  zeros, was 0.040 to 0.064 for TFCE and 0.011 to 0.039 for the cluster test at a nominal 0.05.
  **Who is affected:** runs that stopped with these errors now complete; runs on data with
  constant or all-zero channels now list those channels as excluded where they used to report
  them as tested with no effect. Runs on continuous data are unchanged.

- **Threshold-free cluster enhancement no longer runs out of memory on a very large statistic**
  (`dependencies/ClusterEnhancement.m`). TFCE built a ladder of thresholds up to the largest
  value in the map, so its memory grew with that value. A design with one error degree of
  freedom (for example a repeated-measures table in which only two subjects have both
  conditions) gives statistics with very heavy tails; a permuted value of 9.3 × 10⁶ is a
  legitimate draw, and it asked for about 130 GB. Above 10 000 levels (a statistic above about
  1000 at the default step of 0.1) TFCE now computes the same integral without the ladder,
  rebuilding clusters only where some channel's level count changes and summing each run of
  levels in closed form; its memory and time no longer depend on the size of the statistic, and
  a map holding 10¹² takes about a millisecond. Below 10 000 levels the original computation
  runs unchanged, so every stored reference map is bit-identical. The two computations agree to
  within about 10⁻¹⁵ relative. **Who is affected:** runs that stopped with an out-of-memory error
  in `ClusterEnhancement` now complete. No other result changes.

- **Permutation p-values could be too small where a permuted value tied the observed one**
  (new `dependencies/snpm_count_at_least.m`, used by `snpm_single_threshold_with_TFCE.m`,
  `snpm_perm_correction.m`, `snpm_lmm_TFCE.m`, `snpm_lmm_cluster.m` and
  `scripts/export_parity_fixtures.m`). The corrected p is (b + 1)/(N + 1), where b counts the
  permuted maxima at least as large as the observed value. In the General Linear Model,
  mixed-model and circular tiers the observed map is computed from the data and an equivalent
  permuted labelling from fitted values plus permuted residuals; the two are equal
  mathematically but differ in the last digits, in either direction, so the strict comparison
  dropped about half of the true ties and made p too small. A permuted value within a relative
  10⁻¹² of the observed value now counts as a tie. **What it changed:** in one acceptance case
  (a `mixed2way` simple effect) the corrected p moved from 0.01 to 0.035. **Who is affected:** results where such a tie occurred, which
  is most likely with few distinct labellings (small repeated-measures designs). Only those p
  values go up. The legacy cluster-extent count, which compares integer sizes, is unchanged.

- **Designs with no error degrees of freedom are refused, and an undefined whole-head test is
  reported as undefined** (`core_snpm_glm.m`, `dependencies/generateAnalysisReport.m`,
  `dependencies/func_genSnpmTable.m`). When the number of rows minus the rank of the design
  matrix is below 1, a General Linear Model run stops before any statistic with
  `core_snpm_glm:noErrorDF`, "no error degrees of freedom left after fitting the design"; for
  example `anova1` with one subject per group, or `mixed2way` with one subject per group (the
  whole run is refused). Before, such a design was not detected and a test with no error term
  was run and reported (for a simple effect, as "0 significant"). A
  `mixed2way` simple effect in a group with one subject is now skipped and listed: "simple
  effect: condition within A was not run: only 1 subject in this group; a within-group test
  needs at least 2"; the old check counted rows rather than subjects and let it through. The
  whole-head global test is parametric, so when the model fits the channel average exactly its
  statistic and p are now NaN with `global_reason` "perfect fit, statistic undefined", and the
  global card is still shown with the reason, rather than a number that would read as "no
  evidence" for what is in fact a perfect fit. **Who is affected:** saturated or nearly
  saturated designs only.

### Changed

- **HTML report v2** (`dependencies/generateAnalysisReport.m`, the design in
  `templates/snpm_report_v2_master.html`): every analysis now emits one layout — a sticky
  summary strip with one chip per effect (traffic light, `Clu · TFCE · Unc` counts, global
  statistic and p), a model card, a single amber notice block (excluded channels, Data 2
  positional matching, `normalize` caveat) shown only when non-empty, a per-analysis
  descriptive block (mean topographies for the t-tests, a scatter for the correlations, the
  cell-mean grid and interaction plot for `mixed2way`), one collapsible section per effect
  with a three-way **Cluster / TFCE / Uncorrected** toggle in place of the old two-way
  TFCE / Cluster switch, cluster tables grouped per cluster with signed effects split into
  positive and negative clusters, **the cluster numbers drawn on the Cluster PNG** at each
  cluster's centroid, uncorrected tables capped at 20 rows with a "Show all" button, post-hoc
  pairwise tabs for the omnibus presets, a Global test card, and a **Results CSV** button
  opening the sibling `<base>_report.csv`. The HTML is fully self-contained (inline CSS and
  script, Print / PDF button). The spectral reports from `export_report.m` are unchanged.
  Landed first for `unpairedT`, `correlationP` and `mixed2way`; the remaining analyses follow.
- **Report images are embedded** (`dependencies/generateAnalysisReport.m`): every `<img>` in the
  HTML report is a `data:image/png;base64,` URI by default, so the single file is safe to move,
  email or import into Notion; the PNGs are still written beside it. `params.report_embed_images`
  (default `true`; `false` = relative filenames as before) controls this.
- **Data 2 header guards now apply to every two-file legacy comparison** (`pairedT`,
  `unpairedT`, two-file `onesampleT`, both correlations), not just the correlations. Channel
  selection is matched on Data 1's header and the same positions are applied to Data 2, so a
  metadata column in one file but not the other, or a wider Data 1, used to pair channel *i*
  with a different channel silently (or die on `MATLAB:badsubscript`). Such runs now stop with
  `snpm:correlation:channelOrderMismatch` / `snpm:correlation:widthMismatch`, with a
  `LIKELY CAUSE` clause naming the offending metadata column. **Breaking** for two-file runs
  whose headers differ: Data 1 must be `Subject` + channel columns only. A single-measure Data 2
  under a t-test is refused (`snpm:data2:behaviorNotSupported`).
- **No popup windows during a run.** All figures the pipeline creates are hidden and closed on
  completion; a headless `-batch` run leaves no figure handles behind (pinned by
  `test_legacy_snpm`).

### Removed

- `dependencies/TopoplotSignificant.m` (no callers; `TopoplotSignificant_single.m` is the
  live topoplot path).
- **The internal verification harness and the synthetic-data generators are no longer
  distributed** — the `test_*.m` suites, `test_fixtures/`, and the `gen_synthetic_*.m`
  generators. They remain in use privately; nothing about the statistics changed. What this
  costs a user of a clone is the ability to re-run the checks themselves, so the documentation
  no longer claims otherwise: statements of the form "run this and it comes back green" are
  gone, while statements about what the software *was* checked against (MATLAB's `fitlm`,
  `anovan`, `fitlme`, `ttest`, `corr`) are kept, since those describe the software's provenance
  and are still true. `example_data/` lets you reproduce the planted-cluster-recovery half of
  that verification directly.

## [1.0.0] — 2026-08-08

First public release. Everything below is the content of 1.0.0; the entries are grouped by what
they do rather than by the development commit they arrived in, and the dates mark when the work
landed during private development.

### The toolbox

Permutation-based (TFCE and cluster-extent) statistics on high-density EEG topographies, with
channel-adjacency structure, following Maris & Oostenveld (2007) and Nichols & Holmes (2001).
It ships the core statistics, the GLM and LMM tiers, the App Designer GUI, the reporting stack,
the montage assets and the documentation. Vendored EEGLAB, run outputs, reference PDFs and raw
EEG are deliberately not tracked. (The verification suites named below were part of the private
development tree at 1.0.0 and are not distributed — see Unreleased.)

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
- **Verification coverage for six behaviours that had never been exercised before 1.0.0**: the
  correlation and LMM missing-data invariants, the common-channel rule in the whole-head omnibus,
  the circular tier end to end and its statistics against independent references, the cluster-
  enhancement identity, and the source path. (Development-tree only; not distributed.)

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

Before release, the internal harness ran 16 checks — static analysis, a headless GUI construct,
the statistic-versus-MATLAB-builtin suites, per-analysis output-file assertions, a spectral smoke
run and the edge/failure-mode cases — with 0 skipped and 1 known expected failure (the event
bridge, whose fixture used channel indices that did not match the analysed 178-channel subset;
see the known gaps in the analysis catalog). That harness is not part of the public repository.

[1.0.0]: https://github.com/TancyKao/hdEEG-SnPM/releases/tag/v1.0.0
