# Choosing an analysis: design, missing data, and cost

*Understanding-oriented.* Why the presets are organised the way they are, and how to
pick one for a real sleep-EEG design. For the factual per-analysis table (inputs,
required columns, statistic, fixture, outputs) see the
[Analysis catalog](../reference/ANALYSIS_CATALOG.md). To actually run a chosen
analysis, see the [how-to guides](../how-to/README.md).

## It's one model

The t-tests, ANOVA, ANCOVA and regression are all special cases of the linear model
`Y = Xβ + ε` — the tool runs them through a single OLS engine (`snpm_glm_stat`). Only
the *coding* of the predictors differs, so pick the preset that matches your variable:

- **Categorical** variable of interest (groups) → two-sample t-test / one-way ANOVA / ANCOVA.
- **Continuous** variable of interest → linear regression.

(A two-group t-test = regression on one 0/1 group dummy; k-group ANOVA = k−1 dummies;
ANCOVA = group dummies + a continuous covariate. The presets do the dummy-coding and
choose t vs F for you, so you don't need to hand-build a regression for group comparisons.)

## Unequal n across conditions (e.g. medical dropout)

- **Between-subject** designs — two-sample t-test, one-way ANOVA, ANCOVA, regression,
  correlation — tolerate unequal group sizes **by construction**. Dropout just means a
  smaller group; use freely. These GLM presets assume **one row per subject**: `anova1`,
  `ancova` and `regression` now **error** (`core_snpm_glm:repeatedMeasures`) if a subject
  appears more than once, pointing you to `rmanova` or `mixedmodel` for repeated measures.
- **Within-subject** designs — paired / one-sample t-test, RM-ANOVA, two-way mixed ANOVA —
  need the *same* subjects across conditions. A subject missing a condition is dropped
  (**complete-case** → fewer subjects, less power). Pairing is by `sub-XX` (spectral folder)
  or by row order (CSV).

**For substantial within-subject dropout, prefer the LMM** (`mixedmodel`): it uses all
available observations and is valid under **Missing-At-Random (MAR)**. This matters because
medical dropout is rarely random (sicker patients leave); RM-ANOVA's complete-case analysis
is only unbiased under the stronger **MCAR**.

## But the LMM is computationally heavy

Per-channel *iterative* `fitlme` × permutations is roughly **100–1000× slower** than the
vectorized-OLS GLM presets (minutes-to-an-hour vs seconds). Cheaper routes that usually
suffice:

- **RM-ANOVA GLM** — this tool's `rmanova` models subject as a **fixed effect** (subject
  dummies) by OLS, so it's fast and tolerates *moderate* unbalance; you give up the
  random-effects/partial-pooling and the cleaner MAR theory of a true LMM.
- **Aggregate, then GLM** — if the load is *trial-level*, collapse to subject×condition
  **means** first, then run the fast `rmanova`/`mixed2way`.

Reserve the full LMM for genuinely trial-level data, real random-effects structure, or
heavy non-random dropout. Cost knobs: ≤1000 permutations, keep the `parfor` pool, analyse
fewer channels. The LMM is **script-only** (`scripts/run_lmm_analysis.m`) — not yet in the GUI.

## Sleep question → preset (cheat-sheet)

| Sleep question | Preset |
|---|---|
| Is the per-channel measure ≠ 0 / a fixed reference? | one-sample t (`onesampleT`) |
| Same people, two conditions differ? (baseline vs recovery; pre/post-CPAP) | paired t (`pairedT`) |
| Two independent groups differ? (OSA vs control) | unpaired t (`unpairedT`) |
| 3+ groups differ, then which pair? (control/mild/severe OSA) | `anova1` (+ auto post-hoc) |
| Groups differ after removing age/sex/AHI? | `ancova` |
| Measure scales with a continuous variable (± covariates)? (SWA vs AHI) | `regression` |
| 3+ within-subject conditions/stages differ? (N2/N3/REM) | `rmanova` |
| Does the condition effect differ between groups? (group × stage) | `mixed2way` |
| Monotonic power–behaviour association? | `correlationP` / `correlationS` |
| Trial/event-level effect with many observations per subject? (per-awakening power vs sleep depth, Stephan 2021) | `mixedmodel` |
| Coupling **strength** or event **prevalence** differs? | the ordinary presets — `unpairedT` / `anova1` / `regression`. These are linear quantities and need no circular test. |
| Preferred **phase** differs between two independent groups? | `circ_phase_group` (Hotelling T², primary) or `circ_phase_group_u2` (Watson U², secondary) |
| Does a behavioural measure vary with phase? | `circ_corrAngLinear` |
| Want a **directional** phase claim ("group A couples later")? | the *signed* linearised measure through `unpairedT` / `ancova` — the circular tests are non-directional |

**Caveat locked from the the example study data:** `condition-a` / `condition-b` are within-subject
(same subjects) → compare with `rmanova` / `pairedT`, **not** `anova1` (which assumes
independent groups).

## Method grounding

- **Maris & Oostenveld (2007)** — per-sensor statistic (t for two conditions, F for >2),
  neighbour clustering, cluster-level statistic, max-cluster permutation null for FWER.
- **Nichols & Holmes (2001)** — single-threshold maximum-statistic permutation test;
  sign-flipping for one-sample/within.
- **Winkler et al. (2014)** — Freedman–Lane permutation for the GLM with nuisance regressors.
- **Stephan et al. (2021)** — the event/trial-level LMM tier.

PDFs are in `Ref/`.
