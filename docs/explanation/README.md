# Explanation

*Understanding-oriented.* Background and rationale — read these to understand *why* the tool
works the way it does. They don't give step-by-step instructions (see
[how-to guides](../how-to/README.md)) or exhaustive facts (see
[reference](../reference/ANALYSIS_CATALOG.md)).

## Pages

- **[Detection and statistics are two tools](detection-and-statistics-two-tools.md)** — where
  this toolbox sits relative to its upstream sibling
  [TurtleWave-hdEEG](https://github.com/TancyKao/TurtleWave-hdEEG), why the two halves are
  separate programs joined by a file rather than one, and why TurtleWave is a supported source
  and not a dependency.
- **[Choosing an analysis: design, missing data, and cost](choosing-an-analysis.md)** — why the
  presets are all one linear model, how unequal-n / dropout affects the choice, when the LMM is
  worth its cost, and a sleep-question → preset cheat-sheet.
- **[Two-way mixed ANOVA — how the interaction is computed](mixed2way-interaction.md)** — the
  per-channel GLM, the interaction contrast, the Freedman–Lane within-subject permutation, and
  a fully worked 2×2 example. The deepest single-analysis explanation in the set.
- **[Why this tool exists alongside a PALM-based pipeline](eeg-processor-vs-snpm.md)** — an
  honest comparison against an FSL PALM group-analysis path: what overlaps, and the specific
  reasons this tool exists alongside it (preset-driven UX; the mixed-ANOVA 3-effect report;
  per-channel LMM; circular designs).
- **[Interpreting cluster results](interpreting-cluster-results.md)** — cross-tier interpretation
  points: how each tier scores a cluster (extent on the legacy path; summed or mean Wald on the
  GLM, circular and mixed-model tiers, chosen by `cluster_statistic`, so results aren't
  comparable across tiers or settings), why a significant cluster is a
  regional claim rather than an electrode-by-electrode one, and why a circular map carries no
  direction of effect.
- **[Missing data and the complete-column rule](missing-data-and-excluded-channels.md)** — why a
  correlation analysis drops a whole channel for one missing cell, what pairwise deletion does to
  a permutation null (measured family-wise error up to 0.288 at a nominal 0.05), why the t-tests
  are deliberately not masked, and how to read the excluded-channel counts in your outputs.
- **[Covariates in a correlation: the Freedman–Lane permutation](covariates-in-correlation.md)**
  — why permuting covariate residuals against each other gave a null that was too narrow
  (family-wise error up to 0.86), what the Freedman–Lane re-projection does, why the test is
  approximate, and in which direction at few error degrees of freedom, and what partial Spearman
  controls.
- **[About circular statistics for phase](circular-statistics-for-phase.md)** — why phase needs
  its own machinery, why coupling *strength* and *prevalence* do not, what Hotelling and Watson's
  U² each license you to claim, the honest power figures, the precision confound and why there is
  no Rayleigh mask, and the two places this tier extends the literature rather than following it.

## More background

Deeper design history — the expansion plan and the v2.0 changelog — lives in the
[archive](../archive/README.md). It is frozen and may be stale; the pages above are the
maintained explanation.
