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
  points: why the GLM/LMM cluster statistic is height-weighted (mean Wald) while the legacy path
  is extent-weighted (so results aren't comparable across tiers), why a significant cluster is a
  regional claim rather than an electrode-by-electrode one, and why a circular map carries no
  direction of effect.
- **[Missing data and the complete-column rule](missing-data-and-excluded-channels.md)** — why a
  correlation analysis drops a whole channel for one missing cell, what pairwise deletion does to
  a permutation null (measured family-wise error up to 0.288 at a nominal 0.05), why the t-tests
  are deliberately not masked, and how to read the excluded-channel counts in your outputs.
- **[About circular statistics for phase](circular-statistics-for-phase.md)** — why phase needs
  its own machinery, why coupling *strength* and *prevalence* do not, what Hotelling and Watson's
  U² each license you to claim, the honest power figures, the precision confound and why there is
  no Rayleigh mask, and the two places this tier extends the literature rather than following it.

## More background

Deeper design history — the expansion plan and the v2.0 changelog — lives in the
[archive](../archive/README.md). It is frozen and may be stale; the pages above are the
maintained explanation.
