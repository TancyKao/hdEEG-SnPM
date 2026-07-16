# Explanation

*Understanding-oriented.* Background and rationale — read these to understand *why* the tool
works the way it does. They don't give step-by-step instructions (see
[how-to guides](../how-to/README.md)) or exhaustive facts (see
[reference](../reference/ANALYSIS_CATALOG.md)).

## Pages

- **[Choosing an analysis: design, missing data, and cost](choosing-an-analysis.md)** — why the
  presets are all one linear model, how unequal-n / dropout affects the choice, when the LMM is
  worth its cost, and a sleep-question → preset cheat-sheet.
- **[Two-way mixed ANOVA — how the interaction is computed](mixed2way-interaction.md)** — the
  per-channel GLM, the interaction contrast, the Freedman–Lane within-subject permutation, and
  a fully worked 2×2 example. The deepest single-analysis explanation in the set.
- **[EEG_Processor group analysis vs. SnPM_2025](eeg-processor-vs-snpm.md)** — an honest
  comparison against EEG_Processor's PALM path: what overlaps, and the specific reasons this
  tool exists alongside it (preset-driven UX; the mixed-ANOVA 3-effect report; per-channel LMM;
  circular designs).
- **[Interpreting cluster results](interpreting-cluster-results.md)** — two cross-tier
  interpretation points: why the GLM/LMM cluster statistic is height-weighted (mean Wald) while
  the legacy path is extent-weighted (so results aren't comparable across tiers), and why a
  significant cluster is a regional claim, not an electrode-by-electrode one.

## More background

Deeper design history — the expansion plan, session recaps, and the v2.0 changelog — lives in
the [archive](../archive/README.md). It is frozen and may be stale; the pages above are the
maintained explanation.
