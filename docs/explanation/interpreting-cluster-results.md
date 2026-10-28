# Interpreting cluster results — what a significant cluster does and doesn't say

> *Understanding-oriented.* Two interpretation points that apply to every analysis in the
> tool. For the factual per-analysis table (inputs, statistic, outputs) see the
> [Analysis catalog](../reference/ANALYSIS_CATALOG.md); to run an analysis, see the
> [how-to guides](../how-to/README.md).

Two things routinely get over-read in a topoplot: which cluster statistic a tier uses, and
what a highlighted electrode actually licenses you to claim. Both matter for how you write up
a result.

## How each tier scores a cluster

Every tier forms clusters the same way: neighbouring channels whose per-channel p-value is
below 0.05. They differ in the single number used to score a cluster, and that number is what
the permutation null is built from.

- **The legacy t-test path** (`pairedT`, `onesampleT`, `unpairedT`, and the correlation
  comparisons) uses cluster **extent**: the count of contiguous supra-threshold channels.
  How strong the effect is inside the cluster plays no part.
- **The GLM, circular and mixed-model tiers** (`anova1`, `ancova`, `regression`, `rmanova`,
  `mixed2way`, `circ_phase_group`, `circ_phase_group_u2`, `circ_corrAngLinear`,
  `mixedmodel`) use cluster **mass**, and the option `cluster_statistic` decides which mass:
  - `'sum'` (the default since 2026-10-28) adds up the Wald statistic over the cluster's
    channels: t² for a one-row contrast, F for an omnibus effect (on the circular tier, the
    statistic itself: Hotelling F, Watson U² or the circular-linear F). This is
    a summed-statistic cluster mass in the sense of Maris & Oostenveld (2007), taken here on the Wald scale (t² or F). It rewards both
    height and extent: a broad cluster of moderate channels can outscore a single strong
    channel.
  - `'mean'` divides that sum by the number of channels, which makes it the mean Wald over
    the cluster (the recipe of Stephan et al. 2021, and what every run before 2026-10-28
    used). A broad 20-channel cluster and a focal 3-channel cluster with the same average
    Wald score identically, so spreading an effect over more channels does not raise the
    statistic. It costs power for spatially broad effects: the permutation maximum of a mean
    is set by single strong null channels, and a broad moderate cluster rarely beats them.
    Keep it only to reproduce an older result.

Either choice gives an exact test. The same rule scores the observed clusters and the largest
cluster of every permutation, and the cluster-forming threshold, the Wald scale and the
permutation draws do not depend on it. TFCE is a separate correction and is unaffected: its
maps and p-values are identical under both settings. The setting used is recorded in the
results struct (`cluster_statistic`), in the results text, in a `settings` sheet of the Excel
workbook and in the report's model card and cluster view.

The circular tier's mass differs from Helfrich et al. (2018), who used FieldTrip's `maxsize`
(cluster extent); the results text of every circular run says so.

A practical consequence: **cluster statistics are not directly comparable across tiers, or
across the two `cluster_statistic` settings.** An extent, a summed mass and a mean mass are
different numbers on different scales. When you compare or combine results, compare the
corrected p-values and the per-channel maps, not the raw cluster statistics.

## A significant cluster is a regional claim, not an electrode-by-electrode one

A significant cluster licenses the claim **"there is an effect somewhere in this
region/cluster"** — *not* "these specific electrodes are significant." This is a property of
how permutation cluster and TFCE inference work: they control the family-wise error rate at
the **cluster level**, not at the level of individual channels (Maris & Oostenveld 2007). The
test asks whether a cluster of *this* mass/extent could arise under the null anywhere on the
scalp; a "yes, it's significant" answer attaches to the cluster as a whole.

The topoplots make this easy to forget, because they show highlighted individual electrodes
and invite you to read them one at a time. Resist that. The electrodes inside a significant
cluster mark *where the regional effect is concentrated*, but you cannot single out one
channel and call it individually significant, nor read the exact cluster boundary as a sharp
anatomical edge — a neighbouring channel just outside the outline is not established as
null. Write up the result as a region ("a centro-parietal cluster showed …"), not as a list
of certified channels.

This point holds for **all tiers** — legacy t-tests, GLM presets, the LMM and the circular
analyses — because they all rest on the same neighbour-based permutation cluster / TFCE
machinery.

## A circular significance map carries no direction of effect

There is a further restriction specific to the circular tier, and it is easy to miss because the
map *looks* like every other topoplot.

The Hotelling T², Watson's U² and the circular–linear F are all **non-negative omnibus
statistics**. A significant circular cluster says the two groups' phase distributions *differ*
somewhere in that region — which may be a shift in mean direction, a difference in concentration,
or a mixture of the two. It does not say which, and it does not say which way. There is no sign
to read off the colour bar, which is why `tail` is locked to `both` for these analyses.

So "the patient group's spindles couple **later** in the up-state" is not a claim a circular
cluster supports, however tempting the topography makes it. For a directional statement, run the
**signed** linearised measure through `unpairedT` / `ancova`, which is an ordinary signed t-map
and does carry a direction. The same applies to circular–linear correlation: significance means
the behavioural variable varies with phase, never that it is higher or better at any particular
phase. See
[About circular statistics for phase](circular-statistics-for-phase.md).

## References

- Maris E, Oostenveld R (2007). *Nonparametric statistical testing of EEG- and MEG-data.*
  J Neurosci Methods 164(1), 177–190 — neighbour clustering and cluster-level permutation
  inference for M/EEG; the source of the region-not-channel scope of a cluster.
- Nichols TE, Holmes AP (2001). *Nonparametric permutation tests for functional neuroimaging.*
  Hum Brain Mapp 15(1), 1–25 — max-statistic permutation, family-wise-error control.
- Stephan AM et al. (2021). *Conscious experiences and high-density EEG patterns predicting
  subjective sleep depth.* Curr Biol 31(24), 5487–5500 — the `sum(Wald) / n_channels`
  (mean-Wald) cluster statistic, available as `cluster_statistic = 'mean'` and used by every
  GLM, circular and mixed-model run before 2026-10-28.
