# Interpreting cluster results — what a significant cluster does and doesn't say

> *Understanding-oriented.* Two interpretation points that apply to every analysis in the
> tool. For the factual per-analysis table (inputs, statistic, outputs) see the
> [Analysis catalog](../reference/ANALYSIS_CATALOG.md); to run an analysis, see the
> [how-to guides](../how-to/README.md).

Two things routinely get over-read in a topoplot: which cluster statistic a tier uses, and
what a highlighted electrode actually licenses you to claim. Both matter for how you write up
a result.

## The cluster statistic is height-weighted in the GLM/LMM tier, extent-weighted in the legacy tier

The tool has two families of analyses that summarise a cluster differently, and the two
answer subtly different questions.

- **GLM and LMM presets** (`anova1`, `ancova`, `regression`, `rmanova`, `mixed2way`,
  `mixedmodel`) use cluster **mass normalised by cluster size** — `sum(Wald) / n_channels`,
  the recipe from Stephan et al. (2021). Dividing by the channel count makes this effectively
  a **mean Wald** over the cluster. It removes the reward-for-extent that a plain cluster-mass
  sum carries: a broad 20-channel cluster and a focal 3-channel cluster with the *same average
  Wald* score identically. The consequence is a loss of power for spatially broad effects —
  spreading the same total signal over more channels does not raise the statistic.
- **The legacy t-test path** (`pairedT`, `onesampleT`, `unpairedT`, and the correlation
  comparisons) uses cluster **extent** — cluster *size* only, the count of contiguous
  supra-threshold channels. Here spatial spread is exactly what is rewarded.
- **The circular tier** (`circ_phase_group`, `circ_phase_group_u2`, `circ_corrAngLinear`) uses
  cluster **mass**, like the GLM/LMM tier. Note that this differs from Helfrich et al. (2018),
  who used FieldTrip's `maxsize`; the difference is disclosed in the results text of every
  circular run.

So the two tiers reward different things: the legacy path rewards **spatial extent**, the
GLM/LMM path rewards **per-channel effect height**. A practical consequence: **cluster
statistics are not directly comparable across tiers.** A cluster that looks strong under a
paired t-test (large because it is wide) may look modest under the mean-Wald GLM statistic,
and vice versa, without either result being wrong. When you compare or combine results from
the two families, compare the corrected p-values and the per-channel maps, not the raw cluster
statistics.

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

- Maris & Oostenveld (2007) — neighbour clustering and cluster-level permutation inference for
  M/EEG; the source of the region-not-channel scope of a cluster.
- Nichols & Holmes (2001) — max-statistic permutation, family-wise-error control.
- Stephan et al. (2021) — the `sum(Wald) / n_channels` (mean-Wald) cluster statistic used by
  the GLM/LMM tier.

PDFs are in `Ref/`.
