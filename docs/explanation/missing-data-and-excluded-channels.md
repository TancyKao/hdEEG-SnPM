# Missing data and the complete-column rule

> *Understanding-oriented.* Why a correlation analysis now drops a whole channel for a single
> missing value, why the t-tests drop a much smaller set of channels for a different reason, and
> how to read the excluded-channel counts that appear in your outputs. For
> the factual per-analysis table see the [Analysis catalog](../reference/ANALYSIS_CATALOG.md);
> to run an analysis, see the [how-to guides](../how-to/README.md).

## What changed

Until 2 August 2026 the correlation analyses (`correlationP`, `correlationS`) handled missing
values by *pairwise deletion*: each channel used whatever subject pairs happened to be complete
at that channel, so a channel with two missing subjects was still tested, just on a smaller
sample than its neighbours. That is the standard convenience choice in ordinary correlation
analysis, and outside a permutation test it is defensible.

Since 2 August 2026 those analyses follow the **complete-column rule** instead. A channel is
tested only if *every* analysed subject has a value there. One missing cell and the channel
leaves the analysis entirely, is reported as excluded, and carries no statistic and no p-value
in either the observed map or the null. The GLM presets and the circular tier already worked
this way; the correlation paths have now been brought into line with them.

Nothing changes for complete data. If your file has no missing cells, the results are identical
to what the tool produced before.

## Why pairwise deletion breaks a permutation test

The reasoning is specific to permutation inference, which is why this is not simply a matter of
taste between deletion strategies.

A permutation test builds its null distribution by reordering the rows of one data matrix
against a fixed other matrix, recomputing the whole map, and keeping the largest statistic. Do
that ten thousand times and you have the distribution of "the biggest thing this analysis
produces when there is nothing there", which is what a corrected p-value is scored against.

Now add a partially-missing channel. The correlation at that channel is computed from the
subjects who have a value in *both* matrices — but after a permutation, subject 4's row of the
first matrix sits opposite subject 19's row of the second. A subject who was complete before the
shuffle may now be paired with a missing cell, and vice versa. **The set of complete pairs is
therefore a function of the permutation.** The channel is evaluated at one sample size in the
observed map and at a different, randomly varying sample size in every permuted map.

That matters because the statistic being mapped is the raw correlation *r*, and the spread of
*r* under the null scales as 1/√(n−1). A correlation computed on 22 subjects and a correlation
computed on 15 subjects are simply not on the same scale, so a null assembled from the second is
the wrong ruler for the first. Measured on a 30-subject design with 25% missingness at one
channel: the observed value was computed at n = 22 while the permuted values came from n in the
range 14 to 20. The observed sample size fell outside the entire permuted range — the null was
built at a systematically smaller n, and therefore a systematically wider spread, than the
statistic it was scoring.

The consequence is a real, measured loss of error control rather than a theoretical worry. At a
nominal 0.05, over 2000 replicates, the family-wise false-positive rate of the TFCE-corrected
map reached **0.288** with 40% disjoint missingness and **0.072** at 25%. A second, related
mechanism sat in the covariate path, where missing cells used to be imputed as zeros before
residualization: that erased the missingness pattern altogether and produced a family-wise error
of **0.168** for TFCE and **0.076** for cluster extent. Cluster extent was consistently less
affected than TFCE — it depends on how many channels cross a threshold rather than on the exact
height of each one — but it was still above nominal. Under the complete-column rule the same
measurements sit at nominal (see
[Permutation p-values and sampling](../reference/ANALYSIS_CATALOG.md#permutation-p-values-and-sampling)
in the catalog).

## Why the t-tests do not get the complete-column rule

A paired, one-sample or unpaired t-test still tests a channel that some subjects are missing, and
that is deliberate rather than an oversight.

Those nulls do not reorder one matrix against another. They **sign-flip or relabel whole subject
rows**, so a subject's missing cells travel with the subject. Whatever pattern of missingness a
channel has in the observed map, it has exactly the same pattern under every permuted labelling,
and the effective n at that channel is constant across the null. The measurements confirm it:
the per-channel n was [22, 22] and [14, 14] across sign-flip permutations, and the pooled
unpaired degrees of freedom held at [20, 20] with no NaN permutations at all. Applying the mask
there would cost channels for no inferential gain, so it is not applied.

## The t-tests do have a narrower rule: degeneracy under relabelling

Since 3 August 2026 the t-test paths drop a channel when the relabelling itself stops making
sense there. The argument above holds only while every labelling in the null produces a usable
statistic; a channel with very few remaining subjects breaks that, and it breaks it in the
direction that costs you real effects.

For a paired or one-sample test the threshold is **two usable pairs**. With fewer, the paired t
has zero degrees of freedom and returns ±Inf under every sign flip — and an infinite value does
not merely add a bad entry to the null, it breaks the TFCE threshold ladder outright
(`MATLAB:pmaxsize`), which used to surface as an opaque out-of-memory crash naming no channel.

For an unpaired test the threshold is **`max(n1, n2) + 2` usable subjects pooled**, which is the
point below which some relabelling in the null puts one group at a single subject. MATLAB's
variance of a scalar is 0 rather than NaN, so `ttest2` does not refuse: it returns a *finite* t
computed from a pooled standard deviation that ignores that group, and the value is large. On
pure noise with one subject against ten, the largest |t| over 20 000 draws was 7.4. That single
number then becomes the maximum statistic for that permutation and sets the corrected threshold
for the entire montage.

The cost is a loss of sensitivity everywhere else, not a false positive at the bad channel. On an
8×8 lattice with 800 replicates, a planted, fully-observed nine-channel cluster was TFCE-detected
93.4% of the time with no degenerate channel present and **16.4%** of the time when four channels
were available in only 3 of 20 subjects. Masking those four channels restored detection to 94.1%.
Four unusable channels were costing three quarters of the power of a cluster they had nothing to
do with.

This rule is much narrower than the complete-column rule. A channel missing two subjects out of
twenty is still tested on the eighteen it has; only channels that are near-empty leave. On
complete data it never fires. Its exclusions are reported through the same channels as everything
else described below.

## The whole-head omnibus averages over a common channel set

The omnibus test at the top of every report — the single p-value for "anything at all across the
head" — changed on the same date. It used to average each unit over whatever channels that unit
happened to have; it now averages every unit over the channels that are finite in **every**
analysed row of both arms.

The reason is that scalp power is not spatially flat. NREM delta runs roughly four times higher
frontally than occipitally, so which channels a subject is averaged over sets a per-subject offset
of the same size as the effects this tool looks for. Two subjects averaged over different channel
sets are not measuring the same quantity, and the difference between their means is then partly a
difference between their channel sets.

Measured with no true effect anywhere at a nominal 0.05: an unpaired test in which 4 of 10
group-1 subjects lost an occipital band ran at 0.112. The paired case is worse and more common —
two nights with different bad channels is routine — and when the missing region differed between
condition A and condition B the test rejected **every single time** (1.000), because the
condition difference contained a fixed frontal-minus-occipital offset in every subject. Under the
common-channel rule all of those configurations sit at nominal (0.047 to 0.053), since the offset
is now identical in both arms and cancels in the contrast.

The correlation path is unaffected in practice: its complete-column mask has already made the two
channel sets identical by the time the omnibus runs. The reports state how many channels the
average was taken over, which is the number to check if a global p and a topography disagree.

## What you will see in your outputs

A channel never leaves an analysis silently, whichever rule dropped it. When anything is
excluded, the run produces all of the following:

- A console line and a MATLAB warning (`snpm:channelsExcluded`) naming the count and the
  individual channel labels.
- An **`excludedChannels`** sheet in the `.xlsx`, one row per dropped channel with its index,
  its label, and the reason.
- An amber **Excluded channels** banner near the top of the HTML report, stating how many of how
  many channels were not tested and listing them.
- A `--- Channel Exclusions ---` block in the results text, and
  `results_struct.excluded_channels` in the `.mat`.

Two further warnings fire on top of that. If more than roughly 10% of the montage is excluded,
the run says so (`snpm:manyChannelsExcluded`) and suggests dropping the handful of subjects
driving the missingness instead of losing that much of the head. And because *where* the
exclusions fall matters more than how many there are, the run also profiles their spatial
distribution and warns (`snpm:exclusionsLopsided`) when they concentrate on one side or one
region.

That last warning is the one to take seriously. Losing twelve scattered channels and losing
twelve contiguous frontal channels give the same count and have completely different
consequences: the second quietly removes the region your hypothesis is about, and a
"no frontal effect" conclusion drawn from a map with no frontal channels in it is not a null
result. Read the excluded list, and the topography of the exclusions, before you read the
p-values.

## The trade, stated honestly

The complete-column rule is deliberately conservative. **One missing cell costs the whole
channel**, however many subjects were fine there. A channel that was 29/30 complete is treated
exactly like one that was 15/30.

That is the price of a null that is on the same scale as the statistic it scores, and it is the
same price the GLM tier has always paid — there for a slightly different reason (the
Freedman–Lane nuisance fit works on whole columns, so a single missing cell would turn the
column all-NaN inside each permuted fit while leaving the observed map finite). Having one rule
across the tiers also means the analysed channel set no longer depends on which preset you
picked.

There is a second-order consequence worth naming: dropping incomplete channels is a
complete-case analysis at the channel level, and complete-case analysis is unbiased only when
the missingness is unrelated to the quantity being measured. If channels go missing *because* of
the effect you are studying — an artefact that scales with the measure, say — no deletion rule
saves you, and the exclusion list is telling you something about the data rather than about the
tool.

## Reducing exclusions

Exclusions are a data problem, so the useful fixes are upstream of this tool.

The first question is always *why* those cells are missing. Channels rejected in preprocessing
for a few subjects, an electrode that failed in one recording, a detector that reported nothing
at a channel with too few events, and a spreadsheet with genuinely empty cells all look
identical here and have different remedies. Fixing the cause recovers the channel for everyone.

When that is not possible, the usual arithmetic favours dropping subjects over dropping
channels. If two subjects account for most of the missing cells, removing those two subjects
often returns the whole montage at a small cost in power, whereas keeping them can cost a large
contiguous piece of the head. The per-subject NaN counts printed by the run tell you whether the
missingness is concentrated in a few subjects or spread thin across everyone; only the first case
has this easy fix.

Missing *covariate* values behave differently and more bluntly. A subject with a missing
covariate has no defined design-matrix row at any channel, so that subject is dropped from the
analysis outright rather than costing a channel. Correlation needs at least three subjects to
survive this step.

Finally, if a large share of the montage is excluded no matter what you do, treat that as a
reason to reconsider the analysis rather than to proceed carefully. A map with a third of the
head missing supports far weaker regional claims than its topoplot suggests.

## Further reading

- [Choosing an analysis: design, missing data, and cost](choosing-an-analysis.md) — how
  *subject-level* dropout affects which preset you should pick, and when the LMM's
  missing-at-random handling is worth its cost.
- [Interpreting cluster results](interpreting-cluster-results.md) — why a cluster is a regional
  claim, which is exactly what a lopsided exclusion pattern undermines.
- [Analysis catalog](../reference/ANALYSIS_CATALOG.md) — the per-analysis facts, including the
  excluded-channel outputs and the meaning of the reported per-channel n.
- Primary sources in the code: the comment block in `core_snpm_analysis.m` that computes the
  evaluable mask (both rules are derived there, from the same per-channel availability count),
  and the headers of `dependencies/snpm_corr_columns.m` and `dependencies/global_stat_test.m`
  (the common-channel rule and its measured false-positive table). Each of the three rules —
  correlation, the LMM tier, and the whole-head omnibus — was checked by the same invariance
  argument used above: given the retained channel set, the result must be identical whether or
  not the dropped cells were ever present.
