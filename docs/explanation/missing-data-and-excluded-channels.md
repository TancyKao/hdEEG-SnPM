# Missing data and the complete-column rule

> *Understanding-oriented.* Why a correlation analysis now drops a whole channel for a single
> missing value, why the t-tests drop a much smaller set of channels for a different reason, why
> the two-way mixed design uses one channel set and refuses a table with missing rows, why a
> channel with no variation is left untested, and how to read the excluded-channel counts and
> reasons that appear in your outputs. For
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

One separate case is flagged rather than averaged: a correlation run with `datatype =
'normalize'` z-scores each subject across channels, which makes the channel-averaged signal
analytically zero — the global statistic would be floating-point noise, not a whole-head
correlation. The run marks it as such (`results_struct.global_not_interpretable`, warning
`snpm:globalNotInterpretable`, an amber note on the global card); the channel-wise maps are
unaffected.

## The two-way mixed design: one channel set, and no missing rows

The two-way mixed analysis of variance (`mixed2way`) reports several effects from one table: the
group effect, the condition effect, their interaction, the condition effect within each group
(the simple effects), the pairwise post-hoc comparisons, and a whole-head test for each of the
three main effects. Since 14 October 2026 two rules apply to all of them together.

**Every effect uses one channel set.** A channel is analysed only if it has a value in every row
of the table, and that single set is fixed before any permutation and used for every effect in
the report. The channels excluded for missing data are therefore the same in every section; a
section can list a few more, excluded for having no variation in that test (see
[Channels with no variation](#channels-with-no-variation) below). Before, the group
effect was computed on subject means taken over whatever values a subject had, so it could test
channels that the condition and interaction effects could not, and the three whole-head tests
could rest on three different sets of channels. The rule also holds one level down: a pairwise
post-hoc comparison, in this and in the other General Linear Model presets, tests a channel only
if the omnibus effect tested it.

**Every subject must have every condition.** The table must hold exactly one row per condition
for every subject, and each subject must belong to one group. If a subject is short of a
condition, has a condition twice, or appears in two groups, the run stops before computing
anything with the error `core_snpm_glm:mixed2wayUnbalanced`, which names the subjects and says
what is missing or duplicated.

This is a different kind of missingness from a missing channel value, and it is refused rather
than worked around because of what it does to the group effect. The group effect is tested on one
mean per subject, taken over that subject's conditions. If a subject has no row for one
condition, its mean is taken over the remaining conditions only. Whenever the conditions differ,
that subject's mean is shifted by the condition effect, and the shift is counted as a difference
between subjects. Measured on data with a condition effect and no group effect, with 7 of 20
subjects missing a condition, the group effect rejected 0.61 of the time at a nominal 0.05. With
the requirement in place, every subject mean is a plain mean over the same conditions.

The remedy is the one this page gives throughout. Remove the subjects the error names, or
complete their data, and run again. If a large share of subjects are incomplete, or the data are
at the level of single trials, use the mixed model (`mixedmodel`), which uses every available
row. Unequal group sizes are not a problem and are not affected by this rule.

Within-subject post-hoc comparisons in `rmanova` handle the same situation per pair: a subject
who has only one of the two conditions being compared is left out of that comparison, with a
warning naming the subject, and stays in the other comparisons.

## Channels with no variation

A channel can be complete and still give a permutation test nothing to work with. Event counts
and densities are the common source in sleep EEG: a detector that finds no spindles at a
temporal channel in most subjects leaves a column of zeros. Integer scores, very small samples
and a channel flattened in preprocessing do the same. Since 14 October 2026 every tier handles
these channels with one rule in three parts.

### What happened before

The t-test builtins do not refuse such a channel. If every paired difference is the same
non-zero number, `ttest` returns an infinite t; if two groups are each constant at different
values, `ttest2` returns either an infinite t or, for 0.1 against 0.3, a finite t of about
3 × 10¹⁶, because rounding leaves a tiny non-zero standard deviation. The threshold-free cluster
enhancement (TFCE) step then tried to build a threshold ladder up to that value and stopped with
an out-of-memory or array-size error that named no channel. This is the crash that was found
and deferred in August 2026. It did not need an obviously broken channel: a channel where the
paired differences take only the values +1 and −1 is fine as observed, but one of the sign-flip
relabellings makes every difference equal, and that single relabelling stopped the run. At the
other extreme, a channel constant at a non-zero value (every subject at 0.1, say) gave `ttest2`
a t of exactly 0 and the correlation engine an r of exactly 0, so it was kept, tested and never
mentioned.

### The rule

**First, a channel with nothing to permute is excluded up front and reported with a reason.**
This is decided once, before the permutation loop, on the values the test actually uses: a
paired or one-sample channel whose differences are all zero, an unpaired channel with the same
value in every subject, a correlation channel that is constant on either side (on the ranks for
Spearman, and after removing the covariates when there are covariates), and, in the General
Linear Model presets, a channel whose residuals after the nuisance model are all zero, so that
every permuted map would equal the observed one. In the General Linear Model presets this is
decided separately for each test, on that test's own rows and nuisance model, so a post-hoc pair
can exclude a channel its omnibus tested (a channel constant within two of three groups, for
instance). Each test lists its own exclusions. A constant non-zero paired difference is *not*
excluded here, because sign flips change it.

**Second, at any labelling where the error sum of squares is zero, the statistic is defined as
0** (and its parametric p as 1). This applies to the observed labelling and to every permuted
labelling alike.

**Third, if the observed labelling is one of those, the channel is reported as undefined.** Its
map value is shown as NaN, it has no corrected p, and it is listed as excluded with the reason
"perfect separation, t undefined" (t-tests) or "perfect fit, statistic undefined" (General
Linear Model presets). It stays in every permuted map.

### Why defining the statistic as 0 keeps the test exact

A permutation test is exact for any statistic, as long as the same statistic is computed, by the
same fixed rule, on the observed labelling and on every relabelling. Under the null hypothesis
the labellings are exchangeable, so the observed value is equally likely to sit at any rank among
them, whatever the function is. "Ordinary t, except 0 when the error term is zero" is such a
fixed function. What would break exactness is applying a different rule to the observed map than
to the permuted ones, for example screening only the observed labelling, or dropping a channel
from some permuted maps and not others. The rule avoids both: the channel set is fixed before the
loop, and the zero definition is applied everywhere.

It was also checked by simulation on the data that cause the problem: null count-like data with
about 78% zeros, 40 channels almost always zero, on the EGI 178-channel montage, 199
permutations and 1000 simulated data sets per design. Family-wise error at a nominal 0.05 was
0.040 to 0.064 for TFCE and 0.011 to 0.039 for the cluster test across paired (6, 8 and 12
subjects), unpaired (6 against 6, 8 against 10) and repeated-measures designs. The alternative
considered, screening out up front every channel where any relabelling has a zero denominator,
gave the same error rates and cannot be written for the General Linear Model tier. Without the
rule, up to 36% of those data sets (6 against 6, unpaired) produced a non-finite or absurdly large
statistic somewhere in the 199 permutations, and the share grows with the permutation count.

### Why perfect separation is reported untested

If the observed data separate perfectly (every value identical within each group, or every
paired difference identical and non-zero), the true t is infinite. Scoring it as the strongest
possible effect would need an arbitrary ceiling value fed into TFCE, and that ceiling would set
the corrected threshold for the whole montage. Leaving the channel untested is what the zero
rule gives, it is valid, and it was approved on 14 October 2026. In practice it needs very small
samples or integer data. A channel listed as "perfect separation" may well carry a large real
effect, so look at its raw values rather than reading the NaN as "no effect".

### When a value counts as zero

An error sum of squares counts as zero when it is at most 10⁻²⁰ times the sum of the squared
values the statistic is computed from on that channel, which is a relative standard deviation of
10⁻¹⁰. The squared values are not centred first, so an all-zero channel qualifies. Rounding
leftovers sit near 10⁻¹⁶ in relative standard deviation and real data many orders of magnitude
above, so the line is far from both. A test for exact zero would not do: the 0.1-against-0.3
case above has a pooled sum of squares of about 10⁻³³, not 0. The same tolerance is used in
every tier and in the Python port.

### The reasons you will see

| Analysis | Condition | Reason listed |
|---|---|---|
| `pairedT` | every paired difference is zero | identical values in both conditions |
| `onesampleT` | every value is zero | every value is zero |
| `unpairedT` | the same value in every subject of both groups | same value in every subject |
| `pairedT`, `onesampleT`, `unpairedT` | the observed data separate perfectly | perfect separation, t undefined |
| `correlationP`, `correlationS` | a column is constant across subjects (on the ranks for Spearman) | constant across subjects |
| correlation with covariates | a column is constant once the covariates are removed | fully explained by the covariates |
| `anova1`, `ancova`, `regression`, `rmanova`, `mixed2way`, and each post-hoc and simple effect | no variation left after that test's nuisance model | no variation left after the nuisance model |
| the same | the observed fit is perfect | perfect fit, statistic undefined |

The strings are the same in the Python port. A constant single behavioural score in a
correlation would exclude every channel, so it is refused earlier with its own error
(`snpm:behavior:constant`); if a near-constant score gets past that check, every channel is then
excluded here and the run stops with `snpm:noEvaluableChannels`.

### Designs with no error degrees of freedom

The rule above is about the data at one channel. A design can also leave no error at all,
whatever the data: one subject per group in `anova1`, or a repeated-measures table in which the
subject and condition columns already fit every row. Every channel would then hit the zero rule
and be mislabelled "perfect fit", so the General Linear Model presets refuse such a design before
computing anything, with `core_snpm_glm:noErrorDF` and the message "no error degrees of freedom
left after fitting the design". The check counts rows minus the rank of the design matrix, not
minus its number of columns. For `mixed2way` with one subject per group, all three effects have
zero error degrees of freedom together and the whole run is refused. Inside an otherwise valid
`mixed2way`, the condition effect within a group with a single subject is not run; the report says
"simple effect: condition within A was not run: only 1 subject in this group; a within-group test
needs at least 2". A post-hoc pair with too few rows is skipped with its own message, "too few
rows for these two levels to leave any error degrees of freedom".

### Whole-head tests are parametric, so a perfect fit is reported as undefined

The whole-head global test on the channel average is a parametric F or t test, not a permutation
test. Defining its statistic as 0 with p = 1 would be wrong there: a parametric p claims a position
on the F distribution, and with a perfect fit the true F is infinite, so p = 1 would report the
strongest possible effect as no evidence. A perfectly fitted whole-head average is therefore
reported with statistic and p both NaN and the reason "perfect fit, statistic undefined". The
global card is still shown, with the reason on it, so a missing number is not read as "no
whole-head test was run". In the `.mat` the reason is in `results_struct.effects(k).global_reason`.

### Very large statistics that are not degenerate

A design with one error degree of freedom is valid, and its t or F statistic has very heavy tails:
with one degree of freedom t follows a Cauchy distribution, and a permuted value of 9 × 10⁶ is an
ordinary draw from it, not rounding. Its error sum of squares is small but far above the
tolerance, so the zero rule rightly leaves it alone. The problem was numerical only: TFCE's
threshold ladder grew with the largest value, and that one value asked for about 130 GB. Since
14 October 2026 TFCE computes the same integral without the ladder when the ladder would need
more than 10 000 levels, so its memory and time no longer depend on the size of the statistic;
below that the original computation runs unchanged. Such a channel still dominates the
maximum-statistic null and costs power across the montage, which is valid but worth knowing when
a design has only one or two error degrees of freedom.

## A missing behavioural value removes the subject, not a channel

The correlation analyses also accept Data 2 as a single behavioural score per subject (a
`Subject` column plus one measure column), which is broadcast across every channel of Data 1.
Under the complete-column rule a broadcast `NaN` would make *every* channel incomplete at once —
one missing reaction time would silently empty the whole analysis. So a subject with a missing
behavioural value is instead dropped whole, before anything else happens: the run warns
(`snpm:behaviorNaNRows`) naming the subjects, counts them in
`results_struct.subject_matching.behavior_nan_removed`, and states the drop in the results text.
The remaining subjects are analysed on every channel that is complete for *them*. As with
excluded channels, the removal is never silent — if the analysed n in the report is smaller than
the number of rows in your behavioural file, this is where the difference went.

## What you will see in your outputs

A channel never leaves an analysis silently, whichever rule dropped it. When anything is
excluded, the run produces all of the following:

- A console line and a MATLAB warning (`snpm:channelsExcluded`) naming the count and the
  individual channel labels.
- An **`excludedChannels`** sheet in the `.xlsx`, one row per dropped channel with its index,
  its label, and its own reason (missing data, or one of the no-variation reasons above). In the
  General Linear Model presets each test (omnibus, post-hoc pair, `mixed2way` effect, simple
  effect) has its own rows, and a column says which test a row belongs to.
- An **Excluded channels** entry in the amber notice block under the HTML report's model card,
  stating how many of how many channels were not tested and listing them.
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
survive this step. How the remaining subjects' covariates are then controlled, and why that
permutation is approximate, is explained in
[Covariates in a correlation](covariates-in-correlation.md).

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
  not the dropped cells were ever present. The no-variation rule's tolerance lives in one
  function, `dependencies/snpm_zero_ss.m`, and the t-test form of the zero rule in
  `dependencies/snpm_ttest_novar.m`.
