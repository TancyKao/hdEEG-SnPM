# About circular statistics for phase

> *Understanding-oriented.* Why phase needs its own machinery, what the circular tests do and
> do not license you to claim, and the honest limits of what they can detect at realistic sample
> sizes. To actually run one, see
> [Run a circular (phase) analysis](../how-to/run-a-circular-phase-analysis.md); for inputs and
> outputs, see the [Analysis catalog](../reference/ANALYSIS_CATALOG.md).

## Read this first: you probably do not need a circular test

Phase–amplitude coupling is usually described by three quantities, and only one of them is
circular.

| Quantity | What it is | Kind of number | Analyse it with |
|---|---|---|---|
| **Coupling strength** | how consistently events cluster at *some* phase (mean vector length, modulation index) | ordinary bounded number on [0,1] | `unpairedT`, `anova1`, `regression`, `ancova` … |
| **Prevalence** | how many coupled events there were (count, rate per minute) | ordinary count | the same presets |
| **Preferred phase** | *where on the slow oscillation* the events sit | an angle | the circular tier |

Strength and prevalence are linear quantities. They need no circular analysis at all, they go
through the existing presets today, and they are where the published between-group effects in
this literature predominantly sit. If your question is "is coupling weaker in the patient
group?" or "do patients have fewer coupled spindles?", you are done — use `unpairedT` or
`anova1` and stop reading.

The circular tier exists for the narrower question "do the two groups couple at a *different
point* of the slow oscillation?". That question is real, but it is harder to answer, lower
powered, and easier to over-interpret than the two linear ones. Reach for it deliberately, not
by default.

## Why an angle is not a number

The arithmetic mean of 350° and 10° is 180°. Those two angles are 20° apart and sit either side
of zero, so their true mean direction is 0° — but ordinary averaging puts the answer on the
opposite side of the circle. Every linear operation inherits some version of this: a standard
deviation, a t-test, a z-score, a regression slope all assume the number line runs from −∞ to
+∞ and that 359 is far from 1. On a circle it is not.

This is why the circular analyses refuse `datatype = logscale` and `datatype = normalize`. The
log of a radian is not a quantity (and goes complex for negative angles), and a z-score of an
angle rescales a bounded periodic domain into an unbounded one, which destroys the wrap. Only
`absolute` is accepted, and the tool errors rather than silently ignoring the request.

The circular statistics avoid the problem by embedding each angle as a point on the unit circle,
(cos α, sin α), and working in that two-dimensional space, where ordinary vector arithmetic is
legitimate again.

## The linearised measure: the recommended first analysis

Before reaching for a two-sample circular test, consider turning the angle into a linear
quantity you can feed to the presets you already understand. That is what `snpm_circ_linearise`
does, following Hahn et al. (2020):

- **Signed deviation** — `circ_dist(phase, θ_ref)` on (−π, π], where θ_ref is the *pooled* grand
  mean direction over all subjects in both groups at that channel. Interpretation: this
  subject's phase lead or lag relative to the sample as a whole. Because the reference is pooled,
  it does not move when labels are permuted, which is what keeps the transform exact under
  permutation.
- **Unsigned distance** — `abs(circ_dist(phase, 0))` on [0, π], distance from the up-state peak
  regardless of direction.

Two things make the signed version attractive. First, it runs through `unpairedT`, `anova1`,
`ancova`, `regression`, `rmanova` and `mixed2way` unchanged, so you get covariates, more than
two groups, repeated measures and interactions for free. Second, and more important, **it is
directional**: a positive group difference means one group's phase genuinely leads the other's.
Neither of the two-sample circular tests, nor circular–linear correlation, can tell you that.

The price is an approximation — a bounded circular quantity treated as an unbounded linear one.
That is only defensible when the data are concentrated, so the transform drops any channel where
a subject sits more than 150° from the pooled mean, or where the pooled resultant length falls
below 0.3.

The **unsigned** measure carries an extra hazard and is guarded accordingly. It is a folded,
dispersion-like quantity: a subject whose angle was estimated imprecisely has a noisier angle,
and noise can only push `|deviation|` *up*, never below zero. Any group difference in estimate
precision therefore manufactures a group difference in the measure. The tool stamps the measure's
provenance into the CSV and hard-errors
(`core_snpm:circUnsignedNeedsCountCovariate`) if you point it at a preset with no covariate slot.
Only `ancova` and `regression` — with the precision column actually listed in `covariate_cols` —
are accepted.

## What the two-sample circular tests actually test

Two statistics are available for comparing two groups' phase distributions.

**Hotelling T² on the (cos, sin) embedding** (`circ_phase_group`, the primary) tests whether the
two groups' *mean resultant vectors* differ. Because the vector has both a direction and a
length, this responds to a mean-direction shift, to a concentration difference, or to both. It
takes covariates, which is what makes it the primary: it is the only one of the two that can
adjust for estimate precision. The statistic is F(2, N − q − 3).

**Watson's U²** (`circ_phase_group_u2`, the secondary) is a distance between the two empirical
cumulative distribution functions on the circle — an omnibus test over the whole distribution,
sensitive to shape differences the mean vector misses. It has no covariate slot at all, so it
cannot adjust for precision. Use it as a secondary, corroborating view.

**What a significant cluster licenses you to say.** It says *the two groups' phase distributions
differ somewhere in this region*. That difference may be a shift in mean direction, a difference
in concentration, or a mixture of the two — the test does not decompose them. It does **not**
establish that one group's phase leads the other's, and it does not establish a direction of
effect. If you want a directional claim ("the patient group couples later in the up-state"), the
signed linearised measure is the analysis that supports it, and you should run that.

U² needs one piece of special handling worth knowing about, because it is a good illustration of
why an "exact" workaround can still be wrong. U² spans roughly [0, 0.5], while the TFCE
integration step is tuned for t- and F-scale maps and defaults to 0.1 — about five integration
levels, which flattens most channels to exactly zero. The engine now passes the integration step
`dh = 0.005` straight through to the enhancer. It previously multiplied the map by 20, enhanced
at the default step, and divided back out, which is exactly equivalent in real arithmetic but not
in floating point: at one integration level a channel whose U² was bit-equal to the threshold
0.095 was correctly included by the direct call, but `20 × 0.095` rounds to 1.8999999999999999,
just under 1.9, so the workaround dropped it. One flipped inclusion decision in 23,187, which
then propagated as a small additive offset across the 69-channel cluster it belonged to. No
p-value, cluster membership or significant-channel set changed — but the multiply was introducing
an artefact rather than preventing one, which is why it was removed rather than tidied.

**Circular–linear correlation** (`circ_corrAngLinear`) is non-directional in the same way, and
for the same underlying reason (Mardia 1976). A significant channel means the behavioural
variable varies systematically with phase. It never means the variable is *better* near the
up-state, or worse near the down-state, or higher at any particular angle. Reporting it that way
is a common and serious over-read.

## The honest power figure

At n = 20 per group with about 30 events per subject, a Hahn-magnitude concentration difference
(von Mises κ 14.2 versus 3.94) is detected roughly **22% of the time**. That is an exploratory
sensitivity, not a confirmatory one. A null result at that sample size tells you almost nothing;
a positive result deserves replication before it is believed.

Worse, when precision is confounded with group, adjusting for it removes most of the apparent
effect — *including the part that is real*. In simulation at a point-biserial correlation of
about 0.71 between the log precision covariate and group, a genuine concentration effect went
from a rejection rate of 0.607 unadjusted to **0.097** adjusted, essentially the null rate.
Nothing is broken: the adjusted p-value is valid. The study simply cannot separate the effect
from the confound at that sample size, because at n = 20 the two are near-inseparable.

Design accordingly. See these numbers before you plan a study, not after you run one.

## The precision confound, and why a precision covariate is required

A subject's preferred phase is not measured, it is *estimated* from that subject's detected
events. Estimate it from 50 tightly clustered spindles and it is precise; estimate it from 10
scattered ones and it is noisy. In the (cos, sin) embedding, noise shortens the mean resultant
vector — and a shorter resultant vector **is** a mean shift in that space. So a group whose
per-subject angles were simply estimated less precisely looks different from a group whose were
estimated well, with no difference whatsoever in the underlying phase.

This is not a small effect. In simulation with **identical true phase distributions** and only
the estimate precision differing (here driven by event count, 50 versus 10 or 20 per subject):

| Precision contrast | Type-I rate, unadjusted | Type-I rate, adjusted |
|---|---|---|
| 50 vs 10 events | **0.20 to 0.57** (see below) | 0.031 |
| 50 vs 20 events | 0.205 | 0.039 |

The unadjusted figure at 50-vs-10 is quoted as a **range** on purpose. Two independent
dispersion models gave 0.234 and 0.565 for nominally the same scenario, and the gap turned out
to be informative rather than a contradiction: **the inflation scales inversely with
within-subject coupling strength.** The weaker the within-subject concentration, the noisier each
subject's circular mean and the worse the confound. Since the meta-analytic pooled mean vector
length in this literature is only 0.23 to 0.33 — which is weak — real data sits at the severe end
of that range. The unadjusted false-positive rate can exceed 0.5: a **ten-fold** inflation over
the nominal 0.05.

What matters just as much is that **both dispersion models adjust back to nominal** (0.031 and
0.039). That agreement is the evidence the correction works regardless of which model is right.

### The covariate is log(Rayleigh Z)

The precision covariate is **log(Rayleigh Z)**, where Z = nR² for a subject-channel cell with
*n* events and resultant length *R*.

Z is the right quantity because it subsumes event count entirely: it captures both *how many*
events there were and *how tightly* they clustered, which are the two things that jointly
determine how well that subject's angle is pinned down. Event count alone sees only half of it.
TurtleWave's phase–amplitude coupling module already exports `rayleigh_z` per channel, so using
it costs nothing at the export stage.

**log(event count) is retained only as a fallback** for exports where Z is unavailable. The two
must **never both appear in one model**: they are strongly collinear by construction, and the
second adds nothing but instability.

The covariate is fitted **per channel**, not as a per-subject scalar, because coupling precision
varies strongly across the montage and the confound acts channel by channel. A single
per-subject number under-adjusts exactly at the posterior channels where the problem is worst.

Watson's U² has no covariate slot. When precision is associated with group, the U² map is
unadjusted and the engine warns (`core_snpm:circU2Unadjusted`). Treat it accordingly.

### Conditioning on Z is legitimate, not circular

A careful reader will suspect circularity: Z is computed from the same events as the phase, so
does conditioning on it bias the phase comparison?

It does not, and the reason is a clean classical result. For von Mises data the sample **mean
direction** and the sample **resultant length** are independent, and the distribution of R
depends on the concentration parameter but not on the mean direction (Mardia & Jupp,
*Directional Statistics*, §4.5). Z is therefore **ancillary for the mean direction** — it carries
information about how precisely the direction was estimated, and none about where it points.
Under Fisher's conditionality principle, conditioning on an ancillary statistic is not merely
permissible, it is the textbook-correct move: it makes the inference about the parameter of
interest at the precision actually achieved, rather than averaged over precisions that did not
occur.

### One hard guard: Z is a covariate or an outcome, never both

Z is a **covariate** when the outcome is *phase*. Z is the **outcome** when the analysis is about
*coupling strength* — and coupling-strength analyses are linear analyses that belong in the
ordinary presets, not here. Putting Z on both sides of one model regresses a variable on itself.
The engine rejects that combination.

### When the adjustment has eaten the question

The engine aborts when |point-biserial(log precision, group)| exceeds 0.8
(`circ_collinearity_max`), because past that point the covariate and the group contrast are the
same variable. But the analysis stops being informative well *below* the gate. At |r| ≈ 0.71 —
which passes the gate and produces a report with no error at all — adjustment took a genuine
effect from 0.607 down to 0.097. **A run that completes cleanly at |r| = 0.71 looks fine and is
not.** Read the point-biserial diagnostic the engine prints before you read any p-value.

This diagnostic matters **more** now that the covariate is Z than it did when it was event count.
Z contains coupling strength, and coupling strength is frequently the very thing that differs
between groups — Helfrich found it reduced in ageing, Hahn found it increased with development.
So a high correlation between log(Z) and group is not an accident of data quality; it is often
the biology. Where that correlation is high, phase and strength are **not separable in that
cohort**, and the right response is to report the two side by side rather than presenting one
adjusted for the other. That is what the field does anyway.

## What whole-column deletion does to your montage

A channel is analysed only if **every** subject has an angle there. That is a deliberate choice —
the Freedman–Lane nuisance fit works on whole columns, so a single missing cell would turn the
column all-NaN inside every permutation while leaving the observed map finite, and the channel
would be scored against a null it could never have entered. Dropping the whole column is the
honest alternative.

Note carefully what "missing" means here: an angle that is **genuinely absent** — no
phase–amplitude coupling row for that subject-channel cell, or a channel rejected upstream as
bad. It does **not** mean a weakly coupled cell (see the next section: there is no Rayleigh
mask). On the real data in this project, all 257 channels produce output for every subject, so
this is a structural caution about upstream bad-channel rejection rather than a live problem.

The arithmetic still deserves respect, because survival goes as (1 − q)^N and compounds fast:

| Per-cell absence rate *q* | Subjects *N* | Channels surviving |
|---|---|---|
| 3% | 40 | 0.97⁴⁰ ≈ **30%** (70% of the montage lost) |
| 2% | 40 | 0.98⁴⁰ ≈ **45%** (55% lost) |

A 3% per-cell absence rate sounds negligible and costs you seven channels in ten.

**But the survival number is the lesser problem.** A survival *number* is a power question: if
the dropped channels are a random subset of the montage, you lose sensitivity and nothing else. A
survival *pattern* is a validity question. Upstream channel rejection is not random across the
head — it concentrates where signal quality is poor, which in this literature means posterior
sites. When it does, the surviving set is enriched for exactly the frontal region where a
slow-oscillation/spindle effect is expected, and the analysis **can no longer distinguish "a
frontal effect" from "frontal is all that was measured."**

That distinction is not academic. Helfrich's published topographic result depends precisely on
the parieto-occipital estimates *not* differing between groups — a claim that is unavailable to
you if those channels were never tested.

The excluded-channel list is how you check this. It appears as a banner in the HTML report and as
an `excludedChannels` sheet in the Excel output. **Read it before you look at the topoplot, not
after**, and look at *where* the exclusions are, not just how many.

## Why there is no Rayleigh mask

The intuitive move is to drop any subject-channel cell whose Rayleigh test is not significant, on
the grounds that a non-significant cell has "no real coupling" and so no meaningful phase. That
intuition is wrong, and the simulation is unambiguous, so it is worth setting out.

**The Rayleigh test is not a "do we have enough data" test.** It is a fixed threshold on estimate
precision, and it is essentially invariant to event count. At the α = 0.05 boundary, the
per-subject circular mean has a circular standard deviation of:

| Events in the cell | Circular SD of the mean at the α = 0.05 boundary |
|---|---|
| 50 | 32.6° |
| 200 | 33.1° |
| 418 | 33.1° |

The threshold sits at roughly **33° of estimate error regardless of n**. More events do not buy
you a lower bar; they buy you a better chance of clearing the same bar.

Now compare 33° against the published *between-subject* spreads this analysis is trying to
detect: **15.5°** (Helfrich, young adults), **31.2°** (Helfrich, older adults), **49.7°**
elsewhere. A cell sitting exactly at the significance boundary is measured with an error
comparable to — and in one published case *smaller than* — the between-subject signal itself.
Masking at α = 0.05 therefore deletes usable data.

The general principle is worth stating because it generalises far beyond this analysis:
**inclusion costs variance, deletion costs bias.** A noisy but unbiased per-subject estimate
still carries information about the group mean; it carries it with more variance, symmetrically
and honestly, and the covariate adjustment is there precisely to account for that. Deleting it is
**selective**: it removes weakly-coupled cells non-randomly, and weak coupling is spatially
structured across the head. So masking trades a variance cost you can model for a bias you
cannot.

The Rayleigh p-value is still worth exporting. It populates the descriptive quality panel, and it
supports an optional, explicitly-declared sensitivity analysis if you want to show your result
survives a stricter subset. It does not gate the analysis.

## Two places where we extend the literature rather than follow it

Both of these are defensible. Both are novel enough that a reviewer will ask, so they are stated
plainly here rather than buried.

**No published instance exists of a two-sample circular test statistic feeding a cluster or TFCE
maximum-statistic null.** Helfrich et al. (2018) compared groups' phase at a *single electrode*
(Cz) with a Watson–Williams test; the cluster-corrected topographic analysis in that paper was
circular–**linear** correlation, which we do follow directly. Maris & Oostenveld's framework is
deliberately statistic-agnostic — any per-channel statistic with a permutation null can be
plugged in — so extending it to Hotelling T² or Watson's U² is legitimate. It is nonetheless an
extension, not a citation.

**Our cluster statistic is mass, where Helfrich used FieldTrip's `maxsize`.** Cluster mass sums
the statistic over the cluster; `maxsize` counts supra-threshold channels. Mass rewards height,
size rewards extent. Neither is wrong, but the numbers are not comparable, and the disclosure is
written into the results text of every circular run.

## Why the parametric p is descriptive only

Each circular statistic has an asymptotic parametric p-value, and the engine computes it — but it
is used **only to form clusters** and to populate the descriptive panel. It is never the
inference. All reported significance comes from the permutation null.

The reason to be cautious is small samples. Enumerating the exact permutation null at n = 8 per
group gives rejection rates of 0.1553 / 0.1787 / 0.2412 at nominal α of .10 / .05 / .01, against
the asymptotic 0.152 / 0.187 / 0.268. Agreement is good in the body of the distribution and
conservative in the far tail — which is reassuring, and is exactly why the parametric value is
safe as a cluster-forming threshold. It is not accurate enough to report as a result.

## Why the tails control is locked

The Hotelling F, Watson's U² and the circular–linear F are all **non-negative omnibus
quantities**. There is no left tail and no signed direction: the statistic gets large when the
distributions differ, whichever way they differ. A one-tailed request is not a stricter or
weaker version of these tests, it is undefined. The engine therefore errors
(`core_snpm:circTailNotSupported`) rather than quietly accepting `left` or `right`, and the GUI
Tails control is disabled for these analyses. This is also the reason a circular significance map
carries no direction of effect — see
[Interpreting cluster results](interpreting-cluster-results.md).

## Between-subject now, within-subject later

Everything in this tier compares **two independent groups of subjects**, or correlates an angle
with a per-subject measure. The permutation scheme relabels whole subjects, which is only a valid
null when subjects belong to one group and one group only.

There is a trap here that is easy to fall into. **If your two angle files contain the same
subjects, this is not the right test yet.** Two files of the same people in two conditions is a
within-subject design; relabelling subjects across those files does not generate a null for a
within-subject effect, and the resulting p-values would be meaningless. A within-subject circular
comparison needs a paired scheme that is not yet implemented. Until it is, either restructure the
question as a between-subject one, or use the signed linearised measure with `pairedT` or
`rmanova`, which handles the within-subject case correctly through the existing machinery.

## Further reading

- [Run a circular (phase) analysis](../how-to/run-a-circular-phase-analysis.md) — the task guide,
  including the detector convention table.
- [Interpreting cluster results](interpreting-cluster-results.md) — what a significant cluster
  licenses in general.
- [Choosing an analysis](choosing-an-analysis.md) — where the circular tier sits among the presets.
- [Analysis catalog](../reference/ANALYSIS_CATALOG.md) — inputs, statistics, outputs.
- Hahn et al. (2020) — the linearised phase measure and the concentration effect size.
- Helfrich et al. (2018) — single-electrode group phase comparison; cluster-corrected
  circular–linear correlation.
- Mardia (1976) — circular–linear correlation.
- Maris & Oostenveld (2007), *J Neurosci Methods* 164(1), 177–190; Nichols & Holmes (2001),
  *Hum Brain Mapp* 15(1), 1–25 — the permutation framework this tier plugs into.
