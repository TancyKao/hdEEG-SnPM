# Covariates in a correlation: the Freedman–Lane permutation

This page explains how `correlationP` and `correlationS` control for covariates, why the
permutation scheme changed in October 2026, and what the partial Spearman correlation does and
does not control. For the parameters and the output fields, see the
[analysis catalog](../reference/ANALYSIS_CATALOG.md); for missing covariate values, see
[Missing data and the complete-column rule](missing-data-and-excluded-channels.md).

## What is computed

With covariates, the statistic at each channel is the partial correlation between the EEG value
and the behavioural measure, given the covariates. Both measures are regressed on
`[1, covariates]` and the residuals are correlated. The number in the map is identical to
MATLAB's `partialcorr`, and that has not changed.

What changed is the permutation null against which that number is scored.

## Why the old null was too narrow

The old scheme permuted the EEG residuals against the fixed behavioural residuals and
correlated them. That looks natural, but a reordered residual vector is no longer orthogonal to
the covariates. Its correlation with the behavioural residuals then has a smaller spread than
the observed partial correlation does under the null hypothesis, so the observed value looks
more extreme than it is. The effect grows with the number of covariates and with how much of
the data they explain. In simulation on the 178-channel EGI montage, with covariates that drive
both measures, the family-wise false-positive rate at a nominal 0.05 was:

| Subjects, covariates | Max statistic | TFCE | Cluster |
|---|---|---|---|
| 10, 4 | 0.856 | 0.503 | 0.235 |
| 12, 2 | 0.191 | 0.126 | 0.094 |
| 20, 3 | 0.140 | 0.118 | 0.080 |
| 40, 5 | 0.110 | 0.093 | 0.070 |

## What the toolbox does now

Each permuted map reorders the EEG residuals and then projects them off the covariates again
before correlating them with the behavioural residuals. This is the Freedman–Lane procedure
(Freedman & Lane 1983), the scheme Winkler et al. (2014) recommend for permutation inference
with nuisance variables, and the same one the toolbox's regression and ANCOVA presets use. On
the same simulations, at a nominal 0.05, the family-wise error of the per-channel maximum
statistic was about 0.04 to 0.06, and TFCE and cluster inference were at or below the nominal
rate (cluster rates of 0.035 to 0.038 in several settings, so slightly conservative).

On a single behavioural score, the maximum-statistic corrected p-values are exactly those of the
`regression` preset run on the same permutation orders, because r and t rank every permuted map
in the same order. The TFCE and cluster p-values are not the same as the preset's: TFCE
integrates the r map rather than the t map, and the cluster-forming threshold is set on a
different scale, so those two corrections can give somewhat different answers.

Two properties are worth knowing.

The test is **approximate, not exact.** Residuals are only approximately exchangeable, so even
when all orderings of a small sample are enumerated the p-value is not exact. The report says
"Freedman–Lane, approximate".

**With few error degrees of freedom the approximation is rougher, in both directions.** With
`n − 2 − k` below 5 (n subjects, k covariates), the per-channel maximum statistic may be mildly
liberal: its family-wise error was 0.056 at a nominal 0.05 with 3 error degrees of freedom.
TFCE and cluster inference go the other way and reject less often than the nominal rate (about
0.02 at 4 error degrees of freedom), so a non-significant TFCE or cluster map in that regime is
weak evidence of no effect. The run warns, and borderline results there should be read with
caution. With no error degrees of freedom the run stops.

The covariate count `k` used in every degrees-of-freedom statement is the rank of the covariate
design, so a covariate that duplicates another, or is a combination of others, is not counted
twice.

## What partial Spearman controls

`correlationS` with covariates ranks the data and the covariates, then removes the ranked
covariates by linear regression and correlates the rank residuals. That is what
`partialcorr(..., 'Type', 'Spearman')` computes. The null hypothesis it tests is therefore a
statement about ranks: after the ranked covariates have been removed linearly, the rank
residuals of the EEG value and of the behavioural measure have no linear association. It is not
a test that the two measures are unrelated once the covariates are known, in the general sense,
and a non-significant result does not show that.

In practice it controls a covariate only as far as the covariate's relation to the *ranks* is
linear. A covariate whose effect is linear on the raw scale is not fully removed: with three
such covariates the per-channel maximum-statistic family-wise error was 0.10 to 0.13 in
simulation, whatever permutation scheme is used, because the leftover covariate effect is a
real association, not a permutation artefact. The report states this caveat on every partial
Spearman run. When covariate control is the point of the analysis, use `correlationP` or the
`regression` preset.

## References

- Freedman, D. & Lane, D. (1983). A nonstochastic interpretation of reported significance
  levels. *Journal of Business & Economic Statistics*, 1(4), 292–298.
- Winkler, A. M., Ridgway, G. R., Webster, M. A., Smith, S. M. & Nichols, T. E. (2014).
  Permutation inference for the general linear model. *NeuroImage*, 92, 381–397.
