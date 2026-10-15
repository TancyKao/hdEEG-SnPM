# Two-way mixed ANOVA (`mixed2way`) — how the interaction is computed in hd-EEG

> *Understanding-oriented.* For the factual `mixed2way` row (inputs, roles, outputs) see the
> [Analysis catalog](../reference/ANALYSIS_CATALOG.md); to run it, see the
> [how-to guides](../how-to/README.md).


Scope: the `mixed2way` GLM preset in `core_snpm_glm.m`, with the design built in
`dependencies/snpm_glm_design.m`, the statistic in `dependencies/snpm_glm_stat.m`,
and the permutation in `dependencies/snpm_glm_permute.m`. One GLM is fit **per
channel**; family-wise error is controlled across channels by max-statistic
permutation with TFCE and cluster-mass.

## The model (one GLM per channel)

For a 2×2 mixed design — one between factor (Group), one within factor
(Condition), subjects nested in Group — the design built by `snpm_glm_design.m` for the
interaction is

```
X = [ 1 | Ds | Dc | Dx ]
```

- `Ds` — reference-coded **subject** dummies (n_subj − 1 columns). These
  per-subject intercepts absorb both the subject random effect *and* the
  between-group main effect, because Group is constant within a subject
  (subjects nested in Group).
- `Dc` — reference-coded **condition** dummies (k_cond − 1 columns).
- `Dx` — the **interaction** columns: elementwise products of every Group dummy ×
  every Condition dummy (the `interaction_cols` helper). For 2×2 that is 1 column;
  for g×c levels it is (g−1)(c−1) columns.

The interaction lives entirely in `Dx`; everything else (`1, Ds, Dc`) is nuisance
for that test.

## The interaction statistic

The contrast `C` selects only the `Dx` columns. More than one selected column → **F**; a single interaction column
(2×2) → still routed as F with df1 = 1 (numerically t²).

`snpm_glm_stat.m` fits Y = Xβ + e by OLS across all channels at once (batched
`pinv`), then for the multi-row contrast:

```
CB = C*B                         % interaction betas
M  = C (XᵀX)⁻¹ Cᵀ
F  = ( CBᵀ M⁻¹ CB / q ) / σ̂²     per channel
```

with q = df1 = (g−1)(c−1), df2 = n − rank(X), σ̂² = RSS/df2. This is the standard
reduced-vs-full SS-interaction test, vectorized over channels rather than looped.

## Why this is the interaction, concretely

The interaction β's are the Group×Condition products; testing them jointly asks
"does the condition effect differ across groups." Because subject intercepts
absorb the group main effect, the F is *purely* the interaction — it cannot
accidentally pick up a group difference.

The other side of that: the **between-group main effect cannot be estimated in
this design**, because the subject columns absorb it. The `mixed2way` report
therefore computes the group effect separately, as a one-way test on one mean per
subject (the plain mean of that subject's rows over all conditions), with whole
subjects relabelled in the permutation. Asking `snpm_glm_design` itself for
`effect='group'` is an error that says the same thing.

## The three effects in the report

`mixed2way` always reports all three effects of the design, and the whole-head
tests, simple effects and post-hoc comparisons that go with them.

| Effect | Model | What is permuted |
|---|---|---|
| Group | one-way design on subject means | whole subjects are relabelled |
| Condition | `[1, Ds, Dc, Dxc]`, testing `Dc` | residuals within each subject |
| Interaction | `[1, Ds, Dc, Dx]`, testing `Dx` | residuals within each subject |

The group effect is tested against variation between subjects, with
(subjects − groups) error degrees of freedom. The condition and interaction
effects are both tested against the within-subject residual of the full model,
with (subjects − groups) × (conditions − 1) error degrees of freedom. All effects
use one channel set, the channels that have a value in every row of the table (see
[Missing data](missing-data-and-excluded-channels.md)). On top of that, each test
leaves out the channels that have no variation left in that test, and lists them
with a reason (see
[Channels with no variation](missing-data-and-excluded-channels.md#channels-with-no-variation)),
so one section can list a channel another section tested.

With one subject in each group, all three effects have zero error degrees of
freedom together, and the run is refused before anything is computed
(`core_snpm_glm:noErrorDF`, "no error degrees of freedom left after fitting the
design").

### The condition main effect

The condition effect is tested in the full model, with the group-by-condition
columns present as nuisance terms. Two details make that work.

First, the group-by-condition columns `Dxc` are built from a group code that sums
to zero across groups (the first group's rows are set to −1 in every group
column), not from the 0/1 reference code used for the interaction test. With the
0/1 code and the interaction in the model, the `Dc` coefficients would be the
condition effect in the reference group only. With the sum-to-zero code they are
the condition effect averaged over groups, which is what "main effect" means.

Second, because the interaction columns are in the nuisance model, a real
interaction is fitted and removed before the residuals are permuted. Until 14
October 2026 the condition effect was tested in the reduced model `[1, Ds, Dc]`,
without the interaction columns. That gave the pooled condition effect, but any
real interaction stayed in the residuals that Freedman–Lane permutes, and that
could hide the condition effect completely. In simulation with a crossover interaction, a planted
condition effect was detected 0.00 of the time under the reduced model and 0.93
of the time under the full model.

The condition F from the full model is the F of a standard mixed analysis of
variance. On a check table of 8 + 14 subjects and 3 conditions it reproduces
MATLAB's `ranova`: F(2,40) = 30.9004. When the groups differ in size, the
condition effect is averaged over groups with equal weight per group, not
weighted by the number of subjects in each group. This is the default of standard
software. The two weightings agree when the groups are the same size.

The statistic and its degrees of freedom differ from runs made before 14 October
2026, so a condition main effect reported from an earlier run should be re-run.

### Post-hoc comparisons and simple effects

With three or more conditions (or groups) the report adds pairwise post-hoc maps.
Each pair is a separate two-level test on the rows of those two levels: group
pairs on the subject means of the two groups, condition pairs with the condition
model above on the rows of the two conditions. The other levels are neither
fitted nor permuted, so their effects cannot enter the comparison. The simple
effects are the condition effect within each group, computed on that group's rows
alone, and their post-hoc comparisons follow the same rule. A simple effect needs
at least 2 subjects in its group; a group with one subject has no error term for
it, so that simple effect is not run and the report says why ("simple effect:
condition within A was not run: only 1 subject in this group; a within-group test
needs at least 2"). Post-hoc maps are
corrected across channels and not across pairs. The full rules are in the
[Analysis catalog](../reference/ANALYSIS_CATALOG.md).

## From one channel to the scalp — the correction

Parametric p from the F is used only for cluster *forming*/reporting. Family-wise
error control is by **permutation of the max statistic**, using **Freedman–Lane**
(`snpm_glm_permute.m`, Winkler et al. 2014) because nuisance terms are present:

1. Fit Y on the nuisance Z = `[1, Ds, Dc]`; keep fitted `Zfit` and residuals `R`.
2. Each permutation: shuffle `R`, add back → `Yperm = Zfit + R(perm,:)`.
3. `perm_type = 'within'` with exchangeability block = subject,
   so residuals are shuffled **only within each subject** (`within_block_perm`).
   This respects the repeated-measures structure — it permutes condition labels
   within a subject, the exchangeable unit under H₀ of no interaction, and never
   mixes across subjects or groups.
4. Recompute the F-map per permutation; take the max (for TFCE) and max
   cluster-mass; build the null; get corrected p per channel/cluster as
   (b + 1)/(N + 1), where b counts the permuted maxima at least as large as the
   observed value. A permuted value within a relative 10⁻¹² of the observed one
   counts as a tie (it is the same value up to rounding), which the count used
   to drop about half the time. Neighbour adjacency comes from the montage
   neighbour matrix; TFCE uses E=0.5, H=2.

The three whole-head tests on the channel average are parametric, not permuted.
If the model fits the channel average exactly, that test's statistic and p are
reported as undefined (NaN) with the reason "perfect fit, statistic undefined",
and the global card still appears with that reason.

## Caveats for sleep-EEG use

- **Balance:** every subject must have exactly one row per condition and belong
  to one group. A table that does not meet this stops the run with
  `core_snpm_glm:mixed2wayUnbalanced`, which names the subjects. In the GUI the
  same check runs as soon as the file and the Subject and Condition roles are set:
  Run stays disabled and the "Before you can run" list names the subjects that do
  not have exactly one row per condition, or that appear in more than one group.
  Group sizes may differ; with unequal sizes the condition main effect weights the
  groups equally (see above).
- **Sphericity / covariance:** this is a fixed-effects GLM with subject dummies,
  not a true mixed model. For >2 within levels there is no Greenhouse–Geisser
  correction — but the permutation null sidesteps the parametric sphericity
  assumption for the *inference*, so this matters more for the point estimate
  than for the p.
- If you need random slopes or trial-level nesting, use the `mixedmodel` path
  (`core_snpm_lmm.m`), not this one.
- **Reading a condition main effect next to an interaction.** The condition main
  effect is an average over groups. The test is valid whether or not there is an
  interaction, but when the interaction is strong the average can describe neither
  group well (in a full crossover it can be close to zero while each group has a
  clear condition effect). In that case report the **simple effects** (condition
  within each group) alongside it.

## A worked 2×2 example (one channel)

Four subjects, two per Group (A = control, B = patient), each measured in two
Conditions (Pre, Post) — 8 rows. Power at one channel:

| Subject | Group | Cond | y   | Δ = Post−Pre |
|---------|-------|------|-----|--------------|
| S1      | A     | Pre  | 10  |              |
| S1      | A     | Post | 11  | +1           |
| S2      | A     | Pre  | 12  |              |
| S2      | A     | Post | 15  | +3           |
| S3      | B     | Pre  | 10  |              |
| S3      | B     | Post | 15  | +5           |
| S4      | B     | Pre  | 14  |              |
| S4      | B     | Post | 21  | +7           |

Group-mean within-subject change: Δ̄_A = 2, Δ̄_B = 6.

### The design matrix `X = [1 | Ds | Dc | Dx]`

Reference coding: reference subject S1, reference condition Pre, reference group A
(the group dummy is only used to *form* `Dx`; group itself is absorbed by the
subject dummies). Columns are `Int, S2, S3, S4, Post, Dx`:

| row | Int | S2 | S3 | S4 | Post | Dx |
|-----|-----|----|----|----|------|----|
| S1 Pre  | 1 | 0 | 0 | 0 | 0 | 0 |
| S1 Post | 1 | 0 | 0 | 0 | 1 | 0 |
| S2 Pre  | 1 | 1 | 0 | 0 | 0 | 0 |
| S2 Post | 1 | 1 | 0 | 0 | 1 | 0 |
| S3 Pre  | 1 | 0 | 1 | 0 | 0 | 0 |
| S3 Post | 1 | 0 | 1 | 0 | 1 | 1 |
| S4 Pre  | 1 | 0 | 0 | 1 | 0 | 0 |
| S4 Post | 1 | 0 | 0 | 1 | 1 | 1 |

`Dx = Dg .* Dc` is 1 only for group-B Post rows. p = 6 columns, n = 8, so
error df2 = n − rank(X) = 2.

Interaction contrast: `C = [0 0 0 0 0 1]` → one row → **t/F with df1 = 1**.

### The estimate and F

The interaction β is exactly the **difference-in-differences**:

```
β_Dx = Δ̄_B − Δ̄_A = 6 − 2 = 4
```

The model forces each subject's fitted Post−Pre to equal its group mean Δ̄, so the
residual for each subject is ±(Δ̄_group − Δ_subj)/2:

```
RSS  = Σ (Δ̄_group − Δ_subj)² / 2
     = (2−1)²/2 + (2−3)²/2 + (6−5)²/2 + (6−7)²/2 = 2.0
σ̂²   = RSS / df2 = 2.0 / 2 = 1.0
Var(β_Dx) = 2σ̂² = 2.0        (β_Dx is a difference of two 2-subject means of Δ)
F    = β_Dx² / Var(β_Dx) = 16 / 2 = 8.0     df = (1, 2)
```

Parametric two-sided p ≈ **0.11** — *not* significant, purely because 4 subjects
leave df2 = 2. This is the point: the interaction estimate (Δ-of-Δ = 4) is large
and clean, but the error df is tiny, so a single channel says little on its own.

### What the permutation adds

Across the scalp, the reported p is **not** this parametric 0.11 — it comes from
the within-subject Freedman–Lane null (nuisance Z = `[1, Ds, Post]`, residuals
shuffled within each subject block), taking the max F over channels for TFCE and
max cluster-mass, then TFCE/cluster correction over neighbours. Note the
exchangeability ceiling here: each subject has 2 within-subject rows → 2
arrangements each → only 2⁴ = 16 distinct permutations, so with n = 4 the
smallest achievable corrected p is coarse. Realistic sleep designs (15–30+
subjects) restore both error df and permutation resolution.

## References

- Maris & Oostenveld (2007) — cluster/TFCE permutation for M/EEG.
- Nichols & Holmes (2001) — max-statistic permutation, FWE control.
- Winkler et al. (2014), NeuroImage — Freedman–Lane permutation with nuisance
  and exchangeability blocks.
