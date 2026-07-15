# Two-way mixed ANOVA (`mixed2way`) — how the interaction is computed in hd-EEG

Scope: the `mixed2way` GLM preset in `core_snpm_glm.m`, with the design built in
`dependencies/snpm_glm_design.m`, the statistic in `dependencies/snpm_glm_stat.m`,
and the permutation in `dependencies/snpm_glm_permute.m`. One GLM is fit **per
channel**; family-wise error is controlled across channels by max-statistic
permutation with TFCE and cluster-mass.

## The model (one GLM per channel)

For a 2×2 mixed design — one between factor (Group), one within factor
(Condition), subjects nested in Group — the design (`snpm_glm_design.m:87-116`) is

```
X = [ 1 | Ds | Dc | Dx ]
```

- `Ds` — reference-coded **subject** dummies (n_subj − 1 columns). These
  per-subject intercepts absorb both the subject random effect *and* the
  between-group main effect, because Group is constant within a subject
  (subjects nested in Group).
- `Dc` — reference-coded **condition** dummies (k_cond − 1 columns).
- `Dx` — the **interaction** columns: elementwise products of every Group dummy ×
  every Condition dummy (`interaction_cols`, line 181). For 2×2 that is 1 column;
  for g×c levels it is (g−1)(c−1) columns.

The interaction lives entirely in `Dx`; everything else (`1, Ds, Dc`) is nuisance
for that test.

## The interaction statistic

The contrast `C` selects only the `Dx` columns (`ofint(end-nx+1:end)=true`,
line 101). More than one selected column → **F**; a single interaction column
(2×2) → still routed as F with df1 = 1 (numerically t²).

`snpm_glm_stat.m` fits Y = Xβ + e by OLS across all channels at once (batched
`pinv`), then for the multi-row contrast (lines 76-80):

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

Flip side (a real limitation, enforced in code): the **between-group main effect
is not estimable in this design.** `effect='group'` deliberately errors and
directs you to run `anova1`/`ancova` on subject means (lines 106-110). Only
`'interaction'` (default) and `'condition'` come out of this preset.

## From one channel to the scalp — the correction

Parametric p from the F is used only for cluster *forming*/reporting. Family-wise
error control is by **permutation of the max statistic**, using **Freedman–Lane**
(`snpm_glm_permute.m`, Winkler et al. 2014) because nuisance terms are present:

1. Fit Y on the nuisance Z = `[1, Ds, Dc]`; keep fitted `Zfit` and residuals `R`.
2. Each permutation: shuffle `R`, add back → `Yperm = Zfit + R(perm,:)`.
3. `perm_type = 'within'` with exchangeability block = subject (lines 115-116),
   so residuals are shuffled **only within each subject** (`within_block_perm`).
   This respects the repeated-measures structure — it permutes condition labels
   within a subject, the exchangeable unit under H₀ of no interaction, and never
   mixes across subjects or groups.
4. Recompute the F-map per permutation; take the max (for TFCE) and max
   cluster-mass; build the null; get corrected p per channel/cluster. Neighbour
   adjacency comes from the montage neighbour matrix; TFCE uses E=0.5, H=2.

## Caveats for sleep-EEG use

- **Balance:** the OLS SS for the interaction is clean only if the design is
  reasonably balanced across Group×Condition. Missing cells or very unequal n per
  group make the F less interpretable (Type I vs III SS ambiguity). The engine
  does not warn.
- **Sphericity / covariance:** this is a fixed-effects GLM with subject dummies,
  not a true mixed model. For >2 within levels there is no Greenhouse–Geisser
  correction — but the permutation null sidesteps the parametric sphericity
  assumption for the *inference*, so this matters more for the point estimate
  than for the p.
- If you need random slopes or trial-level nesting, use the `mixedmodel` path
  (`core_snpm_lmm.m`), not this one.

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
