# Why this tool exists alongside a PALM-based pipeline

> *Understanding-oriented.* Why hdEEG-SnPM exists alongside an FSL **PALM** group-analysis
> path — a design-rationale discussion, not instructions. To pick and run an analysis, see the
> [how-to guides](../how-to/README.md).

This is a deliberately unflattering comparison, written in July 2026, of hdEEG-SnPM against
**EEG_Processor** — a separate in-house preprocessing and feature-extraction toolbox that is
*not* distributed with this repository. It is recorded here because the same question applies
to anyone who already has PALM: what does a preset-driven permutation tool add over a general
design-matrix one? The honest answer is "less than you might hope, in a specific direction".

EEG_Processor is not only import/preprocess/feature-extraction. It ships a full group-analysis
module (`funcGroupAnalysis/`: `GroupAnalysis_run.m`, `GroupAnalysis_plot.m`,
`GroupAnalysis_htmlreport.m`) built on a bundled copy of FSL **PALM** (Winkler et al.), which is
statistically in the same family as hdEEG-SnPM.

## What EEG_Processor's group analysis does

Runs **PALM** (permutation GLM). `GroupAnalysis_run.m`:

- Builds a design matrix = intercept + predictors + interactions + random variables.
- Supports **exchangeability blocks** and **variance groups**.
- Runs **T-contrasts and F-tests**.
- Correction options: **TFCE, cluster-mass, cluster-extent, FDR** (configurable TFCE H/E,
  cluster-forming threshold).
- Outputs: FWE p-values, Cohen's d / cope / varcope, cluster tables, topoplots with clusters,
  and a ~1278-line HTML report.

Design/model type is driven by the GUI where the user constructs the design matrix and
contrasts, so it can in principle express any linear model PALM supports (t-test, regression,
ANOVA, ANCOVA, interactions, paired via exchangeability).

## Overlap (the honest part)

Both tools are permutation-based, TFCE/cluster-corrected GLM group analysis of the **same**
power-spectrum features. For a plain t / F / regression contrast with TFCE or cluster
correction, **EEG_Processor already does it**, with a more battle-tested engine (PALM is
peer-reviewed; SnPM reimplements the permutation math). SnPM is partly redundant there.

## Where hdEEG-SnPM genuinely differs / adds value

1. **Audience / UX.** PALM makes the user build the design matrix + contrast vectors. SnPM is
   preset-driven: pick a named analysis, map columns, and statistic/contrast/permutation are
   chosen automatically. Lower friction for non-statisticians. This is its strongest justification.
2. **Sleep-specific designs PALM doesn't package:**
   - Two-way mixed ANOVA rendered as a full **3-effect report** (group / condition / interaction),
     each with the correct per-effect error term + simple effects.
   - Per-channel **LMM** for trial-level / awakening data (Stephan 2021) — a mixed-effects fit
     PALM's fixed-effects GLM cannot do.
   - **Circular-statistics** variants for phase/angle data.
3. **Self-contained** pure MATLAB — but weak, since EEG_Processor bundles `fsl_palm` as MATLAB
   source, so neither needs an external FSL install.

## Where EEG_Processor is arguably better

- PALM is peer-reviewed and heavily validated; SnPM must earn trust through its own test suite.
- Integrated: preprocess + analyze in one app, one BIDS tree, no file bridging.
- General design-matrix builder can express arbitrary models (multi-level exchangeability,
  variance groups, synchronized permutation across contrasts) beyond SnPM's preset menu.

## Bottom line

The current SnPM tool earns its place mainly through (a) preset-driven UX for non-statistician
users and (b) the mixed-ANOVA 3-effect report / per-channel LMM / circular designs that are not
wired into EEG_Processor's PALM GUI. For standard t/F/regression contrasts with TFCE/cluster
correction, EEG_Processor's PALM path already covers it and is more validated — so maintaining a
parallel engine is only justified by the UX and the specialized designs.

This comparison is a snapshot of both tools as they stood in July 2026 and has not been
re-checked since; the exact set of models EEG_Processor's design-builder GUI exposes was never
enumerated, so the redundant-versus-unique split above is a judgement, not an inventory.
