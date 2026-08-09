# Example data — one file (or set) per analysis

Ready-made synthetic data, shipped with the repository. Nothing needs to be generated: point the
GUI or a `scripts/run_*.m` template straight at these files.

All files use the **256-channel EGI montage** (columns `E1`..`E256`), so set **Recording system =
"EGI 256 (HydroCel)"** (the legacy label **"164 channels"** resolves to the same montage) and
**Data Type = "absolute"**. A neighbour cluster of **17 channels** around **E129** carries the
planted effect; everywhere else is noise. Use **Permutations ~ 500-1000** for a quick check
(raise for real runs — final inference wants 5000-10000, because the permutation p is
`(b+1)/(N+1)`). Each analysis should flag a significant cluster near E129.

These files are synthetic. They exist so you can see an analysis run end to end and check that
it recovers an effect that was deliberately planted; they are not sleep data and carry no
scientific meaning of their own.

## GLM presets (single file -> "Data File"; leave Data 2 empty)

| Analysis | File | Meta cols | Roles to set |
|---|---|---|---|
| 3+ groups (ANOVA) | `glm_anova1.csv` | Subject,group | Group=`group` |
| ANCOVA (2 groups -> t) | `glm_ancova.csv` | Subject,group,age,sex | Group=`group`, Covariates=`age,sex` |
| ANCOVA (3 groups -> F + post-hoc) | `glm_ancova3.csv` | Subject,group,age,sex | Group=`group`, Covariates=`age,sex` |
| Regression | `glm_regression.csv` | Subject,cognition,age | Predictor=`cognition`, Covariates=`age` |
| Repeated measures | `glm_rmanova.csv` | Subject,condition | Subject=`Subject`, Condition=`condition` |
| Group x condition (2x2) | `glm_mixed2way.csv` | Subject,group,condition | Group=`group`, Subject=`Subject`, Condition=`condition` |
| Group x condition (2x3) | `glm_mixed2way_2x3.csv` | Subject,group,condition | Group=`group`, Subject=`Subject`, Condition=`condition` |

Levels: `glm_anova1.csv` has three groups `G1/G2/G3`; `glm_rmanova.csv` three conditions
`C1/C2/C3`; `glm_mixed2way.csv` is `G1/G2` x `pre/post` and `glm_mixed2way_2x3.csv` is
`HC/PT` x `pre/mid/post`. `mixed2way` always reports **all three** effects (group, condition,
interaction) in one HTML report, each with its own permutation scheme.

## Legacy / two-file tests (Data 1 + Data 2)

| Analysis | Data 1 | Data 2 |
|---|---|---|
| Paired t-test | `paired_condA.csv` | `paired_condB.csv` |
| One-sample (single condition vs 0) | `onesample_change.csv` | *(leave empty)* |
| One-sample (legacy A-B vs 0) | `onesample_condA.csv` | `onesample_condB.csv` |
| Unpaired t-test | `unpaired_groupA.csv` | `unpaired_groupB.csv` |
| Correlation (Pearson/Spearman) | `corr_eeg.csv` | `corr_behavior.csv` |

## Circular (phase / angle) analyses

All angles are **radians on `[0, 2*pi)` with zero at the slow-oscillation up-state peak**,
so set **Units = `rad`** and **Convention = `literature_uppeak0`**. Tails is locked to
`both` and Data Type to `absolute` (log/z-score of an angle is undefined).

| Analysis (key) | Angles 1 | Angles 2 / measure | Event counts (REQUIRED) |
|---|---|---|---|
| Phase, 2 groups (`circ_phase_group`, Hotelling) | `circ_phase_groupA.csv` | `circ_phase_groupB.csv` | `circ_phase_countsA.csv` + `circ_phase_countsB.csv` |
| Phase, 2 groups (`circ_phase_group_u2`, Watson U2) | same pair | same pair | same pair |
| Concentration demo (either group key) | `circ_conc_groupA.csv` | `circ_conc_groupB.csv` | `circ_conc_countsA.csv` + `circ_conc_countsB.csv` |
| Circular-linear (`circ_corrAngLinear`) | `circ_anglinear_angles.csv` | `circ_anglinear_measure.csv` (`measure` col) | *(not used)* |

**Subject ids differ by design.** The two group fixtures use **disjoint** ids
(`sub001..` vs `sub101..`) because they are independent groups; the circular-linear
fixture uses **matched** ids because it pairs an angle with a measure per subject.
If your own two angle files contain the same subjects, a two-group circular test is
the wrong test.

`circ_phase_group*` plants a **46.3 deg mean-direction shift** (von Mises kappa 4,
16 + 16 subjects) in the cluster and should recover it. `circ_conc_*` plants a
**concentration** difference at the same mean direction (kappa 14.2 vs 3.94, the
Hahn 2020 magnitude, 20 + 20 subjects): it demonstrates what Hotelling and U2 respond
to and the signed linearised measure does not, but at that n it is detected only
about 22% of the time, so do **not** treat a null run on it as a defect.
The count files are the **precision covariate** (`count1_file`/`count2_file`): a
subject's preferred phase is estimated from their detected events, so its precision
depends on how many there were and how tightly they clustered. Counts are balanced
across groups in every fixture, so these fixtures carry no precision confound.
See `docs/explanation/circular-statistics-for-phase.md` for why this input is
required and what it is moving to (`rayleigh_z`).

## LMM (script only: `scripts/run_lmm_analysis.m`)

`lmm_long.csv` — long format, one row per trial: Subject, group, time, sleepDepth (DV) + E1..E256.
Effect: cluster power predicts `sleepDepth`. Set lmm_dv=`sleepDepth`, lmm_fixed=`POWER`, lmm_effect=`POWER`.
