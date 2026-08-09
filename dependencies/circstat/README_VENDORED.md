# Vendored subset of CircStat (2012a)

This folder is a **partial, unmodified copy** of the CircStat toolbox for MATLAB. It is kept in
its own subfolder (not flattened into `dependencies/`) because CircStat ships `Contents.m`,
`license.txt` and `readme.txt`, which would collide with repo files. No path setup is needed:
`core_snpm_analysis.m` calls `addpath(genpath(pwd))`, which recurses into this subfolder.

- **Source:** `CircStat 2012a (upstream distribution)`
- **Version:** CircStat2012a (file dates 2014-02-12; the 2012a release)
- **Vendored on:** 2026-08-02
- **Files are byte-identical copies of the upstream release. Do not edit them.** If a fix is
  needed, write our own function under `dependencies/` instead, so the vendored tree stays a
  clean, verifiable copy of upstream.

## Read this before using anything here for inference

**CircStat contains no Wheeler-Watson uniform-scores test and no two-sample Watson U-squared
test.** Neither statistic exists anywhere in the 2012a release. The similarly-named
`circ_wwtest` is the **Watson-*WILLIAMS*** test: a *parametric* one-way von Mises ANOVA for
equal mean directions that assumes von Mises distributed data with a common concentration
parameter `k`. That is a different test with different assumptions, and it is not a substitute
for either of the two non-parametric statistics this toolbox needs.

Further, `circ_wwtest`'s documented four-argument weighted two-sample form
(`circ_wwtest(alpha1, alpha2, w1, w2)`) is **buggy**. In its `processInput` subfunction:

```matlab
if nargin == 4
  ...
  idx = [ones(size(alpha1)); ones(size(alpha2))];   % <-- both groups labelled 1
```

Both samples are assigned group index `1`, so the group count `s = length(unique(idx))` is `1`
and the F statistic `F = beta * (n-s) * A / (s-1) / B` divides by zero. The two-argument form
(`circ_wwtest(alpha1, alpha2)` with non-integer angles, or `circ_wwtest(alpha, idx)`) labels
groups correctly and does not hit this path.

`circ_corrcl` (circular-linear correlation) returns **`[rho, pval]`** — two outputs, `rho`
first, `pval` computed as `1 - chi2cdf(n*rho^2, 2)`. It is **scalar-only**: it takes a single
angle vector and a single linear vector and returns one scalar `rho`, with no dimension
argument and no way to compute across channels in one call. Wrapping it in a per-channel loop
inside a permutation loop is far too slow for this toolbox's workloads.

**Therefore CircStat is vendored here for descriptive/support use, fixture generation, and as a
validation reference ONLY.** The three core circular statistics for this toolbox are being
written in-repo as vectorised `snpm_circ_*` functions. Nobody should later "simplify" the
engine by swapping our statistic for `circ_wwtest` — it computes a different test, under
different assumptions, and its weighted two-sample form is broken.

## Licence

CircStat is distributed under the **BSD 2-clause licence**. The full text is in `license.txt`
in this folder, retained verbatim as the licence requires. Copyright line, verbatim:

```
Copyright (c) 2011, Philipp Berens
All rights reserved.
```

The licence requires that redistributions of source retain the copyright notice, the list of
conditions and the disclaimer — satisfied by keeping `license.txt` alongside the code here.

## Citation

The upstream `readme.txt` asks: *"Please cite this paper when the provided code is used."*

> Berens P (2009). CircStat: A MATLAB Toolbox for Circular Statistics. *Journal of Statistical
> Software*, 31(10). https://www.jstatsoft.org/v31/i10

This is a **courtesy request in the readme, not a condition of the BSD-2 licence** — the licence
itself imposes no citation requirement. We cite it anyway in any publication or report that uses
these functions.

## What is vendored, and why

| File | Why it is here |
| --- | --- |
| `circ_mean.m` | Mean direction per channel — the descriptive statistic behind circular topoplots. Calls `circ_confmean` when `nargout > 1`. |
| `circ_r.m` | Resultant vector length (phase-locking strength); also the workhorse called by `circ_kappa`, `circ_std`, `circ_confmean`, `circ_wwtest`. |
| `circ_confmean.m` | Confidence limits on the mean direction; required because `circ_mean` calls it whenever more than one output is requested. |
| `circ_kappa.m` | von Mises concentration parameter; needed by `circ_wwtest` and useful for describing phase spread. |
| `circ_std.m` | Circular standard deviation / dispersion around the mean, for descriptive summaries. |
| `circ_dist.m` | Signed angular difference on the circle, wrapped to (-pi, pi] — used wherever phase differences are taken. |
| `circ_rad2ang.m` | Radians to degrees, for reporting and plot axes. |
| `circ_ang2rad.m` | Degrees to radians, for ingesting angle data supplied in degrees. |
| `circ_vmrnd.m` | von Mises random number generator — generates synthetic phase data for test fixtures. Supports a two-entry `n` (e.g. `[10 3]`) for matrix output. |
| `circ_corrcl.m` | **Validation reference only.** Circular-linear correlation; scalar-only, so it is used to check our vectorised implementation, not to run analyses. |
| `circ_wwtest.m` | **Validation reference only.** Parametric Watson-Williams von Mises ANOVA. Kept so the difference from our non-parametric statistics can be demonstrated and cross-checked; see the warning above. Not to be used for inference here. |
| `license.txt` | BSD-2 licence text, retained verbatim as the licence requires. |
| `readme.txt` | Upstream readme — carries the citation request, author contact, and the disclaimer. |

## Dependency closure (verified)

Every `circ_*` function called from a code line (comments excluded) of a vendored file is itself
vendored, so this subset is self-contained:

```
circ_confmean  -> circ_r
circ_kappa     -> circ_r
circ_mean      -> circ_confmean
circ_std       -> circ_r
circ_wwtest    -> circ_kappa, circ_r
circ_ang2rad, circ_corrcl, circ_dist, circ_r, circ_rad2ang, circ_vmrnd  -> (no circ_* calls)
```

Note: `circ_corrcl.m`'s header comment mistakenly names `circ_corrcc` on its usage line
(an upstream typo); it does not call it, so `circ_corrcc.m` is not needed.

Non-CircStat dependencies are MATLAB base plus the Statistics and Machine Learning Toolbox,
which this repo already requires: `chi2inv` (`circ_confmean`), `chi2cdf` and `corr`
(`circ_corrcl`), `fcdf` (`circ_wwtest`), `rand` (`circ_vmrnd`).

## Deliberately NOT vendored

| Excluded | Reason |
| --- | --- |
| `examples/` subfolder | Contains generically-named `parseVarArgs.m` and `formatSubplot.m` that would shadow same-named functions from other toolboxes on the user's MATLAB path once `genpath` picks this folder up. |
| `circ_kuipertest.m`, `kuipertable.mat` | Unused by this toolbox; the lookup-table `.mat` would also land among the montage assets. |
| `Contents.m` | Avoids ambiguity in MATLAB's `help` / `ver` machinery, and would collide conceptually with repo-level contents files. |
| All other `circ_*.m` (e.g. `circ_rtest`, `circ_otest`, `circ_hktest`, `circ_plot`, `circ_clust`, `circ_median`, `circ_var`, ...) | Not needed by the planned circular path. Keeping the vendored surface minimal keeps the closure auditable and reduces path-shadowing risk. |

## Updating this subset

If another CircStat function turns out to be needed, copy it unmodified from the same upstream
release, re-run the closure check (grep the copied files for `circ_*` on non-comment lines and
confirm each callee is present), and add a row to the table above. Do not modify vendored code
in place.
