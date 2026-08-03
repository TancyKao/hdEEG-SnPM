function flctx = snpm_glm_fl_context(Y, D)
% Freedman-Lane permutation context for the GLM tier (Winkler et al. 2014).
% Single source of truth: every GLM caller (core_snpm_glm, export_glm_report,
% the test suite) builds its context here so the evaluable-channel rule below
% cannot drift between them.
%
% The nuisance part of the design is fitted ONCE and held fixed; permutation
% then reshuffles the residuals and adds the nuisance fit back
% (snpm_glm_permute), so the nuisance structure is preserved under the null.
%
% INPUTS
%   Y : nObs x nCh data
%   D : design struct from snpm_glm_design (.X, .nuisance_idx, .eb, .perm_type)
%
% OUTPUT flctx
%   .Zfit      : nObs x nCh nuisance fitted values  Z*(Z\Y)
%   .R         : nObs x nCh residuals               Y - Zfit
%   .eb        : nObs x 1 exchangeability-block id
%   .perm_type : 'free' | 'within'
%   .evaluable : 1 x nCh logical, TRUE for channels that can be evaluated in
%                EVERY permutation (see below)
%
% WHY .evaluable EXISTS
%   Z\Y is a whole-matrix least-squares solve: a SINGLE missing cell makes the
%   entire column of betaZ (hence Zfit and R) NaN, so that channel is NaN in
%   every permuted map. snpm_glm_stat, by contrast, falls back to a per-column
%   solve on the available rows and still returns a finite OBSERVED statistic
%   there. Left alone that mismatch is a null-construction defect, not a
%   cosmetic one: the channel has an observed value but never enters the
%   max-statistic null, so it is scored against a distribution assembled from
%   the other channels only, and its neighbours' observed TFCE is integrated
%   over a channel the permutations never see.
%
%   The fix is to analyse ONE channel set: any channel that cannot be evaluated
%   in every permutation is dropped from the observed map too
%   (snpm_perm_correction applies the mask to both). Column-wise missingness is
%   a property of the data, not of the labels being permuted, and Z is fixed
%   across permutations, so the mask is label-invariant and is computed once,
%   here, before the permutation loop starts.

    Z = D.X(:, D.nuisance_idx);
    if isempty(Z), Z = ones(size(Y, 1), 1); end
    betaZ = Z \ Y;
    Zfit  = Z * betaZ;
    R     = Y - Zfit;

    flctx = struct();
    flctx.Zfit      = Zfit;
    flctx.R         = R;
    flctx.eb        = D.eb;
    flctx.perm_type = D.perm_type;
    % A permuted column is Zfit + R(pe,:); it is finite for every permutation
    % index pe if and only if both Zfit and R are finite throughout the column.
    flctx.evaluable = all(isfinite(Zfit), 1) & all(isfinite(R), 1);
end
