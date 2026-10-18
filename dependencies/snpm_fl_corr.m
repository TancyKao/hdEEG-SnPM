function r = snpm_fl_corr(fl, xr, yr, order, v2)
%SNPM_FL_CORR  One Freedman-Lane partial-correlation map.
%
%   r = SNPM_FL_CORR(fl, xr, yr, order, v2)
%
%   fl    : struct with .Rz (n x n residual-forming matrix I - Q*Q' of the
%           nuisance design [1, covariates]) and .k_eff (its rank minus 1)
%   xr    : n x nCh residuals Rz*x (the EEG side, the side that is permuted)
%   yr    : n x nCh residuals Rz*y (fixed)
%   order : 1 x n row order; 1:n gives the OBSERVED map
%   v2    : 1 x nCh sum(xr.^2, 1), the zero-variance reference
%   r     : 1 x nCh  corr(Rz * xr(order,:), yr), column-wise
%
% The permuted residuals are RE-PROJECTED off the covariates before they are
% correlated (Freedman & Lane 1983; Winkler et al. 2014). Permuting xr
% against yr and stopping there (Kennedy) leaves the permuted residuals
% correlated with the covariates, which narrows the null and inflated the
% family-wise error to 0.87 at n = 10 with 4 covariates (Method Spec approved
% by Tancy 2026-10-18).
%
% The observed map MUST be computed through this function with order = 1:n,
% not through a separate corr()/partialcorr() call, so that the identity row
% of the permutation set reproduces it bit for bit and max-statistic ties are
% exact. It equals partialcorr(x, y, Z) to ~1e-15.
%
% Zero-variance rule 2: a column whose re-projected permuted residual has no
% length by the shared tolerance (||Rz*xr(order,c)||^2 <= 1e-20*||xr_c||^2,
% snpm_zero_ss) scores r = 0 for that order. NaN (masked) columns stay NaN:
% a matrix product never mixes columns.
%
% Used by snpm_single_threshold_with_TFCE and snpm_cluster_analysis
% (correlationPboth), so both engines score the same map for the same order.

    n = size(xr, 1);
    if ~isstruct(fl) || ~isfield(fl, 'Rz') || ~isfield(fl, 'k_eff') || ~isequal(size(fl.Rz), [n n])
        error('snpm:flContext', ['Freedman-Lane context must be struct(''Rz'', n x n, ' ...
            '''k_eff'', k) with n = %d subjects.'], n);
    end
    xp = fl.Rz * xr(order, :);
    r  = snpm_corr_columns(xp, yr);
    z  = snpm_zero_ss(sum(xp.^2, 1), v2) & isfinite(v2);
    r(z) = 0;
end
