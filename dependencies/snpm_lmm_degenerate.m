function [tf, info] = snpm_lmm_degenerate(varargin)
% True when a mixed-model fit has no usable residual variation, so its t/F is
% not defined (the per-labelling backstop of the LMM zero-variance rule).
%
%   tf = snpm_lmm_degenerate(lme)                      % a fitted LinearMixedModel
%   tf = snpm_lmm_degenerate(y, X, Z, one_intercept)   % the same rule on raw arrays
%   [tf, info] = ...   info.rss_x, .rss_xz, .tol, .theta2 (NaN when not computed)
%
% snpm_lmm_fit reports such a channel as stat NaN, wald 0, p NaN, in the
% observed map and in any permuted map; an OBSERVED channel caught here is
% listed by core_snpm_lmm as 'perfect fit, statistic undefined'. Most such
% channels never reach the fit: core_snpm_lmm first excludes, once and for
% every labelling, channels whose value is constant within each subject
% (snpm_lmm_within_ratio <= 1e-6). This is the backstop for what is left.
%
% THE RULE (computed on the data; fitlme's own variance estimates are not
% usable here: on a flat channel it returns SDs of ~1e-16 and an arbitrary
% finite F, and above a subject/residual SD ratio of ~1e4 it is off the optimum):
%   y : response (N x 1), X : fixed-effects design (N x p, full column rank),
%   Z : random-effects design (N x q), tol = 1e-20 * y'y
%   rss_x  = || y - X * (X \ y) ||^2
%   rss_xz = || y - [X Z] * lsqminnorm([X Z], y) ||^2
%   rx = rank(X), rxz = rank([X Z])   (MATLAB rank: SVD, tolerance
%                                        max(size(A)) * eps(norm(A)))
%   degenerate if rss_x <= tol, or (rxz < N and rss_xz <= tol) -- when [X Z]
%   has rank N it fits ANY y exactly (e.g. (POWER|Subject) with 2 trials per
%   subject), so rss_xz = 0 says nothing about the data and is not used;
%   otherwise, when the model has exactly one random-effect term and it is a
%   scalar random intercept (1|G) (one_intercept), with
%     df_b = rxz - rx,  df_w = N - rxz  (both must be > 0, else not degenerate)
%     F_b  = ((rss_x - rss_xz) / df_b) / (rss_xz / df_w)
%     theta2 = (F_b - 1) / (N / q)      (method-of-moments (SD_subject/SD_resid)^2;
%                                        q = number of columns of Z = levels of G)
%   degenerate if theta2 >= 1e3^2.

    ZERO_SS_RTOL = 1e-20;
    THETA_MAX    = 1e3;

    if nargin == 1
        lme = varargin{1};
        y = response(lme);
        X = designMatrix(lme, 'Fixed');
        Z = full(designMatrix(lme, 'Random'));
        psi = covarianceParameters(lme);
        one_intercept = numel(psi) == 1 && size(psi{1}, 1) == 1;
    else
        [y, X, Z, one_intercept] = varargin{1:4};
        Z = full(Z);
    end
    N = numel(y);
    tol = ZERO_SS_RTOL * (y' * y);

    rss_x = sum((y - X * (X \ y)).^2);
    XZ = [X Z];
    rss_xz = sum((y - XZ * lsqminnorm(XZ, y)).^2);
    rxz = rank(XZ); rx = rank(X);
    info = struct('rss_x', rss_x, 'rss_xz', rss_xz, 'tol', tol, 'theta2', NaN, ...
        'rank_x', rx, 'rank_xz', rxz);
    if rss_x <= tol || (rxz < N && rss_xz <= tol)
        tf = true;
        return
    end

    tf = false;
    if one_intercept
        q = size(Z, 2);
        df_b = rxz - rx; df_w = N - rxz;
        if df_b > 0 && df_w > 0
            F_b = ((rss_x - rss_xz) / df_b) / (rss_xz / df_w);
            theta2 = (F_b - 1) / (N / q);
            info.theta2 = theta2;
            tf = theta2 >= THETA_MAX^2;
        end
    end
end
