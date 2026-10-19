function [lme, refit] = snpm_lmm_fitlme(tbl, formula, varargin)
% fitlme with a guard against a false zero-variance solution.
%
%   lme = snpm_lmm_fitlme(tbl, formula, ...) returns fitlme(tbl, formula, ...)
%   unchanged unless some random-effect standard deviation comes back at the
%   boundary (<= BOUNDARY_TOL x the residual SD). In that case the model is
%   refitted once from a random start ('StartMethod','random') drawn from a
%   PRIVATE stream with a fixed seed, and the refit is kept only if its ML (or
%   REML) log-likelihood is higher by more than LL_TOL. [lme, refit] also
%   reports whether the refit was kept.
%
% Why: fitlme optimises the Cholesky factor of the random-effect covariance
% relative to the residual SD (theta). The profiled deviance is symmetric in
% theta, so theta = 0 is always a stationary point (zero gradient). When an
% iterate lands on it (the default start is D = I, theta = 1; the stalled fits
% come back with theta = 0 to within 1e-16), quasinewton stops there even when
% a positive theta has a much higher likelihood. Measured on the lmm_factor parity fixture
% (POWER ~ group + (1|Subject), group_label permutation): 21 channel fits in
% 7 of 400 permuted maps stopped at theta = 0 with a log-likelihood 4-8 lower
% than the optimum at theta ~ 0.35, and an F about 2.6x too large.
%
% Guarantees:
%   * A fit not at the boundary is returned exactly as fitlme returned it (no
%     second fit, bit-identical to the pre-guard behaviour).
%   * The global random stream is never drawn from or advanced: the random
%     start comes from a private RandStream that is swapped in for the one
%     call and the caller's stream object is put back, state untouched, even
%     if the refit throws. Permutation draws are therefore unchanged.
%   * A refit that throws, or that does not improve the likelihood, is
%     discarded and the original fit is returned.

    BOUNDARY_TOL = 1e-3;     % random-effect SD / residual SD regarded as 0
    LL_TOL       = 1e-6;     % log-likelihood gain required to replace the fit
    SEED         = 20261019; % private stream seed (fixed: the refit is deterministic)

    lme = fitlme(tbl, formula, varargin{:});
    refit = false;
    if ~at_boundary(lme, BOUNDARY_TOL)
        return
    end

    % same options, random start (drop any StartMethod the caller passed)
    opts = varargin;
    if ~isempty(opts)
        keep = true(size(opts));
        for k = 1:2:numel(opts) - 1
            if (ischar(opts{k}) || isstring(opts{k})) && strcmpi(opts{k}, 'StartMethod')
                keep([k, k + 1]) = false;
            end
        end
        opts = opts(keep);
    end

    priv = RandStream('mt19937ar', 'Seed', SEED);
    prev = RandStream.setGlobalStream(priv);
    restore = onCleanup(@() RandStream.setGlobalStream(prev));
    try
        alt = fitlme(tbl, formula, opts{:}, 'StartMethod', 'random');
    catch
        alt = [];
    end
    clear restore   % put the caller's global stream back now

    if ~isempty(alt) && alt.LogLikelihood > lme.LogLikelihood + LL_TOL
        lme = alt;
        refit = true;
    end
end

function tf = at_boundary(lme, tol)
% true when any random-effect standard deviation is <= tol x residual SD
    [psi, mse] = covarianceParameters(lme);
    s = sqrt(mse);
    tf = false;
    for g = 1:numel(psi)
        d = sqrt(max(diag(psi{g}), 0));
        if any(d <= tol * s)
            tf = true;
            return
        end
    end
end
