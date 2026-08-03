function flctx = snpm_circ_fl_context(A, logC, usecov, evaluable)
%SNPM_CIRC_FL_CONTEXT  Freedman-Lane context for the circular tier, built with a
%   PER-CHANNEL nuisance covariate.
%
%   flctx = snpm_circ_fl_context(A, logC, usecov, evaluable)
%   Yperm  = snpm_glm_permute(flctx);          % one realisation
%
% The response is the STACKED embedding [cos A, sin A] (nSubj x 2*nCh), so
% snpm_glm_permute applies ONE permutation index to all 2*nCh columns. That is
% the structural reason cos and sin can never desynchronise: there is no code
% path in which they receive different indices.
%
% The nuisance is [1, log(count)_ch] at channels where usecov is true and the
% intercept alone elsewhere, and it is fitted ONCE here, outside the permutation
% loop (Winkler et al. 2014: hold the nuisance fit fixed, permute the residuals,
% add the fit back). The per-channel covariate is why this cannot reuse
% snpm_glm_fl_context, which does a single whole-matrix Z\Y solve with one Z
% shared by every channel.
%
% MISSING DATA. Where usecov is FALSE the nuisance is the intercept alone, the
% column mean is finite everywhere, and a NaN cell simply rides along in the
% residual and is deleted pairwise by the statistic. Where usecov is TRUE the
% caller must have required the column to be COMPLETE, because a covariate fit
% over whole columns cannot tolerate a hole - the same rule
% snpm_glm_fl_context enforces for the whole GLM tier.
%
% INPUTS
%   A         : nSubj x nCh angles in radians (NaN = missing)
%   logC      : nSubj x nCh per-channel covariate (log event count), or []
%   usecov    : 1 x nCh logical, TRUE where logC enters the nuisance
%   evaluable : 1 x nCh logical, the fixed estimability mask
%
% OUTPUT flctx with .Zfit .R .eb .perm_type ('free'), i.e. exactly the struct
% snpm_glm_permute consumes.
%
% See also SNPM_GLM_PERMUTE, SNPM_CIRC_HOTELLING_PERM, SNPM_GLM_FL_CONTEXT.

    [n, nCh] = size(A);
    if isempty(logC), logC = nan(n, nCh); end
    usecov    = reshape(logical(usecov), 1, nCh);
    evaluable = reshape(logical(evaluable), 1, nCh);

    Y    = [cos(A), sin(A)];
    Zfit = nan(n, 2*nCh);
    R    = nan(n, 2*nCh);

    plain = find(evaluable & ~usecov);
    if ~isempty(plain)
        idx = [plain, plain + nCh];
        y   = Y(:, idx);
        mu  = mean(y, 1, 'omitnan');
        Zfit(:, idx) = repmat(mu, n, 1);
        R(:, idx)    = y - Zfit(:, idx);       % NaN cells stay NaN and ride along
    end

    cch = find(evaluable & usecov);
    if ~isempty(cch)
        idx = [cch, cch + nCh];
        c   = [logC(:, cch), logC(:, cch)];    % same covariate for cos and sin
        y   = Y(:, idx);
        cb  = mean(c, 1); yb = mean(y, 1);
        sxx = sum((c - cb) .^ 2, 1);
        b1  = sum((c - cb) .* (y - yb), 1) ./ sxx;
        b0  = yb - b1 .* cb;
        Zfit(:, idx) = b0 + b1 .* c;
        R(:, idx)    = y - Zfit(:, idx);
    end

    % Non-evaluable channels are masked by snpm_perm_correction anyway; zeroing
    % their fitted values keeps a NaN from leaking into a column that IS
    % evaluated (the residual stays NaN, so the column stays NaN).
    Zfit(isnan(Zfit)) = 0;

    flctx = struct('Zfit', Zfit, 'R', R, 'eb', ones(n, 1), 'perm_type', 'free');
end
