function [F, p, df1, df2, rho, extra] = snpm_circ_corr_anglinear(y, A, Cnuis, opts)
% Circular-linear association per channel, as a linear-model F test.
%
%   [F, p, df1, df2, rho, extra] = snpm_circ_corr_anglinear(y, A)
%   [...] = snpm_circ_corr_anglinear(y, A, Cnuis)
%   [...] = snpm_circ_corr_anglinear(y, A, Cnuis, opts)
%
% Tests whether a LINEAR variable y (e.g. a symptom score, an event rate, a
% spectral measure) is associated with the ANGLE at each channel. The angle
% enters the model through its (cos, sin) embedding, so the test is
% rotation-invariant and does not require choosing a reference direction:
%
%   X = [1, Cnuis, cos a, sin a],  contrast selects the two embedding columns
%   F ~ F(2, n - 3 - q)   with q = size(Cnuis, 2)
%
% ALGEBRAIC IDENTITY (removes the guesswork about which quantity to use).
% CircStat's circ_corrcl returns rho with rho^2 exactly equal to the R^2 of y
% regressed on [1, cos a, sin a]. So the covariate-free F here is precisely
% the model F of that regression, and rho is recovered as sqrt(R^2). Both
% identities were validated to 1e-10 against independent references.
%
% rho IS DESCRIPTIVE ONLY, and is reported for the raw (covariate-free) fit
% so that it matches what circ_corrcl would return. Two properties to state
% in any report that prints it:
%   * it is NON-DIRECTIONAL (Mardia 1976) — it has no sign and says nothing
%     about which way the association goes;
%   * it is UPWARD-BIASED at small n, with E[rho^2] ~ 2/(n-1) under the null,
%     so rho ~ 0.32 is the expected null value at n = 20. Never read a raw
%     rho as an effect size without that correction in mind. Inference comes
%     from F / the permutation p, not from rho.
%
% PERMUTATION — STRUCTURAL CONSTRAINT, DO NOT WORK AROUND IT.
% Never permute the angle matrix. Permute the RESPONSE y, through
% dependencies/snpm_glm_permute.m, which applies ONE permutation index to
% every column. Permuting angles invites the failure mode where cos and sin
% are shuffled with different indices (or one channel's cos is paired with
% another's sin), which desynchronises the embedding and silently destroys
% the statistic. Passing y as the permuted quantity makes that class of bug
% unrepresentable. With no covariates, permuting y is exactly as exact as
% permuting the angles; with covariates, Freedman-Lane on y is the correct
% scheme anyway.
%
% INPUTS
%   y     : nSubj x 1 linear variable. NaN = missing.
%   A     : nSubj x nCh angles in RADIANS. NaN = missing at that channel.
%   Cnuis : nSubj x q numeric nuisance covariates, or [] (default). No
%           intercept column — one is added.
%   opts  : struct, optional
%     .min_n           minimum retained subjects per channel. DEFAULT 0 (off);
%                      the gate belongs to the caller and must be held fixed
%                      across permutations (see SNPM_CIRC_HOTELLING).
%     .valid_channels  1 x nCh logical mask; channels outside are NaN.
%     .tol             relative tolerance for the degenerate-design checks.
%                      DEFAULT 1e-12.
%
% MISSING DATA is handled by per-channel (pairwise) deletion, so n varies by
% channel and df2 IS RETURNED AS A PER-CHANNEL VECTOR, with p computed
% channel-wise. Do not assume the scalar df2 that snpm_glm_stat returns.
%
% OUTPUTS (1 x nCh unless noted)
%   F     : model F for the two embedding columns; NaN where not estimable
%   p     : parametric upper-tail p from F(df1, df2), channel-wise
%   df1   : scalar, always 2
%   df2   : 1 x nCh vector, n_ch - 3 - q
%   rho   : 1 x nCh raw circular-linear correlation (covariate-free),
%           descriptive, non-directional, upward-biased at small n
%   extra : .n (retained per channel), .R2 (partial R^2 of the embedding
%           after the covariates), .reason (cellstr)
%
% Reference: Mardia (1976); Zar, Biostatistical Analysis, eq. 27.47.
%
% See also SNPM_GLM_PERMUTE, SNPM_CIRC_HOTELLING, CIRC_CORRCL.

    if nargin < 3, Cnuis = []; end
    if nargin < 4, opts  = struct(); end

    y = y(:);
    [nSubj, nCh] = size(A);
    if numel(y) ~= nSubj
        error('snpm_circ_corr_anglinear:respSize', ...
            'y has %d entries but A has %d rows.', numel(y), nSubj);
    end

    min_n          = getopt(opts, 'min_n', 0);
    valid_channels = getopt(opts, 'valid_channels', []);
    tol            = getopt(opts, 'tol', 1e-12);

    if isempty(Cnuis)
        Cnuis = zeros(nSubj, 0);
    end
    if size(Cnuis, 1) ~= nSubj
        error('snpm_circ_corr_anglinear:covSize', ...
            'Cnuis has %d rows but A has %d rows.', size(Cnuis, 1), nSubj);
    end
    q = size(Cnuis, 2);

    Zfull = [ones(nSubj, 1), Cnuis];
    rowOK = isfinite(y) & all(isfinite(Zfull), 2);

    Zc = cos(A); Zs = sin(A);
    okA = ~isnan(A) & rowOK;
    nn  = sum(okA, 1);

    F   = nan(1, nCh);
    p   = nan(1, nCh);
    df1 = 2;
    df2 = nan(1, nCh);
    rho = nan(1, nCh);
    R2  = nan(1, nCh);
    reason = repmat({''}, 1, nCh);

    % Channels using every covariate-valid row solve in one batch.
    isFull      = all(okA == rowOK, 1) & (nn > 0);
    fullCols    = find(isFull);
    partialCols = find(~isFull);

    if ~isempty(fullCols)
        [Ff, dff, rhof, r2f, whyf] = anglinear_block( ...
            y(rowOK), Zfull(rowOK, :), Zc(rowOK, fullCols), Zs(rowOK, fullCols), q, tol);
        F(fullCols)   = Ff;
        df2(fullCols) = dff;
        rho(fullCols) = rhof;
        R2(fullCols)  = r2f;
        reason(fullCols) = whyf;
    end

    for k = 1:numel(partialCols)
        ch   = partialCols(k);
        rows = okA(:, ch);
        if ~any(rows)
            reason{ch} = 'no observations';
            continue
        end
        [Fc, dfc, rhoc, r2c, whyc] = anglinear_block( ...
            y(rows), Zfull(rows, :), Zc(rows, ch), Zs(rows, ch), q, tol);
        F(ch)   = Fc;
        df2(ch) = dfc;
        rho(ch) = rhoc;
        R2(ch)  = r2c;
        reason(ch) = whyc;
    end

    % ---- estimability gates -------------------------------------------
    if min_n > 0
        bad = nn < min_n;
        for ch = find(bad)
            if isempty(reason{ch})
                reason{ch} = sprintf('below min_n (n=%d)', nn(ch));
            end
        end
        F(bad) = NaN; df2(bad) = NaN;
    end
    if ~isempty(valid_channels)
        if numel(valid_channels) ~= nCh
            error('snpm_circ_corr_anglinear:maskSize', ...
                'opts.valid_channels has %d entries, expected %d.', ...
                numel(valid_channels), nCh);
        end
        masked = ~logical(valid_channels(:)');
        for ch = find(masked)
            if isempty(reason{ch}), reason{ch} = 'excluded by valid_channels'; end
        end
        F(masked) = NaN; df2(masked) = NaN;
    end

    ok = isfinite(F) & isfinite(df2) & df2 > 0;
    p(ok)    = fcdf(F(ok), df1, df2(ok), 'upper');
    F(~ok)   = NaN;
    df2(~ok) = NaN;

    extra = struct('n', nn, 'R2', R2, 'reason', {reason});
end

% ======================================================================
function [F, df2, rho, R2, why] = anglinear_block(y, Z, Zc, Zs, q, tol)
% Partial F for [cos, sin] added to the nuisance model Z, for a set of
% channels sharing the same retained rows. Zc/Zs are nObs x m.
%
% Residualising y and the embedding on Z once, then doing a 2x2 solve per
% channel, is algebraically identical to fitting X = [Z, cos, sin] per channel
% and is fully vectorised across channels.
    m   = size(Zc, 2);
    F   = nan(1, m);
    rho = nan(1, m);
    R2  = nan(1, m);
    why = repmat({''}, 1, m);

    nObs  = numel(y);
    rankZ = rank(Z);
    df2   = nObs - rankZ - 2;

    if rankZ < size(Z, 2)
        why(:) = {'nuisance design rank deficient on retained rows'};
        df2 = NaN;
        return
    end
    if df2 < 1
        why(:) = {sprintf('too few observations (n=%d, q=%d)', nObs, q)};
        df2 = NaN;
        return
    end

    % ---- residualise on the nuisance model (once) ----
    % Z is narrow and full rank here; the explicit inverse keeps the three
    % residualisations below to one factorisation.
    ZtZinv = (Z' * Z) \ eye(size(Z, 2));
    yr = y  - Z * (ZtZinv * (Z' * y));                      % nObs x 1
    cr = Zc - Z * (ZtZinv * (Z' * Zc));                     % nObs x m
    sr = Zs - Z * (ZtZinv * (Z' * Zs));

    RSSred = sum(yr .^ 2);

    a11 = sum(cr .^ 2, 1);
    a22 = sum(sr .^ 2, 1);
    a12 = sum(cr .* sr, 1);
    b1  = yr' * cr;                                         % 1 x m
    b2  = yr' * sr;

    detG = a11 .* a22 - a12 .^ 2;
    scaleG = max(a11 .* a22, eps);
    degen  = ~(detG > tol * scaleG);        % collinear / constant embedding

    extraSS = (b1 .^ 2 .* a22 - 2 * b1 .* b2 .* a12 + b2 .^ 2 .* a11) ./ detG;
    extraSS = max(extraSS, 0);
    RSSfull = RSSred - extraSS;

    flat = ~(RSSred > 0);                   % y constant on retained rows
    sat  = ~(RSSfull > tol * max(RSSred, eps));   % perfect fit -> F undefined

    F = (extraSS / 2) ./ (RSSfull / df2);
    R2 = extraSS ./ RSSred;

    bad = degen | flat | sat;
    F(bad)  = NaN;
    R2(bad) = NaN;
    why(degen) = {'angle embedding degenerate (all angles identical?)'};
    why(flat & ~degen) = {'response constant on retained rows'};
    why(sat & ~degen & ~flat) = {'residual sum of squares numerically zero'};

    % ---- descriptive raw rho: covariate-free R^2 of y on [1, cos, sin] ----
    y0  = y - mean(y);
    c0  = Zc - mean(Zc, 1);
    s0  = Zs - mean(Zs, 1);
    A11 = sum(c0 .^ 2, 1);
    A22 = sum(s0 .^ 2, 1);
    A12 = sum(c0 .* s0, 1);
    B1  = y0' * c0;
    B2  = y0' * s0;
    detG0 = A11 .* A22 - A12 .^ 2;
    TSS   = sum(y0 .^ 2);
    ess0  = (B1 .^ 2 .* A22 - 2 * B1 .* B2 .* A12 + B2 .^ 2 .* A11) ./ detG0;
    r2raw = ess0 / max(TSS, eps);
    r2raw = min(max(r2raw, 0), 1);
    rho   = sqrt(r2raw);
    rho(~(detG0 > tol * max(A11 .* A22, eps)) | TSS <= 0) = NaN;
end

% ======================================================================
function v = getopt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
