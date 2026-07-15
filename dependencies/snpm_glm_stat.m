function [stat, p_model, df1, df2] = snpm_glm_stat(Y, X, C)
% Vectorized General Linear Model statistic across all channels at once.
% Fits Y = X*beta + e by ordinary least squares for every channel
% simultaneously and evaluates the contrast C:
%   - C with one row   -> signed t statistic
%   - C with >1 row    -> F statistic (omnibus, non-negative)
%
% This single computation underlies every GLM preset (t-test, one-way and
% repeated-measures ANOVA, ANCOVA, regression, interaction) — they differ
% only in how X and C are built (see snpm_glm_design).
%
% INPUTS
%   Y : nObs x nCh data (observations x channels)
%   X : nObs x p design matrix (include an intercept column if wanted)
%   C : q x p contrast matrix (q = 1 -> t, q > 1 -> F)
%
% OUTPUTS (each 1 x nCh)
%   stat    : signed t (q==1) or F (q>1); NaN for channels that cannot be fit
%   p_model : parametric two-sided p (used for cluster forming / reporting)
%   df1, df2: contrast and error degrees of freedom (scalars)

    [nObs, nCh] = size(Y);
    q   = size(C, 1);
    df1 = q;

    stat    = nan(1, nCh);
    p_model = nan(1, nCh);

    nanmask    = isnan(Y);
    colHasNan  = any(nanmask, 1);
    colAllNan  = all(nanmask, 1);
    fullCols   = find(~colHasNan);                 % fully present -> batch solve
    partialCols= find(colHasNan & ~colAllNan);     % some NaN      -> per-column

    % ---- batched solve over fully-present channels ----
    if ~isempty(fullCols)
        [s, pm, df2] = local_glm(Y(:, fullCols), X, C, q);
        stat(fullCols)    = s;
        p_model(fullCols) = pm;
    else
        df2 = nObs - rank(X);
    end

    % ---- per-column solve for channels with partial missingness ----
    for k = 1:numel(partialCols)
        ch = partialCols(k);
        rows = ~nanmask(:, ch);
        Xc = X(rows, :);
        if rank(Xc) < size(Xc, 2) || sum(rows) <= size(Xc, 2)
            continue   % under-determined -> leave NaN
        end
        [s, pm] = local_glm(Y(rows, ch), Xc, C, q);
        stat(ch)    = s;
        p_model(ch) = pm;
    end
end

function [stat, p_model, df2] = local_glm(Y, X, C, q)
% OLS for one design X applied to columns of Y (nObs x m), with contrast C.
    [nObs, ~] = size(Y);
    XtXinv = pinv(X' * X);
    B  = XtXinv * (X' * Y);               % p x m
    R  = Y - X * B;                        % residuals
    RSS = sum(R.^2, 1);                    % 1 x m
    rankX = rank(X);
    df2 = nObs - rankX;
    sigma2 = RSS / df2;                    % 1 x m

    CB = C * B;                            % q x m
    M  = C * XtXinv * C';                  % q x q (contrast covariance factor)

    if q == 1
        se   = sqrt(sigma2 * M);           % 1 x m
        stat = CB ./ se;                   % signed t
        p_model = 2 * tcdf(-abs(stat), df2);
    else
        Minv = pinv(M);                    % q x q
        quad = sum((Minv * CB) .* CB, 1);  % 1 x m  = CB' Minv CB per column
        stat = (quad / q) ./ sigma2;       % F
        p_model = fcdf(stat, q, df2, 'upper');
    end
end
