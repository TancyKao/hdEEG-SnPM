function [stat, p_model, df1, df2, zero_err] = snpm_glm_stat(Y, X, C)
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
%   zero_err: 1 x nCh logical, TRUE where the full-model residual sum of
%             squares counts as zero (zero-variance rule 2, below)
%
% ZERO-VARIANCE RULE 2 (Method Spec 2026-10-14, GLM tier). Where the full-model
% residual sum of squares meets the tolerance against sum(Y.^2) of the data
% being fitted (snpm_zero_ss: RSS <= 1e-20 * sum(Y.^2)), the statistic is
% DEFINED as 0 and p as 1. Without it a zero error term gives F = Inf/NaN or a
% 1e30-scale F that breaks ClusterEnhancement's threshold ladder. Because this
% function computes the observed AND every permuted statistic, the rule
% applies to every labelling alike, so it is one fixed statistic and the
% permutation test stays exact. The caller uses zero_err on the OBSERVED fit
% to report those channels as undefined (rule 3). On data with a non-zero
% error term everywhere this changes nothing.

    [nObs, nCh] = size(Y);
    q   = size(C, 1);
    df1 = q;

    stat    = nan(1, nCh);
    p_model = nan(1, nCh);
    zero_err = false(1, nCh);

    nanmask    = isnan(Y);
    colHasNan  = any(nanmask, 1);
    colAllNan  = all(nanmask, 1);
    fullCols   = find(~colHasNan);                 % fully present -> batch solve
    partialCols= find(colHasNan & ~colAllNan);     % some NaN      -> per-column

    % ---- batched solve over fully-present channels ----
    if ~isempty(fullCols)
        [s, pm, df2, z] = local_glm(Y(:, fullCols), X, C, q);
        stat(fullCols)    = s;
        p_model(fullCols) = pm;
        zero_err(fullCols) = z;
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
        [s, pm, ~, z] = local_glm(Y(rows, ch), Xc, C, q);
        stat(ch)    = s;
        p_model(ch) = pm;
        zero_err(ch) = z;
    end
end

function [stat, p_model, df2, zero_err] = local_glm(Y, X, C, q)
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

    % Zero-variance rule 2: a zero error term defines the statistic as 0, p 1.
    zero_err = snpm_zero_ss(RSS, sum(Y.^2, 1));
    stat(zero_err)    = 0;
    p_model(zero_err) = 1;
end
