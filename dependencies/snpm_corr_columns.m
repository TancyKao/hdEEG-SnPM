function r = snpm_corr_columns(X, Y)
%SNPM_CORR_COLUMNS  Column-wise Pearson correlation of two matched matrices.
%
%   r = SNPM_CORR_COLUMNS(X, Y) returns a 1 x nCol vector with
%   r(i) = correlation between X(:,i) and Y(:,i), computed for every column at
%   once (no per-channel loop, no per-channel call into corr).
%
%   X and Y must have the same size: rows are the analysed units (subjects),
%   columns are channels.
%
% COMPLETE-COLUMN RULE (this is the point of the function)
%   A column that contains ANY NON-FINITE value (NaN, +Inf or -Inf) in X or in Y
%   returns NaN. There is no pairwise (per-channel) deletion here, deliberately.
%
%   Non-finite, not just NaN: log10(0) = -Inf, which datatype 'logscale'
%   produces from a zero-power cell. Centring such a column gives Inf and NaN
%   entries and the r below comes out non-finite, which the clamp at the end
%   turns into NaN. Callers must therefore build their evaluable-channel mask
%   with isfinite rather than ~isnan, or the channel set they report will not be
%   the channel set this function actually evaluates (see core_snpm_analysis).
%
%   Pairwise deletion inside a permutation loop is an inference bug, not a
%   convenience: the permutation reorders the rows of X against a fixed Y, so
%   the set of complete pairs at a partially-missing channel is a FUNCTION OF
%   THE PERMUTATION. The observed map is then computed on a different sample
%   size from the permuted maps, and since the map statistic is raw r, whose
%   null spread scales as 1/sqrt(n-1), the null is mis-scaled. Measured on a
%   30-subject design with 25% missingness at one channel: observed n = 22
%   while the permuted n ranged [14, 20] -- the observed value sat outside the
%   entire permuted range. Family-wise error at nominal 0.05 reached 0.288
%   (TFCE, 40% disjoint missingness, 2000 replicates).
%
%   Returning NaN for an incomplete column makes the analysed row set the same
%   for the observed labelling and for every permutation, which is what the
%   permutation null assumes. Callers are expected to decide the evaluable
%   channel set ONCE, before the permutation loop (see core_snpm_analysis), and
%   this rule is the engine-level backstop that makes a direct caller safe too.
%
% NUMERICAL NOTE
%   Follows MATLAB's own corrPearson: centre, then divide the inner product by
%   the two column norms, then clamp |r| <= 1. Agreement with a per-channel
%   corr(x,y) call is to ~2e-16 (one unit in the last place), not bitwise --
%   no vectorised form reproduces corr() bitwise, because corr() on an n x 1
%   pair and on an n x m block take different BLAS paths. Measured on
%   randn(30,178): max |difference| = 1.665e-16.
%
%   A constant column gives 0/0 = NaN, matching corr().
%
%   Spearman: rank-transform the columns ONCE with tiedrank and pass the ranks
%   in. Row permutation commutes with ranking (permuting rows then ranking ==
%   ranking then permuting rows), so hoisting the rank transform out of the
%   permutation loop is exact, not an approximation.

    if ~isequal(size(X), size(Y))
        error('snpm_corr_columns:sizeMismatch', ...
            'X is %s and Y is %s; column-wise correlation needs identical sizes.', ...
            mat2str(size(X)), mat2str(size(Y)));
    end

    n = size(X, 1);
    if n < 3
        error('snpm_corr_columns:tooFewRows', ...
            'Column-wise correlation needs at least 3 rows; got %d.', n);
    end

    Xc = X - sum(X, 1) / n;          % centre (NaN propagates down the column)
    Yc = Y - sum(Y, 1) / n;

    dx = sqrt(sum(Xc .^ 2, 1));
    dy = sqrt(sum(Yc .^ 2, 1));

    r = sum(Xc .* Yc, 1) ./ (dx .* dy);

    % Constant / degenerate columns -> NaN (corr() behaviour); rounding can put
    % |r| a hair above 1, clamp exactly as corr() does.
    r(~isfinite(r)) = NaN;
    over = abs(r) > 1;
    r(over) = sign(r(over));
end
