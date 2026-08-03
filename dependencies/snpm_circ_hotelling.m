function [F, p, df1, df2, extra] = snpm_circ_hotelling(A, g, Cnuis, opts)
% Two-sample Hotelling T^2 on the (cos, sin) embedding, per channel.
%
%   [F, p, df1, df2, extra] = snpm_circ_hotelling(A, g)
%   [...] = snpm_circ_hotelling(A, g, Cnuis)
%   [...] = snpm_circ_hotelling(A, g, Cnuis, opts)
%
% This is the primary Tier-2 circular statistic: it asks whether the mean
% RESULTANT VECTOR (direction and concentration jointly) differs between two
% groups of subjects at each channel. Angles are embedded as the 2-vector
% (cos a, sin a) and the group difference in that 2-vector is tested with a
% multivariate contrast, which is what makes the test rotation-invariant: an
% arbitrary rotation of the reference direction is an orthogonal rotation of
% the embedding, and T^2 is invariant under it.
%
% COVARIATE-ADJUSTED GENERAL FORM (Method Spec Addendum, Amendment 2). The
% design is X = [1, Cnuis, g] with g the 0/1 group indicator; the response is
% Z = [cos a, sin a] (N x 2); the contrast c selects the g column:
%
%   nu = N - rank(X)                      % = N - q - 2
%   B  = (X'X)^-1 X' Z                    % p x 2
%   S  = (Z - X B)' (Z - X B) / nu        % 2 x 2 residual covariance
%   d  = (c' B)'                          % 2 x 1
%   v  = c' (X'X)^-1 c                    % scalar
%   T2 = d' S^-1 d / v
%   F  = T2 * (nu - 1) / (2 * nu)   ~  F(2, nu - 1)
%
% At q = 0 (no covariates) this reduces exactly to the classical unadjusted
% two-sample form, T2 = (n1 n2 / N) d' S^-1 d with S the pooled covariance on
% nu = N - 2 df, and F = T2 (N - 3) / (2 (N - 2)) on F(2, N - 3). That
% identity is asserted to machine precision in test_circ_stats.
%
% INPUTS
%   A     : nSubj x nCh angles in RADIANS. NaN = missing at that channel.
%   g     : nSubj x 1 two-level group indicator. Accepts logical, numeric
%           (any two distinct values), categorical, string or cellstr. The
%           contrast is (second level) - (first level), where the levels are
%           taken in sorted/category order; extra.levels records which is
%           which. The statistic is invariant to that choice (see the group
%           symmetry check in test_circ_stats), so it only affects labelling.
%   Cnuis : nSubj x q numeric nuisance covariates, or [] (default). Must NOT
%           contain an intercept column — one is added. Rows with NaN in
%           Cnuis are dropped at every channel.
%   opts  : struct, all fields optional
%     .min_per_group   minimum non-missing subjects per group per channel
%                      before the channel is usable. DEFAULT 0 (off) — see
%                      the estimability note below.
%     .valid_channels  1 x nCh logical mask applied to the output. Channels
%                      outside the mask are returned as NaN.
%     .det_floor       determinant floor for the 2x2 residual covariance S.
%                      DEFAULT 1e-12. cos/sin are bounded in [-1, 1] so the
%                      entries of S are bounded by 1 and an absolute floor is
%                      meaningful.
%
% ESTIMABILITY GATE — READ THIS BEFORE WIRING A PERMUTATION LOOP.
% min_per_group defaults to OFF on purpose. The gate that decides which
% channels enter the analysis must be computed ONCE from the observed labels
% and then held FIXED across every permutation, otherwise the permuted
% max-statistic null is taken over a channel set that changes from
% permutation to permutation and the family-wise correction is no longer
% calibrated. The intended usage is therefore:
%
%   [~, ~, ~, ~, ex] = snpm_circ_hotelling(A, g, Cnuis);
%   ok = (ex.n1 >= 8) & (ex.n2 >= 8);          % ONCE, from observed labels
%   o.valid_channels = ok;
%   ... then pass o to every permuted call ...
%
% Inside the function only hard non-estimability is NaN'd (too few
% observations for the design, rank-deficient X on the retained rows, or a
% singular S), which cannot be avoided.
%
% MISSING DATA is handled by per-channel (pairwise) deletion, so N varies by
% channel. CONSEQUENCE: df2 IS RETURNED AS A PER-CHANNEL VECTOR, not the
% scalar that snpm_glm_stat returns, and p is computed channel-wise from it.
% Any caller that assumes a scalar df2 will be wrong here.
%
% SINGULAR S. S is 2x2 and is singular only when every retained angle at a
% channel is identical (resultant length exactly 1), which the estimability
% gate makes unreachable in practice. It is nevertheless checked explicitly:
% the channel is NaN'd with a recorded reason rather than pushed through a
% pseudo-inverse, which would return a finite but meaningless statistic.
%
% OUTPUTS (all 1 x nCh unless noted)
%   F     : Hotelling F statistic, non-negative; NaN where not estimable
%   p     : parametric upper-tail p from F(df1, df2), channel-wise
%   df1   : scalar, always 2
%   df2   : 1 x nCh vector, nu - 1 = N_ch - q - 3
%   extra : struct with
%           .T2       1 x nCh Hotelling T^2
%           .n1,.n2   1 x nCh retained subjects per group
%           .n        1 x nCh retained subjects total
%           .dcos,.dsin  1 x nCh adjusted group difference in the embedding
%           .levels   1 x 2 cell, the two group levels as text
%           .reason   1 x nCh cellstr, '' when the channel is usable
%
% Reference: Mardia & Jupp, Directional Statistics, §10.4 (two-sample tests
% on the embedded resultant); Method Spec Addendum Amendment 2 for the
% covariate-adjusted form.
%
% See also SNPM_CIRC_WATSONS_U2, SNPM_CIRC_CORR_ANGLINEAR, SNPM_CIRC_LINEARISE.

    if nargin < 3, Cnuis = []; end
    if nargin < 4, opts  = struct(); end

    [nSubj, nCh] = size(A);
    if numel(g) ~= nSubj
        error('snpm_circ_hotelling:groupSize', ...
            'g has %d entries but A has %d rows.', numel(g), nSubj);
    end

    min_per_group = getopt(opts, 'min_per_group', 0);
    det_floor     = getopt(opts, 'det_floor', 1e-12);
    valid_channels = getopt(opts, 'valid_channels', []);

    [gi, levels] = two_level_indicator(g);        % nSubj x 1 in {0,1}

    if isempty(Cnuis)
        Cnuis = zeros(nSubj, 0);
    end
    if size(Cnuis, 1) ~= nSubj
        error('snpm_circ_hotelling:covSize', ...
            'Cnuis has %d rows but A has %d rows.', size(Cnuis, 1), nSubj);
    end
    X    = [ones(nSubj, 1), Cnuis, gi];           % nSubj x (q+2)
    pcol = size(X, 2);
    c    = [zeros(pcol - 1, 1); 1];               % contrast selects g

    F   = nan(1, nCh);
    p   = nan(1, nCh);
    df1 = 2;
    df2 = nan(1, nCh);

    T2   = nan(1, nCh);
    dcos = nan(1, nCh);
    dsin = nan(1, nCh);
    reason = repmat({''}, 1, nCh);

    rowOK = all(isfinite(X), 2);                  % covariate/group validity
    okA   = ~isnan(A) & rowOK;                    % nSubj x nCh usable cells

    n1 = sum(okA & (gi == 1), 1);
    n2 = sum(okA & (gi == 0), 1);
    nn = n1 + n2;

    Zc = cos(A); Zs = sin(A);

    % Channels that use every covariate-valid row can be solved in one batch.
    fullCols    = find(all(okA == rowOK, 1) & nn > 0);
    partialCols = find(~(all(okA == rowOK, 1) & nn > 0));

    if ~isempty(fullCols)
        rows = rowOK;
        Xf   = X(rows, :);
        [Ff, T2f, nuf, dcf, dsf, whyf] = hotelling_block( ...
            Zc(rows, fullCols), Zs(rows, fullCols), Xf, c, det_floor);
        F(fullCols)    = Ff;
        T2(fullCols)   = T2f;
        df2(fullCols)  = nuf - 1;
        dcos(fullCols) = dcf;
        dsin(fullCols) = dsf;
        reason(fullCols) = whyf;
    end

    for k = 1:numel(partialCols)
        ch   = partialCols(k);
        rows = okA(:, ch);
        if ~any(rows)
            reason{ch} = 'no observations';
            continue
        end
        Xc = X(rows, :);
        [Fc, T2c, nuc, dcc, dsc, whyc] = hotelling_block( ...
            Zc(rows, ch), Zs(rows, ch), Xc, c, det_floor);
        F(ch)    = Fc;
        T2(ch)   = T2c;
        df2(ch)  = nuc - 1;
        dcos(ch) = dcc;
        dsin(ch) = dsc;
        reason(ch) = whyc;
    end

    % ---- estimability gates -------------------------------------------
    if min_per_group > 0
        bad = (n1 < min_per_group) | (n2 < min_per_group);
        for ch = find(bad)
            if isempty(reason{ch})
                reason{ch} = sprintf('below min_per_group (%d/%d)', n1(ch), n2(ch));
            end
        end
        F(bad) = NaN; T2(bad) = NaN; df2(bad) = NaN;
    end
    if ~isempty(valid_channels)
        if numel(valid_channels) ~= nCh
            error('snpm_circ_hotelling:maskSize', ...
                'opts.valid_channels has %d entries, expected %d.', ...
                numel(valid_channels), nCh);
        end
        masked = ~logical(valid_channels(:)');
        for ch = find(masked)
            if isempty(reason{ch}), reason{ch} = 'excluded by valid_channels'; end
        end
        F(masked) = NaN; T2(masked) = NaN; df2(masked) = NaN;
    end

    ok = isfinite(F) & isfinite(df2) & df2 > 0;
    p(ok)   = fcdf(F(ok), df1, df2(ok), 'upper');
    F(~ok)  = NaN;
    T2(~ok) = NaN;
    df2(~ok) = NaN;

    extra = struct('T2', T2, 'n1', n1, 'n2', n2, 'n', nn, ...
        'dcos', dcos, 'dsin', dsin, 'levels', {levels}, 'reason', {reason});
end

% ======================================================================
function [F, T2, nu, dcos, dsin, why] = hotelling_block(Zc, Zs, X, c, det_floor)
% Covariate-adjusted two-sample Hotelling on a set of channels that share the
% same design X (i.e. the same retained rows). Zc/Zs are nObs x m.
    m   = size(Zc, 2);
    F    = nan(1, m);
    T2   = nan(1, m);
    dcos = nan(1, m);
    dsin = nan(1, m);
    why  = repmat({''}, 1, m);

    nObs  = size(X, 1);
    rankX = rank(X);
    nu    = nObs - rankX;

    if rankX < size(X, 2)
        why(:) = {'design rank deficient on retained rows'};
        return
    end
    if nu < 2
        why(:) = {sprintf('too few observations (N=%d, rank(X)=%d)', nObs, rankX)};
        return
    end

    % (X'X)^-1 is formed explicitly because the contrast variance factor
    % v = c'(X'X)^-1 c needs the inverse itself, not just a solve. X is
    % (q+2)-wide and full rank here, so this is cheap and well conditioned.
    XtXinv = (X' * X) \ eye(size(X, 2));
    Bc = XtXinv * (X' * Zc);                                % p x m
    Bs = XtXinv * (X' * Zs);
    Rc = Zc - X * Bc;
    Rs = Zs - X * Bs;

    Scc = sum(Rc .^ 2, 1) / nu;                             % 1 x m
    Sss = sum(Rs .^ 2, 1) / nu;
    Scs = sum(Rc .* Rs, 1) / nu;

    dc = c' * Bc;                                           % 1 x m
    ds = c' * Bs;
    v  = c' * XtXinv * c;                                   % scalar

    detS = Scc .* Sss - Scs .^ 2;
    sing = ~(detS > det_floor);                             % catches NaN too

    quad = (dc .^ 2 .* Sss - 2 * dc .* ds .* Scs + ds .^ 2 .* Scc) ./ detS;
    t2   = quad / v;
    f    = t2 * (nu - 1) / (2 * nu);

    f(sing)  = NaN;
    t2(sing) = NaN;
    why(sing) = {'residual covariance singular (all angles identical?)'};

    F    = f;
    T2   = t2;
    dcos = dc;
    dsin = ds;
end

% ======================================================================
function [gi, levels] = two_level_indicator(g)
% Map an arbitrary two-level grouping variable to a 0/1 column.
    g = g(:);
    if islogical(g)
        gi = double(g);
        levels = {'false', 'true'};
        return
    end
    if isnumeric(g)
        u = unique(g(~isnan(g)));
        if numel(u) ~= 2
            error('snpm_circ_hotelling:twoLevels', ...
                'g must have exactly two distinct levels (found %d).', numel(u));
        end
        gi = nan(size(g));
        gi(g == u(1)) = 0;
        gi(g == u(2)) = 1;
        levels = {num2str(u(1)), num2str(u(2))};
        return
    end
    if iscellstr(g) || isstring(g) || ischar(g)
        gc = categorical(cellstr(g));
    elseif iscategorical(g)
        gc = removecats(g);
    else
        error('snpm_circ_hotelling:groupType', ...
            'Unsupported group variable class ''%s''.', class(g));
    end
    cats = categories(gc);
    if numel(cats) ~= 2
        error('snpm_circ_hotelling:twoLevels', ...
            'g must have exactly two distinct levels (found %d).', numel(cats));
    end
    gi = double(gc == cats{2});
    gi(isundefined(gc)) = NaN;
    levels = cats(:)';
end

% ======================================================================
function v = getopt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
