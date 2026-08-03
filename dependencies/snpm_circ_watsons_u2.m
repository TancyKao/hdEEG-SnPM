function [U2, p, extra] = snpm_circ_watsons_u2(A, g, opts)
% Watson's two-sample U^2 statistic, per channel, in ECDF form.
%
%   [U2, p, extra] = snpm_circ_watsons_u2(A, g)
%   [...] = snpm_circ_watsons_u2(A, g, opts)
%
% This is the secondary Tier-2 circular statistic. Where the Hotelling test
% asks whether the mean resultant vectors differ, U^2 is an omnibus test of
% whether the two angular DISTRIBUTIONS differ at all — it will pick up a
% difference in spread or shape that leaves the mean direction untouched.
% It is rotation-invariant and it handles ties natively.
%
% STATISTIC. U^2 is the integral of the squared, mean-centred difference of
% the two empirical distribution functions against the POOLED empirical
% measure. Over the pooled sample sorted by angle, let i_k and j_k be the
% counts of group-1 and group-2 members at or below the k-th DISTINCT pooled
% value v_k, let m_k be that value's multiplicity, and d_k = i_k/n1 - j_k/n2:
%
%   U^2 = (n1 n2 / N^2) * [ sum_k m_k d_k^2 - (sum_k m_k d_k)^2 / N ]
%
% With no ties every m_k = 1 and this is the familiar running-count form
%   U^2 = (n1 n2 / N^2) * [ sum d^2 - (sum d)^2 / N ]
% evaluated at every observation.
%
% TIES ARE HANDLED BY MULTIPLICITY, NOT BY TIE-BREAKING — this matters.
% The textbook running-count form evaluates d after every individual
% observation, so inside a block of tied angles its value depends on the
% arbitrary order in which the tied observations happen to be listed. On
% whole-degree data (Luna and TurtleWave both export whole degrees in at
% least one column) that ambiguity reaches ~2e-2 in U^2, which is over 10% of
% the 0.187 critical value at alpha = .05 — not acceptable. Evaluating the
% distribution functions only at the DISTINCT values and weighting each by
% its multiplicity is the correct reading of the defining integral (the
% pooled empirical measure puts mass m_k/N at v_k, not m_k separate unit
% masses at m_k distinct places) and it removes the ambiguity entirely: the
% result no longer depends on how tied observations are ordered.
%
% OTHER IMPLEMENTATION NOTES
%  * The pooled sort order and the tie-block structure are LABEL-FREE, so they
%    are computed once and permutation only reshuffles labels along that fixed
%    order. That is both a large speed win and a structural guarantee that the
%    permuted statistic is the same function of the data with a different
%    label vector, i.e. the permutation test stays exact.
%  * MATLAB sorts NaN last. Missing subjects therefore sit at the tail of the
%    order, where the running counts have already saturated and d_k = 0, so
%    they contribute nothing to either sum. No extra masking is needed —
%    PROVIDED the division by N uses the per-channel finite count, which it
%    does (n1 and n2 are the retained counts, not size(A,1)). NaN ~= NaN, so
%    each missing entry is its own tie block and cannot absorb a real one.
%  * NOT batched over permutations on purpose: a naive P x N x nCh array is
%    about 1 GB at P = 5000, while the serial loop already runs the whole
%    5000-permutation correction in a few seconds.
%
% INPUTS
%   A    : nSubj x nCh angles in RADIANS. NaN = missing at that channel.
%          Angles outside [0, 2*pi) are wrapped internally, so the statistic
%          is invariant to adding any constant rotation.
%   g    : nSubj x 1 two-level group indicator (logical / numeric / string /
%          cellstr / categorical). Level order only affects labelling; U^2 is
%          symmetric under swapping the two groups.
%   opts : struct, optional
%     .min_per_group   minimum retained subjects per group per channel.
%                      DEFAULT 0 (off) — see the estimability note in
%                      SNPM_CIRC_HOTELLING; the gate belongs to the caller and
%                      must be fixed across permutations.
%     .valid_channels  1 x nCh logical mask; channels outside are NaN.
%
% OUTPUTS (1 x nCh)
%   U2    : Watson's U^2, non-negative; NaN where not estimable
%   p     : asymptotic upper-tail p from WATSON_U2_PVALUE
%   extra : .n1, .n2, .n (retained counts), .reason (cellstr)
%
% Reference: Watson (1962), Biometrika; Zar, Biostatistical Analysis, Ch. 27.
%
% See also WATSON_U2_PVALUE, SNPM_CIRC_HOTELLING.

    if nargin < 3, opts = struct(); end

    [nSubj, nCh] = size(A);
    if numel(g) ~= nSubj
        error('snpm_circ_watsons_u2:groupSize', ...
            'g has %d entries but A has %d rows.', numel(g), nSubj);
    end

    min_per_group  = getopt(opts, 'min_per_group', 0);
    valid_channels = getopt(opts, 'valid_channels', []);

    gi = two_level_indicator(g);                 % nSubj x 1 in {0,1}, NaN ok

    % Wrap to [0, 2*pi) so the sort is over a well-defined circular order.
    Aw = mod(A, 2 * pi);
    Aw(~isfinite(A)) = NaN;
    Aw(~isfinite(gi), :) = NaN;                  % subjects with no group label

    % Pooled sort per channel. NaN goes last (MATLAB default). The sort order
    % and the tie-block structure below depend only on the angles, never on
    % the labels.
    [As, ord] = sort(Aw, 1, 'ascend');           % nSubj x nCh

    isv = ~isnan(As);                            % nSubj x nCh, trailing false

    % Last position of each tie block. NaN ~= NaN is true, so every missing
    % entry ends its own block and cannot merge with a real value.
    if nSubj > 1
        isEnd = [As(1:end-1, :) ~= As(2:end, :); true(1, nCh)];
    else
        isEnd = true(nSubj, nCh);
    end

    G1 = (gi == 1);
    G2 = (gi == 0);
    in1 = G1(ord) & isv;
    in2 = G2(ord) & isv;

    ik = cumsum(in1, 1);
    jk = cumsum(in2, 1);

    n1 = ik(end, :);
    n2 = jk(end, :);
    nn = n1 + n2;

    usable = (n1 > 0) & (n2 > 0) & (nn > 2);

    U2 = nan(1, nCh);
    reason = repmat({''}, 1, nCh);
    for ch = find(~usable)
        reason{ch} = sprintf('degenerate split (n1=%d, n2=%d)', n1(ch), n2(ch));
    end

    if any(usable)
        u = usable;
        d = ik(:, u) ./ n1(u) - jk(:, u) ./ n2(u);   % trailing rows are 0

        % Multiplicity weighting: replace d inside a tie block by its value at
        % the END of the block. Summing the filled vector over all positions
        % is exactly sum_k m_k * d_k over distinct values, with no dependence
        % on how the tied observations happen to be ordered. The scan is
        % nSubj (~tens) vector operations across all channels at once.
        eu = isEnd(:, u);
        for k = nSubj-1:-1:1
            carry = ~eu(k, :);
            d(k, carry) = d(k+1, carry);
        end

        sd  = sum(d, 1);
        sd2 = sum(d .^ 2, 1);
        U2(u) = (n1(u) .* n2(u) ./ nn(u) .^ 2) .* (sd2 - sd .^ 2 ./ nn(u));
    end

    % ---- estimability gates -------------------------------------------
    if min_per_group > 0
        bad = (n1 < min_per_group) | (n2 < min_per_group);
        for ch = find(bad)
            if isempty(reason{ch})
                reason{ch} = sprintf('below min_per_group (%d/%d)', n1(ch), n2(ch));
            end
        end
        U2(bad) = NaN;
    end
    if ~isempty(valid_channels)
        if numel(valid_channels) ~= nCh
            error('snpm_circ_watsons_u2:maskSize', ...
                'opts.valid_channels has %d entries, expected %d.', ...
                numel(valid_channels), nCh);
        end
        masked = ~logical(valid_channels(:)');
        for ch = find(masked)
            if isempty(reason{ch}), reason{ch} = 'excluded by valid_channels'; end
        end
        U2(masked) = NaN;
    end

    p = watson_u2_pvalue(U2);

    extra = struct('n1', n1, 'n2', n2, 'n', nn, 'reason', {reason});
end

% ======================================================================
function gi = two_level_indicator(g)
% Map an arbitrary two-level grouping variable to a 0/1 column (NaN allowed).
    g = g(:);
    if islogical(g)
        gi = double(g);
        return
    end
    if isnumeric(g)
        u = unique(g(~isnan(g)));
        if numel(u) ~= 2
            error('snpm_circ_watsons_u2:twoLevels', ...
                'g must have exactly two distinct levels (found %d).', numel(u));
        end
        gi = nan(size(g));
        gi(g == u(1)) = 0;
        gi(g == u(2)) = 1;
        return
    end
    if iscellstr(g) || isstring(g) || ischar(g)
        gc = categorical(cellstr(g));
    elseif iscategorical(g)
        gc = removecats(g);
    else
        error('snpm_circ_watsons_u2:groupType', ...
            'Unsupported group variable class ''%s''.', class(g));
    end
    cats = categories(gc);
    if numel(cats) ~= 2
        error('snpm_circ_watsons_u2:twoLevels', ...
            'g must have exactly two distinct levels (found %d).', numel(cats));
    end
    gi = double(gc == cats{2});
    gi(isundefined(gc)) = NaN;
end

% ======================================================================
function v = getopt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f))
        v = s.(f);
    else
        v = d;
    end
end
