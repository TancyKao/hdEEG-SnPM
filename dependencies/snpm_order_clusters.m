function [order, sgn] = snpm_order_clusters(Clusters, T, alpha)
%SNPM_ORDER_CLUSTERS  Display order + sign of the SIGNIFICANT clusters.
%   [ORDER, SGN] = SNPM_ORDER_CLUSTERS(CLUSTERS, T) returns the indices of the
%   clusters that the report lists -- those with a cluster-level p <= ALPHA and
%   at least one channel -- in the order they are numbered, plus each one's
%   sign. ORDER(j) is the index into CLUSTERS of the cluster labelled "j"; SGN(j)
%   is +1 / -1 (the sign of the mean statistic over the cluster's channels).
%
%   ONE numbering rule, used by BOTH the HTML cluster table
%   (generateAnalysisReport) and the numbers drawn on the Cluster PNG
%   (TopoplotSignificant_single), so a number in the table always means the same
%   cluster as that number on the map:
%
%     * SIGNED statistic (any negative value anywhere in T -- t, r, LMM
%       continuous, post-hoc contrast): positive clusters first, ordered by
%       ascending cluster p, then negative clusters, ordered by ascending
%       cluster p. This matches the "Positive clusters (t > 0)" /
%       "Negative clusters (t < 0)" dividers in the table.
%     * UNSIGNED statistic (F, Hotelling T^2, Watson U^2 -- T has no negative
%       entry): ascending cluster p only. Every cluster then counts as positive,
%       so the two branches agree and no sign divider is drawn.
%
%   Ties in p keep the engine's own cluster order (MATLAB's sort is stable).
%
%   INPUTS
%     Clusters : struct array with fields .channels (channel indices) and .p
%     T        : per-channel statistic vector (Clusters index into it)
%     alpha    : significance level for the cluster p (default 0.05)
%
%   An empty / malformed CLUSTERS gives empty outputs.

    if nargin < 3 || isempty(alpha), alpha = 0.05; end
    order = []; sgn = [];
    if isempty(Clusters) || ~isstruct(Clusters), return; end
    if ~isfield(Clusters, 'channels') || ~isfield(Clusters, 'p'), return; end

    T = double(T(:).');
    signed = any(T < 0);

    n = numel(Clusters);
    keep = false(1, n); pv = nan(1, n); sv = ones(1, n);
    for i = 1:n
        ch = Clusters(i).channels;
        p  = Clusters(i).p;
        if isempty(ch) || isempty(p), continue; end
        p = p(1);
        if ~isfinite(p) || p > alpha, continue; end
        ch = double(ch(:).');
        ch = ch(ch >= 1 & ch <= numel(T));
        if isempty(ch), continue; end
        keep(i) = true; pv(i) = p;
        if signed
            m = mean(T(ch), 'omitnan');
            if isfinite(m) && m < 0, sv(i) = -1; end
        end
    end

    idx = find(keep);
    if isempty(idx), return; end
    if signed
        pos = idx(sv(idx) > 0); neg = idx(sv(idx) < 0);
        [~, op] = sort(pv(pos)); [~, on] = sort(pv(neg));
        order = [pos(op), neg(on)];
    else
        [~, o] = sort(pv(idx));
        order = idx(o);
    end
    sgn = sv(order);
end
