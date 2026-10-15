function TFCEdata = ClusterEnhancement(unenhanced_data,sparse_channel_adjacency_matrix,E,H,dh,branch)
% Threshold-Free Cluster Enhancement over a channel-adjacency graph.
%   TFCEdata = ClusterEnhancement(map, sparse_adj)            % E=.5, H=2, dh=.1
%   TFCEdata = ClusterEnhancement(map, sparse_adj, E, H)
%   TFCEdata = ClusterEnhancement(map, sparse_adj, E, H, dh)
%   TFCEdata = ClusterEnhancement(map, sparse_adj, E, H, dh, branch)
%
% Integrates cluster_size^E * h^H * dh over thresholds h, separately for the
% positive and the negative side of the map, and returns pos - neg.
%
% dh (optional, 5th argument) is the integration step. It defaults to 0.1,
% which is the value every pre-2026-08 caller used and which is appropriate
% for t- and F-scale maps. It is load-bearing for **Watson's U^2 only**: U^2
% lives on roughly [0, 0.5], so dh=0.1 leaves about 5 integration levels and
% most channels enhanced to exactly 0 — pass dh=0.005 (~91 levels) for U^2
% maps. Do NOT rescale dh for F maps: an F-quantile-derived step is coarser
% than the 0.1 default, not finer.
%
% The negative-enhancement half is skipped automatically when the map has no
% negative entries (every non-NaN value >= 0), which is the case for every
% omnibus statistic in this toolbox (F, Hotelling T^2 -> F, Wald, U^2). On a
% non-negative map the negative half contributes exactly zero, so skipping it
% is numerically identical and about 70% faster. This is deliberately an
% auto-detect and not a flag, so existing callers get the saving without a
% signature change.
%
% TWO BRANCHES, ONE INTEGRAL (round-4 spec, approved 2026-10-14).
% The dense branch builds the threshold ladder dh:dh:max|map|+dh and finds the
% clusters at every level; memory is (levels x channels) doubles, so a single
% legitimate but enormous statistic (a t or F with 1 error df is Cauchy-like;
% 9.3e6 was observed) asked for 130+ GB. The dense branch is used, unchanged
% and bit-identical, whenever the ladder has at most MAX_DENSE_LEVELS = 1e4
% levels (|stat| below about 1000 at dh = 0.1; the golden maps need 3,121).
% Above that the EVENT branch computes the same discretised integral with
% levels i*dh, i = 1..K (K = the dense ladder's level count): clusters are
% rebuilt only at the distinct per-channel level counts (at most one cluster
% search per channel and side), and sum(i^H) over each run of levels is taken
% in closed form. Its cost does not depend on the size of the statistic.
% It differs from the dense branch only by rounding (dense uses MATLAB's
% colon-built level values, the event branch i*dh, which the Python port
% reproduces exactly). The event branch accepts integer H in 0..3 only
% (H = 2 everywhere in this toolbox).
%
% branch (optional, 6th argument, for tests): 'auto' (default), 'dense' or
% 'event'. Production callers never pass it.
%
% Validated by isequaln against stored golden maps (test_fixtures/
% ce_golden_178.mat). 204 of the 206 were captured before the 2026-08 edit;
% entries 199 (all-zero) and 203 (zero+NaN) were replaced on 2026-10-13 after
% the sum-dimension fix below (sum(dhtot,1,...)), which makes a map whose
% every |value| is under the first threshold return 1 x numch instead of a
% scalar or a short vector.

%%default
if nargin < 3 % default E H
    E  = .5;
    H  = 2;
end
if nargin < 5 || isempty(dh)
    dh = 0.1;
end
if nargin < 6 || isempty(branch)
    branch = 'auto';
end

MAX_DENSE_LEVELS = 1e4;

top = max(abs(unenhanced_data));
switch branch
    case 'auto'
        use_event = isfinite(top) && ladder_levels(top, dh) > MAX_DENSE_LEVELS;
    case 'dense'
        use_event = false;
    case 'event'
        use_event = true;
    otherwise
        error('ClusterEnhancement:badBranch', 'branch must be ''auto'', ''dense'' or ''event''.');
end
if use_event
    TFCEdata = enhance_event(unenhanced_data, sparse_channel_adjacency_matrix, E, H, dh, top);
    return
end

numch=length(unenhanced_data);
hrange=dh:dh:max(abs(unenhanced_data))+dh;

nani=isnan(unenhanced_data);

% A map with no negative values has an identically-zero negative half.
skip_negative_half = all(unenhanced_data >= 0 | nani);

%%% positive enhancement
temp_data=unenhanced_data;
temp_data(temp_data<0)=0;
dhtot=NaN(length(hrange),numch);
for hindex = 1:length(hrange)
    clusters_above_threshold = snpm_find_clusters_graphalgs(temp_data,hrange(hindex),sparse_channel_adjacency_matrix);
    for c = 1:length(clusters_above_threshold)
        cluster_size=length(clusters_above_threshold{c});
        dhtot(hindex,clusters_above_threshold{c}) = ...
            cluster_size^E*hrange(hindex)^H*dh;
    end
end

enhanced_pos=sum(dhtot,1,'omitnan');

%%% negative enhancement
if skip_negative_half
    enhanced_neg=0;
else
    temp_data=-unenhanced_data;
    temp_data(temp_data<0)=0;
    dhtot=NaN(length(hrange),numch);

    for hindex = 1:length(hrange)
        clusters_above_threshold = snpm_find_clusters_graphalgs(temp_data,hrange(hindex),sparse_channel_adjacency_matrix);
        for c = 1:length(clusters_above_threshold)
            cluster_size=length(clusters_above_threshold{c});
            dhtot(hindex,clusters_above_threshold{c}) = ...
                cluster_size^E*hrange(hindex)^H*dh;
        end
    end

    enhanced_neg=sum(dhtot,1,'omitnan');
end

TFCEdata=enhanced_pos-enhanced_neg;
TFCEdata(nani)=NaN;
end

% ---------------------------------------------------------------------------
function K = ladder_levels(top, dh)
% numel(dh:dh:top+dh) without building the vector: MATLAB's colon count rule
% (round to the nearest step, then drop the last if it overshoots the end by
% more than 2*eps of the larger end point).
    a = dh; b = top + dh;
    tol = 2 * eps * max(abs(a), abs(b));
    nn = floor((b - a) / dh + 0.5);
    if a + nn * dh - b > tol, nn = nn - 1; end
    K = nn + 1;
end

function out = enhance_event(v, adj, E, H, dh, top)
% Event branch: same integral as the dense loop, levels i*dh (i = 1..K).
    if ~(H == round(H) && H >= 0 && H <= 3)
        error('ClusterEnhancement:eventH', ['The event branch of ClusterEnhancement needs an ' ...
            'integer H in 0..3 (got %g); the dense ladder would need more than 1e4 levels here.'], H);
    end
    v = reshape(v, 1, []);
    nani = isnan(v);
    if isnan(top)                                       % all-NaN map, as the dense branch
        out = NaN(1, numel(v));
        return
    end
    if ~isfinite(top)
        error('ClusterEnhancement:nonFinite', ['ClusterEnhancement: the map holds a non-finite ' ...
            'value (%g); TFCE is undefined.'], top);
    end
    out = side(max(v, 0));                              % NaN -> max(NaN,0) = 0, never active
    if ~all(v >= 0 | nani)
        out = out - side(max(-v, 0));
    end
    out(nani) = NaN;

    function s = side(x)
        x(isnan(x)) = 0;
        % g(c) = number of levels i*dh (i = 1..K) strictly below x(c), i.e. the
        % levels at which channel c is supra-threshold (the cluster finder
        % keeps x > h). ceil(x/dh)-1 is exact up to rounding; correct it
        % against the same i*dh expression the cluster search uses below.
        % Not clipped to K (as the Python port): x <= top < K*dh, so g < K.
        % Two fixed passes, up then down on the updated g, exactly as the
        % Python port (tfce._enhance_side_event). No convergence error: once
        % x/dh passes 2^53 (|stat| above ~9e14 at dh 0.1) (g+1)*dh cannot
        % move and g stays at ceil(x/dh)-1, which is then the count to
        % within the resolution of a double. Such values are legitimate
        % (zero-variance rule 2 admits F up to ~1e20).
        g = max(ceil(x / dh) - 1, 0);
        for it = 1:2
            up = (g + 1) * dh < x;
            g(up) = g(up) + 1;
            dn = g > 0 & ~(g * dh < x);
            g(dn) = g(dn) - 1;
        end
        s = zeros(1, numel(x));
        br = sort(unique(g(g > 0)), 'descend');
        for j = 1:numel(br)
            hi = br(j);
            if j < numel(br), lo = br(j + 1) + 1; else, lo = 1; end
            % Levels lo..hi share one supra-threshold set, {c : g(c) >= hi}.
            % Clusters are found on that indicator (as the Python port), not
            % by re-thresholding x at hi*dh, so the set is the one g defines
            % even where x/dh is beyond double resolution.
            w = dh^(H + 1) * power_run(lo, hi, H);       % sum_{i=lo}^{hi} (i*dh)^H * dh
            cl = snpm_find_clusters_graphalgs(double(g >= hi), 0.5, adj);
            for k = 1:numel(cl)
                s(cl{k}) = s(cl{k}) + numel(cl{k})^E * w;
            end
        end
    end
end

function y = power_run(lo, hi, H)
% sum_{i=lo}^{hi} i^H for integer H in 0..3, as sum_{k=0}^{m} (lo+k)^H with
% m = hi - lo, expanded binomially over Faulhaber sums S_j(m) = sum_{k=0}^{m}
% k^j. Every term is non-negative, so a short run far up the ladder keeps full
% precision (the textbook S(hi) - S(lo-1) loses about eps*hi relative: 1e-6
% near 1e11 levels). Same operations in the same order as the Python port's
% tfce._power_run, so the two agree bit for bit.
    m = hi - lo;
    faulhaber = [m + 1, m * (m + 1) / 2, m * (m + 1) * (2 * m + 1) / 6, (m * (m + 1) / 2)^2];
    binom = {1, [1 1], [1 2 1], [1 3 3 1]};
    b = binom{H + 1};
    y = 0;
    for j = 0:H
        y = y + b(j + 1) * lo^(H - j) * faulhaber(j + 1);
    end
end
