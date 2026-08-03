function TFCEdata = ClusterEnhancement(unenhanced_data,sparse_channel_adjacency_matrix,E,H,dh)
% Threshold-Free Cluster Enhancement over a channel-adjacency graph.
%   TFCEdata = ClusterEnhancement(map, sparse_adj)            % E=.5, H=2, dh=.1
%   TFCEdata = ClusterEnhancement(map, sparse_adj, E, H)
%   TFCEdata = ClusterEnhancement(map, sparse_adj, E, H, dh)
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
% Pinned by test_cluster_enhancement_identity.m (isequaln against goldens
% captured before the 2026-08 edit).

%%default
if nargin < 3 % default E H
    E  = .5;
    H  = 2;
end
if nargin < 5 || isempty(dh)
    dh = 0.1;
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

enhanced_pos=sum(dhtot,'omitnan');

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

    enhanced_neg=sum(dhtot,'omitnan');
end

TFCEdata=enhanced_pos-enhanced_neg;
TFCEdata(nani)=NaN;
