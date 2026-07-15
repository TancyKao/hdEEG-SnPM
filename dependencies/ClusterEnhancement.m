function TFCEdata = ClusterEnhancement(unenhanced_data,sparse_channel_adjacency_matrix,E,H)
%%default
if nargin < 3 % default E H
    E  = .5;
    H  = 2;
end

dh = 0.1;
numch=length(unenhanced_data);
hrange=dh:dh:max(abs(unenhanced_data))+dh;

nani=isnan(unenhanced_data);
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

TFCEdata=enhanced_pos-enhanced_neg;
TFCEdata(nani)=NaN;