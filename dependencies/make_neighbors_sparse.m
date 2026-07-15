function sparse_channel_adjacency_matrix = make_neighbors_sparse(neighbors,num_channels);

channel_adjacency_matrix = zeros(num_channels, num_channels);
 
for channel_num = 1:num_channels
    channel_neighbors = neighbors(channel_num, ~isnan(neighbors(channel_num, :)));
   
    channel_adjacency_matrix(channel_num, channel_neighbors) = 1;
    channel_adjacency_matrix(channel_neighbors, channel_num) = 1;
end
 
sparse_channel_adjacency_matrix = sparse(channel_adjacency_matrix);