function clusters = snpm_find_clusters_graphalgs(t_values, threshold, sparse_adjacency_matrix)
    if threshold > 0
        significant_vertices_logical = t_values > threshold;
    else
        significant_vertices_logical = t_values < threshold;
    end
    
    significant_vertices                = find(significant_vertices_logical);
    significant_sparse_adjacency_matrix = sparse_adjacency_matrix(significant_vertices_logical, significant_vertices_logical);
    
    % strongly connected component - every vertexd is reachable from every other vertex
    %[num_clusters cluster_nums] = graphalgs('scc', 0, 0, significant_sparse_adjacency_matrix);
    
    % using Bioinformatics Toolbox 4.9 
    % added by tancy 
    % replace-code of graphconncomp 28112022

    %function [S,C] = conncomp(G)
      % CONNCOMP Drop in replacement for graphconncomp.m from the bioinformatics
      % toobox. G is an n by n adjacency matrix, then this identifies the S
      % connected components C. This is also an order of magnitude faster.
      %
      % [S,C] = conncomp(G)
      %
      % Inputs:
      %   G  n by n adjacency matrix
      % Outputs:
      %   S  scalar number of connected components
      %   C  
   %   [p,~,r] = dmperm(G+speye(size(G)));
   %   S = numel(r)-1;
   %   C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
   %   C(p) = C;
    %end


    %[num_clusters cluster_nums] = graphconncomp(significant_sparse_adjacency_matrix);
    G = significant_sparse_adjacency_matrix+speye(size(significant_sparse_adjacency_matrix));

    [p,~,r] = dmperm(G);
    S = numel(r)-1;
    C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
    C(p) = C;
    num_clusters = S;
    cluster_nums = C;


    % view graph 
    %h = view(biograph(significant_sparse_adjacency_matrix));
    %colors = jet(num_clusters);
    %for i = 1:numel(h.nodes)
    %  h.Nodes(i).Color = colors(cluster_nums(i),:);
    %end

    clusters = cell(1, num_clusters);
    
    for cluster_num = 1:num_clusters
        clusters{cluster_num} = significant_vertices(cluster_nums == cluster_num);
    end
end
