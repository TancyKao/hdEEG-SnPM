function [Clusters] = snpm_lmm_cluster(power, meta, spec, neighbors, alpha, permutations)
% Cluster-based permutation correction for the per-channel LMM, following
% Stephan et al. 2021: neighbouring channels with model p < alpha form a
% cluster; the cluster statistic is the mean Wald inside the cluster
% (sum(Wald)/n_channels); the null is the distribution of the maximum
% cluster statistic over permutations; clusters above the 95th percentile
% (i.e. permutation p <= alpha) are significant.
%
% Mirrors the output contract of snpm_cluster_analysis: a struct array with
% .channels and .p so the downstream code (SnPMsigch derivation, Excel and
% HTML reports) works unchanged.
%
% INPUTS
%   power        : trials x nCh channel matrix
%   meta         : trial-level table of model variables
%   spec         : LMM spec (see snpm_lmm_fit); spec.perm selects the
%                  permutation scheme ('within_subject' | 'group_label')
%   neighbors    : channels x maxNeighbors adjacency (NaN padded)
%   alpha        : cluster-forming AND significance level (default 0.05)
%   permutations : number of permutations
%
% OUTPUT
%   Clusters : struct array with fields .channels, .mass, .p, .threshold,
%              .permutations

    if nargin < 5 || isempty(alpha), alpha = 0.05; end

    nCh = size(power, 2); %#ok<NASGU>
    sparse_adj = make_neighbors_sparse(neighbors, size(neighbors, 1));
    possible_permutations = permutations;

    % ---- real clusters ----
    [~, real_wald, real_p] = snpm_lmm_fit(power, meta, spec);
    real_clusters = local_find_pclusters(real_p, real_wald, alpha, sparse_adj);

    % ---- null distribution of the max cluster mass ----
    max_cluster_mass = zeros(possible_permutations, 1);
    for permIndex = 1:possible_permutations
        if mod(permIndex, 100) == 0
            disp([num2str(permIndex), ' out of ', num2str(possible_permutations), ' LMM cluster permutations completed...']);
        end
        meta_perm = snpm_lmm_permute_meta(meta, spec, spec.perm);
        [~, pwald, pp] = snpm_lmm_fit(power, meta_perm, spec);
        pc = local_find_pclusters(pp, pwald, alpha, sparse_adj);
        if isempty(pc)
            max_cluster_mass(permIndex) = 0;
        else
            max_cluster_mass(permIndex) = max([pc.mass]);
        end
    end
    max_cluster_mass = sort(max_cluster_mass, 'descend');

    % ---- assemble Clusters with permutation p-values ----
    if isempty(real_clusters)
        Clusters = struct('channels', {[]}, 'mass', {0}, 'p', {1}, ...
            'threshold', {alpha}, 'permutations', {possible_permutations});
        return
    end

    Clusters = struct([]);
    for cli = 1:numel(real_clusters)
        Clusters(cli, 1).channels = real_clusters(cli).channels;
        Clusters(cli, 1).mass     = real_clusters(cli).mass;
        rank = find(max_cluster_mass < real_clusters(cli).mass, 1, 'first');
        if isempty(rank)
            Clusters(cli, 1).p = 1;
        else
            Clusters(cli, 1).p = rank / possible_permutations;
        end
        Clusters(cli, 1).threshold    = alpha;
        Clusters(cli, 1).permutations = possible_permutations;
    end
end

function clusters = local_find_pclusters(p_model, wald, alpha, sparse_adj)
% Connected components of channels with model p < alpha, with each cluster's
% mass = sum(Wald)/n_channels. Reuses snpm_find_clusters_graphalgs by feeding
% it a 0/1 significance mask thresholded at 0.5.
    mask = double(p_model < alpha);
    mask(isnan(p_model)) = 0;
    cc = snpm_find_clusters_graphalgs(mask, 0.5, sparse_adj);
    clusters = struct([]);
    k = 0;
    for c = 1:numel(cc)
        ch = cc{c};
        if isempty(ch)
            continue
        end
        k = k + 1;
        clusters(k, 1).channels = ch;
        clusters(k, 1).mass     = sum(wald(ch), 'omitnan') / numel(ch);
    end
end
