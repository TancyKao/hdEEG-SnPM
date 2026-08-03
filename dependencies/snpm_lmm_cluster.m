function [Clusters] = snpm_lmm_cluster(power, meta, spec, neighbors, alpha, permutations, evaluable)
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
%   evaluable    : (optional) 1 x nCh logical, TRUE for channels that can be
%                  evaluated in EVERY permutation. Defaults to
%                  all(isfinite(power), 1). Channels outside the mask form no
%                  cluster and contribute no mass, in the observed map and in
%                  every permuted map alike, so the observed cluster statistic
%                  and its null are built on one channel set. See
%                  snpm_lmm_TFCE for why per-channel row deletion otherwise
%                  breaks exchangeability.
%
% OUTPUT
%   Clusters : struct array with fields .channels, .mass, .p, .threshold,
%              .permutations

    if nargin < 5 || isempty(alpha), alpha = 0.05; end

    nCh = size(power, 2);
    sparse_adj = make_neighbors_sparse(neighbors, size(neighbors, 1));
    possible_permutations = permutations;

    if nargin < 7 || isempty(evaluable)
        evaluable = all(isfinite(power), 1);
    end
    evaluable = reshape(logical(evaluable), 1, []);
    if numel(evaluable) ~= nCh
        error('snpm_lmm_cluster:evaluableSize', ...
            'evaluable mask has %d entries but power has %d channels.', ...
            numel(evaluable), nCh);
    end
    if ~any(evaluable)
        error('snpm_lmm_cluster:noEvaluableChannels', ...
            ['No channel is complete across all %d analysed rows, so no channel can be ' ...
             'tested against a permutation null. Drop the incomplete subjects/trials ' ...
             'instead of the channels.'], size(power, 1));
    end

    % ---- real clusters ----
    [~, real_wald, real_p] = snpm_lmm_fit(power, meta, spec, evaluable);
    [real_p, real_wald] = mask_map(real_p, real_wald, evaluable);
    real_clusters = local_find_pclusters(real_p, real_wald, alpha, sparse_adj);

    % ---- null distribution of the max cluster mass ----
    max_cluster_mass = zeros(possible_permutations, 1);
    for permIndex = 1:possible_permutations
        if mod(permIndex, 100) == 0
            disp([num2str(permIndex), ' out of ', num2str(possible_permutations), ' LMM cluster permutations completed...']);
        end
        meta_perm = snpm_lmm_permute_meta(meta, spec, spec.perm);
        [~, pwald, pp] = snpm_lmm_fit(power, meta_perm, spec, evaluable);
        [pp, pwald] = mask_map(pp, pwald, evaluable);
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
        % Minimum-bias FWE p (Nichols & Holmes 2001; Phipson & Smyth 2010):
        %   p = (#{max-null cluster mass >= observed} + 1)/(N + 1), never zero.
        Clusters(cli, 1).p = (sum(max_cluster_mass >= real_clusters(cli).mass) + 1) / (possible_permutations + 1);
        Clusters(cli, 1).threshold    = alpha;
        Clusters(cli, 1).permutations = possible_permutations;
    end
end

function [p_model, wald] = mask_map(p_model, wald, evaluable)
% One channel set for the observed and the permuted maps: a channel outside
% the evaluable mask has no p (so it can never join a cluster) and no mass.
    p_model(~evaluable) = NaN;
    wald(~evaluable)    = 0;
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
