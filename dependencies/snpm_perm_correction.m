function [T, p, Clusters] = snpm_perm_correction(real_stat, real_p, perm_stat_fn, ...
        neighbors, E, H, alpha, permutations, contrast_type)
% Shared permutation correction driver for the GLM tier. Given the real
% per-channel statistic map and a function that returns one permuted
% statistic map, it builds BOTH corrections in a single permutation loop:
%   - TFCE (ClusterEnhancement) with a max-statistic null -> p.correctedTFCE
%   - cluster mass (channels with model p<alpha, mass = sum(wald)/n) with a
%     max-cluster-mass null -> Clusters(i).p
% Emits the same T/p/Clusters field shapes the output chain expects
% (func_genSnpmTable, generateAnalysisReport, TopoplotSignificant_single),
% so it is interchangeable with snpm_single_threshold_with_TFCE + cluster.
%
% INPUTS
%   real_stat    : 1 x nCh observed statistic (signed t, or F)
%   real_p       : 1 x nCh observed parametric p (for cluster forming)
%   perm_stat_fn : handle, [stat, p] = perm_stat_fn() returns one permutation
%   neighbors    : channels x maxNeighbors adjacency (NaN padded)
%   E, H         : TFCE exponents
%   alpha        : significance / cluster-forming level
%   permutations : number of permutations
%   contrast_type: 't' (signed -> two-sided abs) | 'F' (non-negative)
%
% OUTPUTS: T (real_T, real_TFCE, tMaxTFCE, critical_T_TFCE),
%          p (real, correctedTFCE), Clusters (channels, mass, p, ...)

    sparse_adj = make_neighbors_sparse(neighbors, size(neighbors, 1));
    nCh = numel(real_stat);
    is_F = strcmpi(contrast_type, 'F');
    waldfun = @(s) wald_of(s, is_F);

    % ---- real maps ----
    T.real_T    = real_stat;
    T.real_TFCE = ClusterEnhancement(real_stat, sparse_adj, E, H);
    p.real      = real_p;
    real_clusters = find_pclusters(real_p, waldfun(real_stat), alpha, sparse_adj);

    % ---- permutation null (TFCE max + cluster-mass max in one loop) ----
    maxTFCE = zeros(permutations, 1);
    maxMass = zeros(permutations, 1);
    for i = 1:permutations
        if mod(i, 100) == 0
            disp([num2str(i), ' out of ', num2str(permutations), ' GLM permutations completed...']);
        end
        [s, pp] = perm_stat_fn();
        tfce = ClusterEnhancement(s, sparse_adj, E, H);
        if is_F
            maxTFCE(i) = max(tfce, [], 'omitnan');
        else
            maxTFCE(i) = max(abs(tfce), [], 'omitnan');
        end
        pc = find_pclusters(pp, waldfun(s), alpha, sparse_adj);
        if ~isempty(pc), maxMass(i) = max([pc.mass]); end
    end
    T.tMaxTFCE = sort(maxTFCE, 'descend');
    maxMass    = sort(maxMass, 'descend');

    % ---- FWE-corrected TFCE p per channel ----
    critical_T_indx   = floor(alpha * permutations) + 1;
    T.critical_T_TFCE = T.tMaxTFCE(critical_T_indx);
    p.correctedTFCE = ones(1, nCh);
    p.correctedTFCE(isnan(T.real_TFCE)) = NaN;
    for i = 1:nCh
        if is_F
            if T.real_TFCE(i) > T.tMaxTFCE(end)
                p.correctedTFCE(i) = find(T.tMaxTFCE < T.real_TFCE(i), 1, 'first') / permutations;
            end
        else
            if abs(T.real_TFCE(i)) > T.tMaxTFCE(end)
                p.correctedTFCE(i) = find(T.tMaxTFCE < abs(T.real_TFCE(i)), 1, 'first') / permutations;
            end
        end
    end

    % ---- cluster p-values from the max-mass null ----
    if isempty(real_clusters)
        Clusters = struct('channels', {[]}, 'mass', {0}, 'p', {1}, ...
            'threshold', {alpha}, 'permutations', {permutations});
    else
        Clusters = struct([]);
        for cli = 1:numel(real_clusters)
            Clusters(cli, 1).channels = real_clusters(cli).channels;
            Clusters(cli, 1).mass     = real_clusters(cli).mass;
            rank = find(maxMass < real_clusters(cli).mass, 1, 'first');
            if isempty(rank), Clusters(cli, 1).p = 1; else, Clusters(cli, 1).p = rank / permutations; end
            Clusters(cli, 1).threshold    = alpha;
            Clusters(cli, 1).permutations = permutations;
        end
    end
end

function w = wald_of(s, is_F)
    if is_F, w = s; else, w = s.^2; end   % F is already a ratio; t -> Wald = t^2
end

function clusters = find_pclusters(p_model, wald, alpha, sparse_adj)
% Connected components of channels with model p < alpha; mass = sum(Wald)/n.
    mask = double(p_model < alpha);
    mask(isnan(p_model)) = 0;
    cc = snpm_find_clusters_graphalgs(mask, 0.5, sparse_adj);
    clusters = struct([]);
    k = 0;
    for c = 1:numel(cc)
        ch = cc{c};
        if isempty(ch), continue; end
        k = k + 1;
        clusters(k, 1).channels = ch;
        clusters(k, 1).mass     = sum(wald(ch), 'omitnan') / numel(ch);
    end
end
