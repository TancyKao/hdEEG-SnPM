function [T, p] = snpm_lmm_TFCE(power, meta, spec, neighbors, E, H, alpha, tail, permutations)
% Per-channel linear mixed model with TFCE multiple-comparison correction.
% Mirrors the output contract of snpm_single_threshold_with_TFCE so the
% downstream output chain (func_genSnpmTable, TopoplotSignificant_single,
% generateAnalysisReport) works unchanged.
%
% The per-channel statistic is the LMM effect of interest (signed t for a
% continuous predictor, F for a categorical factor) from snpm_lmm_fit, and
% the null distribution is built by refitting the model on permuted data
% (snpm_lmm_permute_meta), enhancing each statistic map with the existing
% ClusterEnhancement TFCE, and keeping the max across channels (FWE).
%
% INPUTS
%   power        : trials x nCh channel matrix
%   meta         : trial-level table of model variables
%   spec         : LMM spec (see snpm_lmm_fit); must also carry the
%                  permutation scheme in spec.perm ('within_subject' |
%                  'group_label')
%   neighbors    : channels x maxNeighbors adjacency (NaN padded)
%   E, H         : TFCE exponents (defaults 0.5, 2)
%   alpha        : significance level (default 0.05)
%   tail         : 'both' (default), 'left', 'right'
%   permutations : number of permutations
%
% OUTPUTS
%   T : struct with real_T, real_TFCE, tMaxTFCE, critical_T_TFCE
%   p : struct with real (uncorrected model p) and correctedTFCE (FWE)

    if nargin < 5 || isempty(E), E = 0.5; end
    if nargin < 6 || isempty(H), H = 2;   end
    if nargin < 7 || isempty(alpha), alpha = 0.05; end
    if nargin < 8 || isempty(tail), tail = 'both'; end

    nCh = size(power, 2);
    % Match the production convention: adjacency is sized to the neighbour
    % matrix, and ClusterEnhancement requires the statistic vector length to
    % equal that size, i.e. size(power,2) == size(neighbors,1).
    sparse_adj = make_neighbors_sparse(neighbors, size(neighbors, 1));

    % ---- real fit ----
    [real_stat, ~, real_p] = snpm_lmm_fit(power, meta, spec);
    T.real_T    = real_stat;
    T.real_TFCE = ClusterEnhancement(real_stat, sparse_adj, E, H);
    p.real      = real_p;

    % ---- permutation null of max TFCE ----
    possible_permutations = permutations;
    maxTFCE = zeros(possible_permutations, 1);
    for permIndex = 1:possible_permutations
        if mod(permIndex, 100) == 0
            disp([num2str(permIndex), ' out of ', num2str(possible_permutations), ' LMM permutations completed...']);
        end
        meta_perm = snpm_lmm_permute_meta(meta, spec, spec.perm);
        perm_stat = snpm_lmm_fit(power, meta_perm, spec);
        tfce = ClusterEnhancement(perm_stat, sparse_adj, E, H);
        switch tail
            case 'both'
                maxTFCE(permIndex) = max(abs(tfce), [], 'omitnan');
            otherwise
                maxTFCE(permIndex) = max(tfce, [], 'omitnan');
        end
    end
    T.tMaxTFCE = sort(maxTFCE, 'descend');

    % ---- FWE-corrected p per channel ----
    % Minimum-bias FWE p (Nichols & Holmes 2001; Phipson & Smyth 2010):
    %   p = (#{max-null >= obs} + 1)/(N + 1), never zero.
    critical_T_indx     = floor(alpha * possible_permutations) + 1;
    T.critical_T_TFCE   = T.tMaxTFCE(critical_T_indx);

    p.correctedTFCE = ones(1, nCh);
    p.correctedTFCE(isnan(T.real_TFCE)) = NaN;
    for i = 1:length(T.real_TFCE)
        if isnan(T.real_TFCE(i)), continue; end
        switch tail
            case 'both'
                obs = abs(T.real_TFCE(i));
            otherwise
                obs = T.real_TFCE(i);
        end
        p.correctedTFCE(i) = (sum(T.tMaxTFCE >= obs) + 1) / (possible_permutations + 1);
    end
end
