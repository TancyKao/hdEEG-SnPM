function [T, p] = snpm_lmm_TFCE(power, meta, spec, neighbors, E, H, alpha, tail, permutations, evaluable)
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
%   evaluable    : (optional) 1 x nCh logical, TRUE for channels that can be
%                  evaluated in EVERY permutation. Defaults to
%                  all(isfinite(power), 1) — see below. Channels that are
%                  FALSE are NaN in the observed map as well as in every
%                  permuted map, so both are defined on one channel set.
%
% OUTPUTS
%   T : struct with real_T, real_TFCE, tMaxTFCE, critical_T_TFCE,
%       excluded_channels
%   p : struct with real (uncorrected model p) and correctedTFCE (FWE); NaN
%       (not 1) on channels outside the evaluable set — they were not tested.
%
% WHY THE EVALUABLE MASK EXISTS (null-construction, not cosmetics)
%   snpm_lmm_fit deletes rows per channel (valid = ~isnan(col)), so a channel
%   present for only a subset of the rows is fitted on that subset. Two things
%   then break exchangeability:
%     1. snpm_lmm_permute_meta permutes labels across ALL rows/subjects, so a
%        relabelling can make the group factor constant among that channel's
%        valid rows. fitlme does not throw and anova(lme) still returns the
%        term — with FStat = NaN. Measured on a channel present for 3 of 12
%        subjects: the observed F was finite while 58 of 300 permuted fits
%        came back NaN. The channel therefore has an observed statistic but is
%        missing from a fifth of the null it is scored against.
%     2. The within_subject scheme shuffles the DV across all of a subject's
%        rows, but the fit uses only the channel's valid rows. When
%        missingness is associated with the DV the permuted fits see a
%        different DV distribution than the observed fit (measured gap: 5.5
%        permutation SE), so observed and permuted statistics are not computed
%        on comparable data.
%   Both vanish if a channel is analysed only when it is COMPLETE across all
%   analysed rows: then `valid` is all-true, no relabelling can make the
%   factor constant on a subset, and the shuffle acts on exactly the rows that
%   enter the fit. The mask is a property of the data, not of the labels, so
%   it is computed once, here, before the permutation loop.

    if nargin < 5 || isempty(E), E = 0.5; end
    if nargin < 6 || isempty(H), H = 2;   end
    if nargin < 7 || isempty(alpha), alpha = 0.05; end
    if nargin < 8 || isempty(tail), tail = 'both'; end

    nCh = size(power, 2);

    if nargin < 10 || isempty(evaluable)
        evaluable = all(isfinite(power), 1);
    end
    evaluable = reshape(logical(evaluable), 1, []);
    if numel(evaluable) ~= nCh
        error('snpm_lmm_TFCE:evaluableSize', ...
            'evaluable mask has %d entries but power has %d channels.', ...
            numel(evaluable), nCh);
    end
    if ~any(evaluable)
        error('snpm_lmm_TFCE:noEvaluableChannels', ...
            ['No channel is complete across all %d analysed rows, so no channel can be ' ...
             'tested against a permutation null. Drop the incomplete subjects/trials ' ...
             'instead of the channels.'], size(power, 1));
    end
    % Match the production convention: adjacency is sized to the neighbour
    % matrix, and ClusterEnhancement requires the statistic vector length to
    % equal that size, i.e. size(power,2) == size(neighbors,1).
    sparse_adj = make_neighbors_sparse(neighbors, size(neighbors, 1));

    % ---- real fit ----
    % Restrict the observed map to the evaluable set BEFORE enhancement, so the
    % observed TFCE integral sees exactly the channels the permuted ones do.
    [real_stat, ~, real_p] = snpm_lmm_fit(power, meta, spec, evaluable);
    real_stat(~evaluable) = NaN;
    real_p(~evaluable)    = NaN;
    T.real_T    = real_stat;
    T.real_TFCE = enhance_map(real_stat, sparse_adj, E, H, nCh);
    p.real      = real_p;

    % ---- permutation null of max TFCE ----
    possible_permutations = permutations;
    maxTFCE = zeros(possible_permutations, 1);
    n_degenerate = 0;
    for permIndex = 1:possible_permutations
        if mod(permIndex, 100) == 0
            disp([num2str(permIndex), ' out of ', num2str(possible_permutations), ' LMM permutations completed...']);
        end
        meta_perm = snpm_lmm_permute_meta(meta, spec, spec.perm);
        perm_stat = snpm_lmm_fit(power, meta_perm, spec, evaluable);
        perm_stat(~evaluable) = NaN;   % identical channel set to the observed map
        tfce = enhance_map(perm_stat, sparse_adj, E, H, nCh);
        switch tail
            case 'both'
                cand = abs(tfce(evaluable));
            otherwise
                cand = tfce(evaluable);
        end
        % The max is taken over the MASKED channel set, never over the whole
        % map. If a permutation yields no statistic at all (every evaluable
        % channel non-finite) the entry is set to +Inf rather than left NaN:
        % NaN >= obs is FALSE, so a NaN entry would silently drop out of
        % sum(tMaxTFCE >= obs) and make the p-value anti-conservative.
        mx = max(cand, [], 'omitnan');
        if isempty(mx) || ~isfinite(mx)
            mx = Inf;
            n_degenerate = n_degenerate + 1;
        end
        maxTFCE(permIndex) = mx;
    end
    if n_degenerate > 0
        warning('snpm_lmm_TFCE:degeneratePermutations', ...
            ['%d of %d permutations produced no finite statistic on any evaluable ' ...
             'channel; those null entries were set to +Inf (conservative) rather ' ...
             'than dropped. Check model convergence.'], n_degenerate, possible_permutations);
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
    p.correctedTFCE(~evaluable) = NaN;   % explicit: not tested, not "p = 1"

    % Carry the exclusion with the result so a caller cannot report a channel
    % set without also being able to report what left it.
    T.excluded_channels = find(~evaluable);
end

function m = enhance_map(stat, sparse_adj, E, H, nCh)
% ClusterEnhancement with a shape guard.
%
% When the WHOLE statistic map lies below the TFCE integration step (dh=0.1),
% ClusterEnhancement's threshold range collapses to a single level, its
% per-channel sum degenerates to a scalar 0, and the trailing
% TFCEdata(isnan(...)) = NaN assignment then grows that scalar only as far as
% the LAST NaN channel — so the returned vector can be SHORTER than nCh. The
% integral in that situation is genuinely zero at every channel, so the map is
% right-padded with zeros and the NaN channels are re-marked. Without this the
% map cannot be indexed by the evaluable mask (and, pre-existing, the
% corrected-p loop would silently leave the tail channels at p = 1).
    m = reshape(ClusterEnhancement(stat, sparse_adj, E, H), 1, []);
    if numel(m) < nCh
        m(numel(m)+1:nCh) = 0;
    elseif numel(m) > nCh
        error('snpm_lmm_TFCE:enhancedSize', ...
            'ClusterEnhancement returned %d values for %d channels.', numel(m), nCh);
    end
    m(isnan(stat)) = NaN;
end
