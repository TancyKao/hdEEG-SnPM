function excl = snpm_excluded_channel_info(evaluable, chanlocs, context_label, reason, why)
%SNPM_EXCLUDED_CHANNEL_INFO  Report the channels that left an analysis.
%
%   excl = SNPM_EXCLUDED_CHANNEL_INFO(evaluable, chanlocs, context_label)
%   excl = SNPM_EXCLUDED_CHANNEL_INFO(evaluable, chanlocs, context_label, reason, why)
%
%   evaluable      1 x nCh logical, TRUE for channels that can be evaluated with
%                  the SAME row set in the observed map and in every permutation.
%   chanlocs       EEGLAB chanlocs struct array (for labels); may be shorter.
%   context_label  free text put in the report ('correlationP', preset name...).
%   reason         OPTIONAL short reason string stored in excl.reason.
%   why            OPTIONAL sentence explaining the consequence, spliced into the
%                  printed line and the warning.
%
%   TWO DIFFERENT RULES USE THIS REPORTER AND THEY MUST NOT SHARE A SENTENCE.
%   The correlation / GLM tiers drop a channel that is INCOMPLETE (the
%   complete-column rule -- a validity fix, because permutation re-pairs the
%   rows and the analysed set would otherwise change with the permutation). The
%   t-test tier drops only the channels that are DEGENERATE under some labelling
%   -- a power fix, since its analysed row set is already permutation-invariant
%   and it is deliberately NOT masked for ordinary missingness. Telling a t-test
%   user their channel was "not complete across all analysed subjects" would
%   describe a rule that tier does not apply, so the default wording below is
%   only used when the caller supplies none.
%
%   Returns the struct the output chain already knows how to render:
%     .n .n_channels .index .labels .context .reason
%   consumed by dependencies/func_genSnpmTable.m (sheet 'excludedChannels') and
%   dependencies/generateAnalysisReport.m (amber banner). Shape is identical to
%   the local helpers in core_snpm_glm.m / core_snpm_circ.m -- this is the
%   shared version; those two can be pointed at it in a later cleanup.
%
%   A channel must never leave an analysis silently, so this both prints and
%   raises a warning when anything is excluded.

    DEFAULT_REASON = 'missing data: channel not complete across all analysed subjects';
    DEFAULT_WHY = ['They are not complete across all analysed subjects, so the number of ' ...
        'subjects entering them would depend on the permutation. They are NaN in the ' ...
        'observed map as well as in the null (a value scored against a null built at a ' ...
        'different sample size is invalid).'];
    if nargin < 4 || isempty(reason), reason = DEFAULT_REASON; end
    if nargin < 5 || isempty(why),    why    = DEFAULT_WHY;    end

    evaluable = reshape(logical(evaluable), 1, []);
    idx = find(~evaluable);
    labels = cell(1, numel(idx));
    for k = 1:numel(idx)
        lab = '';
        if idx(k) <= numel(chanlocs) && isfield(chanlocs, 'labels')
            lab = chanlocs(idx(k)).labels;
        end
        if isempty(lab), lab = sprintf('#%d', idx(k)); end
        labels{k} = lab;
    end

    excl = struct('n', numel(idx), 'n_channels', numel(evaluable), ...
        'index', idx, 'labels', {labels}, 'context', context_label, ...
        'reason', reason);

    if isempty(idx)
        return;
    end

    msg = sprintf('%s: EXCLUDED %d of %d channels. %s Channels: %s', ...
        context_label, numel(idx), numel(evaluable), why, strjoin(labels, ', '));
    fprintf('%s\n', msg);
    warning('snpm:channelsExcluded', '%s', msg);
end
