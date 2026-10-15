function excl = snpm_excluded_add(excl, drop, reason, chanlocs)
%SNPM_EXCLUDED_ADD  Add channels with their own reason to an excluded-channel report.
%
%   excl = SNPM_EXCLUDED_ADD(excl, drop, reason, chanlocs)
%
%   excl     : struct from snpm_excluded_channel_info (or the core_snpm_glm
%              local of the same shape): .n .n_channels .index .labels
%              .context .reason, optionally .reasons
%   drop     : 1 x nCh logical, TRUE for the channels to add. Channels already
%              listed keep their first reason.
%   reason   : the reason for THESE channels (one string)
%   chanlocs : EEGLAB chanlocs (for labels)
%
% A loop can now exclude channels for more than one reason (missing data, no
% variation left after the nuisance model, an observed fit that is perfect), so
% the struct carries a per-channel .reasons cellstr parallel to .index. The
% scalar .reason stays for every existing reader: it is the single reason when
% there is one, else the distinct reasons joined with '; '. Entries are kept in
% channel-index order.
%
% Prints and warns, like snpm_excluded_channel_info: a channel must never leave
% an analysis silently.

    drop = reshape(logical(drop), 1, []);
    if ~isfield(excl, 'reasons') || numel(excl.reasons) ~= numel(excl.index)
        excl.reasons = repmat({char(excl.reason)}, 1, numel(excl.index));
    end
    new = setdiff(find(drop), excl.index);
    if isempty(new), return; end

    labs = cell(1, numel(new));
    for k = 1:numel(new)
        lab = '';
        if new(k) <= numel(chanlocs) && isfield(chanlocs, 'labels')
            lab = chanlocs(new(k)).labels;
        end
        if isempty(lab), lab = sprintf('#%d', new(k)); end
        labs{k} = lab;
    end

    idx  = [reshape(excl.index, 1, []), new];
    lab  = [reshape(excl.labels, 1, []), labs];
    rsn  = [reshape(excl.reasons, 1, []), repmat({char(reason)}, 1, numel(new))];
    [idx, o] = sort(idx);
    excl.index   = idx;
    excl.labels  = lab(o);
    excl.reasons = rsn(o);
    excl.n       = numel(idx);
    u = unique(excl.reasons, 'stable');
    excl.reason  = strjoin(u, '; ');

    msg = sprintf('%s: EXCLUDED %d of %d channels: %s. Channels: %s', ...
        char(excl.context), numel(new), excl.n_channels, char(reason), strjoin(labs, ', '));
    fprintf('%s\n', msg);
    warning('snpm:channelsExcluded', '%s', msg);
end
