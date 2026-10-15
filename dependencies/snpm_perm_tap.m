function out = snpm_perm_tap(cmd, tag, value)
%SNPM_PERM_TAP  Inert recorder for the permutation indices an analysis drew.
%
%   snpm_perm_tap('on')                  clear the log and start recording
%   snpm_perm_tap('off')                 stop recording (the log is kept)
%   snpm_perm_tap('clear')               empty the log, leave on/off unchanged
%   tf  = snpm_perm_tap('ison')          true while recording
%   snpm_perm_tap('record', tag, value)  append one entry; returns at once when
%                                        recording is off
%   log = snpm_perm_tap('get')           K x 1 struct array, in call order:
%           .tag     char, what was drawn (see TAGS below)
%           .value   the array that was drawn, copied as-is ([] for a marker)
%           .site    'function:line' of the function that called 'record'
%           .caller  'function:line' of THAT function's caller ('' at top level)
%
% WHY IT EXISTS
%   The Python port of this toolbox is trusted only when it reproduces the
%   MATLAB result on shared inputs using the SAME permutation indices. No
%   engine returns the indices it drew, so this function is called at every
%   draw site and copies them out. scripts/export_parity_fixtures.m is the
%   consumer.
%
% WHAT IT MUST NEVER DO
%   Off (the default): no random number is consumed, no output changes and the
%   call costs one persistent-variable test. On: still no random number is
%   consumed and no output changes -- it only copies an array the caller has
%   already drawn (dbstack, used for .site/.caller, does not touch the
%   generator). Nothing in the toolbox reads the log back, so a statistic
%   cannot come to depend on it.
%
% TAGS
%   'signflip:exact' | 'signflip:random'   snpm_signflip_patterns
%       value = Nperm x nSubj logical, true = pair kept, false = sign-flipped
%   'relabel:exact'  | 'relabel:random'    snpm_relabel_assignments
%       value = Nperm x nSubj, row = subject indices, first nGrp are group 1
%   'rowperm:exact'  | 'rowperm:random'    snpm_row_permutations
%       value = Nperm x nSubj, row = the order applied to data_x's rows
%   'glm_permute:free' | 'glm_permute:within'   snpm_glm_permute, ONE draw
%       value = the row index pe, so that Yperm = Zfit + R(pe,:)
%   'lmm:within_subject' | 'lmm:group_label'    snpm_lmm_permute_meta, ONE draw
%       value = row index into the trial table (within_subject) or the
%               subject-level index perm (group_label)
%   'circ_u2:labels'                       core_snpm_circ Watson U2, ONE draw
%       value = the index applied to the group-label vector
%   'legacy:engine_inputs'                 core_snpm_analysis, once per run
%       value = struct: data_x, data_y exactly as handed to global_stat_test and
%               both engines (after matching, channel selection, transform,
%               masking, residualization), comparison (the key the engines get),
%               tail, ncov, alpha, threshold (cluster-forming)
%   'legacy:cluster_null'                  snpm_cluster_analysis, once per run
%       value = struct: null_sorted (the sorted max-cluster-size null the
%               cluster p is computed against), thresholds, threshold
%   'glm:cluster_mass_null'                snpm_perm_correction, once per call
%       value = the sorted max cluster-mass null the cluster p is computed
%               against (GLM and circular tiers share this driver)
%   'loop:...' | 'scope:...'               markers, value = []. Written by the
%       orchestrators immediately before a permutation loop so the per-draw
%       entries that follow can be attributed to it.
%   All indices are 1-based.
%
% The legacy generators cannot know which loop asked for them; .caller says
% so (snpm_single_threshold_with_TFCE vs snpm_cluster_analysis).

    persistent recording entries n
    if isempty(recording)
        recording = false;
        entries = cell(0, 4);
        n = 0;
    end

    % Hot path first: this is what every draw site executes.
    if strcmp(cmd, 'record')
        if recording
            st = dbstack(1);
            site = ''; caller = '';
            if ~isempty(st)
                site = sprintf('%s:%d', st(1).name, st(1).line);
            end
            if numel(st) >= 2
                caller = sprintf('%s:%d', st(2).name, st(2).line);
            end
            n = n + 1;
            if n > size(entries, 1)
                entries{max(256, 2 * n), 4} = [];   % amortised growth
            end
            entries(n, :) = {tag, value, site, caller};
        end
        return
    end

    switch cmd
        case 'on'
            entries = cell(0, 4); n = 0;
            recording = true;
        case 'off'
            recording = false;
        case 'clear'
            entries = cell(0, 4); n = 0;
        case 'ison'
            out = recording;
        case 'get'
            out = cell2struct(entries(1:n, :), {'tag', 'value', 'site', 'caller'}, 2);
        otherwise
            error('snpm_perm_tap:badCommand', ...
                'Unknown command ''%s'' (use on | off | clear | ison | record | get).', cmd);
    end
end
