function [meta, power, dropped] = snpm_lmm_drop_missing_meta(meta, power, spec)
% Drop trials with a missing MODEL variable before anything else, and check
% that the design left is one the permutation scheme can be valid for.
%
%   [meta, power, dropped] = snpm_lmm_drop_missing_meta(meta, power, spec)
%
% Called by core_snpm_lmm immediately after the metadata and channel matrix
% are split from the input table, before the datatype transform, the channel
% setup and the evaluable-channel mask. meta and power lose the same rows.
%
% THE VARIABLE SET (M). Only what enters the model or the permutation:
%   spec.dv (when it is a meta column; a DV of POWER is the channel, whose
%   missingness is the evaluable-channel rule's job), every meta column named
%   in spec.fixed and spec.random, spec.subject, and spec.group when
%   spec.perm = 'group_label'. Other metadata columns (e.g. a time column the
%   model does not use) never drop a row.
% "Missing" = ismissing: NaN, '' / <missing> text, <undefined> categorical.
%
% WHY. fitlme silently deletes rows with a missing model variable, so before
% this rule every channel was fitted on fewer rows than n_trials reported, and
%   - under group_label, snpm_lmm_permute_meta takes each subject's group from
%     its FIRST trial (subj_groups = grp(ia)); a group label missing on only
%     some trials (or varying within a subject) meant even the identity
%     permutation did not reproduce the observed table -- the null was built
%     on a different design from the observed statistic (invalid);
%   - under within_subject, a missing DV was shuffled onto other trials, so
%     the fitted row set changed from permutation to permutation.
% Dropping the rows once fixes one row set for the observed fit and every
% permutation. It assumes missingness is unrelated to the effect tested
% (e.g. awakenings with no report), as any complete-case analysis does.
%
% CHECKS after the drop (run whether or not anything was dropped):
%   fewer than 3 subjects                     -> error core_snpm_lmm:tooFewSubjects
%   group_label: group varies within subject  -> error core_snpm_lmm:groupVariesWithinSubject
%   group_label: < 2 levels, or a level with
%                fewer than 2 subjects         -> error core_snpm_lmm:groupLevelTooFewSubjects
%   within_subject: a subject with 1 trial    -> warning snpm:lmmSingleTrialSubject
%                   (its DV cannot move under the shuffle)
%
% OUTPUT dropped (plain types, saved in results_struct.lmm.dropped_rows):
%   .n               total rows dropped
%   .n_input         rows before the drop
%   .columns         1 x k cellstr, the variable set M
%   .per_column      1 x k double, missing cells per column of M (a row can
%                    count under several columns)
%   .rows            1 x n double, dropped row indices into the input table
%   .subjects_lost   1 x s cellstr, subjects that lost every trial

    if isfield(spec, 'subject') && ~isempty(spec.subject), subj_col = spec.subject; else, subj_col = 'Subject'; end
    if isfield(spec, 'group')   && ~isempty(spec.group),   grp_col  = spec.group;   else, grp_col  = 'group';   end
    vn = meta.Properties.VariableNames;
    if ~ismember(subj_col, vn)
        error('core_snpm_lmm:noSubjectColumn', ...
            'The subject column ''%s'' is not among the metadata columns (lmm_meta_cols).', subj_col);
    end
    is_grp_perm = strcmpi(spec.perm, 'group_label');
    if is_grp_perm && ~ismember(grp_col, vn)
        error('core_snpm_lmm:noGroupColumn', ...
            'group_label permutation needs the group column ''%s'' among the metadata columns.', grp_col);
    end

    % ---- the variable set M ----
    toks = regexp([char(spec.fixed) ' ' char(spec.random)], '[A-Za-z]\w*', 'match');
    M = [{char(spec.dv)}, toks, {subj_col}];
    if is_grp_perm, M{end+1} = grp_col; end
    M = M(~strcmp(M, 'POWER'));
    M = intersect(vn, M, 'stable');          % meta columns only, table order

    miss = ismissing(meta(:, M));            % n x k
    drop = any(miss, 2);
    n_in = height(meta);
    subj_before = unique_ids(meta.(subj_col)(~ismissing(meta.(subj_col))));

    meta  = meta(~drop, :);
    power = power(~drop, :);

    subj_after = unique_ids(meta.(subj_col));
    dropped = struct();
    dropped.n             = sum(drop);
    dropped.n_input       = n_in;
    dropped.columns       = M;
    dropped.per_column    = double(sum(miss, 1));
    dropped.rows          = find(drop)';
    dropped.subjects_lost = setdiff(subj_before, subj_after, 'stable');

    if dropped.n > 0
        parts = cell(1, 0);
        for k = 1:numel(M)
            if dropped.per_column(k) > 0
                parts{end+1} = sprintf('%s %d', M{k}, dropped.per_column(k)); %#ok<AGROW>
            end
        end
        lost = 'none';
        if ~isempty(dropped.subjects_lost), lost = strjoin(dropped.subjects_lost, ', '); end
        warning('snpm:lmmMetaNaNRows', ...
            ['%d of %d rows (trials) dropped for a missing model variable (%s); ' ...
             'subjects that lost every trial: %s. The dropped rows are left out of the ' ...
             'observed fit and every permutation alike; this assumes the missingness is ' ...
             'unrelated to the effect tested.'], dropped.n, n_in, strjoin(parts, ', '), lost);
    end

    % ---- design checks on what is left ----
    [~, ~, sidx] = unique(id_strings(meta.(subj_col)), 'stable');
    nsub = max([sidx; 0]);
    if nsub < 3
        error('core_snpm_lmm:tooFewSubjects', ...
            'Only %d subject(s) left after dropping rows with missing model variables; at least 3 are needed.', nsub);
    end

    if is_grp_perm
        [glev, ~, gidx] = unique(id_strings(meta.(grp_col)));
        per_subj = accumarray(sidx, gidx, [], @(v) numel(unique(v)));
        if any(per_subj > 1)
            us = unique(id_strings(meta.(subj_col)), 'stable');
            bad = us(per_subj > 1);
            error('core_snpm_lmm:groupVariesWithinSubject', ...
                ['group_label permutation relabels whole subjects, so the group must be constant ' ...
                 'within each subject; it varies within %d subject(s): %s.'], numel(bad), ...
                strjoin(bad(1:min(10, end))', ', '));
        end
        sub_grp = accumarray(sidx, gidx, [], @(v) v(1));
        n_per_level = accumarray(sub_grp, 1, [numel(glev), 1]);
        if numel(glev) < 2 || any(n_per_level < 2)
            error('core_snpm_lmm:groupLevelTooFewSubjects', ...
                ['group_label permutation needs at least 2 groups with at least 2 subjects each; ' ...
                 'subjects per group after the drop: %s.'], ...
                strjoin(cellfun(@(l, n) sprintf('%s %d', l, n), glev(:)', num2cell(n_per_level(:)'), ...
                'UniformOutput', false), ', '));
        end
    elseif strcmpi(spec.perm, 'within_subject')
        ntr = accumarray(sidx, 1);
        if any(ntr == 1)
            us = unique(id_strings(meta.(subj_col)), 'stable');
            warning('snpm:lmmSingleTrialSubject', ...
                ['%d subject(s) have a single trial (%s): the within-subject shuffle cannot move ' ...
                 'their DV, so they add no permutation variability.'], sum(ntr == 1), ...
                strjoin(us(ntr == 1)', ', '));
        end
    end
end

function s = id_strings(x)
% Any ID column (numeric, cellstr, string, categorical) as a cellstr column.
    s = cellstr(string(x));
    s = s(:);
end

function u = unique_ids(x)
    u = unique(id_strings(x), 'stable')';
end
