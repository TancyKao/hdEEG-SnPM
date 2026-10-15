function pr = snpm_glm_pair(preset, meta, opts, factor_col, a, b, label)
%SNPM_GLM_PAIR  One pairwise post-hoc comparison as its own two-level test.
%
%   pr = snpm_glm_pair(preset, meta, opts, factor_col, a, b, label)
%
% Method (Method Spec A, approved 2026-10-14). The comparison "level a vs
% level b" of a factor with three or more levels is tested on the ROWS OF a
% AND b ONLY, with the preset's design rebuilt on those rows. Under the
% hypothesis that a and b do not differ, only the a and b rows are
% exchangeable; the other levels take no part, so they are neither fitted nor
% permuted.
%
% This replaces the earlier scheme, which tested the pair as a contrast in the
% omnibus model under the omnibus nuisance fit. That left the OTHER levels'
% effects in the residuals being permuted: with a = b and a third level raised
% head-wide the cluster test rejected 0.17-0.59 of the time (nominal 0.05),
% and 0.19-0.33 when the third level was merely less variable.
%
% INPUTS
%   preset     : 'anova1' | 'ancova' | 'rmanova' | 'mixed2way'
%                ('mixed2way' builds the CONDITION effect model, Spec C)
%   meta       : table the omnibus design was built from (one row per row of Y)
%   opts       : the omnibus column mapping (group_col, condition_col,
%                subject_col, covariate_cols)
%   factor_col : the factor whose levels are compared (D.factor_col)
%   a, b       : level indices in order of first appearance in meta.(factor_col)
%                (D.posthoc(k).a / .b); a comes before b
%   label      : 'A vs B', used in warnings
%
% OUTPUT pr
%   .ok        : false when the pair cannot be tested
%   .reason    : '' or, when ~ok, why the pair was not run (plain words; this
%                string is what the results text, the workbook and the report
%                print, so the reader of an output file sees it, not only the
%                console)
%   .rows      : indices into meta (and Y) of the rows used, in table order
%   .D         : design on those rows (snpm_glm_design; .X .nuisance_idx .eb .perm_type)
%   .C         : 1 x p contrast, signed so the statistic is mean(a) - mean(b)
%   .n_rows    : number of rows used
%   .n_subjects: within designs, subjects that have both levels (the n of the
%                test); between designs, n_a + n_b
%   .n_subjects_total : within designs, subjects in the effect's table
%   .group_sizes      : between designs, [n_a n_b]; [] for within designs
%   .dropped   : cellstr, subjects left out because they have only one of the
%                two levels (within designs)
%   .run_on    : what the pair was run on, in words, e.g. '23 of 25 subjects'
%                or '8 and 9 subjects'
%   .n_small   : complete subjects (within) or the smaller group's size
%
% RULES
%   Within-subject designs (rmanova, mixed2way): a subject with only one of
%   the two levels contributes nothing to a within-subject contrast (its row
%   is fitted exactly by its own intercept), so it is dropped, with a warning
%   naming it. Fewer than 2 complete subjects: skipped. Fewer than 6: warned,
%   because the smallest two-sided permutation p is 2/2^n.
%   Between-subject designs (anova1, ancova): warned under 4 in the smaller
%   group. Covariates are re-centred on the subset, categorical covariates are
%   coded from the levels present in the subset, and a covariate constant on
%   the subset is dropped (all three inside snpm_glm_design).
%   The design must have full column rank and at least 1 error degree of
%   freedom; otherwise the pair is skipped with a warning.
%   Permutation: Freedman-Lane with every column except the two-level dummy
%   as nuisance; 'free' over the subset rows (between) or 'within' subject
%   blocks, as returned in pr.D.

    pr = struct('ok', false, 'reason', '', 'rows', [], 'D', [], 'C', [], 'n_rows', 0, ...
        'n_subjects', 0, 'n_subjects_total', 0, 'group_sizes', [], 'dropped', {{}}, ...
        'run_on', '', 'n_small', 0);
    preset = lower(preset);
    % Level order: the same one the omnibus design used (snpm_glm_levels), so
    % a and b here are the a and b of D.posthoc.
    order = [];
    if isfield(opts, 'condition_levels') && ~isempty(opts.condition_levels) ...
            && isfield(opts, 'condition_col') && strcmp(factor_col, opts.condition_col)
        order = opts.condition_levels;
    end
    [~, gi] = snpm_glm_levels(meta.(factor_col), order);
    keep = (gi == a) | (gi == b);
    is_within = ismember(preset, {'rmanova', 'mixed2way'});

    if is_within
        [subj, ~, si] = unique(meta.(opts.subject_col), 'stable');
        ns = numel(subj);
        hasA = accumarray(si(gi == a), 1, [ns 1]) > 0;
        hasB = accumarray(si(gi == b), 1, [ns 1]) > 0;
        both = hasA & hasB;
        single = find((hasA | hasB) & ~both);
        pr.dropped = arrayfun(@(k) tostr(subj, k), single(:)', 'UniformOutput', false);
        if ~isempty(single)
            warning('snpm_glm_pair:singleLevelSubjects', ...
                ['Post-hoc %s: subject(s) %s have only one of the two levels and are left ' ...
                 'out of this comparison.'], label, strjoin(pr.dropped, ', '));
        end
        keep = keep & both(si);
        pr.n_small = sum(both);
        pr.n_subjects = sum(both);
        pr.n_subjects_total = ns;
        if pr.n_small < 2
            pr = skip(pr, label, 'fewer than 2 subjects have both levels');
            return
        end
        pr.run_on = sprintf('%d of %d subjects', pr.n_subjects, ns);
        if pr.n_small < 6
            pr.run_on = sprintf('%s (smallest possible two-sided p-value %.3g)', pr.run_on, 2 / 2^pr.n_small);
            warning('snpm_glm_pair:fewSubjects', ...
                ['Post-hoc %s: only %d subjects have both levels, so the smallest possible ' ...
                 'two-sided permutation p-value is %.3g.'], label, pr.n_small, 2 / 2^pr.n_small);
        end
    else
        na = sum(gi == a); nb = sum(gi == b);
        pr.group_sizes = [na nb];
        pr.n_subjects = na + nb;
        pr.n_subjects_total = numel(gi);
        pr.n_small = min(na, nb);
        pr.run_on = sprintf('%d and %d subjects', na, nb);
        if pr.n_small < 4
            warning('snpm_glm_pair:smallGroup', ...
                'Post-hoc %s: only %d in the smaller group.', label, pr.n_small);
        end
    end

    rows = find(keep);
    sub  = meta(rows, :);
    o = opts;
    if strcmp(preset, 'mixed2way'), o.effect = 'condition'; end
    D = snpm_glm_design(preset, sub, o);
    if size(D.C, 1) ~= 1 || rank(D.X) < size(D.X, 2)
        pr = skip(pr, label, ['the design for these two levels cannot be estimated (a covariate ' ...
            'may coincide with the two groups)']);
        return
    end
    if numel(rows) - size(D.X, 2) < 1
        pr = skip(pr, label, 'too few rows for these two levels to leave any error degrees of freedom');
        return
    end

    % The builder's own contrast is (its second level) - (its reference level).
    % Its reference is the first level of the SUBSET in the builder's level
    % order: first appearance in the subset, or opts.condition_levels when that
    % fixes the order. Ask the same function the builder asked, then report
    % mean(a) - mean(b).
    [~, gs] = snpm_glm_levels(sub.(factor_col), order);
    ref_is_a = (gi(rows(find(gs == 1, 1))) == a);
    if ref_is_a, sgn = -1; else, sgn = 1; end

    pr.ok     = true;
    pr.rows   = rows;
    pr.n_rows = numel(rows);
    pr.D      = D;
    pr.C      = sgn * D.C;
end

function pr = skip(pr, label, reason)
    pr.ok = false;
    pr.reason = reason;
    warning('snpm_glm_pair:skipped', 'Post-hoc %s was not run: %s.', label, reason);
end

function s = tostr(levels, k)
    if iscell(levels), v = levels{k}; else, v = levels(k); end
    if isnumeric(v) || islogical(v), s = num2str(v); else, s = char(string(v)); end
end
