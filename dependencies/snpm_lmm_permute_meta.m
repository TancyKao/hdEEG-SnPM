function meta_perm = snpm_lmm_permute_meta(meta, spec, scheme)
% Produce one permuted copy of the trial-level metadata table for the LMM
% permutation null. The scheme is COUPLED to the effect being tested:
%
%   'within_subject' : shuffle the dependent variable within each subject's
%                      own trials. Breaks the trial-level predictor->DV link
%                      while preserving subject structure. Correct null for a
%                      within-subject / continuous predictor effect (Stephan
%                      et al. 2021).
%
%   'group_label'    : permute the group label across subjects (each subject
%                      keeps all its trials but may receive a different
%                      group). Correct null for a BETWEEN-group effect; a
%                      within-subject DV shuffle leaves group means unchanged
%                      and therefore provides no null for that effect.
%
% INPUTS
%   meta   : trial-level table (one row per trial)
%   spec   : struct with .dv, and .subject (default 'Subject'),
%            .group (default 'group')
%   scheme : 'within_subject' | 'group_label'
%
% OUTPUT
%   meta_perm : copy of meta with either the DV or the group column permuted.

    if isfield(spec, 'subject') && ~isempty(spec.subject)
        subj_col = spec.subject;
    else
        subj_col = 'Subject';
    end
    if isfield(spec, 'group') && ~isempty(spec.group)
        grp_col = spec.group;
    else
        grp_col = 'group';
    end

    meta_perm = meta;

    switch lower(scheme)
        case 'within_subject'
            dv = meta.(spec.dv);
            subj = meta.(subj_col);
            [~, ~, gidx] = unique(subj, 'stable');
            for g = 1:max(gidx)
                rows = find(gidx == g);
                if numel(rows) > 1
                    dv(rows) = dv(rows(randperm(numel(rows))));
                end
            end
            meta_perm.(spec.dv) = dv;

        case 'group_label'
            subj = meta.(subj_col);
            grp  = meta.(grp_col);
            % sidx maps each trial to its subject index; broadcasting the
            % permuted subject-level labels back with sidx is type-agnostic
            % (works for string / cellstr / categorical / numeric IDs, unlike
            % a `subj == usubj(s)` comparison which errors on cellstr).
            [usubj, ia, sidx] = unique(subj, 'stable');
            subj_groups = grp(ia);                      % each subject's group label
            perm = randperm(numel(usubj));
            meta_perm.(grp_col) = subj_groups(perm(sidx)); % permuted labels per trial

        otherwise
            error('snpm_lmm_permute_meta:badScheme', ...
                'scheme must be ''within_subject'' or ''group_label''');
    end
end
