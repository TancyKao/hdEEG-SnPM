function info = snpm_glm_estimability(D, mode)
% Is the effect of a GLM design estimable apart from its Freedman-Lane nuisance?
%
%   info = snpm_glm_estimability(D)             % report only
%   info = snpm_glm_estimability(D, 'enforce')  % refuse / warn as below
%
% D is a snpm_glm_design struct (X, nuisance_idx, effect_label, col_names).
% Z = X(:, nuisance_idx), Xe = X(:, ~nuisance_idx), and
% tol = max(size(X)) * 2.220446049250313e-16 * s_max(X), computed once on X.
%
% EFFECT ALIASED. The identity contrast on Xe is estimable exactly when
% rank([Z Xe], tol) == rank(Z, tol) + cols(Xe). When it is not, snpm_glm_stat
% (pinv) still returns a finite t/F -- the minimum-norm split of one coefficient
% over two aliased columns -- and every p-value is arbitrary. 'enforce' refuses
% with core_snpm_glm:effectAliased, even when only one effect column is
% aliased (a partly aliased F would test the estimable remainder with the wrong
% df1 under a label that is false).
%
% NUISANCE ONLY. When Z alone is rank-deficient but the condition above holds,
% the contrast is estimable, the Freedman-Lane fitted values Z*(Z\Y) are the
% unique projection and df2 = n - rank(X) is already right, so the analysis
% runs: 'enforce' warns snpm:glmNuisanceRankDeficient.
%
% NAMING. Columns are walked left to right over Z first, then Xe: a column is
% kept when it raises rank(., tol) over the kept set, else it is aliased and
% listed with the kept columns k for which |w_k| * norm(x_k) / norm(x_j) > 1e-6
% in the least-squares fit of x_j on them ("is zero on every row" when none).
% Effect items name only aliased Xe columns, nuisance items only aliased Z columns.
%
% info fields: estimable (logical), nuisance_deficient (logical),
% effect_items, nuisance_items (cellstr).
%
% Message texts are shared with the Python port; do not reword.

    if nargin < 2, mode = ''; end
    X = D.X;
    nuis = logical(D.nuisance_idx(:)');
    names = D.col_names;
    if numel(names) ~= size(X, 2)
        names = arrayfun(@(j) sprintf('column %d', j), 1:size(X, 2), 'UniformOutput', false);
    end
    info = struct('estimable', true, 'nuisance_deficient', false, ...
        'effect_items', {{}}, 'nuisance_items', {{}});
    if isempty(X), return, end

    s = svd(X);
    tol = max(size(X)) * 2.220446049250313e-16 * s(1);
    Z = X(:, nuis);
    nXe = sum(~nuis);
    rZ = 0;
    if ~isempty(Z), rZ = rank(Z, tol); end
    info.estimable = rank(X(:, [find(nuis) find(~nuis)]), tol) == rZ + nXe;
    info.nuisance_deficient = rZ < size(Z, 2);

    if ~info.estimable || info.nuisance_deficient
        kept = [];
        for j = [find(nuis) find(~nuis)]
            if rank(X(:, [kept j]), tol) > numel(kept)
                kept(end+1) = j; %#ok<AGROW>
                continue
            end
            xj = X(:, j);
            of = {};
            if ~isempty(kept) && norm(xj) > 0
                w = X(:, kept) \ xj;
                rel = abs(w(:))' .* vecnorm(X(:, kept)) / norm(xj);
                of = names(kept(rel > 1e-6));
            end
            if isempty(of)
                item = sprintf('%s is zero on every row', names{j});
            else
                item = sprintf('%s is a linear combination of %s', names{j}, strjoin(of, ', '));
            end
            if nuis(j)
                info.nuisance_items{end+1} = item;
            else
                info.effect_items{end+1} = item;
            end
        end
    end

    if ~strcmpi(mode, 'enforce'), return, end
    if ~info.estimable
        error('core_snpm_glm:effectAliased', ...
            ['The effect tested (%s) cannot be separated from the nuisance columns of the design: %s. ' ...
             'Its contrast is not estimable, so the statistic and every p-value would be arbitrary. ' ...
             'Remove the covariate that repeats the effect, or test a different effect.'], ...
            D.effect_label, strjoin(info.effect_items, '; '));
    end
    if info.nuisance_deficient
        warning('snpm:glmNuisanceRankDeficient', ...
            ['The nuisance columns of the design are linearly dependent: %s. The effect is still ' ...
             'estimable and the analysis runs; the adjustment is for the space these columns span, ' ...
             'which the redundant column does not enlarge. Check the covariate list.'], ...
            strjoin(info.nuisance_items, '; '));
    end
end
