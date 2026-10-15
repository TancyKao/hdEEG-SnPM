function [levels, gi] = snpm_glm_levels(x, order)
%SNPM_GLM_LEVELS  Levels of a design column and the level index of each row.
%
%   [levels, gi] = snpm_glm_levels(x)         levels in order of first appearance
%   [levels, gi] = snpm_glm_levels(x, order)  levels in the order given by ORDER
%
% ORDER is a list of level values (same type as x). The levels present in x are
% returned in ORDER's sequence, followed by any level of x that ORDER does not
% name (in order of first appearance); a level in ORDER that x does not contain
% is ignored. With ORDER empty or absent this is unique(x,'stable').
%
% WHY ORDER EXISTS. Reference coding makes the first level the reference, and
% the pairwise post-hocs are labelled and signed "a vs b" with a before b in
% level order. A mixed2way report runs the condition effect once on the whole
% table and once inside each group; the level order inside one group's rows is
% first appearance WITHIN that group, which need not be the whole table's. One
% group could then show "mid vs pre" where the other shows "pre vs mid", with
% opposite sign. Passing the whole table's condition order here makes every
% part of the report list and sign its pairs identically.
%
% The design builder (snpm_glm_design) and the pairwise builder
% (snpm_glm_pair) both take their level order from this one function, so the
% level indices a / b in D.posthoc mean the same thing in both.

    [levels, ~, gi] = unique(x, 'stable');
    gi = gi(:);
    if nargin < 2 || isempty(order), return; end
    order = order(:);
    [inOrder, pos] = ismember(levels, order);
    rankv = zeros(numel(levels), 1);
    rankv(inOrder)  = pos(inOrder);                              % position in ORDER
    rankv(~inOrder) = numel(order) + (1:sum(~inOrder))';         % then the rest, as they appeared
    [~, perm] = sort(rankv);                                     % perm(k) = old index of new level k
    newidx = zeros(numel(levels), 1); newidx(perm) = 1:numel(levels);
    levels = levels(perm);
    gi = newidx(gi);
end
