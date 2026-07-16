function orders = snpm_row_permutations(nSubj, requested_permutations)
%SNPM_ROW_PERMUTATIONS  Row reorderings for a correlation permutation null.
%
%   orders = SNPM_ROW_PERMUTATIONS(nSubj, requested_permutations)
%
%   Returns an (Nperm x nSubj) matrix; each row is a permutation of 1:nSubj used
%   to shuffle one variable's subject order, breaking the pairing for a
%   permutation test of association. Row 1 is ALWAYS the identity 1:nSubj, so the
%   observed alignment is one of the Nperm permutations (Nichols & Holmes 2001).
%
%   The permutation group has F = nSubj! distinct orderings.
%     * If requested_permutations >= F AND nSubj <= 10, ALL F orderings are
%       enumerated exactly via perms (an exact test); Nperm = F. The nSubj <= 10
%       ceiling bounds perms() memory (10! rows); nSubj >= 11 always falls to the
%       Monte-Carlo branch regardless of the request.
%     * Otherwise, requested_permutations random orderings are drawn WITH
%       REPLACEMENT (randperm), with NO uniqueness loop, and row 1 is forced to
%       the identity. With-replacement sampling removes the exhaustible pool that
%       caused the historical infinite-hang when requested_permutations >= F (the
%       bug this function replaces).
%
%   Mirrors snpm_signflip_patterns / snpm_relabel_assignments: exact when the
%   group fits the request, otherwise a Monte-Carlo sample with the observed
%   ordering forced into row 1.

    if nSubj <= 10 && requested_permutations >= factorial(nSubj)
        % Exact test: enumerate every ordering once, identity forced to row 1
        % (perms lists the identity last).
        orders = perms(1:nSubj);
        idRow = find(all(orders == (1:nSubj), 2), 1);
        orders([1 idRow], :) = orders([idRow 1], :);
    else
        % Monte-Carlo WITH REPLACEMENT (no uniqueness loop -> cannot hang).
        orders = zeros(requested_permutations, nSubj);
        orders(1, :) = 1:nSubj;                       % observed alignment
        for r = 2:requested_permutations
            orders(r, :) = randperm(nSubj);
        end
    end
end
