function assignments = snpm_relabel_assignments(nSubj, nGrp, requested_permutations)
%SNPM_RELABEL_ASSIGNMENTS  Two-group relabelings for an unpaired permutation null.
%
%   assignments = SNPM_RELABEL_ASSIGNMENTS(nSubj, nGrp, requested_permutations)
%
%   Returns an (Nperm x nSubj) index matrix. Columns 1:nGrp are the group-1
%   subject indices for each relabeling; columns nGrp+1:end are the group-2
%   indices. Row 1 is ALWAYS the observed labelling [1:nGrp, nGrp+1:nSubj], so
%   the observed statistic is one of the Nperm permutations (Nichols & Holmes
%   2001).
%
%   The relabeling group for a two-group split has G = nchoosek(nSubj,nGrp)
%   DISTINCT group assignments (NOT nSubj!): the statistic is invariant to the
%   order of subjects within each group, so only the split matters.
%     * If requested_permutations >= G, ALL G assignments are enumerated exactly
%       (an exact test); Nperm = G. enumerate_combinations lists them with the
%       observed [1:nGrp | nGrp+1:nSubj] as its first row.
%     * Otherwise, requested_permutations random splits are drawn WITH
%       REPLACEMENT (randperm), with NO uniqueness loop, and row 1 is forced to
%       the observed labelling. Sampling with replacement is what makes the run
%       terminate even when requested_permutations >= G, which would exhaust a
%       without-replacement pool and hang forever (the historical bug this
%       function replaces).
%
%   Mirrors snpm_signflip_patterns (the paired/one-sample sign-flip null): exact
%   when the group fits the request, otherwise a Monte-Carlo sample with the
%   observed labelling forced into row 1.

    G = nchoosek(nSubj, nGrp);

    if requested_permutations >= G
        % Exact test: every distinct two-group assignment exactly once. Row 1 of
        % enumerate_combinations is [1:nGrp | nGrp+1:nSubj] (the observed split).
        assignments = enumerate_combinations(nSubj, nGrp);
    else
        % Monte-Carlo WITH REPLACEMENT (no uniqueness loop -> cannot hang).
        assignments = zeros(requested_permutations, nSubj);
        assignments(1, :) = 1:nSubj;                 % observed labelling
        for r = 2:requested_permutations
            assignments(r, :) = randperm(nSubj);
        end
    end
end
