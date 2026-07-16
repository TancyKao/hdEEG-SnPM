function patterns = snpm_signflip_patterns(nSubj, requested_permutations)
%SNPM_SIGNFLIP_PATTERNS  Sign-flip labellings for a paired/one-sample null.
%
%   patterns = SNPM_SIGNFLIP_PATTERNS(nSubj, requested_permutations)
%
%   Returns an (Nperm x nSubj) LOGICAL matrix of sign-flip labellings for the
%   Nichols & Holmes (2001) sign-flip permutation group used by the paired and
%   one-sample tests. Rows are permutations, columns are subjects.
%
%   Semantics match the paired/one-sample engines: patterns(r,i) == true means
%   subject i's paired difference is KEPT for permutation r; false means it is
%   sign-flipped (data_x/data_y swapped for that subject). The observed /
%   identity labelling is the all-true row (no flips), which reproduces the
%   observed statistic exactly.
%
%   Sampling of the sign-flip group (2^nSubj patterns):
%     * If 2^nSubj <= requested_permutations, ALL 2^nSubj patterns are
%       enumerated exactly (an exact test); Nperm = 2^nSubj. The all-true
%       identity row is guaranteed to be present.
%     * Otherwise, requested_permutations RANDOM sign vectors are drawn, each
%       subject independently flipped with probability 0.5 so every subject has
%       flip frequency ~0.5 across the null (no subject is ever frozen). Row 1
%       is forced to the all-true identity so the observed labelling is
%       included, per Nichols & Holmes.
%
%   This replaces the previous first-K-integers enumeration
%   (dec2binvec(1:K, nSubj)), which is NOT a uniform sample of the sign-flip
%   group: for K < 2^nSubj it froze the high-index subjects (flip frequency 0)
%   and miscalibrated family-wise error.

exact_count = 2^nSubj;

if exact_count <= requested_permutations
    % Exact test: enumerate every sign pattern once.
    patterns = false(exact_count, nSubj);
    for k = 1:(exact_count - 1)
        patterns(k+1, :) = logical(bitget(k, 1:nSubj));
    end
    % k = 0 leaves row 1 all-false; the all-true identity is at k = exact_count-1.
else
    % Random sampling of the sign-flip group. Each entry is an independent
    % fair coin; row 1 is the observed (all-keep) labelling.
    patterns = rand(requested_permutations, nSubj) > 0.5;
    patterns(1, :) = true;
end
end
