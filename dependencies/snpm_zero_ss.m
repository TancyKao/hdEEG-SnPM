function tf = snpm_zero_ss(ss, v2)
%SNPM_ZERO_SS  The zero-variance tolerance, shared by every tier.
%
%   tf = SNPM_ZERO_SS(ss, v2)
%
%   ss : sum of squares that would be the error term (or the centred sum of
%        squares of a column), any shape
%   v2 : sum(v.^2) of the UNCENTRED values the statistic is computed from on
%        that channel, same shape as ss
%   tf : true where ss counts as zero, i.e. ss <= 1e-20 * v2
%
% Zero-variance rule (Method Spec, biostatistician round 3, approved by Tancy
% 2026-10-14). 1e-20 on a sum of squares is a relative SD of 1e-10: rounding
% leftovers sit near 1e-16 relative SD, real data many orders above, so MATLAB
% and the Python port make the same call. The comparison is against the
% UNCENTRED sum of squares so all-zero data (ss = 0, v2 = 0) also qualifies.
% An exact-zero test is not enough: two groups at 0.1 and 0.3 give a pooled SS
% of about 1e-33, not 0.
%
% Used by: snpm_glm_stat (GLM rule 2), snpm_glm_fl_context (GLM rule 1),
% snpm_ttest_novar (t-test rule 2), core_snpm_analysis (t-test and
% correlation rule 1). Change the constant here and nowhere else.

    TOL2 = 1e-20;
    tf = ss <= TOL2 .* v2;
end
