% LMM_BOUNDARY_STALL_RATE  How often fitlme stops at a false zero subject variance.
%
% Reproduces the counts quoted in CHANGELOG.md ("Mixed model: a fit that
% stopped at zero subject variance"). Synthetic data only, fixed seeds.
%
% Design: 20 subjects x 15 trials, 3 groups (7/7/6 subjects), no group effect,
%   y = u_subject + e,  u ~ N(0, sd_u^2),  e ~ N(0, 1)
% fitted with POWER ~ group + (1|Subject) (ML, fitlme defaults, as
% snpm_lmm_fit uses), for 200 seeds (mt19937ar, seed = 1..200) at each sd_u.
%
% Columns printed:
%   at_boundary : plain fitlme returned subject SD <= 1e-3 x residual SD
%   stalled     : of those, snpm_lmm_fitlme's random-start refit found a
%                 log-likelihood higher by > 1e-6, i.e. the zero was NOT the
%                 ML optimum (the rest are genuine boundary optima)
%
% Run from the repo root (about 1 minute):
%   matlab -batch "run('scripts/lmm_boundary_stall_rate.m')"

ROOT = fileparts(fileparts(mfilename('fullpath')));
if isempty(ROOT), ROOT = pwd; end
addpath(fullfile(ROOT, 'dependencies'));

nS = 20; nT = 15; NSEED = 200;
subj = repelem((1:nS)', nT);
g = repelem([ones(7, 1); 2 * ones(7, 1); 3 * ones(6, 1)], nT);
formula = 'POWER ~ group + (1|Subject)';

fprintf('%-8s %-14s %-10s\n', 'sd_u', 'at_boundary', 'stalled');
for sdu = [0.2 0.4 1 2]
    nbd = 0; nst = 0;
    for seed = 1:NSEED
        rs = RandStream('mt19937ar', 'Seed', seed);
        u = sdu * randn(rs, nS, 1);
        y = u(subj) + randn(rs, nS * nT, 1);
        tbl = table(categorical(subj), categorical(g), y, 'VariableNames', {'Subject', 'group', 'POWER'});
        lme = fitlme(tbl, formula);
        [psi, mse] = covarianceParameters(lme);
        nbd = nbd + (sqrt(psi{1}) <= 1e-3 * sqrt(mse));
        [~, refit] = snpm_lmm_fitlme(tbl, formula);
        nst = nst + refit;
    end
    fprintf('%-8.1f %3d/%-10d %3d/%d\n', sdu, nbd, NSEED, nst, NSEED);
end
