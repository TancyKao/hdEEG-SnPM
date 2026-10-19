% EXPORT_LMM_BACKSTOP_SWEEP  Keep/drop sweep of the LMM per-labelling backstop.
%
% Writes, for the Python port to replay (synthetic data, fixed seed):
%   test_fixtures/parity/lmm_backstop_sweep_data.csv      Subject, group, y_<k>
%   test_fixtures/parity/lmm_backstop_sweep_expected.csv  one row per y_<k>
% Model: y ~ group + (1|Subject), 12 subjects x 5 trials, 2 groups (6/6).
% y_k = 10 + 0.3*group + a_k * b(Subject) + e, with b and e fixed and a_k solved
% so that the method-of-moments subject/residual SD ratio of
% dependencies/snpm_lmm_degenerate.m equals the target exactly: 1e2, 5e2,
% 0.999e3, 1.001e3, 2e3, 1e4. The backstop flags theta >= 1e3, so the expected
% decisions are 0 0 0 1 1 1. Each column is also fitted with fitlme and the
% fitted-model call must give the same decision as the raw-array call.
%
%   matlab -batch "run('scripts/export_lmm_backstop_sweep.m')"

ROOT = fileparts(fileparts(mfilename('fullpath')));
if isempty(ROOT), ROOT = pwd; end
addpath(fullfile(ROOT, 'dependencies'));
OUT = fullfile(ROOT, 'test_fixtures', 'parity');
if ~exist(OUT, 'dir'), mkdir(OUT); end

[T, X, Z, targets, a] = lmm_backstop_sweep_data();
names = arrayfun(@(k) sprintf('y_%d', k), 1:numel(targets), 'UniformOutput', false);
rows = cell(numel(targets), 8);
for k = 1:numel(targets)
    y = T.(names{k});
    [tf, info] = snpm_lmm_degenerate(y, X, Z, true);
    tbl = table(categorical(T.Subject), categorical(T.group), y, 'VariableNames', {'Subject', 'group', 'POWER'});
    lme = fitlme(tbl, 'POWER ~ group + (1|Subject)');
    tf_lme = snpm_lmm_degenerate(lme);
    assert(tf == tf_lme, 'raw-array and fitted-model decisions differ for %s', names{k});
    rows(k, :) = {names{k}, sprintf('%.17g', targets(k)), sprintf('%.17g', a(k)), ...
        sprintf('%.17g', sqrt(info.theta2)), sprintf('%.17g', info.rss_x), sprintf('%.17g', info.rss_xz), ...
        sprintf('%d', tf), sprintf('%.17g', snpm_lmm_within_ratio(y, T.Subject))};
end
fid = fopen(fullfile(OUT, 'lmm_backstop_sweep_data.csv'), 'w');
fprintf(fid, 'Subject,group,%s\n', strjoin(names, ','));
M = T{:, names};
for r = 1:height(T)
    fprintf(fid, '%d,%s%s\n', T.Subject(r), T.group{r}, sprintf(',%.17g', M(r, :)));
end
fclose(fid);
fid = fopen(fullfile(OUT, 'lmm_backstop_sweep_expected.csv'), 'w');
fprintf(fid, 'column,target_theta,scale_a,theta_mom,rss_x,rss_xz,degenerate,within_ratio\n');
for k = 1:size(rows, 1), fprintf(fid, '%s\n', strjoin(rows(k, :), ',')); end
fclose(fid);
disp(cell2table(rows, 'VariableNames', {'column','target','a','theta_mom','rss_x','rss_xz','degenerate','w'}));

function [T, X, Z, targets, a] = lmm_backstop_sweep_data()
    rs = RandStream('mt19937ar', 'Seed', 20261020);
    q = 12; n = 5; N = q * n;
    S = repelem((1:q)', n);
    g = repelem([zeros(q/2, 1); ones(q/2, 1)], n);
    b = randn(rs, q, 1); e = randn(rs, N, 1);
    X = [ones(N, 1), g];
    Z = double(S == 1:q);
    targets = [1e2 5e2 0.999e3 1.001e3 2e3 1e4];
    a = zeros(size(targets));
    T = table(S, repmat({'G1'}, N, 1), 'VariableNames', {'Subject', 'group'});
    T.group(g == 1) = {'G2'};
    for k = 1:numel(targets)
        f = @(la) log(theta_of(10 + 0.3 * g + exp(la) * b(S) + e, X, Z)) - log(targets(k));
        a(k) = exp(fzero(f, log(targets(k) / 2), optimset('TolX', 1e-14)));
        T.(sprintf('y_%d', k)) = 10 + 0.3 * g + a(k) * b(S) + e;
    end
end

function th = theta_of(y, X, Z)
    [~, info] = snpm_lmm_degenerate(y, X, Z, true);
    th = sqrt(info.theta2);
end
