%% run_glm_report.m  —  headless 3+ group ANOVA spectral report (no GUI; for NCI)
% Between-subject one-way ANOVA (omnibus F + post-hoc pairwise t), per stage x
% band, on the EGI montage, TFCE + cluster correction, group-mean / F / post-hoc
% topoplots, k-group periodograms, and a ready-to-open HTML report.
% Wraps export_glm_report.m.
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_glm_report.m')"

% ======================= CONFIG (edit me) =======================
ROOT = fileparts(fileparts(mfilename('fullpath')));   % scripts/ -> repo root
if isempty(ROOT), ROOT = pwd; end
DATA = fullfile(ROOT,'test_data','glm_report_demo');   % subjects.csv + <stage>_<power>_<band>.csv + gmPSD_*
OUT  = fullfile(ROOT,'test_data','glm_report_out');    % filled HTML + PNGs land here

opts = struct();
% preset picks the analysis; the report auto-routes by statistic:
%   'anova1'     group_col           -> F (3+ groups) or t (2 groups)
%   'ancova'     group_col + covariate_cols
%   'regression' predictor_col + covariate_cols  -> signed t
% F presets render the ANOVA template (omnibus F + post-hoc); t presets the
% t template (signed map, no post-hoc). ANCOVA reuses whichever by group count.
opts.preset         = 'anova1';
opts.group_col      = 'group';
% opts.predictor_col = 'cognition';        % regression
% opts.covariate_cols = {'age','sex'};     % ancova / regression nuisance
opts.powers       = {'absolute','normalised'};
opts.stages       = {'n2','n3'};
opts.bands        = {'delta','sigma','beta'};
opts.permutations = 1000;
opts.alpha        = 0.05;
% opts.title = 'Spectral power: 3-group ANOVA (omnibus F)';   % auto-derived if omitted
% opts.html  = ...                                            % override template path
% ================================================================

addpath(genpath(ROOT));
export_glm_report(DATA, OUT, opts);
fprintf('Done: ANOVA report -> %s\n', fullfile(OUT,'anova_report.html'));
