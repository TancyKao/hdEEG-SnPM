%% run_spectral_analysis.m  —  multi-folder spectral SnPM (no GUI; for NCI)
% Each BIDS spectral-power folder (EEG_processor output: a folder of
% sub-*_..._desc-<stage>_powerspect.mat) is ONE level of the design factor.
% List the folders + labels, pick a band + power type, and the analysis compares
% them per stage. A band x stage SWEEP writes a one-row-per-cell results grid.
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_spectral_analysis.m')"
%
% Supported comparisons (folder = factor level):
%   pairedT / onesampleT  - 2 folders, same subjects (within; matched by sub-XX)
%   unpairedT             - 2 folders, independent subjects (between)
%   anova1                - >=2 folders, between-group
%   rmanova               - >=2 folders, same subjects across all (within)
% (ANCOVA / regression / correlation / two-way mixed / LMM need per-subject
%  covariates or trial-level data -> use the CSV scripts, not this one.)

% ======================= CONFIG (edit me) =======================
ROOT = '/path/to/hdEEG-SnPM';
OUT  = fullfile(ROOT,'test_data','spectral_analysis_out');

% One folder per factor level; LABELS is parallel to FOLDERS.
FOLDERS = { ...
    fullfile(ROOT,'test_data','hdeeg_analysis_all_sub','condition-a'), ...
    fullfile(ROOT,'test_data','hdeeg_analysis_all_sub','condition-b') };
LABELS  = {'a','b'};

comparison = 'pairedT';        % pairedT | onesampleT | unpairedT | anova1 | rmanova
type       = 'absolute';       % 'absolute' | 'normalized'  (band POWER type)
sweep_bands  = {'sigma'};      % low-delta delta theta alpha sigma beta gamma
sweep_stages = {'n2'};         % n1 n2 n3 rem  (band x stage grid)

channels   = '178 channels';
permutations = 1000;
tail       = 'both';           % both | left | right
datatype   = '';               % '' = auto (absolute->logscale, normalized->none); or absolute|logscale|normalize
level_A    = '';               % optional: pick/order the 2 levels for 2-folder tests
level_B    = '';
MAKE_DASHBOARD = true;         % also build ONE faceted dashboard HTML (abs+relative x all stages x
                               % Uncorrected/TFCE/Cluster) via export_report -- 2-level contrasts only
                               % (pairedT/onesampleT/unpairedT). Ignores sweep_bands/stages/type (does all).
% ================================================================

cd(ROOT); addpath(genpath(ROOT));
if ~exist(OUT,'dir'), mkdir(OUT); end
if isempty(datatype), datatype = ternary(strcmpi(type,'absolute'),'logscale','absolute'); end

grid = {};
for bi = 1:numel(sweep_bands)
    for si = 1:numel(sweep_stages)
        band = sweep_bands{bi}; stg = sweep_stages{si};
        celldir = fullfile(OUT, sprintf('%s_%s_%s', comparison, band, stg));
        if ~exist(celldir,'dir'), mkdir(celldir); end

        o = struct('band',band, 'type',type, 'stages',{{stg}}, 'output_path',celldir, ...
            'channels',channels, 'datatype',datatype, 'permutations',permutations, 'tail',tail);
        o.folders = FOLDERS; o.labels = LABELS;
        if ~isempty(level_A) && ~isempty(level_B), o.level_A = level_A; o.level_B = level_B; end

        [params, ns] = spectral_to_snpm_params('', comparison, o);
        fprintf('\n==== %s | band=%s type=%s stage=%s | n=%d ====\n', comparison, band, type, stg, ns);
        rs = core_snpm_analysis(params);
        grid(end+1,:) = {comparison, band, type, stg, ns, ...
            numel(rs.uncorrsigch), numel(rs.correctTFCEsigch), numel(rs.SnPMsigch), min_clusterp(rs)}; %#ok<SAGROW>
        fprintf('   uncorr=%d  TFCE=%d  cluster=%d  minClusterP=%.3f\n', grid{end,6:9});
    end
end

G = cell2table(grid, 'VariableNames', ...
    {'comparison','band','type','stage','n','nUncorr','nTFCE','nCluster','minClusterP'});
gridfile = fullfile(OUT, sprintf('SWEEP_grid_%s_%s.csv', comparison, type));
writetable(G, gridfile);
fprintf('\n===== SWEEP GRID (%d cells) -> %s =====\n', height(G), gridfile);
disp(G);

% ---- faceted dashboard (one HTML: abs+relative x all stages x Uncorrected/TFCE/Cluster) ----
if MAKE_DASHBOARD
    if ~ismember(lower(comparison), {'pairedt','onesamplet','unpairedt'})
        fprintf('Dashboard skipped: only 2-level contrasts (pairedT/onesampleT/unpairedT) are supported; got %s.\n', comparison);
    else
        dash = fullfile(OUT, 'dashboard'); if ~exist(dash,'dir'), mkdir(dash); end
        copyfile(fullfile(ROOT,'templates','sleep_eeg_report.html'), fullfile(dash,'sleep_eeg_report.html'));
        dopts = struct('comparison',comparison, 'permutations',permutations, ...
            'powers',{{'absolute','normalised'}}, 'condA',LABELS{1}, 'condB',LABELS{2});
        dopts.folders = FOLDERS; dopts.labels = LABELS;
        export_report('', dash, dopts);
        fprintf('Dashboard -> %s\n', fullfile(dash,'sleep_eeg_report_filled.html'));
    end
end
fprintf('Done.\n');

% ============================================================ local functions
function p = min_clusterp(rs)
    p = NaN;
    if isfield(rs,'Clusters') && ~isempty(rs.Clusters) && isfield(rs.Clusters,'p')
        ps = [rs.Clusters.p]; if ~isempty(ps), p = min(ps); end
    end
end

function y = ternary(c, a, b), if c, y=a; else, y=b; end; end
