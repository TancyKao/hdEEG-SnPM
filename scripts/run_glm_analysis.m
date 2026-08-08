%% run_glm_analysis.m  —  headless GLM/legacy SnPM runner (no GUI; for NCI batch)
% Runs a per-channel SnPM analysis (anova1 / ancova / regression / rmanova /
% mixed2way, or legacy pairedT / unpairedT / onesampleT / correlationP /
% correlationS).
%
% Two input shapes (DATA_MODE):
%   'single' — ONE CSV that already holds the metadata columns + one column per
%              channel (E1..E256 / Cz). This is the LONG format required by
%              rmanova / mixed2way (one row per Subject x condition) and also
%              works for any preset. Set DATA_FILE + META_COLS.
%   'join'   — merge a subject x channel MATRIX_CSV (rows = Subject, cols =
%              channels) with a SUBJECTS_CSV of per-subject metadata on the
%              'Subject' key. One row per subject -> between-subjects designs
%              (anova1 / ancova / regression / 2-group). NOT for mixed2way.
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_glm_analysis.m')"
% As shipped this runs the two-way mixed ANOVA example end to end. Edit the
% CONFIG block for your own data.

% ======================= CONFIG (edit me) =======================
ROOT        = fileparts(fileparts(mfilename('fullpath')));   % scripts/ -> repo root
if isempty(ROOT), ROOT = pwd; end

comparison  = 'mixed2way';     % anova1|ancova|regression|rmanova|mixed2way|pairedT|unpairedT|onesampleT|correlationP|correlationS
DATA_MODE   = 'single';        % 'single' (one long/wide CSV) | 'join' (matrix + subjects.csv)

channels    = 'egi';           % recording system: 'egi' (EGI 256/HydroCel) | 'compu' (Compumedics 257)
                               %   legacy aliases '164 channels'/'178 channels' still resolve to egi
datatype    = 'absolute';      % absolute|logscale|normalize
permutations= 500;             % production: 10000+ ; quick check: 500-1000

% ---- 'single' mode: one CSV with meta columns + channel columns ----
% Long format for mixed2way/rmanova: one row per Subject x condition.
DATA_FILE   = fullfile(ROOT,'test_data','examples','mixed2way_2x3','glm_mixed2way_2x3.csv');
META_COLS   = {'Subject','group','condition'};   % non-channel columns; everything else = channels

% ---- 'join' mode: subject x channel matrix + per-subject metadata ----
MATRIX_CSV  = fullfile(ROOT,'test_data','synthetic_events','eventStat_density_spindle_12-15.csv'); % rows=Subject, cols=channels
SUBJECTS_CSV= fullfile(ROOT,'test_data','synthetic_events','subjects.csv');   % Subject + group/cognition/...
DROP_CZ     = true;            % keep E1..E256, drop a trailing Cz column if present

output_dir  = fullfile(ROOT,'test_data','examples','mixed2way_2x3');   % as shipped: reproduce the example in place

% column roles (set only those the chosen comparison needs; leave '' otherwise)
group_col     = 'group';       % anova1/ancova/mixed2way (group main effect)
predictor_col = '';            % regression
condition_col = 'condition';   % rmanova/mixed2way (within factor)
subject_col   = 'Subject';     % rmanova/mixed2way
covariate_cols= {};            % ancova/regression, e.g. {'age','sex'}
% ================================================================

addpath(genpath(ROOT));
if ~exist(output_dir,'dir'), mkdir(output_dir); end

switch lower(DATA_MODE)
    case 'single'
        % Use the CSV as-is; META_COLS names the non-channel columns.
        glm_csv  = DATA_FILE;
        metacols = META_COLS;

    case 'join'
        M    = readtable(MATRIX_CSV,   'TextType','string','VariableNamingRule','preserve');
        subj = readtable(SUBJECTS_CSV, 'TextType','string','VariableNamingRule','preserve');
        chanCols = setdiff(M.Properties.VariableNames,{'Subject'},'stable');
        if DROP_CZ
            keep = chanCols(~strcmpi(chanCols,'Cz'));
            M = M(:, ['Subject', keep]);
        end
        T = innerjoin(subj, M, 'Keys','Subject');       % metadata + channels
        glm_csv  = fullfile(output_dir,'glm_input.csv'); writetable(T, glm_csv);
        metacols = subj.Properties.VariableNames;        % all subjects.csv columns are meta

    otherwise
        error('run_glm_analysis:mode', 'DATA_MODE must be ''single'' or ''join'' (got ''%s'').', DATA_MODE);
end

params = struct('comparison',comparison,'data_file',glm_csv,'data_sheet','CSV File', ...
    'output_path',output_dir,'channels',channels,'datatype',datatype,'permutations',permutations, ...
    'meta_cols',{metacols});
if ~isempty(group_col),      params.group_col=group_col; end
if ~isempty(predictor_col),  params.predictor_col=predictor_col; end
if ~isempty(condition_col),  params.condition_col=condition_col; end
if ~isempty(subject_col),    params.subject_col=subject_col; end
if ~isempty(covariate_cols), params.covariate_cols=covariate_cols; end

[results_struct, results_text] = core_snpm_analysis(params); %#ok<ASGLU>
fprintf('\n%s\n', strjoin(results_text, newline));
fprintf('Done: %s analysis -> %s\n', comparison, output_dir);
