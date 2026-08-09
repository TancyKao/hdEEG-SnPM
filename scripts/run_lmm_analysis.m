%% run_lmm_analysis.m  —  headless per-channel LMM SnPM runner (no GUI; for NCI)
% Per-channel linear mixed model (fitlme) with permutation FWE (Stephan 2021).
% Input is ONE long-format CSV: one row per observation, meta columns +
% per-channel value columns (E1.., or the montage labels). The fixed/random
% formula and the tested effect are set below.
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_lmm_analysis.m')"

% ======================= CONFIG (edit me) =======================
ROOT        = fileparts(fileparts(mfilename('fullpath')));   % scripts/ -> repo root
if isempty(ROOT), ROOT = pwd; end
data_file   = fullfile(ROOT,'example_data','lmm_long.csv'); % long format: meta cols + channel cols
output_dir  = fullfile(ROOT,'results_lmm');

channels    = '178 channels';   % '164 channels' | '178 channels'
datatype    = 'absolute';       % absolute|logscale|normalize
permutations= 100;
tail        = 'both';           % both|right|left

% model spec (the per-channel value enters the formula as POWER)
lmm_dv          = 'sleepDepth';     % dependent variable column (behaviour), OR 'POWER' if modelling power
lmm_fixed       = 'POWER*group';          % fixed-effect terms, e.g. 'POWER' or 'group*POWER + age'
lmm_random      = '(1|Subject)';    % random-effect spec
lmm_effect      = 'POWER';          % term whose effect is permutation-tested
lmm_effect_type = 'continuous';     % continuous (within-subject shuffle) | group (label shuffle)
lmm_meta_cols   = {'Subject','group','time','sleepDepth'};   % all non-channel columns
% ================================================================

addpath(genpath(ROOT));
if ~exist(output_dir,'dir'), mkdir(output_dir); end

params = struct('comparison','mixedmodel','data_file',data_file,'data_sheet','CSV File', ...
    'output_path',output_dir,'channels',channels,'datatype',datatype, ...
    'permutations',permutations,'tail',tail, ...
    'lmm_dv',lmm_dv,'lmm_fixed',lmm_fixed,'lmm_random',lmm_random, ...
    'lmm_effect',lmm_effect,'lmm_effect_type',lmm_effect_type, ...
    'lmm_meta_cols',{lmm_meta_cols});

[results_struct, results_text] = core_snpm_analysis(params); %#ok<ASGLU>
fprintf('\n%s\n', strjoin(results_text, newline));
fprintf('Done: LMM (%s ~ %s) -> %s\n', lmm_dv, lmm_fixed, output_dir);
