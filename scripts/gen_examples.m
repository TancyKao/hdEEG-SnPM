%% gen_examples.m  —  build test_data/examples/<analysis>/ : input CSV(s) + report
% One folder per analysis, each holding its synthetic_gui import CSV(s) and a
% freshly generated combined HTML report (.mat/.xlsx/topo PNGs too). Mirrors the
% GUI analyses; pairs with the Notion "SnPMAnalysisGUI-v2" docs.
%
% Run headless from repo root:
%   matlab -batch "run('scripts/gen_examples.m')"

ROOT = fileparts(fileparts(mfilename('fullpath')));   % scripts/ -> repo root
if isempty(ROOT), ROOT = pwd; end
SRC  = fullfile(ROOT,'test_data','synthetic_gui');
DST  = fullfile(ROOT,'test_data','examples');
addpath(genpath(ROOT));

channels = 'egi'; datatype = 'absolute'; nperm = 500;
if ~exist(DST,'dir'), mkdir(DST); end

% ---- GLM presets: single file -> core_snpm_glm ----
% columns: folder | comparison | csv | meta_cols | roles-struct
% (folder and comparison differ only for mixed2way_2x3, a 2nd mixed2way shape)
glm = {
  'anova1'        'anova1'    'glm_anova1.csv'    {'Subject','group'}              struct('group_col','group');
  'ancova'        'ancova'    'glm_ancova.csv'    {'Subject','group','age','sex'}  struct('group_col','group','covariate_cols',{{'age','sex'}});
  'ancova3'       'ancova'    'glm_ancova3.csv'   {'Subject','group','age','sex'}  struct('group_col','group','covariate_cols',{{'age','sex'}});
  'regression'    'regression' 'glm_regression.csv' {'Subject','cognition','age'}  struct('predictor_col','cognition','covariate_cols',{{'age'}});
  'rmanova'       'rmanova'   'glm_rmanova.csv'   {'Subject','condition'}          struct('subject_col','Subject','condition_col','condition');
  'mixed2way'     'mixed2way' 'glm_mixed2way.csv' {'Subject','group','condition'}  struct('group_col','group','subject_col','Subject','condition_col','condition');
  'mixed2way_2x3' 'mixed2way' 'glm_mixed2way_2x3.csv' {'Subject','group','condition'} struct('group_col','group','subject_col','Subject','condition_col','condition');
};
for i = 1:size(glm,1)
    key = glm{i,1}; comp = glm{i,2}; csv = glm{i,3}; meta = glm{i,4}; roles = glm{i,5};
    outdir = fullfile(DST,key); if ~exist(outdir,'dir'), mkdir(outdir); end
    copyfile(fullfile(SRC,csv), fullfile(outdir,csv));
    p = struct('comparison',comp,'data_file',fullfile(outdir,csv),'data_sheet','CSV File', ...
        'output_path',outdir,'channels',channels,'datatype',datatype,'permutations',nperm, ...
        'snpm_path',ROOT,'meta_cols',{meta});
    fn = fieldnames(roles); for k=1:numel(fn), p.(fn{k}) = roles.(fn{k}); end
    core_snpm_analysis(p);
    fprintf('Done: %s\n', key);
end

% ---- LMM: long-format single file -> core_snpm_lmm ----
outdir = fullfile(DST,'lmm'); if ~exist(outdir,'dir'), mkdir(outdir); end
copyfile(fullfile(SRC,'lmm_long.csv'), fullfile(outdir,'lmm_long.csv'));
p = struct('comparison','mixedmodel','data_file',fullfile(outdir,'lmm_long.csv'),'data_sheet','CSV File', ...
    'output_path',outdir,'channels',channels,'datatype',datatype,'permutations',nperm,'tail','both', ...
    'snpm_path',ROOT,'lmm_dv','sleepDepth','lmm_fixed','POWER','lmm_random','(1|Subject)', ...
    'lmm_effect','POWER','lmm_effect_type','continuous','lmm_meta_cols',{{'Subject','group','time','sleepDepth'}});
core_snpm_analysis(p);
fprintf('Done: lmm\n');

% ---- One-sample: single condition tested vs 0 (no Data 2) ----
outdir = fullfile(DST,'onesampleT'); if ~exist(outdir,'dir'), mkdir(outdir); end
copyfile(fullfile(SRC,'onesample_change.csv'), fullfile(outdir,'onesample_change.csv'));
p = struct('comparison','onesampleT','data1_file',fullfile(outdir,'onesample_change.csv'), ...
    'data1_sheet','CSV File','data2_file','','output_path',outdir,'channels',channels, ...
    'datatype',datatype,'permutations',nperm,'tail','both','use_covariates',false,'snpm_path',ROOT);
core_snpm_analysis(p);
fprintf('Done: onesampleT (single condition vs 0)\n');

% ---- Legacy two-file tests -> core_snpm_analysis (compstring) ----
leg = {
  'pairedT'     'paired_condA.csv'    'paired_condB.csv'    'pairedT';
  'unpairedT'   'unpaired_groupA.csv' 'unpaired_groupB.csv' 'unpairedT';
  'correlation' 'corr_eeg.csv'        'corr_behavior.csv'   'correlationP';
};
for i = 1:size(leg,1)
    key = leg{i,1}; a = leg{i,2}; b = leg{i,3}; comp = leg{i,4};
    outdir = fullfile(DST,key); if ~exist(outdir,'dir'), mkdir(outdir); end
    copyfile(fullfile(SRC,a), fullfile(outdir,a));
    copyfile(fullfile(SRC,b), fullfile(outdir,b));
    p = struct('comparison',comp, ...
        'data1_file',fullfile(outdir,a),'data1_sheet','CSV File', ...
        'data2_file',fullfile(outdir,b),'data2_sheet','CSV File', ...
        'output_path',outdir,'channels',channels,'datatype',datatype, ...
        'permutations',nperm,'tail','both','use_covariates',false,'snpm_path',ROOT);
    core_snpm_analysis(p);
    fprintf('Done: %s (%s)\n', key, comp);
end

% ---- Circular: CSV-only (CircStat not bundled -> no HTML) ----
outdir = fullfile(DST,'circular'); if ~exist(outdir,'dir'), mkdir(outdir); end
copyfile(fullfile(SRC,'circ_condA.csv'), fullfile(outdir,'circ_condA.csv'));
copyfile(fullfile(SRC,'circ_condB.csv'), fullfile(outdir,'circ_condB.csv'));
fid = fopen(fullfile(outdir,'README.md'),'w');
fprintf(fid, ['# Circular (Wheeler-Watson / Watson U2)\n\n', ...
    'Input CSVs only. These analyses are **not yet available** \x2014 the CircStat\n', ...
    'toolbox is not bundled in `dependencies/`, so no HTML report is generated.\n']);
fclose(fid);
fprintf('Done: circular (CSV only)\n');

% ---- examples/README.md : folder -> analysis -> GUI settings ----
fid = fopen(fullfile(DST,'README.md'),'w');
fprintf(fid, ['# Example analyses \x2014 one folder per analysis\n\n', ...
 'Each folder holds the import CSV(s) and the report SnPM produces from them.\n', ...
 'All use the **EGI 256 (HydroCel)** system (E1..E256), Data Type **absolute**, with a\n', ...
 'planted 17-channel cluster near E129. Generated by `scripts/gen_examples.m` (500 perms).\n\n', ...
 '| folder | analysis | input CSV(s) | roles |\n|---|---|---|---|\n', ...
 '| `anova1` | One-way ANOVA (3+ groups) | glm_anova1.csv | Group=group |\n', ...
 '| `ancova` | One-way ANCOVA (2 groups \x2192 t) | glm_ancova.csv | Group=group, Covariates=age,sex |\n', ...
 '| `ancova3` | One-way ANCOVA (3 groups \x2192 F + post-hoc) | glm_ancova3.csv | Group=group, Covariates=age,sex |\n', ...
 '| `regression` | Multiple linear regression | glm_regression.csv | Predictor=cognition, Covariates=age |\n', ...
 '| `rmanova` | Repeated-measures ANOVA | glm_rmanova.csv | Subject=Subject, Condition=condition |\n', ...
 '| `mixed2way` | Two-way mixed ANOVA (2\xD72: pre/post) | glm_mixed2way.csv | Group=group, Subject=Subject, Condition=condition |\n', ...
 '| `mixed2way_2x3` | Two-way mixed ANOVA (2\xD73: pre/mid/post) | glm_mixed2way_2x3.csv | Group=group, Subject=Subject, Condition=condition |\n', ...
 '| `pairedT` | Paired t-test | paired_condA.csv + paired_condB.csv | \x2014 |\n', ...
 '| `onesampleT` | One-sample t-test (single condition vs 0) | onesample_change.csv | \x2014 |\n', ...
 '| `unpairedT` | Two independent-samples t-test | unpaired_groupA.csv + unpaired_groupB.csv | \x2014 |\n', ...
 '| `correlation` | Pearson correlation | corr_eeg.csv + corr_behavior.csv | \x2014 |\n', ...
 '| `lmm` | Linear Mixed Model | lmm_long.csv | DV=sleepDepth, fixed=POWER |\n', ...
 '| `circular` | Wheeler-Watson / Watson U2 | circ_condA.csv + circ_condB.csv | not available (no CircStat) |\n']);
fclose(fid);
fprintf('\nAll examples written to %s\n', DST);
