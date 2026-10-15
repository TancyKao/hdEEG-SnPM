%% export_parity_fixtures.m  —  MATLAB reference fixtures for the Python port
% MATLAB is the reference implementation: a Python analysis is trusted only when
% it reproduces the MATLAB result on the same inputs using the SAME permutation
% indices. For every case below this script
%   1. seeds the generator            rng(SEED, 'twister')
%   2. switches the recorder on       snpm_perm_tap('on')
%   3. runs the production pipeline   core_snpm_analysis(params)
%   4. collects the indices it drew   snpm_perm_tap('get')
%   5. saves ONE fixture              <OUT_DIR>/<case_name>.mat   (-v7)
% and then writes <OUT_DIR>/README.md (one row per fixture + field dictionary).
%
% Fixtures hold only numeric arrays, logicals, char row vectors, cellstr and
% plain structs / struct arrays of those, so scipy.io.loadmat reads them:
%   scipy.io.loadmat(path, simplify_cells=True)
%
% Datatype: every case is 'absolute' except pairedT_logscale and
% correlationP_normalize (legacy tier). The GLM replay supports 'absolute' only.
%
% Scope: phase 1 = pairedT, unpairedT, onesampleT, correlationP, correlationS
% (with and without covariates), anova1, ancova, regression, rmanova, mixed2way.
% The circular and mixed-model tiers are recorded by the same hook but are not
% exported here.
%
% Nothing in the production results_struct is changed for this. What the
% pipeline computes but does not keep comes from two places. The hook records
% the legacy engine inputs, the legacy cluster-size null and the GLM
% cluster-mass null where the pipeline computes them. The GLM design matrices
% are rebuilt HERE with the same dependency functions; each GLM loop is then
% replayed from the recorded indices and must reproduce the pipeline's own
% output, including the recorded mass null, before the fixture is saved
% (glm.loops(k).replay_ok).
%
% Run headless:
%   matlab -batch "run('scripts/export_parity_fixtures.m')"

% ======================= CONFIG (edit me) =======================
ROOT         = fileparts(fileparts(mfilename('fullpath')));   % scripts/ -> repo root
if isempty(ROOT), ROOT = pwd; end

OUT_DIR      = fullfile(ROOT, 'test_fixtures', 'parity');     % git-ignored
SEED         = 20261013;      % rng(SEED,'twister') immediately before every case
PERMUTATIONS = 500;           % requested permutations per loop
INPUT_SEED   = 4242;          % seeds the inputs this script BUILDS (compu, NaN, covariates)
ONLY_CASES   = {};            % {} = all cases; or a cellstr of case names to (re)export
RUN_DIR      = fullfile(tempdir, 'snpm_parity_runs');   % each run's own .mat/.xlsx/.html/PNGs
                              % (kept OUTSIDE the repo: core_snpm_analysis does
                              % addpath(genpath(pwd)), so a run folder under ROOT
                              % would land on the MATLAB path)
RESUME       = false;         % true = keep a fixture that already exists with this SEED and
                              % PERMUTATIONS and only run the missing cases (use after an
                              % interrupted export; false regenerates everything)
KEEP_RUN_OUTPUTS = false;     % false = delete each run's outputs once its fixture is saved
% ================================================================

cd(ROOT);                     % core_snpm_analysis resolves assets from pwd
addpath(genpath(ROOT));
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
IN_DIR  = fullfile(OUT_DIR, 'inputs');
if ~exist(IN_DIR, 'dir'),  mkdir(IN_DIR);  end
if ~exist(RUN_DIR, 'dir'), mkdir(RUN_DIR); end
EX = fullfile(ROOT, 'example_data');

built = build_inputs(EX, IN_DIR, INPUT_SEED);
cases = case_table(EX, built);
if ~isempty(ONLY_CASES)
    unknown = setdiff(ONLY_CASES, {cases.name});
    if ~isempty(unknown)
        error('export_parity:unknownCase', 'Unknown case name(s): %s', strjoin(unknown, ', '));
    end
    cases = cases(ismember({cases.name}, ONLY_CASES));
end

fprintf('export_parity_fixtures: %d cases, %d permutations, seed %d -> %s\n', ...
    numel(cases), PERMUTATIONS, SEED, OUT_DIR);

summary = cell(numel(cases), 1);
for ci = 1:numel(cases)
    c = cases(ci);
    out_file = fullfile(OUT_DIR, [c.name '.mat']);
    if RESUME && exist(out_file, 'file')
        old = load(out_file); %#ok<UNRCH>
        if isequal(old.seed, SEED) && isequal(old.permutations, PERMUTATIONS)
            summary{ci} = readme_row(c, old, SEED, PERMUTATIONS);
            fprintf('[%2d/%2d] %-34s kept (RESUME)\n', ci, numel(cases), c.name);
            continue
        end
    end
    P = c.params;
    P.permutations = PERMUTATIONS;
    P.output_path  = fullfile(RUN_DIR, c.name);
    if exist(P.output_path, 'dir'), rmdir(P.output_path, 's'); end
    mkdir(P.output_path);

    t0 = tic;
    [results_struct, results_text, raw_log] = run_case(P, SEED);
    % Three kinds of entry are not permutation draws: the matrices the legacy
    % engines were handed and the legacy cluster-size null (-> fx.engine), and
    % the GLM cluster-mass null of each loop (-> fx.glm.loops).
    is_aux  = ismember({raw_log.tag}, {'legacy:engine_inputs', 'legacy:cluster_null', 'glm:cluster_mass_null'});
    aux_log = raw_log(is_aux);
    raw_log = raw_log(~is_aux);
    close all force

    fx = struct();
    fx.case_name      = c.name;
    fx.tier           = c.tier;
    fx.description    = c.note;
    fx.seed           = SEED;
    fx.rng_algorithm  = 'twister';
    fx.permutations   = PERMUTATIONS;
    fx.matlab_version = version;
    fx.created        = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    fx.params         = flat(P);
    fx.source_files   = source_files(P, ROOT);
    fx.perm_log       = group_log(raw_log);
    fx.results        = flat(results_struct);
    fx.results_text   = flat(results_text(:));
    if strcmp(c.tier, 'glm')
        fx.inputs  = glm_inputs(P);
        [fx.montage, fx.glm] = glm_capture(P, results_struct, fx.perm_log, aux_log);
        replay_txt = sprintf('replay %d/%d loops ok', ...
            sum([fx.glm.loops.replay_ok]), numel(fx.glm.loops));
        if ~all([fx.glm.loops.replay_ok])
            error('export_parity:replayMismatch', ...
                ['%s: the recorded indices did not reproduce the pipeline''s own null ' ...
                 '(%s). The fixture would not be a valid reference; not saved.'], ...
                c.name, replay_txt);
        end
    else
        fx.inputs  = legacy_inputs(P);
        fx.montage = legacy_montage(P, fx.inputs, results_struct);
        fx.engine  = legacy_engine(aux_log, results_struct, c.name);
        replay_txt = 'replay n/a (legacy)';
    end

    assert_plain(fx, c.name);
    save(out_file, '-struct', 'fx', '-v7');
    chk = load(out_file);                      % every fixture must load back
    assert(isequaln(chk.results, fx.results) && isequaln(chk.perm_log, fx.perm_log), ...
        'export_parity:roundTrip', '%s: fixture did not round-trip through save/load.', c.name);

    if ~KEEP_RUN_OUTPUTS, rmdir(P.output_path, 's'); end
    summary{ci} = readme_row(c, fx, SEED, PERMUTATIONS);
    fprintf('[%2d/%2d] %-34s %2d perm loops, %s, %.0f s\n', ci, numel(cases), c.name, ...
        numel(fx.perm_log), replay_txt, toc(t0));
end
if ~KEEP_RUN_OUTPUTS && exist(RUN_DIR, 'dir') && numel(dir(RUN_DIR)) <= 2
    rmdir(RUN_DIR);
end

if isempty(ONLY_CASES)
    write_readme(fullfile(OUT_DIR, 'README.md'), cases, summary, SEED, PERMUTATIONS, INPUT_SEED, built);
    fprintf('README written: %s\n', fullfile(OUT_DIR, 'README.md'));
else
    fprintf('ONLY_CASES set: README.md left as it was (run with ONLY_CASES = {} to rewrite it).\n');
end
fprintf('export_parity_fixtures: DONE (%d fixtures)\n', numel(cases));

% =======================================================================
% ONE RUN
% =======================================================================
function [results_struct, results_text, raw_log] = run_case(P, SEED) %#ok<INUSD>
% Seed, switch the recorder on, run the pipeline, collect the log. The
% onCleanup lives in this function's workspace, so the recorder is switched
% off however the run ends -- normal return, error, or Ctrl-C -- and an
% interrupted export cannot leave it recording.
    snpm_perm_tap('off'); snpm_perm_tap('clear');
    guard = onCleanup(@() snpm_perm_tap('off'));
    rng(SEED, 'twister');
    snpm_perm_tap('on');
    [~, results_struct, results_text] = evalc('core_snpm_analysis(P)');
    snpm_perm_tap('off');
    raw_log = snpm_perm_tap('get');
    snpm_perm_tap('clear');
end

% =======================================================================
% CASES
% =======================================================================
function cases = case_table(EX, B)
    cases = struct('name', {}, 'tier', {}, 'analysis', {}, 'params', {}, 'note', {});
    ex = @(f) fullfile(EX, f);

    % ---- 1. every phase-1 analysis on example_data, EGI, tail 'both' ----
    cases(end+1) = legacy('pairedT_both', 'pairedT', 'both', ex('paired_condA.csv'), ex('paired_condB.csv'), 'egi', '', ...
        'Paired t-test, 24 subjects, sign-flip null (random branch).');
    cases(end+1) = legacy('unpairedT_both', 'unpairedT', 'both', ex('unpaired_groupA.csv'), ex('unpaired_groupB.csv'), 'egi', '', ...
        'Unpaired t-test, 16 vs 16, group-relabelling null (random branch).');
    cases(end+1) = legacy('onesampleT_single', 'onesampleT', 'both', ex('onesample_change.csv'), '', 'egi', '', ...
        'One-sample t-test of one file against 0 (no Data 2; data_y is a zero matrix).');
    cases(end+1) = legacy('onesampleT_pair', 'onesampleT', 'both', ex('onesample_condA.csv'), ex('onesample_condB.csv'), 'egi', '', ...
        'One-sample t-test, legacy two-file form (A - B against 0).');
    cases(end+1) = legacy('correlationP_wide', 'correlationP', 'both', ex('corr_eeg.csv'), ex('corr_behavior.csv'), 'egi', '', ...
        'Pearson, channel i vs channel i, 30 subjects, row-permutation null.');
    cases(end+1) = legacy('correlationS_wide', 'correlationS', 'both', ex('corr_eeg.csv'), ex('corr_behavior.csv'), 'egi', '', ...
        'Spearman, channel i vs channel i, 30 subjects.');
    cases(end+1) = legacy('correlationP_single', 'correlationP', 'both', ex('corr_eeg.csv'), ex('corr_behavior_single.csv'), 'egi', '', ...
        'Pearson against ONE behavioural column (MRRT), broadcast across channels.');
    cases(end+1) = legacy('correlationS_single', 'correlationS', 'both', ex('corr_eeg.csv'), ex('corr_behavior_single.csv'), 'egi', '', ...
        'Spearman against ONE behavioural column (MRRT), broadcast across channels.');
    cases(end+1) = legacy('correlationP_wide_cov', 'correlationP', 'both', ex('corr_eeg.csv'), ex('corr_behavior.csv'), 'egi', B.cov, ...
        'Partial Pearson: both matrices residualized on [1, age, bmi] before permutation.');
    cases(end+1) = legacy('correlationS_wide_cov', 'correlationS', 'both', ex('corr_eeg.csv'), ex('corr_behavior.csv'), 'egi', B.cov, ...
        'Partial Spearman: rank first, residualize the ranks, then Pearson engines.');
    cases(end+1) = legacy('correlationP_single_cov', 'correlationP', 'both', ex('corr_eeg.csv'), ex('corr_behavior_single.csv'), 'egi', B.cov, ...
        'Partial Pearson against one behavioural column.');
    cases(end+1) = legacy('correlationS_single_cov', 'correlationS', 'both', ex('corr_eeg.csv'), ex('corr_behavior_single.csv'), 'egi', B.cov, ...
        'Partial Spearman against one behavioural column.');

    cases(end+1) = glm('anova1', 'anova1', ex('glm_anova1.csv'), 'egi', {'Subject','group'}, ...
        struct('group_col', 'group'), 'One-way ANOVA, 3 groups: omnibus F + 3 pairwise post-hoc t.');
    cases(end+1) = glm('ancova_2group', 'ancova', ex('glm_ancova.csv'), 'egi', {'Subject','group','age','sex'}, ...
        struct('group_col', 'group', 'covariate_cols', {{'age','sex'}}), ...
        'ANCOVA, 2 groups -> t; nuisance = intercept + age + sex (sex dummy-coded).');
    cases(end+1) = glm('ancova_3group', 'ancova', ex('glm_ancova3.csv'), 'egi', {'Subject','group','age','sex'}, ...
        struct('group_col', 'group', 'covariate_cols', {{'age','sex'}}), ...
        'ANCOVA, 3 groups -> F + pairwise post-hoc; nuisance = intercept + age + sex.');
    cases(end+1) = glm('regression', 'regression', ex('glm_regression.csv'), 'egi', {'Subject','cognition','age'}, ...
        struct('predictor_col', 'cognition', 'covariate_cols', {{'age'}}), ...
        'Regression slope of cognition (t), nuisance = intercept + age.');
    cases(end+1) = glm('rmanova', 'rmanova', ex('glm_rmanova.csv'), 'egi', {'Subject','condition'}, ...
        struct('subject_col', 'Subject', 'condition_col', 'condition'), ...
        'Repeated-measures ANOVA, 3 conditions: F + post-hoc, within-subject permutation.');
    cases(end+1) = glm('mixed2way_2x2', 'mixed2way', ex('glm_mixed2way.csv'), 'egi', {'Subject','group','condition'}, ...
        struct('group_col', 'group', 'condition_col', 'condition', 'subject_col', 'Subject'), ...
        'Two-way mixed ANOVA 2x2: group, condition, interaction + simple effects.');
    cases(end+1) = glm('mixed2way_2x3', 'mixed2way', ex('glm_mixed2way_2x3.csv'), 'egi', {'Subject','group','condition'}, ...
        struct('group_col', 'group', 'condition_col', 'condition', 'subject_col', 'Subject'), ...
        'Two-way mixed ANOVA 2x3: condition F + post-hoc, interaction F, simple effects with post-hoc.');

    % ---- 2. one-sided tails ----
    cases(end+1) = legacy('pairedT_left', 'pairedT', 'left', ex('paired_condA.csv'), ex('paired_condB.csv'), 'egi', '', ...
        'Paired t, tail left (the engines swap data_x and data_y and test right).');
    cases(end+1) = legacy('pairedT_right', 'pairedT', 'right', ex('paired_condA.csv'), ex('paired_condB.csv'), 'egi', '', ...
        'Paired t, tail right.');
    cases(end+1) = legacy('unpairedT_left', 'unpairedT', 'left', ex('unpaired_groupA.csv'), ex('unpaired_groupB.csv'), 'egi', '', ...
        'Unpaired t, tail left (the engines stack [data_y; data_x] and test right).');
    cases(end+1) = legacy('unpairedT_right', 'unpairedT', 'right', ex('unpaired_groupA.csv'), ex('unpaired_groupB.csv'), 'egi', '', ...
        'Unpaired t, tail right.');

    % ---- 3. exact-enumeration branches (pin the enumeration row order) ----
    cases(end+1) = legacy('exact_pairedT_n6', 'pairedT', 'both', B.small_pairA, B.small_pairB, 'egi', '', ...
        'First 6 subjects of the paired example: 2^6 = 64 sign patterns enumerated.');
    cases(end+1) = legacy('exact_unpairedT_4v4', 'unpairedT', 'both', B.small_unpA, B.small_unpB, 'egi', '', ...
        'First 4 subjects of each unpaired example file: nchoosek(8,4) = 70 splits enumerated.');
    cases(end+1) = legacy('exact_correlationP_n5', 'correlationP', 'both', B.small_corrX, B.small_corrY, 'egi', '', ...
        'First 5 subjects of the correlation example: 5! = 120 orderings enumerated.');

    % ---- 4. Compumedics 257 recording system ----
    cases(end+1) = legacy('compu_pairedT', 'pairedT', 'both', B.compu_pairA, B.compu_pairB, 'compu', '', ...
        'Compumedics montage, paired t, 12 subjects, planted neighbour cluster.');
    cases(end+1) = legacy('compu_unpairedT', 'unpairedT', 'both', B.compu_unpA, B.compu_unpB, 'compu', '', ...
        'Compumedics montage, unpaired t, 10 vs 10, planted neighbour cluster.');
    cases(end+1) = legacy('compu_correlationP', 'correlationP', 'both', B.compu_corrX, B.compu_corrY, 'compu', '', ...
        'Compumedics montage, Pearson channel-vs-channel, 14 subjects.');
    cases(end+1) = glm('compu_anova1', 'anova1', B.compu_anova, 'compu', {'Subject','group'}, ...
        struct('group_col', 'group'), 'Compumedics montage, one-way ANOVA, 3 groups of 8.');

    % ---- 5. missing data ----
    cases(end+1) = legacy('miss_pairedT_scattered', 'pairedT', 'both', B.miss_pairA, B.miss_pairB, 'egi', '', ...
        'Scattered NaN cells in both conditions; the t-map is NOT complete-column masked.');
    cases(end+1) = legacy('miss_correlationP_scattered', 'correlationP', 'both', B.miss_corrX, B.miss_corrY, 'egi', '', ...
        'Scattered NaN cells; complete-column rule excludes every channel with any NaN.');
    cases(end+1) = legacy('miss_unpairedT_group_channel', 'unpairedT', 'both', B.miss_unpA, B.miss_unpB, 'egi', '', ...
        ['One channel all-NaN in group B only, one channel with too few usable subjects ' ...
         '(both removed by the degeneracy mask), one channel with 2 NaN that stays in.']);
    cases(end+1) = legacy('miss_correlationP_cov_nan', 'correlationP', 'both', ex('corr_eeg.csv'), ex('corr_behavior.csv'), 'egi', B.cov_nan, ...
        'Partial Pearson where one subject has a missing covariate value (subject dropped whole).');

    % ---- 6. branches of the uncorrected Spearman p, and the datatype transforms ----
    cases(end+1) = legacy('exact_correlationS_n7', 'correlationS', 'both', B.s7_corrX, B.s7_corrY, 'egi', '', ...
        ['First 7 subjects, Spearman, no ties: corr() uses its exact permutation p (n < 10). ' ...
         '7! = 5040 > 500, so the permutation loops are on the random branch.']);
    cases(end+1) = legacy('exact_correlationS_n7_ties', 'correlationS', 'both', B.s7_corrX_ties, B.s7_corrY, 'egi', '', ...
        'Same 7 subjects with one tied pair planted in 5 EEG channels: exact-with-ties p on those channels.');
    cases(end+1) = legacy('correlationS_ties_n30', 'correlationS', 'both', B.ties_corrX, ex('corr_behavior.csv'), 'egi', '', ...
        '30 subjects, Spearman, one tied pair planted in 5 EEG channels: large-sample p with and without ties.');
    cases(end+1) = legacy('pairedT_logscale', 'pairedT', 'both', ex('paired_condA.csv'), ex('paired_condB.csv'), 'egi', '', ...
        'Paired t with datatype logscale (log10 of every cell; the example data are all positive).');
    cases(end).params.datatype = 'logscale';
    cases(end+1) = legacy('correlationP_normalize', 'correlationP', 'both', ex('corr_eeg.csv'), ex('corr_behavior.csv'), 'egi', '', ...
        ['Pearson with datatype normalize (each subject z-scored across the analysed channels). ' ...
         'The global test is flagged not interpretable by the pipeline in this configuration.']);
    cases(end).params.datatype = 'normalize';
end

function c = legacy(name, comparison, tail, f1, f2, channels, covfile, note)
    P = struct('comparison', comparison, 'tail', tail, ...
        'data1_file', f1, 'data1_sheet', 'CSV File', ...
        'data2_file', f2, 'data2_sheet', 'CSV File', ...
        'channels', channels, 'datatype', 'absolute', ...
        'use_covariates', ~isempty(covfile), 'covariate_file', covfile);
    if isempty(f2), P.data2_sheet = ''; end
    c = struct('name', name, 'tier', 'legacy', 'analysis', [comparison ' / ' tail], ...
        'params', P, 'note', note);
end

function c = glm(name, preset, f, channels, meta_cols, roles, note)
    P = struct('comparison', preset, 'data_file', f, 'data_sheet', 'CSV File', ...
        'channels', channels, 'datatype', 'absolute');
    P.meta_cols = meta_cols;
    fn = fieldnames(roles);
    for k = 1:numel(fn), P.(fn{k}) = roles.(fn{k}); end
    c = struct('name', name, 'tier', 'glm', 'analysis', preset, 'params', P, 'note', note);
end

% =======================================================================
% INPUTS THIS SCRIPT BUILDS (seeded; written under OUT_DIR/inputs)
% =======================================================================
function B = build_inputs(EX, IN_DIR, seed)
% A private stream, so building inputs never touches the global generator the
% analyses are seeded on.
    rs = RandStream('mt19937ar', 'Seed', seed);
    in = @(f) fullfile(IN_DIR, f);
    rd = @(f) readtable(fullfile(EX, f));
    B = struct();

    % ---- covariates for the correlation example (example_data ships none) ----
    X = rd('corr_eeg.csv');
    n = height(X);
    C = table(X.Subject, round(40 + 12 * randn(rs, n, 1)), round(26 + 3.5 * randn(rs, n, 1), 1), ...
        'VariableNames', {'Subject', 'age', 'bmi'});
    B.cov = in('corr_covariates.csv');                 writetable(C, B.cov);
    Cn = C; Cn.age(7) = NaN;                           % subject 7 has no age
    B.cov_nan = in('corr_covariates_missing.csv');     writetable(Cn, B.cov_nan);

    % ---- exact-enumeration inputs: the first rows of the example files ----
    A = rd('paired_condA.csv');   Bt = rd('paired_condB.csv');
    B.small_pairA = in('small_paired_condA_n6.csv');   writetable(A(1:6, :),  B.small_pairA);
    B.small_pairB = in('small_paired_condB_n6.csv');   writetable(Bt(1:6, :), B.small_pairB);
    Ua = rd('unpaired_groupA.csv'); Ub = rd('unpaired_groupB.csv');
    B.small_unpA = in('small_unpaired_groupA_n4.csv'); writetable(Ua(1:4, :), B.small_unpA);
    B.small_unpB = in('small_unpaired_groupB_n4.csv'); writetable(Ub(1:4, :), B.small_unpB);
    Y = rd('corr_behavior.csv');
    B.small_corrX = in('small_corr_eeg_n5.csv');       writetable(X(1:5, :), B.small_corrX);
    B.small_corrY = in('small_corr_behavior_n5.csv');  writetable(Y(1:5, :), B.small_corrY);

    % ---- missing data ----
    scalp = analysed_labels('egi', X.Properties.VariableNames);   % only scalp channels are analysed
    % (a) scattered cells
    B.miss_pairA = in('miss_paired_condA.csv');   writetable(scatter_nan(A,  scalp, 40, rs), B.miss_pairA);
    B.miss_pairB = in('miss_paired_condB.csv');   writetable(scatter_nan(Bt, scalp, 40, rs), B.miss_pairB);
    B.miss_corrX = in('miss_corr_eeg.csv');       writetable(scatter_nan(X, scalp, 12, rs), B.miss_corrX);
    B.miss_corrY = in('miss_corr_behavior.csv');  writetable(scatter_nan(Y, scalp, 8,  rs), B.miss_corrY);
    % (b) unpaired: group-confined and degenerate channels
    pick = scalp(randperm(rs, numel(scalp), 3));
    Ma = Ua; Mb = Ub;
    Mb.(pick{1})(:) = NaN;                      % all of group B missing -> 16 usable < 18
    Ma.(pick{2})(4:end) = NaN;                  % 3 usable in A ...
    Mb.(pick{2})(5:end) = NaN;                  % ... 4 usable in B -> 7 usable < 18
    Ma.(pick{3})([2 9]) = NaN;                  % 30 usable >= 18 -> stays in the map
    B.miss_unpA = in('miss_unpaired_groupA.csv'); writetable(Ma, B.miss_unpA);
    B.miss_unpB = in('miss_unpaired_groupB.csv'); writetable(Mb, B.miss_unpB);
    B.miss_unp_channels = pick;

    % ---- Compumedics 257: synthetic, planted effect on one neighbour cluster ----
    m = snpm_montage_registry('compu');
    labels = m.labels(:)';
    assert(all(cellfun(@isvarname, labels)), 'compu labels must be valid column names');
    S = snpm_setup_channels('compu', {zeros(1, numel(labels))}, labels);
    lab = {S.chanlocs.labels};
    % cluster = the best-connected channel, its neighbours and their neighbours
    [~, centre] = max(sum(~isnan(S.neighbors), 2));
    ring1 = S.neighbors(centre, :);       ring1 = ring1(~isnan(ring1));
    ring2 = S.neighbors(ring1, :);        ring2 = ring2(~isnan(ring2));
    cl = ismember(labels, lab(unique([centre, ring1, ring2(:)'])));
    nCh = numel(labels);
    subj = @(k) arrayfun(@(i) sprintf('sub%03d', i), (k(:)), 'UniformOutput', false);
    mk = @(ids, M) [table(subj(ids), 'VariableNames', {'Subject'}), array2table(M, 'VariableNames', labels)];

    Pa = randn(rs, 12, nCh); Pb = Pa * 0.5 + randn(rs, 12, nCh); Pa(:, cl) = Pa(:, cl) + 1.5;
    B.compu_pairA = in('compu_paired_condA.csv');   writetable(mk(1:12, Pa), B.compu_pairA);
    B.compu_pairB = in('compu_paired_condB.csv');   writetable(mk(1:12, Pb), B.compu_pairB);

    Ga = randn(rs, 10, nCh); Gb = randn(rs, 10, nCh); Ga(:, cl) = Ga(:, cl) + 1.8;
    B.compu_unpA = in('compu_unpaired_groupA.csv'); writetable(mk(1:10, Ga),    B.compu_unpA);
    B.compu_unpB = in('compu_unpaired_groupB.csv'); writetable(mk(101:110, Gb), B.compu_unpB);

    z = randn(rs, 14, 1);
    Cx = randn(rs, 14, nCh); Cy = randn(rs, 14, nCh);
    Cx(:, cl) = 0.5 * Cx(:, cl) + z; Cy(:, cl) = 0.5 * Cy(:, cl) + z;
    B.compu_corrX = in('compu_corr_eeg.csv');       writetable(mk(1:14, Cx), B.compu_corrX);
    B.compu_corrY = in('compu_corr_behavior.csv');  writetable(mk(1:14, Cy), B.compu_corrY);

    An = randn(rs, 24, nCh); g = repmat({'G1'; 'G2'; 'G3'}, 8, 1);
    An(strcmp(g, 'G3'), cl) = An(strcmp(g, 'G3'), cl) + 1.8;
    T = mk(1:24, An); T = [T(:, 1), table(g, 'VariableNames', {'group'}), T(:, 2:end)];
    B.compu_anova = in('compu_glm_anova1.csv');     writetable(T, B.compu_anova);
    B.compu_cluster = labels(cl);

    % ---- Spearman p-value branches (drawn LAST so every input above keeps the
    %      values it had before these were added) ----
    tie_ch = scalp(randperm(rs, numel(scalp), 5));
    Xt = X;
    for i = 1:numel(tie_ch), Xt.(tie_ch{i})(2) = Xt.(tie_ch{i})(1); end   % subjects 1 and 2 share a value
    B.s7_corrX      = in('small_corr_eeg_n7.csv');        writetable(X(1:7, :),  B.s7_corrX);
    B.s7_corrY      = in('small_corr_behavior_n7.csv');   writetable(Y(1:7, :),  B.s7_corrY);
    B.s7_corrX_ties = in('small_corr_eeg_n7_ties.csv');   writetable(Xt(1:7, :), B.s7_corrX_ties);
    B.ties_corrX    = in('ties_corr_eeg.csv');            writetable(Xt,         B.ties_corrX);
    B.tie_channels  = tie_ch;
end

function T = scatter_nan(T, scalp, k, rs)
% Set k cells to NaN, each in a different scalp channel and a random row.
    cols = scalp(randperm(rs, numel(scalp), k));
    for i = 1:k
        T.(cols{i})(randi(rs, height(T))) = NaN;
    end
end

function lab = analysed_labels(system, varnames)
    cols = varnames(~contains(lower(varnames), 'subject'));
    S = snpm_setup_channels(system, {zeros(1, numel(cols))}, cols);
    lab = {S.chanlocs.labels};
end

% =======================================================================
% CAPTURE: inputs as loaded
% =======================================================================
function I = legacy_inputs(P)
    I = struct();
    I.data1 = table_block(P.data1_file, {});
    if isfield(P, 'data2_file') && ~isempty(P.data2_file)
        I.data2 = table_block(P.data2_file, {});
    else
        I.data2 = [];                      % one-sample against 0: no second file
    end
    if isfield(P, 'use_covariates') && P.use_covariates && ~isempty(P.covariate_file)
        % Same reading rule as core_snpm_analysis: column 1 = subject id, the
        % rest are covariates.
        T = readtable(P.covariate_file);
        I.covariates = struct('file', P.covariate_file, ...
            'subject_ids', {id_column(T{:, 1})}, ...
            'names', {T.Properties.VariableNames(2:end)}, ...
            'values', T{:, 2:end});
    else
        I.covariates = [];
    end
end

function I = glm_inputs(P)
    I = struct();
    I.data = table_block(P.data_file, P.meta_cols);
end

function blk = table_block(file, meta_cols)
% One input file exactly as readtable returns it (the call the core makes).
%   meta_cols {} -> legacy rule: the FIRST column whose name contains 'subject'
%                   is the subject id; every other column is data.
    T = readtable(file);
    vn = T.Properties.VariableNames;
    if isempty(meta_cols)
        k = find(contains(lower(vn), 'subject'), 1);
        meta_cols = vn(k);
    end
    data_cols = setdiff(vn, meta_cols, 'stable');
    blk = struct();
    blk.file          = file;
    blk.all_columns   = vn;                         % header order in the file
    blk.column_labels = data_cols;                  % columns of .values
    blk.values        = T{:, data_cols};            % rows x numel(column_labels)
    blk.meta_names    = meta_cols;
    blk.meta          = struct();
    for i = 1:numel(meta_cols)
        blk.meta.(meta_cols{i}) = id_column(T.(meta_cols{i}));
    end
    k = find(contains(lower(meta_cols), 'subject'), 1);
    if isempty(k), blk.subject_ids = {}; else, blk.subject_ids = blk.meta.(meta_cols{k}); end
end

function v = id_column(v)
% A metadata column as a numeric column vector or an n x 1 cellstr.
    if isnumeric(v) || islogical(v)
        v = double(v(:));
    else
        v = cellstr(string(v(:)));
    end
end

function S = source_files(P, ROOT)
    S = struct();
    for f = {'data1_file', 'data2_file', 'data_file', 'covariate_file'}
        if isfield(P, f{1}) && ~isempty(P.(f{1}))
            a = P.(f{1});
            r = a;
            if strncmp(a, [ROOT filesep], numel(ROOT) + 1), r = a(numel(ROOT) + 2:end); end
            S.(f{1}) = struct('absolute', a, 'relative_to_repo', strrep(r, filesep, '/'));
        end
    end
end

% =======================================================================
% CAPTURE: montage (what the pipeline analysed, in its order)
% =======================================================================
function M = legacy_montage(P, I, rs)
    S = snpm_setup_channels(P.channels, {zeros(1, numel(I.data1.column_labels))}, I.data1.column_labels);
    M = montage_struct(P.channels, S);
    assert(isequal(M.labels, {rs.chanlocs.labels}), 'export_parity:montageOrder', ...
        'Re-derived channel order differs from results_struct.chanlocs.');
end

function M = montage_struct(key, S)
    reg = snpm_montage_registry(key);
    M = struct();
    M.system_key   = reg.key;
    M.display      = reg.display;
    M.labels       = {S.chanlocs.labels};                               % 1 x nCh, analysed order
    M.data_columns = S.channel_mapping.original_data_columns(:)';       % into inputs.*.column_labels
    M.neighbors    = S.neighbors;                                       % nCh x maxNb, NaN padded
    xyz = nan(numel(S.chanlocs), 3);
    if all(isfield(S.chanlocs, {'X', 'Y', 'Z'}))
        for i = 1:numel(S.chanlocs)
            if ~isempty(S.chanlocs(i).X)
                xyz(i, :) = [S.chanlocs(i).X, S.chanlocs(i).Y, S.chanlocs(i).Z];
            end
        end
    end
    M.xyz = xyz;
    M.tfce_E = 0.5; M.tfce_H = 2; M.tfce_dh = 0.1; M.alpha = 0.05;
end

% =======================================================================
% CAPTURE: what the legacy engines were handed, and what they do not return
% =======================================================================
function E = legacy_engine(aux, rs, name)
% data_x / data_y are the matrices core_snpm_analysis passed to
% global_stat_test, snpm_single_threshold_with_TFCE and snpm_cluster_analysis
% (recorded immediately before the first of those calls; nothing modifies them
% in between). The cluster-size null is recorded inside snpm_cluster_analysis
% after its sort. channels_kept is obtained by calling global_stat_test again
% on the recorded matrices, and that call must reproduce the run's global.
    ei = aux(strcmp({aux.tag}, 'legacy:engine_inputs'));
    cn = aux(strcmp({aux.tag}, 'legacy:cluster_null'));
    assert(isscalar(ei) && isscalar(cn), 'export_parity:auxCount', ...
        '%s: expected one engine_inputs and one cluster_null record, got %d and %d.', ...
        name, numel(ei), numel(cn));
    v = ei.value;
    E = struct();
    E.data_x     = v.data_x;
    E.data_y     = v.data_y;
    E.comparison = v.comparison;
    E.tail       = v.tail;
    E.ncov       = v.ncov;
    E.alpha      = v.alpha;
    E.cluster_threshold   = cn.value.threshold;          % as passed to the engine
    E.cluster_thresholds  = cn.value.thresholds(:)';      % as used: [-thr thr] for tail both
    E.cluster_null_sorted = cn.value.null_sorted;
    assert(isequal(v.threshold, cn.value.threshold), 'export_parity:thresholdMismatch', ...
        '%s: the core and the cluster engine disagree on the cluster-forming threshold.', name);
    [~, G] = evalc('global_stat_test(v.data_x, v.data_y, v.alpha, v.comparison, v.tail, v.ncov)');
    assert(isequaln(G.stat, rs.global_stat) && isequaln(G.pval, rs.global_pval), ...
        'export_parity:globalMismatch', ...
        '%s: global_stat_test on the recorded engine inputs did not reproduce the run''s global.', name);
    E.global_channels_kept = find(G.channels_kept(:)');   % 1-based, analysed order
    E.global_stat = G.stat;
    E.global_pval = G.pval;
    % The recorded null must be the one the run's cluster p-values came from.
    N = size(E.cluster_null_sorted, 1);
    for k = 1:numel(rs.Clusters)
        col = 1;
        if size(E.cluster_null_sorted, 2) > 1
            col = find(rs.Clusters(k).threshold == E.cluster_thresholds);
        end
        pk = (sum(E.cluster_null_sorted(:, col) >= numel(rs.Clusters(k).channels)) + 1) / (N + 1);
        assert(isequal(pk, rs.Clusters(k).p), 'export_parity:clusterNullMismatch', ...
            '%s: recorded cluster null does not reproduce Clusters(%d).p.', name, k);
    end
end

% =======================================================================
% CAPTURE: permutation log, one entry per permutation LOOP, in call order
% =======================================================================
function L = group_log(raw)
% Legacy generators record one whole matrix per loop. The GLM tier records one
% row index per draw, preceded by 'scope:' / 'loop:' markers; those are stacked
% into one matrix per loop here.
    L = struct('name', {}, 'tag', {}, 'kind', {}, 'branch', {}, 'site', {}, 'caller', {}, ...
        'n_draws', {}, 'n_units', {}, 'indices', {});
    scope = ''; open = 0; buf = {};
    for i = 1:numel(raw)
        tag = raw(i).tag;
        if strncmp(tag, 'scope:', 6)
            L = close_loop(L, open, buf); open = 0; buf = {};
            scope = tag(7:end);
        elseif strncmp(tag, 'loop:', 5)
            L = close_loop(L, open, buf); buf = {};
            if isempty(scope), nm = tag(6:end); else, nm = [scope '/' tag(6:end)]; end
            L(end+1) = new_loop(nm, '', 'row_index', '', raw(i).site, ''); %#ok<AGROW>
            open = numel(L);
        elseif any(strncmp(tag, {'signflip:', 'relabel:', 'rowperm:'}, 8))
            L = close_loop(L, open, buf); open = 0; buf = {};
            parts = strsplit(tag, ':');
            engine = strtok(raw(i).caller, ':');
            switch engine
                case 'snpm_single_threshold_with_TFCE', nm = 'legacy/TFCE';
                case 'snpm_cluster_analysis',           nm = 'legacy/cluster';
                otherwise,                              nm = ['legacy/' engine];
            end
            e = new_loop(nm, tag, parts{1}, parts{2}, raw(i).site, raw(i).caller);
            e.indices = double(raw(i).value);
            e.n_draws = size(e.indices, 1); e.n_units = size(e.indices, 2);
            L(end+1) = e; %#ok<AGROW>
        else
            if open == 0
                error('export_parity:orphanDraw', ...
                    'Draw ''%s'' from %s is not inside a marked loop.', tag, raw(i).site);
            end
            if isempty(L(open).tag)
                L(open).tag = tag; L(open).caller = raw(i).caller;
                parts = strsplit(tag, ':'); L(open).branch = parts{end};
            elseif ~strcmp(L(open).tag, tag)
                error('export_parity:mixedLoop', 'Loop %s mixes %s and %s.', L(open).name, L(open).tag, tag);
            end
            buf{end+1} = double(raw(i).value(:)'); %#ok<AGROW>
        end
    end
    L = close_loop(L, open, buf);
    L = L(:);
end

function e = new_loop(name, tag, kind, branch, site, caller)
    e = struct('name', name, 'tag', tag, 'kind', kind, 'branch', branch, 'site', site, ...
        'caller', caller, 'n_draws', 0, 'n_units', 0, 'indices', []);
end

function L = close_loop(L, open, buf)
    if open == 0, return; end
    if ~isempty(buf)
        L(open).indices = vertcat(buf{:});
    end
    L(open).n_draws = size(L(open).indices, 1);
    L(open).n_units = size(L(open).indices, 2);
end

% =======================================================================
% CAPTURE: GLM designs + replay of every permutation loop
% =======================================================================
function [M, G] = glm_capture(P, rs, plog, aux)
% Rebuilds, with the same dependency functions core_snpm_glm calls, the design
% of every permutation loop the run executed, replays the loop from the
% RECORDED indices, and checks the replay against the pipeline's own output.
% What the replay adds to the fixture is the part the pipeline does not keep:
% X, C, nuisance columns and blocks. The sorted cluster-mass null saved is the
% one snpm_perm_correction recorded; the replay's own is only a cross-check.
    assert(strcmp(P.datatype, 'absolute'), 'glm_capture replays datatype ''absolute'' only.');
    alpha = 0.05; E = 0.5; H = 2;
    T = readtable(P.data_file);
    channel_cols = setdiff(T.Properties.VariableNames, P.meta_cols, 'stable');
    meta  = T(:, P.meta_cols);
    raw   = T{:, channel_cols};
    S = snpm_setup_channels(P.channels, {raw}, channel_cols);
    power = S.data{1};
    M = montage_struct(P.channels, S);
    assert(isequal(M.labels, {rs.chanlocs.labels}), 'export_parity:montageOrder', ...
        'Re-derived channel order differs from results_struct.chanlocs.');

    opts = struct();
    for f = {'group_col', 'condition_col', 'subject_col', 'predictor_col', 'covariate_cols', 'effect'}
        if isfield(P, f{1}), opts.(f{1}) = P.(f{1}); end
    end

    % ---- the plan: one entry per permutation loop, in the pipeline's order ----
    % Method of 2026-10-14: every pairwise post-hoc is its own two-level design
    % on the rows of its two levels (snpm_glm_pair, the function the core
    % calls), and every mixed2way loop shares one channel mask.
    plan = struct('role', {}, 'Y', {}, 'rows', {}, 'D', {}, 'C', {}, 'ctype', {}, 'mask', {}, ...
        'T', {}, 'p', {}, 'Clusters', {});
    every = (1:size(power, 1))';
    if strcmpi(P.comparison, 'mixed2way')
        % one mask: channels finite in every row of the long table
        mask = all(isfinite(power), 1);
        % (1) group main effect on subject means -- mirrors compute_group_effect
        %     (plain means: the pipeline has already required a balanced table)
        [sid, ~, si] = unique(meta.(opts.subject_col), 'stable');
        nS = numel(sid);
        subj_power = nan(nS, size(power, 2)); first = zeros(nS, 1);
        for k = 1:nS
            r = (si == k);
            subj_power(k, :) = mean(power(r, :), 1);
            first(k) = find(r, 1);
        end
        gv = meta.(opts.group_col);
        subj_meta = table(gv(first), 'VariableNames', {opts.group_col});
        gopts = struct('group_col', opts.group_col);
        plan = add_effect(plan, 'group', subj_power, first, 'anova1', subj_meta, gopts, mask, rs.effects(1), true);
        % (2) condition main effect: full model, group-by-condition columns in the nuisance
        oc = opts; oc.effect = 'condition';
        plan = add_effect(plan, 'condition', power, every, 'mixed2way', meta, oc, mask, rs.effects(2), true);
        % (3) interaction (no post-hoc)
        ox = opts; ox.effect = 'interaction';
        plan = add_effect(plan, 'interaction', power, every, 'mixed2way', meta, ox, mask, rs.effects(3), false);
        % simple effects: condition within each group
        [glev, ~, gi] = unique(meta.(opts.group_col), 'stable');
        se = rs.effects(3).simple_effects;
        assert(numel(se) == numel(glev), 'simple-effect count differs from the number of groups');
        sopts = struct('subject_col', opts.subject_col, 'condition_col', opts.condition_col);
        sopts.condition_levels = snpm_glm_levels(meta.(opts.condition_col));   % whole-table order, as the core passes
        for g = 1:numel(glev)
            rows = find(gi == g);
            plan = add_effect(plan, se(g).key, power(rows, :), rows, 'rmanova', meta(rows, :), sopts, mask, se(g), true);
        end
    else
        eff = struct('T', rs.T, 'p', rs.p, 'Clusters', rs.Clusters, 'posthoc', rs.posthoc);
        plan = add_effect(plan, 'main', power, every, P.comparison, meta, opts, true(1, size(power, 2)), eff, true);
    end

    assert(numel(plan) == numel(plog), 'export_parity:loopCount', ...
        'Rebuilt %d GLM loops but the run recorded %d.', numel(plan), numel(plog));
    % One recorded cluster-mass null per snpm_perm_correction call, i.e. per
    % loop, in the same order as the loops.
    mass_rec = aux(strcmp({aux.tag}, 'glm:cluster_mass_null'));
    assert(numel(mass_rec) == numel(plog), 'export_parity:massNullCount', ...
        'The run recorded %d cluster-mass nulls for %d loops.', numel(mass_rec), numel(plog));

    G = struct();
    G.preset = P.comparison;
    G.scheme = '2026-10-14';      % pairwise post-hoc on its own rows; mixed2way one mask + full-model condition effect
    G.alpha = alpha; G.tfce_E = E; G.tfce_H = H;
    G.loops = struct('name', {}, 'role', {}, 'Y', {}, 'rows', {}, 'X', {}, 'C', {}, ...
        'contrast_type', {}, 'nuisance_idx', {}, 'nuisance_cols', {}, 'eb', {}, 'perm_type', {}, ...
        'effect_label', {}, 'evaluable', {}, 'stat', {}, 'p_uncorrected', {}, 'tfce', {}, ...
        'tfce_null_sorted', {}, 'cluster_mass_null_sorted', {}, 'cluster_mass_null_source', {}, ...
        'replay_ok', {}, 'replay_detail', {});
    sparse_adj = make_neighbors_sparse(S.neighbors, size(S.neighbors, 1));
    for k = 1:numel(plan)
        q = plan(k);
        idx = plog(k).indices;
        flctx = snpm_glm_fl_context(q.Y, q.D);
        ev = flctx.evaluable & q.mask;
        isF = strcmpi(q.ctype, 'F');
        [stat, pm] = snpm_glm_stat(q.Y, q.D.X, q.C);
        stat(~ev) = NaN; pm(~ev) = NaN;
        tfce = ClusterEnhancement(stat, sparse_adj, E, H);
        nP = size(idx, 1);
        maxT = zeros(nP, 1); maxM = zeros(nP, 1);
        for i = 1:nP
            Yp = flctx.Zfit + flctx.R(idx(i, :), :);          % snpm_glm_permute, recorded pe
            [s, pp] = snpm_glm_stat(Yp, q.D.X, q.C);
            s(~ev) = NaN; pp(~ev) = NaN;
            e = ClusterEnhancement(s, sparse_adj, E, H);
            if isF, maxT(i) = max(e, [], 'omitnan'); else, maxT(i) = max(abs(e), [], 'omitnan'); end
            maxM(i) = max_cluster_mass(pp, s, isF, alpha, sparse_adj);
        end
        maxT = sort(maxT, 'descend'); maxM = sort(maxM, 'descend');

        % ---- the replay must reproduce the pipeline's own numbers exactly ----
        % The cluster-mass null SAVED is the one the pipeline recorded; the
        % local recomputation (max_cluster_mass, a second implementation of
        % the cluster finder) is only a cross-check and must equal it exactly.
        recM = mass_rec(k).value;
        okv = [isequaln(stat, q.T.real_T), isequaln(pm, q.p.real), ...
               isequaln(tfce, q.T.real_TFCE), isequaln(maxT, q.T.tMaxTFCE), ...
               cluster_p_match(q.Clusters, recM, nP), isequaln(maxM, recM)];
        L = struct();
        L.name = plog(k).name;  L.role = q.role;
        L.Y = q.Y;              L.rows = q.rows(:);
        L.X = q.D.X;            L.C = q.C;
        L.contrast_type = q.ctype;
        L.nuisance_idx  = logical(q.D.nuisance_idx(:)');
        L.nuisance_cols = find(q.D.nuisance_idx(:)');
        L.eb = q.D.eb(:);       L.perm_type = q.D.perm_type;
        L.effect_label = q.D.effect_label;
        L.evaluable = ev;
        L.stat = stat;          L.p_uncorrected = pm;   L.tfce = tfce;
        L.tfce_null_sorted = maxT;
        L.cluster_mass_null_sorted = recM;
        L.cluster_mass_null_source = 'recorded in snpm_perm_correction (tap glm:cluster_mass_null)';
        L.replay_ok = all(okv);
        L.replay_detail = double(okv);       % [stat, p, tfce, tfce-null, cluster-p, mass-null]
        G.loops(end+1) = L;
        assert(strcmp(plog(k).tag, ['glm_permute:' q.D.perm_type]), 'export_parity:schemeMismatch', ...
            'Loop %s recorded %s but its design says %s.', plog(k).name, plog(k).tag, q.D.perm_type);
    end
    G.loops = G.loops(:);
end

function plan = add_effect(plan, role, Y, rows, preset, meta, opts, mask, eff, with_posthoc)
% The omnibus loop of one effect, then (pipeline order) one loop per pairwise
% post-hoc when the effect is an F test with three or more levels.
%   Y, rows : the effect's data and the input-table row behind each row of Y
%   preset, meta, opts : what the effect's design is built from
%   mask    : channels the effect may test (fixed before any permutation)
% Each post-hoc is the two-level design snpm_glm_pair returns for that pair --
% the same call core_snpm_glm/compute_posthoc makes -- on the pair's rows, and
% is masked with the omnibus effect's evaluable set. A pair snpm_glm_pair
% declines is skipped here exactly as it is in the pipeline.
    D = snpm_glm_design(preset, meta, opts);
    fl = snpm_glm_fl_context(Y, D);
    evaluable = fl.evaluable & mask;
    plan(end+1) = struct('role', role, 'Y', Y, 'rows', rows, 'D', D, 'C', D.C, ...
        'ctype', D.contrast_type, 'mask', evaluable, 'T', eff.T, 'p', eff.p, 'Clusters', {eff.Clusters});
    if with_posthoc && strcmpi(D.contrast_type, 'F') && ~isempty(D.posthoc)
        ws = warning('off', 'snpm_glm_pair:skipped');          % the pipeline already said so
        warning('off', 'snpm_glm_pair:singleLevelSubjects');
        warning('off', 'snpm_glm_pair:fewSubjects');
        warning('off', 'snpm_glm_pair:smallGroup');
        restore = onCleanup(@() warning(ws));
        j = 0;
        for ph = 1:numel(D.posthoc)
            pr = snpm_glm_pair(preset, meta, opts, D.factor_col, D.posthoc(ph).a, D.posthoc(ph).b, ...
                D.posthoc(ph).label);
            if ~pr.ok, continue; end
            j = j + 1;
            assert(j <= numel(eff.posthoc) && strcmp(eff.posthoc(j).label, D.posthoc(ph).label), ...
                'export_parity:posthocOrder', 'Post-hoc %s is not where the pipeline put it.', D.posthoc(ph).label);
            plan(end+1) = struct('role', [role '/posthoc:' D.posthoc(ph).label], 'Y', Y(pr.rows, :), ...
                'rows', rows(pr.rows), 'D', pr.D, 'C', pr.C, 'ctype', 't', 'mask', evaluable, ...
                'T', eff.posthoc(j).T, 'p', eff.posthoc(j).p, 'Clusters', {eff.posthoc(j).Clusters}); %#ok<AGROW>
        end
        assert(j == numel(eff.posthoc), 'export_parity:posthocCount', ...
            'Rebuilt %d post-hoc loops for %s but the pipeline ran %d.', j, role, numel(eff.posthoc));
    end
end

function m = max_cluster_mass(p_model, s, isF, alpha, sparse_adj)
% snpm_perm_correction/find_pclusters: connected components of p < alpha,
% mass = sum(Wald)/n_channels, Wald = F or t^2. Returns the largest (0 if none).
    if isF, w = s; else, w = s.^2; end
    mask = double(p_model < alpha); mask(isnan(p_model)) = 0;
    cc = snpm_find_clusters_graphalgs(mask, 0.5, sparse_adj);
    m = 0;
    for c = 1:numel(cc)
        ch = cc{c};
        if isempty(ch), continue; end
        m = max(m, sum(w(ch), 'omitnan') / numel(ch));
    end
end

function ok = cluster_p_match(Clusters, maxM, nP)
    ok = true;
    for c = 1:numel(Clusters)
        if isempty(Clusters(c).channels), continue; end
        ok = ok && isequal(Clusters(c).p, (snpm_count_at_least(maxM, Clusters(c).mass) + 1) / (nP + 1));
    end
end

% =======================================================================
% PLAIN-TYPE CONVERSION (scipy.io.loadmat-safe)
% =======================================================================
function v = flat(v)
    if istable(v)
        s = struct();
        vn = v.Properties.VariableNames;
        for i = 1:numel(vn), s.(vn{i}) = flat(v.(vn{i})); end
        v = s;
    elseif isstring(v)
        if isscalar(v), v = char(v); else, v = cellstr(v); end
    elseif iscategorical(v) || isdatetime(v) || isduration(v)
        v = cellstr(string(v));
    elseif isa(v, 'function_handle')
        v = func2str(v);
    elseif iscell(v)
        for i = 1:numel(v), v{i} = flat(v{i}); end
    elseif isstruct(v)
        fn = fieldnames(v);
        for i = 1:numel(v)
            for j = 1:numel(fn), v(i).(fn{j}) = flat(v(i).(fn{j})); end
        end
    elseif issparse(v)
        v = full(v);
    elseif ~(isnumeric(v) || islogical(v) || ischar(v))
        v = ['<' class(v) '>'];
    end
end

function assert_plain(v, where)
    if iscell(v)
        for i = 1:numel(v), assert_plain(v{i}, where); end
    elseif isstruct(v)
        fn = fieldnames(v);
        for i = 1:numel(v)
            for j = 1:numel(fn), assert_plain(v(i).(fn{j}), [where '.' fn{j}]); end
        end
    elseif ~(isnumeric(v) || islogical(v) || ischar(v)) || issparse(v)
        error('export_parity:notPlain', '%s holds a %s, which scipy cannot read.', where, class(v));
    end
end

% =======================================================================
% README
% =======================================================================
function s = readme_row(c, fx, SEED, NPERM)
    s = sprintf('| `%s` | %s | %s | %d | %d | %s |', c.name, c.analysis, ...
        input_list(fx.source_files), SEED, NPERM, loop_list(fx.perm_log));
end

function s = input_list(S)
    fn = fieldnames(S); parts = cell(1, numel(fn));
    for i = 1:numel(fn), parts{i} = ['`' S.(fn{i}).relative_to_repo '`']; end
    s = strjoin(parts, '<br>');
end

function s = loop_list(L)
    parts = cell(1, numel(L));
    for i = 1:numel(L)
        parts{i} = sprintf('%d. %s [%s, %dx%d]', i, L(i).name, L(i).tag, L(i).n_draws, L(i).n_units);
    end
    s = strjoin(parts, '<br>');
end

function write_readme(path, cases, rows, SEED, NPERM, INPUT_SEED, B)
    fid = fopen(path, 'w');
    assert(fid > 0, 'Cannot write %s', path);
    c = onCleanup(@() fclose(fid));
    w = @(varargin) fprintf(fid, [varargin{1} '\n'], varargin{2:end});
    w('# Parity fixtures: MATLAB reference results for the Python port');
    w('');
    w('Written by `scripts/export_parity_fixtures.m` (MATLAB %s) on %s. Do not edit by hand;', ...
        version, char(datetime('now', 'Format', 'yyyy-MM-dd')));
    w('re-run the script. This folder is git-ignored.');
    w('');
    w('MATLAB is the reference. A Python analysis passes when, given the inputs in a fixture and');
    w('the permutation indices in `perm_log`, it reproduces `results`. Each fixture is one run of');
    w('`core_snpm_analysis(params)` after `rng(%d,''twister'')`, %d requested permutations.', SEED, NPERM);
    w('Inputs this script builds (under `inputs/`) come from a private stream seeded %d.', INPUT_SEED);
    w('');
    w('Load with `scipy.io.loadmat(path, simplify_cells=True)`. MATLAB struct -> dict, struct array ->');
    w('list of dicts (a 1-element struct array collapses to ONE dict, not a list of one), cellstr ->');
    w('array of str (a 1-element cellstr collapses to ONE str), logical -> uint8 0/1, 1xN and Nx1 numeric');
    w('-> 1-D array, 1x1 -> scalar, empty -> empty array. Whole-number arrays (indices, counts, 0/1 designs)');
    w('load as uint8/uint16 because MATLAB v7 stores them compactly: cast to int64 or float64 before doing');
    w('arithmetic on them (uint8 - 1 wraps, uint8 products overflow). Because `simplify_cells` drops singleton dimensions, the orientation given below is');
    w('the MATLAB one; a 1-D array has lost it.');
    w('');
    w('## Conventions');
    w('');
    w('- Data matrices are **subjects (rows) x channels (columns)**; GLM matrices are observations x channels.');
    w('- Every index is **1-based** (channel indices, row indices, permutation indices, neighbours). Subtract 1 in Python.');
    w('- Channel indices in `results` (`uncorrsigch`, `Clusters.channels`, `excluded_channels.index`, ...) index the');
    w('  **analysed** channel order `montage.labels`, not the columns of the input file.');
    w('- `montage.data_columns` maps analysed channel k to column `data_columns(k)` of `inputs.*.values`.');
    w('- A channel that was not tested is `NaN` in the statistic and p maps, never 0 or 1.');
    w('- Corrected p-values are `(b + 1) / (N + 1)`, b = number of null maxima >= the observed value, N = the number of');
    w('  permutations actually run (`n_draws`), which is smaller than requested on an exact-enumeration branch.');
    w('- TFCE: E = 0.5, H = 2, step dh = 0.1 (`ClusterEnhancement`); alpha = 0.05.');
    w('- Text read from CSV by Python must be parsed round-trip exact (`pandas.read_csv(..., float_precision="round_trip")`)');
    w('  to equal `inputs.*.values` bit for bit.');
    w('');
    w('## Fixtures');
    w('');
    w('| Fixture (`.mat`) | Analysis | Source files | Seed | Requested permutations | Permutation loops, in call order [tag, draws x units] |');
    w('|---|---|---|---|---|---|');
    for i = 1:numel(rows), fprintf(fid, '%s\n', rows{i}); end
    w('');
    w('What each case is for:');
    w('');
    for i = 1:numel(cases), w('- `%s` (%s tier): %s', cases(i).name, cases(i).tier, cases(i).note); end
    w('');
    w('`example_data/` ships no covariate file, so the covariate cases use `inputs/corr_covariates.csv`');
    w('(Subject, age, bmi), built by this script. `inputs/corr_covariates_missing.csv` is the same file with');
    w('subject 7''s age blank. The Compumedics inputs are synthetic noise with an effect planted on one channel');
    w('and its first- and second-ring neighbours (%d channels: %s). The missing-data inputs are copies of the', ...
        numel(B.compu_cluster), strjoin(B.compu_cluster, ', '));
    w('example files with cells set to NaN. In `miss_unpaired_group*.csv`: `%s` is all-NaN in group B only,', B.miss_unp_channels{1});
    w('`%s` keeps 3 subjects in A and 4 in B (both fall below the 18 usable subjects the unpaired test needs', B.miss_unp_channels{2});
    w('and are excluded), and `%s` has 2 NaN in group A and stays in the map on 30 subjects.', B.miss_unp_channels{3});
    w('The tie inputs (`ties_corr_eeg.csv`, `small_corr_eeg_n7_ties.csv`) are `corr_eeg.csv` with subject 2 given');
    w('subject 1''s value in %s, so each of those channels has exactly one tied pair in the EEG matrix; the', strjoin(B.tie_channels, ', '));
    w('behaviour matrix has no ties. `results.p.real` for correlationS is MATLAB `corr(..., ''Type'', ''Spearman'')`,');
    w('whose p-value method depends on n and on ties; the statistic `T.real_T` does not.');
    w('');
    w('## Order of permutation loops per analysis');
    w('');
    w('The generator is seeded once per fixture and never re-seeded, so the loops below consume one random');
    w('stream in this order. Reproducing the draws needs this order; reproducing the results needs only');
    w('`perm_log`, which already holds every index.');
    w('');
    w('**Legacy tier** (pairedT, unpairedT, onesampleT, correlationP, correlationS): two loops, always in this order.');
    w('');
    w('1. `legacy/TFCE` -- `snpm_single_threshold_with_TFCE` asks its generator for one whole matrix.');
    w('2. `legacy/cluster` -- `snpm_cluster_analysis` asks the same generator again.');
    w('');
    w('On a random branch these are two DIFFERENT matrices (the cluster test does not reuse the TFCE test''s');
    w('permutations). On an exact branch both are the same full enumeration. The global test draws nothing.');
    w('The tap also records two non-draw entries per legacy run (`legacy:engine_inputs`, `legacy:cluster_null`); they are');
    w('stored under `engine`, not in `perm_log`.');
    w('Covariate residualization and the missing-data masks happen before either loop and draw nothing.');
    w('');
    w('| Analysis | Generator | `indices` (n_draws x n_units) | How a row is applied |');
    w('|---|---|---|---|');
    w('| pairedT, onesampleT | `snpm_signflip_patterns` | 0/1, units = subjects | 1 = pair kept, 0 = the subject''s two conditions are swapped (difference sign-flipped). The engine builds `x_perm = [x(keep,:); y(~keep,:)]`, `y_perm = [y(keep,:); x(~keep,:)]`, i.e. rows are also regrouped, kept subjects first. |');
    w('| unpairedT | `snpm_relabel_assignments` | subject indices into the stacked matrix, units = n1 + n2 | first nGrp columns = group 1, the rest = group 2. Stacking: tail both/right `data = [data_x; data_y]`, nGrp = n1. Tail left `data = [data_y; data_x]`, nGrp = n2, tested right. |');
    w('| correlationP, correlationS | `snpm_row_permutations` | row orders, units = subjects | `r = corr(data_x(order,:), data_y)` column by column; data_y is not moved. Spearman ranks both matrices once, before the loop. |');
    w('');
    w('Random branches: `signflip:random` = `rand(N, n) > 0.5` (one call, filled column by column); `relabel:random`');
    w('and `rowperm:random` = one `randperm(n)` per row for rows 2..N.');
    w('Row 1 of every random matrix is the observed labelling (all ones / `1:n`). Exact branches, and the row order');
    w('they pin:');
    w('');
    w('- `signflip:exact` when `2^n <= requested`: row k+1 = bits of k, least-significant bit = subject 1, k = 0..2^n-1.');
    w('  Row 1 is therefore all zeros (every pair flipped) and the observed labelling is the LAST row.');
    w('- `relabel:exact` when `nchoosek(n1+n2, nGrp) <= requested`: `enumerate_combinations`, lexicographic in the');
    w('  group-1 index set, row 1 = observed split; within each row group 1 then group 2 are each ascending.');
    w('- `rowperm:exact` when `n <= 10` and `n! <= requested`: MATLAB `perms(1:n)` (reverse lexicographic) with the');
    w('  identity row swapped into row 1 (and the old row 1 moved to where the identity was, the last row).');
    w('');
    w('For pairedT tail left the engines swap data_x and data_y before the loop and run a right-tailed test.');
    w('');
    w('**GLM tier** (anova1, ancova, regression, rmanova, mixed2way): one row index `pe` per permutation, drawn by');
    w('`snpm_glm_permute` inside `snpm_perm_correction`; TFCE and cluster-mass nulls share the same draws.');
    w('`perm_log(k).indices` row i is `pe` for permutation i: `Y_perm = Zfit + R(pe,:)` (Freedman-Lane), where');
    w('`Z = X(:, nuisance_idx)`, `Zfit = Z*(Z\\Y)`, `R = Y - Zfit`. The observed labelling is NOT forced into row 1 here.');
    w('`glm_permute:free` = `randperm(nObs)`; `glm_permute:within` = for each block of `eb`, in order of first');
    w('appearance, `randperm(block size)` applied to that block''s rows (blocks of one row draw nothing).');
    w('');
    w('| Preset | Loops, in order |');
    w('|---|---|');
    w('| regression, ancova with 2 groups | 1 loop: `glm:main:<effect>` (t) |');
    w('| anova1, ancova with 3+ groups, rmanova | `glm:main:<effect>` (F), then one `glm:posthoc:<A vs B>` (t) per pair, in design order. Each post-hoc is its OWN two-level test on the rows of A and B only (see below) and draws its own fresh indices over those rows. |');
    w('| mixed2way | `mixed2way:group/glm:effect:...` on SUBJECT MEANS (anova1 design, free; post-hoc only with 3+ groups) -> `mixed2way:condition/glm:effect:...` (within; post-hoc when 3+ conditions) -> `mixed2way:interaction/glm:effect:...` (within; never post-hoc) -> for each group in order of first appearance, `mixed2way:simple:<group>/glm:effect:...` (rmanova design on that group''s rows, within) followed by its post-hoc loops when 3+ conditions. The global tests draw nothing. |');
    w('');
    w('**Method of 2026-10-14** (`glm.scheme = ''2026-10-14''` marks a fixture that carries it).');
    w('');
    w('- **Pairwise post-hoc.** "A vs B" is tested on the rows whose level is A or B, in table order, with the preset''s design');
    w('  rebuilt on those rows: anova1 `[1, D]`; ancova `[1, D, covariates]` with covariates re-centred on the subset, categorical');
    w('  covariates coded from the levels present, and covariates constant on the subset dropped; rmanova and simple effects');
    w('  `[1, Ds, Dc]` with dummies for the subjects present; mixed2way condition post-hocs the condition model below on the');
    w('  two-condition rows. The contrast is the single two-level dummy, signed so the statistic is mean(A) - mean(B) (the');
    w('  builder''s reference level is the first to appear in the subset). Nuisance = every other column; `free` over the');
    w('  subset rows for between-subject pairs, `within` subject blocks otherwise. So `perm_log(k).indices` of a post-hoc loop');
    w('  has the SUBSET''s row count and indexes rows of that loop''s `Y`, not rows of the input table. Within-subject pairs');
    w('  drop a subject that has only one of the two levels. A pair that cannot be estimated is skipped and has no loop.');
    w('  Channels: the omnibus effect''s evaluable set, and finite in every row of the pair.');
    w('  `results.posthoc(j)` (and each effect''s `.posthoc(j)`) records what the pair ran on: `.n_rows`, `.n_subjects`,');
    w('  `.n_subjects_total`, `.group_sizes` (between-subject pairs), `.dropped_subjects`, `.run_on` (the same in words).');
    w('  Pairs that were not run are listed in `results.posthoc_skipped` (`.label`, `.reason`), per effect for mixed2way.');
    w('  Level order: first appearance in the table. The simple effects inside each group use the WHOLE table''s condition');
    w('  order, so every part of a mixed2way report lists and signs its pairs identically.');
    w('- **mixed2way balance and one mask.** Every subject must have exactly one row per condition and one group (otherwise the');
    w('  run errors). One channel set, the channels finite in every row of the long table, is used for all effects, all');
    w('  post-hocs and the three whole-head averages. Subject means for the group effect are plain means.');
    w('- **mixed2way condition main effect.** `X = [1, Ds, Dc, Dxc]`: `Ds` subject dummies and `Dc` condition dummies, both');
    w('  reference-coded on the first level to appear; `Dg` reference-coded group dummies; `E = Dg` with the first group''s rows');
    w('  set to -1 in every column (sum-to-zero); `Dxc = [E(:,1).*Dc, E(:,2).*Dc, ...]`, group outer, condition inner.');
    w('  Contrast = the `Dc` columns, nuisance = everything else, within-subject permutation. The interaction model is');
    w('  unchanged: `[1, Ds, Dc, Dx]` with `Dx` from the reference-coded `Dg`, contrast = the `Dx` columns.');
    w('');
    w('`glm.loops(k)` describes the same loop as `perm_log(k)`.');
    w('');
    w('## Field dictionary');
    w('');
    w('Top level of every fixture:');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `case_name` | char | fixture name (= file name without `.mat`) |');
    w('| `tier` | char | `legacy` or `glm` |');
    w('| `description` | char | what the case is for |');
    w('| `seed`, `rng_algorithm` | scalar, char | `rng(seed, rng_algorithm)` was called immediately before the run |');
    w('| `permutations` | scalar | permutations REQUESTED (`params.permutations`) |');
    w('| `matlab_version`, `created` | char | provenance |');
    w('| `params` | struct | the params struct passed to `core_snpm_analysis` (paths are absolute, this machine) |');
    w('| `source_files` | struct | per input (`data1_file`, `data2_file`, `data_file`, `covariate_file`): `.absolute` and `.relative_to_repo` |');
    w('| `inputs` | struct | the input files as `readtable` returned them, before any processing (below) |');
    w('| `montage` | struct | the channels the run analysed and their adjacency (below) |');
    w('| `perm_log` | struct array, K x 1 | one entry per permutation loop, in call order (below) |');
    w('| `results` | struct | the complete `results_struct` returned by the pipeline, plain types (below) |');
    w('| `results_text` | cellstr | the text summary the pipeline returned |');
    w('| `engine` | struct | legacy tier only: the matrices the engines were handed, the cluster-size null, the global channel set (below) |');
    w('| `glm` | struct | GLM tier only: designs and replayed nulls (below) |');
    w('');
    w('`inputs` -- legacy tier has `data1`, `data2` (empty for `onesampleT_single`), `covariates` (empty when unused);');
    w('GLM tier has `data`. Each data block:');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `file` | char | absolute path that was read |');
    w('| `all_columns` | cellstr 1 x nCols | header, file order |');
    w('| `column_labels` | cellstr 1 x nData | the data columns (header minus metadata), file order. For `corr_behavior_single.csv` this is the single column `MRRT`. |');
    w('| `values` | double nRows x nData | the data, rows in file order, NaN = missing. NOT yet channel-selected. |');
    w('| `meta_names` | cellstr | metadata columns: legacy = the first column whose name contains `subject`; GLM = `params.meta_cols` |');
    w('| `meta.<name>` | double nRows x 1 or cellstr nRows x 1 | each metadata column (group, condition, age, sex, cognition, ...) |');
    w('| `subject_ids` | cellstr nRows x 1 | the subject column again, for convenience |');
    w('');
    w('`inputs.covariates`: `file`, `subject_ids` (cellstr n x 1, first column of the file), `names` (cellstr), `values` (n x nCov, NaN = missing).');
    w('');
    w('`montage`:');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `system_key`, `display` | char | recording system (`egi` / `compu`) |');
    w('| `labels` | cellstr 1 x nCh | analysed channels, in analysed order (= `results.chanlocs(k).labels`). This is the order of every per-channel map. It is the data file''s column order restricted to the system''s scalp labels. |');
    w('| `data_columns` | 1 x nCh, 1-based | column of `inputs.data1.values` (or `inputs.data.values`) holding analysed channel k |');
    w('| `neighbors` | nCh x maxNb, 1-based, NaN padded | row k lists the analysed-order indices adjacent to channel k |');
    w('| `xyz` | nCh x 3 | channel coordinates from the chanlocs asset |');
    w('| `tfce_E`, `tfce_H`, `tfce_dh`, `alpha` | scalars | 0.5, 2, 0.1, 0.05 |');
    w('');
    w('`perm_log(k)`:');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `name` | char | `legacy/TFCE`, `legacy/cluster`, or the GLM loop name built from the markers |');
    w('| `tag` | char | generator and branch: `signflip:exact|random`, `relabel:exact|random`, `rowperm:exact|random`, `glm_permute:free|within` |');
    w('| `kind` | char | `signflip`, `relabel`, `rowperm` or `row_index` |');
    w('| `branch` | char | `exact`, `random`, `free` or `within` |');
    w('| `site`, `caller` | char | `function:line` that recorded the draw, and its caller |');
    w('| `n_draws`, `n_units` | scalars | size of `indices` |');
    w('| `indices` | double n_draws x n_units | row i = permutation i, in the order the loop used them. NOT sorted. See the tables above for how a row is applied. |');
    w('');
    w('`engine` (legacy tier only). `data_x` / `data_y` were recorded in `core_snpm_analysis` immediately before the');
    w('`global_stat_test` call; the same two variables, unmodified, are then passed to both engines. Subject matching,');
    w('channel selection, the datatype transform, the missing-data masks and covariate residualization are all upstream');
    w('of that point, so this is the stage to compare against before any statistic is computed.');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `data_x`, `data_y` | nSubj_x x nCh, nSubj_y x nCh | engine inputs, analysed channel order, rows = analysed subjects in analysed order. A masked channel is an all-NaN column in BOTH; other NaN are original missing cells (t-tests only). With covariates these are residuals (of ranks for correlationS). One-sample against 0: `data_y` is all zeros. Single behavioural column: `data_y` is that column repeated across channels. |');
    w('| `comparison` | char | the key the engines were given. Equals `params.comparison` except partial Spearman, which is handed on as `correlationP` (ranking and residualizing are already done). |');
    w('| `tail` | char | `both`, `left` or `right` |');
    w('| `ncov` | scalar | number of covariates (reduces the df of the global correlation p only) |');
    w('| `alpha` | scalar | 0.05 |');
    w('| `cluster_threshold` | scalar | cluster-forming threshold passed to `snpm_cluster_analysis` (t: `tinv(1-alpha/2, df)`; r: step function of n) |');
    w('| `cluster_thresholds` | 1 x 2 or scalar | thresholds actually applied: `[-thr, +thr]` for tail both (clusters of stat < -thr and of stat > +thr are formed separately), `thr` for a one-sided tail |');
    w('| `cluster_null_sorted` | n_draws x 1, **sorted descending** | the vector every `Clusters(c).p` is computed against: per permutation the largest cluster SIZE (channel count; for tail both the larger of the negative-side and positive-side maxima), 0 when the permutation has no cluster. Uses the `legacy/cluster` permutations, `perm_log(2)`. |');
    w('| `global_channels_kept` | 1 x K, 1-based, analysed order | channels the global test averaged over (finite in every analysed row of both matrices); K = `results.global_n_channels` |');
    w('| `global_stat`, `global_pval` | scalars | `global_stat_test` re-run on `data_x`/`data_y`; asserted equal to `results.global_stat`/`global_pval` |');
    w('');
    w('`results`, legacy tier (fields absent in MATLAB are absent here):');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `T.real_T` | 1 x nCh | observed statistic map: t, or r for correlations (partial r after residualization) |');
    w('| `T.real_TFCE` | 1 x nCh | TFCE-enhanced observed map (signed) |');
    w('| `T.tMax` | n_draws x 1, **sorted descending** | max-statistic null (max over channels of |stat|, or of stat for a one-sided tail) |');
    w('| `T.tMaxTFCE` | n_draws x 1, **sorted descending** | max-TFCE null, same rule |');
    w('| `T.chk_T` | 1 x nCh | mean over permutations of |stat| (stat for one-sided) |');
    w('| `T.critical_T`, `T.critical_T_TFCE` | scalars | element `floor(alpha*n_draws)+1` of the sorted nulls |');
    w('| `p.real` | 1 x nCh | uncorrected parametric p (`ttest`/`ttest2`/`corr`) |');
    w('| `p.corrected` | 1 x nCh | max-statistic corrected p |');
    w('| `p.correctedTFCE` | 1 x nCh | TFCE corrected p |');
    w('| `Clusters` | struct array | cluster-EXTENT test. Tail both: clusters of stat < -threshold are listed first, then stat > +threshold; a threshold with no cluster contributes one entry with empty `channels` and `p` = 1. |');
    w('| `Clusters(c).channels` | 1 x k, 1-based | members, analysed order |');
    w('| `Clusters(c).threshold` | scalar | signed cluster-forming threshold (t: `tinv(1-alpha/2, df)`; r: 0.50/0.40/0.35/... by n) |');
    w('| `Clusters(c).p` | scalar | `(#{null max cluster SIZE >= size} + 1)/(N + 1)`; the cluster statistic is the number of channels |');
    w('| `Clusters(c).permutations` | scalar | N |');
    w('| `uncorrsigch`, `correctTFCEsigch`, `SnPMsigch` | 1-based index vectors | channels with p.real <= alpha, p.correctedTFCE <= alpha, and members of clusters with p <= alpha |');
    w('| `alpha` | scalar | 0.05 |');
    w('| `chanlocs` | struct array 1 x nCh | EEGLAB channel locations, analysed order |');
    w('| `per_channel_n` | 1 x nCh | units the statistic was computed on (0 on an excluded channel) |');
    w('| `n_available` | 1 x nCh | units available before masking (pairs; pooled subjects for unpairedT) |');
    w('| `n_analysed_constant` | scalar or empty | correlation: the one n used on every retained channel |');
    w('| `excluded_channels` | struct | `.n`, `.n_channels`, `.index` (1-based, analysed order), `.labels`, `.context`, `.reason`; absent or `.n` = 0 when nothing was excluded |');
    w('| `exclusion_spatial` | struct | where the exclusions fell (present only with exclusions) |');
    w('| `global_stat`, `global_pval` | scalars | omnibus test on the channel average |');
    w('| `global_n_channels`, `global_n_channels_total` | scalars | channels the average was taken over, and the montage size |');
    w('| `subject_matching`, `covariates` | structs | correlation only: matched ids; covariate names and the matched covariate matrix |');
    w('| `behavior_col`, `data2_matched`, `hide_condition_b`, `is_source`, `data_summary`, `descriptive` | mixed | reporting metadata |');
    w('');
    w('`results`, GLM tier: `T.real_T` (t or F), `T.real_TFCE`, `T.tMaxTFCE` (**sorted descending**, max of |TFCE| for t,');
    w('of TFCE for F), `T.critical_T_TFCE`, `T.excluded_channels` (1-based), `p.real`, `p.correctedTFCE`, `Clusters(c)`');
    w('with `.channels`, `.mass` (= sum(Wald)/n_channels, Wald = F or t^2; cluster-MASS, unlike the legacy tier),');
    w('`.p`, `.threshold` (= alpha: clusters are formed from p.real < alpha), `.permutations`; `posthoc(j)` with `.label`,');
    w('`.T`, `.p`, `.Clusters`, `.uncorrsigch`, `.correctTFCEsigch`, `.SnPMsigch`; `glm.preset`, `.effect_label`,');
    w('`.contrast_type` (`t`/`F`), `.perm_type`, `.n_obs`, `.n_subjects`, `.group_sizes`, `.group_size_labels`.');
    w('GLM t contrasts are always two-sided. For mixed2way, `effects(1..3)` = group, condition, interaction, each with');
    w('the fields above plus `.key`, `.df1`, `.df2`, `.global_stat`, `.global_pval`, `.global_df1`, `.global_df2`,');
    w('`.global_n_channels`, `.excluded_channels`; `effects(3).simple_effects(g)` has the same layout (`.key` =');
    w('`condition@<group>`); the top-level `T`/`p`/`Clusters` repeat the interaction.');
    w('');
    w('`glm` (GLM tier only; recomputed by the export script from the recorded indices and checked against `results`):');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `preset`, `alpha`, `tfce_E`, `tfce_H` | char, scalars | |');
    w('| `scheme` | char | `2026-10-14`: the fixture was produced with the method described under "Method of 2026-10-14" above. A GLM fixture without this field predates it and must not be used. |');
    w('| `loops(k).name`, `.role` | char | same loop as `perm_log(k)`; role = `main`, `group`, `condition`, `interaction`, `condition@<group>`, with `/posthoc:<label>` appended for post-hoc loops |');
    w('| `loops(k).Y` | nObs x nCh | the data matrix that loop tested, channel-selected (analysed order): the rows of the two levels for a post-hoc loop, subject means for the mixed2way group effect. `perm_log(k).indices` index ITS rows. |');
    w('| `loops(k).rows` | nObs x 1, 1-based | rows of `inputs.data.values` behind each row of `Y` (for subject means: each subject''s first row; for a post-hoc: the rows of its two levels, in table order) |');
    w('| `loops(k).X` | nObs x p | design matrix |');
    w('| `loops(k).C` | q x p | contrast (q = 1 -> t, q > 1 -> F). A post-hoc contrast is +1 or -1 on the two-level dummy of ITS design, so that the statistic is mean(A) - mean(B) |');
    w('| `loops(k).contrast_type` | char | `t` or `F` |');
    w('| `loops(k).nuisance_idx` / `.nuisance_cols` | logical 1 x p / 1-based indices | columns of X held fixed under Freedman-Lane |');
    w('| `loops(k).eb` | nObs x 1 | exchangeability block id (all ones for `free`) |');
    w('| `loops(k).perm_type` | char | `free` or `within` |');
    w('| `loops(k).evaluable` | logical 1 x nCh | channels this loop tested; the rest are NaN everywhere. For a post-hoc: the omnibus effect''s set AND finite in every row of the pair. For mixed2way: the one shared mask (finite in every row of the long table), identical in every loop. |');
    w('| `loops(k).stat`, `.p_uncorrected`, `.tfce` | 1 x nCh | observed statistic, parametric p, TFCE map (= `results` values) |');
    w('| `loops(k).tfce_null_sorted` | n_draws x 1, **sorted descending** | max-TFCE null (= `T.tMaxTFCE`) |');
    w('| `loops(k).cluster_mass_null_sorted` | n_draws x 1, **sorted descending** | max cluster-mass null (0 when a permutation has no cluster), **recorded inside `snpm_perm_correction`** after its sort: the vector every `Clusters(c).p` of that loop is computed against |');
    w('| `loops(k).cluster_mass_null_source` | char | provenance of the field above |');
    w('| `loops(k).replay_ok`, `.replay_detail` | logical, 1 x 6 | the export script''s replay reproduced [stat, p, tfce, tfce-null, cluster p from the recorded mass null, the recorded mass null itself] exactly |');
    w('');
    w('## Not in these fixtures');
    w('');
    w('- **A global test for single-effect GLM presets.** The pipeline computes none (mixed2way has one per effect).');
    w('- Circular and mixed-model analyses (later phases).');
end
