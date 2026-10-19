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
% Phase 2 = the circular tier (circ_phase_group, circ_phase_group_u2,
% circ_corrAngLinear), replayed like the GLM tier (fx.circ.replay_ok). The
% mixed-model tier is recorded by the same hook but is not exported here.
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
    % Four kinds of entry are not permutation draws: the matrices the legacy
    % engines were handed and the legacy cluster-size null (-> fx.engine), the
    % cluster-mass null of each GLM / circular loop (-> fx.glm.loops, fx.circ),
    % and what the circular engine was handed (-> fx.circ).
    is_aux  = ismember({raw_log.tag}, {'legacy:engine_inputs', 'legacy:cluster_null', ...
        'glm:cluster_mass_null', 'circ:engine_inputs'});
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
    elseif strcmp(c.tier, 'circ')
        fx.inputs = circ_inputs(P);
        [fx.montage, fx.circ] = circ_capture(P, results_struct, fx.perm_log, aux_log, c.name);
        replay_txt = sprintf('replay %d/%d checks ok', sum(fx.circ.replay_detail), ...
            numel(fx.circ.replay_detail));
        if ~fx.circ.replay_ok
            error('export_parity:replayMismatch', ...
                ['%s: the recorded indices did not reproduce the pipeline''s own null ' ...
                 '(%s; failed: %s). The fixture would not be a valid reference; not saved.'], ...
                c.name, replay_txt, strjoin(fx.circ.replay_checks(~fx.circ.replay_detail), ', '));
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

    % ---- 7. correlation with covariates: backlog of 2026-10-18 ----
    cases(end+1) = legacy('miss_correlationP_cov_eeg_nan', 'correlationP', 'both', B.miss_corrX_cov, ex('corr_behavior.csv'), 'egi', B.cov, ...
        ['Partial Pearson (Freedman-Lane) with missing EEG cells, not missing covariates: the complete-column ' ...
         'rule masks those channels and every subject is kept.']);
    cases(end+1) = legacy('exact_correlationP_n5_cov', 'correlationP', 'both', B.small_corrX, B.small_corrY, 'egi', B.small_cov, ...
        ['First 5 subjects, partial Pearson on ONE covariate (age): 5! = 120 orderings enumerated under ' ...
         'Freedman-Lane; df = 5 - 2 - 1 = 2.']);

    % ---- 8. circular tier (phase 2) ----
    cases(end+1) = circ('circ_phase_group', 'circ_phase_group', ex('circ_phase_groupA.csv'), ex('circ_phase_groupB.csv'), ...
        ex('circ_phase_countsA.csv'), ex('circ_phase_countsB.csv'), 'rad', 'literature_uppeak0', ...
        ['Hotelling T^2 on (cos, sin), 16 + 16 subjects, per-channel log(event count) covariate, ' ...
         'Freedman-Lane on the stacked response. Planted 46 deg mean-direction shift.']);
    cases(end+1) = circ('circ_phase_group_u2', 'circ_phase_group_u2', ex('circ_phase_groupA.csv'), ex('circ_phase_groupB.csv'), ...
        ex('circ_phase_countsA.csv'), ex('circ_phase_countsB.csv'), 'rad', 'literature_uppeak0', ...
        'Watson U^2 on the same pair (no covariate slot; the count files still feed the gate). TFCE dh = 0.005.');
    cases(end+1) = circ('circ_corrAngLinear', 'circ_corrAngLinear', ex('circ_anglinear_angles.csv'), ex('circ_anglinear_measure.csv'), ...
        '', '', 'rad', 'literature_uppeak0', ...
        'Circular-linear F (y on [1, cos, sin]), 30 matched subjects, the measure y permuted (free).');
    cases(end).params.measure_col = 'measure';
    cases(end+1) = circ('circ_phase_group_gate', 'circ_phase_group', B.circ_gateA, B.circ_gateB, ...
        B.circ_gate_countsA, B.circ_gate_countsB, 'rad', 'literature_uppeak0', ...
        ['Hotelling on a copy of the example pair edited so that every exclusion route fires: one channel ' ...
         'below 8 per group, one incomplete channel under the covariate, one channel with log(count) ' ...
         'disjoint by group (collinear), and one constant-count channel with a missing angle and a zero ' ...
         'count that stays in on the covariate-free (pairwise-deletion) path.']);
    cases(end+1) = circ('circ_phase_group_u2_gate', 'circ_phase_group_u2', B.circ_gateA, B.circ_gateB, ...
        B.circ_gate_countsA, B.circ_gate_countsB, 'rad', 'literature_uppeak0', ...
        'Watson U^2 on the gate inputs: same evaluable set, missing angles deleted per channel.');
    cases(end+1) = circ('circ_phase_group_deg_luna', 'circ_phase_group', B.circ_lunaA, B.circ_lunaB, ...
        ex('circ_phase_countsA.csv'), ex('circ_phase_countsB.csv'), 'deg', 'luna_zerocross0', ...
        ['The example pair exported in DEGREES with the Luna convention (up-state peak at 270 deg): ' ...
         'import converts units and rotates by +90 deg.']);
    cases(end+1) = circ('circ_phase_group_u2_turtlewave', 'circ_phase_group_u2', B.circ_twA, B.circ_twB, ...
        ex('circ_phase_countsA.csv'), ex('circ_phase_countsB.csv'), 'rad', 'turtlewave_pre_v4', ...
        ['The example pair inverted by 180 deg (TurtleWave pre-v4.0) and declared so: import rotates by ' ...
         '180 deg and the inversion check passes. circ.inversion_probe records that the same files ' ...
         'declared literature_uppeak0 are refused.']);
    cases(end+1) = circ('circ_corrAngLinear_miss_y', 'circ_corrAngLinear', B.circ_missy_angles, B.circ_missy_measure, ...
        '', '', 'rad', 'literature_uppeak0', ...
        ['Circular-linear F with 6 of 30 measure values missing (one Inf) and 3 missing angles in each of ' ...
         '10 channels: the 6 subjects are dropped whole before anything is computed (results.circ.' ...
         'dropped_subjects), the missing angles are deleted per channel on the fixed 24-subject row set.']);
    cases(end).params.measure_col = 'measure';
end

function c = circ(name, comparison, f1, f2, k1, k2, units, convention, note)
    P = struct('comparison', comparison, 'tail', 'both', ...
        'data1_file', f1, 'data1_sheet', 'CSV File', ...
        'data2_file', f2, 'data2_sheet', 'CSV File', ...
        'channels', 'egi', 'datatype', 'absolute', ...
        'circ_units', units, 'circ_convention', convention);
    if ~isempty(k1)
        P.count1_file = k1; P.count1_sheet = 'CSV File';
        P.count2_file = k2; P.count2_sheet = 'CSV File';
    end
    c = struct('name', name, 'tier', 'circ', 'analysis', comparison, 'params', P, 'note', note);
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

    % ---- inputs added with the circular tier: a SECOND private stream, so
    %      nothing above draws a different number because these exist ----
    rs2 = RandStream('mt19937ar', 'Seed', seed + 1);

    % correlation with covariates, missing EEG cells (not covariates)
    nan_ch = scalp(randperm(rs2, numel(scalp), 2));
    Xm = X;
    Xm.(nan_ch{1})(randi(rs2, n)) = NaN;
    Xm.(nan_ch{2})(randi(rs2, n)) = NaN;
    B.miss_corrX_cov = in('miss_corr_eeg_cov.csv');       writetable(Xm, B.miss_corrX_cov);
    B.miss_cov_channels = nan_ch;
    % exact enumeration with one covariate: the first 5 subjects, age only
    B.small_cov = in('small_corr_covariates_n5.csv');     writetable(C(1:5, {'Subject', 'age'}), B.small_cov);

    % circular: the example pair edited so that every exclusion route fires
    Ga = rd('circ_phase_groupA.csv'); Gb = rd('circ_phase_groupB.csv');
    Ka = rd('circ_phase_countsA.csv'); Kb = rd('circ_phase_countsB.csv');
    cscalp = analysed_labels('egi', Ga.Properties.VariableNames);
    pk = cscalp(randperm(rs2, numel(cscalp), 4));
    na = height(Ga); nb = height(Gb);
    Ea = Ga; Eb = Gb; Qa = Ka; Qb = Kb;
    Ea.(pk{1})(randperm(rs2, na, na - 7)) = NaN;          % (1) 7 angles left in A  -> below the gate
    Qa.(pk{1})(:) = 20; Qb.(pk{1})(:) = 20;               %     (constant count, so the gate is its ONLY reason)
    Eb.(pk{2})(randi(rs2, nb)) = NaN;                     % (2) one missing angle under the covariate
    Qa.(pk{3}) = max(Kb.(pk{3})) + (1:na)';               % (3) every A count above every B count
    Qa.(pk{4})(:) = 20; Qb.(pk{4})(:) = 20;               % (4) constant count -> covariate dropped ...
    Qa.(pk{4})(3) = 0;                                    %     ... A row 3 has no usable count (angle dropped)
    Eb.(pk{4})(5) = NaN;                                  %     ... B row 5 has no angle; the channel stays in
    B.circ_gateA = in('circ_gate_groupA.csv');          writetable(Ea, B.circ_gateA);
    B.circ_gateB = in('circ_gate_groupB.csv');          writetable(Eb, B.circ_gateB);
    B.circ_gate_countsA = in('circ_gate_countsA.csv');  writetable(Qa, B.circ_gate_countsA);
    B.circ_gate_countsB = in('circ_gate_countsB.csv');  writetable(Qb, B.circ_gate_countsB);
    B.circ_gate_channels = pk;

    % circular: the example pair in other export conventions
    ach = setdiff(Ga.Properties.VariableNames, {'Subject'}, 'stable');
    La = Ga; Lb = Gb; Wa = Ga; Wb = Gb;
    La{:, ach} = mod(Ga{:, ach} * 180 / pi - 90, 360);   % degrees, up-state peak at 270
    Lb{:, ach} = mod(Gb{:, ach} * 180 / pi - 90, 360);
    Wa{:, ach} = mod(Ga{:, ach} + pi, 2 * pi);           % radians, inverted
    Wb{:, ach} = mod(Gb{:, ach} + pi, 2 * pi);
    B.circ_lunaA = in('circ_luna_deg_groupA.csv');      writetable(La, B.circ_lunaA);
    B.circ_lunaB = in('circ_luna_deg_groupB.csv');      writetable(Lb, B.circ_lunaB);
    B.circ_twA   = in('circ_tw_inverted_groupA.csv');   writetable(Wa, B.circ_twA);
    B.circ_twB   = in('circ_tw_inverted_groupB.csv');   writetable(Wb, B.circ_twB);

    % ---- circular-linear with a missing measure: a THIRD private stream, so
    %      every input above keeps the values it had before this was added ----
    rs3 = RandStream('mt19937ar', 'Seed', seed + 2);
    Ya = rd('circ_anglinear_angles.csv'); Ym = rd('circ_anglinear_measure.csv');
    ny = height(Ym);
    ymiss = sort(randperm(rs3, ny, round(0.2 * ny)));     % 20% of subjects lose y ...
    Ym.measure(ymiss) = NaN;
    Ym.measure(ymiss(end)) = Inf;                         % ... one of them as Inf (also missing)
    lscalp = analysed_labels('egi', Ya.Properties.VariableNames);
    ach_nan = lscalp(randperm(rs3, numel(lscalp), 10));   % 10 channels, 3 missing angles each
    for i = 1:numel(ach_nan)
        Ya.(ach_nan{i})(randperm(rs3, height(Ya), 3)) = NaN;
    end
    B.circ_missy_angles  = in('circ_anglinear_angles_miss.csv');   writetable(Ya, B.circ_missy_angles);
    B.circ_missy_measure = in('circ_anglinear_measure_miss_y.csv'); writetable(Ym, B.circ_missy_measure);
    B.circ_missy_subjects = cellstr(string(Ym.Subject(ymiss)))';
    B.circ_missy_channels = ach_nan;
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
    for f = {'data1_file', 'data2_file', 'data_file', 'covariate_file', 'count1_file', 'count2_file'}
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
    if isfield(v, 'fl_Rz')                 % partial correlation (Freedman-Lane) only
        E.fl_Rz = v.fl_Rz;
        E.k_eff = v.k_eff;
    end
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
% CAPTURE: circular tier
% =======================================================================
function I = circ_inputs(P)
% The input files as readtable returns them (raw units and convention).
    I = struct();
    I.angles1 = table_block(P.data1_file, {});
    if strcmp(P.comparison, 'circ_corrAngLinear')
        I.angles2 = [];
        I.measure = table_block(P.data2_file, {});
        I.counts1 = [];
        I.counts2 = [];
    else
        I.angles2 = table_block(P.data2_file, {});
        I.measure = [];
        I.counts1 = table_block(P.count1_file, {});
        I.counts2 = table_block(P.count2_file, {});
    end
end

function [M, C] = circ_capture(P, rs, plog, aux, name)
% Rebuilds every stage of core_snpm_circ from the input FILES with the same
% dependency functions (import, montage, count rule, gate, Freedman-Lane
% context, observed statistic), checks each stage against what the engine
% recorded (tap 'circ:engine_inputs'), then replays the one permutation loop
% from the RECORDED indices and checks it against the pipeline's own nulls,
% corrected p-values and clusters.
    alpha = 0.05; E = 0.5; H = 2;
    ei = aux(strcmp({aux.tag}, 'circ:engine_inputs'));
    mn = aux(strcmp({aux.tag}, 'glm:cluster_mass_null'));
    assert(isscalar(ei) && isscalar(mn) && isscalar(plog), 'export_parity:auxCount', ...
        ['%s: expected one circ:engine_inputs record, one cluster-mass null and one ' ...
         'permutation loop; got %d, %d and %d.'], name, numel(ei), numel(mn), numel(plog));
    v = ei.value;
    is_u2   = strcmp(P.comparison, 'circ_phase_group_u2');
    is_corr = strcmp(P.comparison, 'circ_corrAngLinear');
    is_grp  = ~is_corr;
    ck = struct();                         % named replay checks, in this order

    % ---- montage: the angle file's header against the registry ----
    W1 = circ_wide(P.data1_file);
    W2 = circ_wide(P.data2_file);
    S = snpm_setup_channels(P.channels, {zeros(1, numel(W1.cols))}, W1.cols);
    M = montage_struct(P.channels, S);
    nCh = numel(M.labels);
    kept = M.data_columns(:)';
    ck.montage = isequal(M.labels, {rs.chanlocs.labels}) && isequal(double(kept), double(v.keptcols(:)'));

    % ---- import (units, convention, rotation) ----
    gate_mode = 'warn'; if is_u2, gate_mode = 'error'; end
    iopt = struct('units', P.circ_units, 'convention', P.circ_convention, ...
        'zero_offset_deg', pget(P, 'circ_zero_offset_deg', []), 'resolution_gate', gate_mode, ...
        'label', 'Angles file 1');
    [A1, imp] = snpm_circ_import(W1.M, iopt);
    if is_grp
        [tf, loc] = ismember(W1.cols, W2.cols);
        assert(all(tf), 'export_parity:circColumns', '%s: group B lacks channels of group A.', name);
        iopt.label = 'Angles file 2';
        A2 = snpm_circ_import(W2.M(:, loc), iopt);
        n1 = size(A1, 1); n2 = size(A2, 1);
        Araw = [A1; A2];
        gi = [ones(n1, 1); zeros(n2, 1)];
        y = [];
        subj = [W1.subj; W2.subj];
        src_file = [ones(n1, 1); 2 * ones(n2, 1)];
        src_row  = [(1:n1)'; (1:n2)'];
        measure_row = [];
        measure_col = '';
    else
        measure_col = pget(P, 'measure_col', '');
        if isempty(measure_col), measure_col = W2.cols{1}; end
        [subj, i1, i2] = intersect(W1.subj, W2.subj, 'stable');
        y = W2.M(i2, strcmp(W2.cols, measure_col));
        % Subjects without a finite measure are dropped whole, before anything
        % else (core_snpm_circ, core_snpm:circMeasureNaNRows).
        keep = isfinite(y);
        ck_dropped = isequal(subj(~keep)', rs.circ.dropped_subjects(:)') && ...
            rs.circ.measure_nan_removed == sum(~keep);
        subj = subj(keep); i1 = i1(keep); i2 = i2(keep); y = y(keep);
        Araw = A1(i1, :);
        gi = [];
        src_file = ones(numel(i1), 1);
        src_row  = i1(:);
        measure_row = i2(:);
    end
    A_imported = Araw(:, kept);
    ck.import = isequaln(A_imported, v.A_imported) && isequaln(gi, v.gi) && ...
        isequaln(y, v.y) && isequaln(imp, v.import);
    if is_corr, ck.import = ck.import && ck_dropped; end

    % ---- event counts -> NaN rule -> per-channel log covariate ----
    A = A_imported;
    Cnt = []; logC = [];
    if is_grp
        K1 = circ_wide(P.count1_file); K2 = circ_wide(P.count2_file);
        Cnt = [circ_align(K1, W1.cols, W1.subj); circ_align(K2, W1.cols, W2.subj)];
        Cnt = Cnt(:, kept);
        no_count = isfinite(A) & ~(isfinite(Cnt) & Cnt > 0);
        A(no_count) = NaN;
        assert(strcmp(v.covariate.mode, 'per_channel'), 'export_parity:circCovMode', ...
            '%s: the replay implements the per_channel count mode only.', name);
        Cu = Cnt; Cu(~(isfinite(Cu) & Cu > 0)) = NaN;
        logC = log(Cu);
    end
    ck.count_rule = isequaln(A, v.A) && isequaln(logC, v.logC);

    % ---- covariate rules, rebuilt from logC and the observed labels ----
    % (snpm_circ_count_covariate: dropped where sd(log count) < 0.05 or no
    % finite value; collinear where log count is disjoint by group or
    % abs(point-biserial) > circ_collinearity_max, default 0.8)
    if is_grp
        pb = nan(1, nCh);
        for k = 1:nCh
            okk = isfinite(logC(:, k));
            if sum(okk) < 3 || numel(unique(gi(okk))) < 2, continue; end
            cc = corrcoef(logC(okk, k), gi(okk));
            pb(k) = cc(1, 2);
        end
        sdl = std(logC, 0, 1, 'omitnan');
        usecov = ~(sdl < 0.05) & any(isfinite(logC), 1);
        disjoint = false(1, nCh);
        for k = 1:nCh
            a1 = logC(gi == 1, k); a1 = a1(isfinite(a1));
            a0 = logC(gi == 0, k); a0 = a0(isfinite(a0));
            if isempty(a1) || isempty(a0), disjoint(k) = true; continue; end
            disjoint(k) = (min(a1) > max(a0)) || (min(a0) > max(a1));
        end
        collinear = usecov & (disjoint | abs(pb) > pget(P, 'circ_collinearity_max', 0.8));
        usecov(collinear) = false;
        ck.covariate_rules = isequaln(pb, v.covariate.pb) && ...
            isequal(usecov, logical(v.usecov(:)')) && isequal(collinear, logical(v.collinear(:)'));
    else
        usecov = false(1, nCh); collinear = false(1, nCh);
        ck.covariate_rules = ~any(v.usecov) && ~any(v.collinear);
    end

    % ---- gate masks, rebuilt from the observed labels ----
    if is_grp
        n_g1 = sum(isfinite(A) & gi == 1, 1);
        n_g2 = sum(isfinite(A) & gi == 0, 1);
        below = ~(n_g1 >= 8 & n_g2 >= 8);
        n_fin = n_g1 + n_g2;
        gate0 = ~below & ~collinear;
    else
        n_g1 = []; n_g2 = [];
        n_fin = sum(isfinite(A) & isfinite(y), 1);
        below = n_fin < 8;
        gate0 = (sum(isfinite(A), 1) >= 8) & ~collinear;
    end
    gate = gate0 & ~below;
    complete = all(isfinite(A), 1);
    ev = gate & (~usecov | complete);
    incomplete = usecov & ~complete;
    ck.gate = isequal(gate0, logical(v.gate0(:)')) && isequal(gate, logical(v.gate(:)')) && ...
        isequal(complete, logical(v.complete(:)')) && isequal(ev, logical(v.evaluable(:)')) && ...
        isequal(~ev, below | collinear | incomplete);
    reason = repmat({''}, 1, nCh);
    for k = find(~ev)
        r = {};
        if below(k),      r{end+1} = 'fewer than 8 subjects per group'; end %#ok<AGROW>
        if collinear(k),  r{end+1} = 'log(count) collinear with group'; end %#ok<AGROW>
        if incomplete(k), r{end+1} = 'incomplete column under the per-channel covariate'; end %#ok<AGROW>
        reason{k} = strjoin(r, '; ');
    end
    ck.excluded = isequal(find(~ev), double(rs.excluded_channels.index(:)'));

    % ---- inversion check on the pooled sample ----
    [~, inv] = evalc('snpm_circ_check_inversion(A, gate0, S.chanlocs, imp)');
    ck.inversion = isequaln(inv, v.inversion);
    probe = [];
    if strcmpi(P.circ_convention, 'turtlewave_pre_v4')
        % Same files declared as if they were NOT inverted: the guard must refuse.
        Q = P; Q.circ_convention = 'literature_uppeak0';
        probe = struct('convention', Q.circ_convention, 'error_id', '', 'message', '');
        try
            evalc('core_snpm_analysis(Q)');
        catch ME
            probe.error_id = ME.identifier;
            probe.message  = ME.message;
        end
        ck.inversion = ck.inversion && strcmp(probe.error_id, 'core_snpm:circPhaseConventionInverted');
    end

    % ---- observed statistic and the permuted-statistic function ----
    sopt = struct('valid_channels', ev);
    fl = [];
    if is_u2
        [stat, pm] = snpm_circ_watsons_u2(A, gi, sopt);
        permfn = @(ix) snpm_circ_watsons_u2(A, gi(ix), sopt);
        tag = 'circ_u2:labels'; dh = 0.005;
    elseif is_grp
        hctx = snpm_circ_hotelling_perm(gi, logC, usecov, ev);
        [stat, pm] = hctx.eval(cos(A), sin(A));
        fl = snpm_circ_fl_context(A, logC, usecov, ev);
        permfn = @(ix) hotelling_eval(hctx, fl.Zfit + fl.R(ix, :), nCh);
        tag = 'glm_permute:free'; dh = [];
    else
        [stat, pm] = snpm_circ_corr_anglinear(y, A, [], sopt);
        n = numel(y);
        fl = struct('Zfit', repmat(mean(y), n, 1), 'R', y - mean(y), ...
            'eb', ones(n, 1), 'perm_type', 'free');
        permfn = @(ix) snpm_circ_corr_anglinear(fl.Zfit + fl.R(ix, :), A, [], sopt);
        tag = 'glm_permute:free'; dh = [];
    end
    ck.fl_context = isempty(fl) || (isequaln(fl.Zfit, v.fl_Zfit) && isequaln(fl.R, v.fl_R));
    ck.observed_engine = isequaln(stat, v.stat) && isequaln(pm, v.p) && isequaln(dh, v.tfce_dh);
    stat(~ev) = NaN; pm(~ev) = NaN;
    adj = make_neighbors_sparse(S.neighbors, size(S.neighbors, 1));
    if isempty(dh), enh = @(x) ClusterEnhancement(x, adj, E, H);
    else,           enh = @(x) ClusterEnhancement(x, adj, E, H, dh); end
    tfce = enh(stat);
    ck.observed_maps = isequaln(stat, rs.T.real_T) && isequaln(pm, rs.p.real) && ...
        isequaln(tfce, rs.T.real_TFCE);

    % ---- replay the loop from the recorded indices ----
    idx = plog.indices;
    nP = size(idx, 1);
    ck.loop_tag = strcmp(plog.tag, tag) && nP == P.permutations;
    maxT = zeros(nP, 1); maxM = zeros(nP, 1);
    for i = 1:nP
        [s, pp] = permfn(idx(i, :)');
        s(~ev) = NaN; pp(~ev) = NaN;
        maxT(i) = max(enh(s), [], 'omitnan');
        maxM(i) = max_cluster_mass(pp, s, true, alpha, adj);
    end
    maxT = sort(maxT, 'descend'); maxM = sort(maxM, 'descend');
    ck.tfce_null = isequaln(maxT, rs.T.tMaxTFCE);
    ck.mass_null = isequaln(maxM, mn.value);
    ck.cluster_p = cluster_p_match(rs.Clusters, mn.value, nP);
    pT = nan(1, nCh);
    for k = find(isfinite(tfce))
        pT(k) = (snpm_count_at_least(maxT, tfce(k)) + 1) / (nP + 1);
    end
    ck.p_tfce = isequaln(pT, rs.p.correctedTFCE);

    % ---- descriptive panel, recomputed from the analysed angles ----
    if is_grp, sets = {gi == 1, gi == 0}; else, sets = {true(size(A, 1), 1)}; end
    ok = numel(rs.circ_descriptive.group) == numel(sets);
    for k = 1:numel(sets)
        a = A(sets{k}, ev); a = a(isfinite(a)); nn = numel(a);
        r  = abs(sum(exp(1i * a))) / nn;
        mu = mod(rad2deg(angle(sum(exp(1i * a)))), 360);
        sd = rad2deg(sqrt(-2 * log(max(r, realmin))));
        pv = min(max(exp(sqrt(1 + 4*nn + 4*(nn^2 - (nn*r)^2)) - (1 + 2*nn)), 0), 1);
        g = rs.circ_descriptive.group(k);
        ok = ok && isequal(g.n, nn) && isequal(g.mean_deg, mu) && isequal(g.sd_deg, sd) && ...
            isequal(g.R, r) && isequal(g.rayleigh_z, nn * r^2) && isequal(g.rayleigh_p, pv);
    end
    ck.descriptive = ok;

    % ---- assemble ----
    C = struct();
    C.comparison = P.comparison;
    if is_u2
        C.statistic = 'watson_U2';
    elseif is_grp
        C.statistic = 'hotelling_F';
    else
        C.statistic = 'circlinear_F';
    end
    C.alpha = alpha; C.tfce_E = E; C.tfce_H = H;
    if isempty(dh), C.tfce_dh = 0.1; else, C.tfce_dh = dh; end
    C.cluster_statistic = 'mass = sum(stat)/n_channels over channels with p_uncorrected < alpha';
    C.import = imp;
    C.subject_ids = subj(:);
    C.source_file = src_file;
    C.source_row  = src_row;
    C.measure_col = measure_col;
    C.measure_row = measure_row;
    C.angles_imported = A_imported;
    C.angles = A;
    C.group  = gi;
    C.y      = y;
    C.counts = Cnt;
    C.log_count = logC;
    C.covariate = v.covariate;
    C.usecov    = usecov;
    C.collinear = collinear;
    C.n_group1 = n_g1; C.n_group2 = n_g2; C.n_finite = n_fin;
    C.gate0 = gate0; C.gate = gate;
    C.complete = complete; C.evaluable = ev;
    C.excluded_below_gate  = below & ~ev;
    C.excluded_collinear   = collinear & ~ev;
    C.excluded_incomplete  = incomplete & ~ev;
    C.excluded_reason      = reason;
    C.inversion       = v.inversion;
    C.inversion_probe = probe;
    if isempty(fl)
        C.fl_Zfit = []; C.fl_R = []; C.perm_type = 'labels';
    else
        C.fl_Zfit = fl.Zfit; C.fl_R = fl.R; C.perm_type = fl.perm_type;
    end
    C.stat = stat; C.p_uncorrected = pm; C.tfce = tfce;
    C.tfce_null_sorted = maxT;
    C.cluster_mass_null_sorted = mn.value;
    C.cluster_mass_null_source = 'recorded in snpm_perm_correction (tap glm:cluster_mass_null)';
    C.p_correctedTFCE = pT;
    C.replay_checks = fieldnames(ck)';
    C.replay_detail = cellfun(@(f) logical(ck.(f)), C.replay_checks);
    C.replay_ok = all(C.replay_detail);
end

function [s, pp] = hotelling_eval(hctx, Yp, nCh)
% The stacked response: columns 1..nCh = cos, nCh+1..2nCh = sin.
    [s, pp] = hctx.eval(Yp(:, 1:nCh), Yp(:, nCh+1:end));
end

function S = circ_wide(file)
% core_snpm_circ/read_wide: subject column = first name containing 'subject',
% ids as sub%03d when numeric.
    T = readtable(file);
    vn = T.Properties.VariableNames;
    is_s = contains(lower(vn), 'subject');
    raw = T{:, find(is_s, 1)};
    if isnumeric(raw)
        S.subj = arrayfun(@(x) sprintf('sub%03d', x), raw, 'UniformOutput', false);
    else
        S.subj = cellstr(string(raw));
    end
    S.subj = S.subj(:);
    S.cols = vn(~is_s);
    S.M = T{:, ~is_s};
end

function M = circ_align(S, cols, subj)
% Count matrix in the angle file's channel and subject order (exact set match).
    [tf, cl] = ismember(cols, S.cols);
    [ts, rl] = ismember(subj, S.subj);
    assert(all(tf) && all(ts) && numel(S.cols) == numel(cols) && numel(S.subj) == numel(subj), ...
        'export_parity:circCounts', 'count file does not match its angle file');
    M = S.M(rl, cl);
end

function v = pget(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
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
    w('`inputs/miss_corr_eeg_cov.csv` is `corr_eeg.csv` with one EEG cell blank in each of %s (covariates complete),', ...
        strjoin(B.miss_cov_channels, ' and '));
    w('so those two channels are masked by the complete-column rule and the Freedman-Lane fit uses all 30 subjects.');
    w('`inputs/small_corr_covariates_n5.csv` is the first 5 rows of `corr_covariates.csv`, `Subject` and `age` only.');
    w('');
    w('Circular inputs. `example_data/circ_*` are radians on [0, 2*pi) with zero at the up-state peak. Built here:');
    w('`inputs/circ_gate_*.csv` are copies of the `circ_phase_*` pair with four channels edited. `%s`: 9 of 16 group-A', B.circ_gate_channels{1});
    w('angles blank (7 left, below the 8-per-group gate) and its counts set to 20 everywhere, so the gate is its only reason.');
    w('`%s`: one group-B angle blank while the count covariate is in use, so the column is incomplete and excluded.', B.circ_gate_channels{2});
    w('`%s`: every group-A count set above every group-B count, so log(count) is disjoint by group (collinear) and', B.circ_gate_channels{3});
    w('the channel is excluded. `%s`: counts 20 everywhere (constant, so the covariate is dropped and the channel', B.circ_gate_channels{4});
    w('uses the covariate-free Hotelling with pairwise deletion), group-A row 3 count 0 (its angle is dropped by the');
    w('count rule) and group-B row 5 angle blank; it stays in with 15 + 15 angles.');
    w('`inputs/circ_luna_deg_group*.csv` = `mod(rad2deg(a) - 90, 360)` of the example angles (degrees, Luna convention,');
    w('declared `deg` + `luna_zerocross0`; import adds 90 deg). `inputs/circ_tw_inverted_group*.csv` = `mod(a + pi, 2*pi)`');
    w('(declared `rad` + `turtlewave_pre_v4`; import adds 180 deg). Both reuse the example count files. Text round-trips');
    w('through CSV, so the imported angles equal the example angles only to ~1e-14 rad, not bit for bit.');
    w('`inputs/circ_anglinear_measure_miss_y.csv` is `circ_anglinear_measure.csv` with `measure` blank for %s and', ...
        strjoin(B.circ_missy_subjects(1:end-1), ', '));
    w('`Inf` for %s (Inf counts as missing). `inputs/circ_anglinear_angles_miss.csv` is `circ_anglinear_angles.csv` with', ...
        B.circ_missy_subjects{end});
    w('3 angles blank in each of %s.', strjoin(B.circ_missy_channels, ', '));
    w('circ_corrAngLinear drops every subject whose measure is not finite right after subject matching: they are absent');
    w('from `circ.subject_ids`, `circ.y`, `circ.angles` and every map, and listed in `results.circ.dropped_subjects`');
    w('(count `results.circ.measure_nan_removed`). Missing angles stay with their row and are deleted per channel, so');
    w('`results.per_channel_n` is the same in the observed map and in every permuted map.');
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
    w('| correlationP/S **with covariates** (`*_cov`, `miss_correlationP_cov_nan`) | `snpm_row_permutations` | row orders, units = subjects | Freedman-Lane (since 2026-10-18): `r = corr(engine.fl_Rz * data_x(order,:), data_y)`, i.e. the permuted residuals are re-projected off the covariates before correlating; the observed map is the same expression with `order = 1:n`. A column with `sum((fl_Rz*data_x(order,c)).^2) <= 1e-20*sum(data_x(:,c).^2)` scores 0. Partial Spearman is handed on as `correlationP` on rank residuals. Uncorrected p: `t = r*sqrt(df/(1-r^2))`, `df = n - 2 - engine.k_eff`. |');
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
    w('**Circular tier** (circ_phase_group, circ_phase_group_u2, circ_corrAngLinear): ONE loop, `circ:<comparison>`,');
    w('inside `snpm_perm_correction`; the TFCE and cluster-mass nulls share its draws. There is no global test and no');
    w('post-hoc. The observed labelling is NOT forced into row 1. Each row of `perm_log(1).indices` is one `randperm(n)`:');
    w('');
    w('| Analysis | Tag | How a row `ix` is applied |');
    w('|---|---|---|');
    w('| circ_phase_group (Hotelling) | `glm_permute:free` | Freedman-Lane on the STACKED response `Y = [cos(A), sin(A)]`, n x 2*nCh (columns 1..nCh = cos of channel 1..nCh, nCh+1..2nCh = sin of the same channels). `Y_perm = circ.fl_Zfit + circ.fl_R(ix, :)`: ONE index moves every column, so cos and sin of a subject always move together. The permuted statistic is evaluated on `Y_perm(:, 1:nCh)`, `Y_perm(:, nCh+1:end)` directly (it is off the unit circle; there is no angle). `circ.group` and the count covariate stay fixed. |');
    w('| circ_phase_group_u2 (Watson U^2) | `circ_u2:labels` | label permutation: the permuted map is U^2 of the FIXED angles with group vector `circ.group(ix)`. Nothing is re-centred. |');
    w('| circ_corrAngLinear | `glm_permute:free` | Freedman-Lane with the intercept as the only nuisance: `y_perm = circ.fl_Zfit + circ.fl_R(ix)` (n x 1; `fl_Zfit` = mean of y; subjects with a non-finite y were dropped before this, so y is complete), against the FIXED angles. |');
    w('');
    w('**Event-count covariate (Hotelling only).** Per channel ch, the nuisance is `[1, log(count(:, ch))]` where');
    w('`circ.usecov(ch)` is true and the intercept alone elsewhere; the design for the statistic is `X = [1, logC_ch, g]`');
    w('(or `[1, g]`), `g` = `circ.group` (1 = group A, 0 = group B). The same `logC_ch` column is used for the cos and the sin');
    w('column of that channel. `fl_Zfit`/`fl_R` are the per-channel least-squares fit and residual of each stacked column on');
    w('its nuisance, fitted once before the loop (`snpm_circ_fl_context`). Covariate channels: `F = T2 (nu-1)/(2 nu)` with');
    w('`nu = n - 3`, `T2 = (c''B) S^-1 (c''B)'' / (c''(X''X)^-1 c)`, `S` the 2 x 2 residual covariance of the full model on `nu`, `F ~');
    w('F(2, nu - 1)`. Covariate-free channels: classic pooled two-sample Hotelling with per-channel pairwise deletion,');
    w('`nu = n1 + n2 - 2`. A permuted residual row carries its NaN with it. `usecov` is false where sd(log count) < 0.05 or');
    w('the channel is collinear. A cell with an angle but count <= 0 or missing has its ANGLE set to NaN before anything else.');
    w('Watson U^2 never uses the covariate, but the count files still decide the evaluable set (same gate as Hotelling).');
    w('');
    w('**Gate (computed once, from the observed labels, never per permutation).** A channel is evaluable when it has >= 8');
    w('finite angles in each group (circ_corrAngLinear: >= 8 subjects with an angle, y being complete), is not collinear, and, if');
    w('`usecov`, is complete. Non-evaluable channels are NaN in the observed map and in every permuted map, before TFCE.');
    w('TFCE: `ClusterEnhancement(stat, adj, 0.5, 2, dh)` with dh = 0.1 (Hotelling, circular-linear) or 0.005 (U^2) for');
    w('the observed AND every permuted map. Cluster forming: channels with parametric p < 0.05, connected through');
    w('`montage.neighbors`; cluster mass = sum(stat)/n_channels (stat = F or U^2 itself, NOT squared). TFCE null = max over');
    w('channels of the TFCE map (non-negative statistics). p = (b + 1)/(N + 1) with ties at relative 1e-12 counted.');
    w('');
    w('## Field dictionary');
    w('');
    w('Top level of every fixture:');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `case_name` | char | fixture name (= file name without `.mat`) |');
    w('| `tier` | char | `legacy`, `glm` or `circ` |');
    w('| `description` | char | what the case is for |');
    w('| `seed`, `rng_algorithm` | scalar, char | `rng(seed, rng_algorithm)` was called immediately before the run |');
    w('| `permutations` | scalar | permutations REQUESTED (`params.permutations`) |');
    w('| `matlab_version`, `created` | char | provenance |');
    w('| `params` | struct | the params struct passed to `core_snpm_analysis` (paths are absolute, this machine) |');
    w('| `source_files` | struct | per input (`data1_file`, `data2_file`, `data_file`, `covariate_file`, `count1_file`, `count2_file`): `.absolute` and `.relative_to_repo` |');
    w('| `inputs` | struct | the input files as `readtable` returned them, before any processing (below) |');
    w('| `montage` | struct | the channels the run analysed and their adjacency (below) |');
    w('| `perm_log` | struct array, K x 1 | one entry per permutation loop, in call order (below) |');
    w('| `results` | struct | the complete `results_struct` returned by the pipeline, plain types (below) |');
    w('| `results_text` | cellstr | the text summary the pipeline returned |');
    w('| `engine` | struct | legacy tier only: the matrices the engines were handed, the cluster-size null, the global channel set (below) |');
    w('| `glm` | struct | GLM tier only: designs and replayed nulls (below) |');
    w('| `circ` | struct | circular tier only: every stage of the circular pipeline and the replayed null (below) |');
    w('');
    w('`inputs` -- legacy tier has `data1`, `data2` (empty for `onesampleT_single`), `covariates` (empty when unused);');
    w('GLM tier has `data`; circular tier has `angles1`, `angles2`, `counts1`, `counts2` (two-group analyses) or `angles1`,');
    w('`measure` (circ_corrAngLinear), the others empty, all in the file''s own units and convention. Each data block:');
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
    w('| `ncov` | scalar | covariate rank `k_eff = rank([1, covariates]) - 1` (0 without covariates); reduces the df of the global correlation p and, with covariates, of the uncorrected per-channel p |');
    w('| `fl_Rz`, `k_eff` | n x n, scalar | covariate cases only (absent otherwise): the Freedman-Lane residual-forming matrix `I - Q*Q''` of `[1, covariates]` on the analysed rows (rank-revealing QR) and its rank minus 1. `data_x = fl_Rz * x`, `data_y = fl_Rz * y`. |');
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
    w('`results`, circular tier: `T.real_T` (F or U^2), `T.real_TFCE`, `T.tMaxTFCE` (sorted descending), `T.critical_T_TFCE`,');
    w('`T.excluded_channels`, `p.real`, `p.correctedTFCE`, `Clusters(c)` (`.channels`, `.mass`, `.p`, `.threshold` = alpha,');
    w('`.permutations`), the significant-channel sets, `alpha`, `chanlocs`, `excluded_channels` (one combined reason string),');
    w('`per_channel_n` (finite angles per channel), `hide_condition_b`, `circ` (`.import`, `.rotation_deg`, `.convention`,');
    w('`.units`, `.covariate`, `.gate_min_per_group` = 8, `.evaluable`, `.n_group1`, `.n_group2`, and `.rho` for');
    w('circ_corrAngLinear), `circ_descriptive` (`.group(k)`: `.label`, `.n` angles pooled over evaluable channels, `.mean_deg`,');
    w('`.sd_deg` = circular SD in degrees, `.R`, `.rayleigh_z` = nR^2, `.rayleigh_p`; `.n_channels`, `.note`), `glm.effect_label`,');
    w('`data_summary`. No `global_stat`/`global_pval`: the tier has no whole-head test, deliberately.');
    w('');
    w('`circ` (circular tier only; rebuilt by the export script from the input files with the pipeline''s own functions, each');
    w('stage checked against what `core_snpm_circ` recorded, then the loop replayed from `perm_log(1).indices`):');
    w('');
    w('| Field | Type / shape | Meaning |');
    w('|---|---|---|');
    w('| `comparison`, `statistic` | char | `hotelling_F`, `watson_U2` or `circlinear_F` |');
    w('| `alpha`, `tfce_E`, `tfce_H`, `tfce_dh` | scalars | 0.05, 0.5, 2, and 0.1 or 0.005 (U^2) |');
    w('| `cluster_statistic` | char | how cluster mass is computed |');
    w('| `import` | struct | `snpm_circ_import` record of angle file 1: `.units`, `.convention`, `.rotation_deg` (added before wrapping), `.resolution_deg`, `.quantised`, `.whole_degree`, `.n_finite`, `.label`. Import = `mod((deg + rotation_deg) * pi/180, 2*pi)`, `deg` = raw (deg) or raw*180/pi (rad); NaN stays NaN. File 2 gets the same transform. |');
    w('| `subject_ids`, `source_file`, `source_row` | n x 1 | analysed row r is row `source_row(r)` of angle file `source_file(r)` (1 = `inputs.angles1`, 2 = `inputs.angles2`). Two-group: group A rows then group B rows, file order. |');
    w('| `measure_col`, `measure_row` | char, n x 1 | circ_corrAngLinear: the y column and the row of `inputs.measure` matched to each analysed row (subjects intersected in angle-file order) |');
    w('| `angles_imported` | n x nCh, radians | after import and channel selection (analysed order), BEFORE the count rule |');
    w('| `angles` | n x nCh, radians | as analysed: `angles_imported` with the angle of every cell whose count is missing or <= 0 set to NaN |');
    w('| `group` | n x 1 | 1 = group A (file 1), 0 = group B; empty for circ_corrAngLinear |');
    w('| `y` | n x 1 | circ_corrAngLinear measure; empty otherwise |');
    w('| `counts` | n x nCh | event counts aligned to `angles` (by subject id and channel name); empty for circ_corrAngLinear |');
    w('| `log_count` | n x nCh | `log(counts)`, NaN where the count is missing or <= 0: the per-channel covariate |');
    w('| `covariate` | struct | `snpm_circ_count_covariate` report: `.mode` (`per_channel`), `.pb` (1 x nCh point-biserial r of log count vs group), `.pb_max`, `.pb_median`, `.n_dropped_constant`, `.n_collinear`, `.n_cells_no_count`, `.frac_cells_no_count`, `.collinearity_max` (0.8); circ_corrAngLinear: `.mode` = `none` |');
    w('| `usecov` | logical 1 x nCh | the count covariate enters the channel''s design (false where constant, sd < 0.05, or collinear) |');
    w('| `collinear` | logical 1 x nCh | log count disjoint by group or abs(pb) > 0.8 |');
    w('| `n_group1`, `n_group2`, `n_finite` | 1 x nCh | finite angles per group / in total (circ_corrAngLinear: subjects with angle and y) |');
    w('| `gate0` | logical 1 x nCh | >= 8 per group and not collinear: the channel set the inversion check pools over |');
    w('| `gate`, `complete`, `evaluable` | logical 1 x nCh | gate after the statistic''s own per-group count; column complete; `evaluable = gate & (~usecov | complete)`, the ONE channel set of the observed and every permuted map |');
    w('| `excluded_below_gate`, `excluded_collinear`, `excluded_incomplete` | logical 1 x nCh | why each non-evaluable channel left (a channel can have several reasons) |');
    w('| `excluded_reason` | cellstr 1 x nCh | the same in words, `''''` for an analysed channel |');
    w('| `inversion` | struct | the convention check: `.channels` (frontal third by chanlocs X among `gate0`, 1-based), `.n_angles` pooled, `.mean_deg` (pooled mean direction), `.deviation_deg` (wrapped to (-180, 180]), `.verdict` (`pass`, `suspect` = warning at abs(dev) > 70, `skipped`). abs(dev) >= 110 is an error and never reaches a fixture. |');
    w('| `inversion_probe` | struct or empty | turtlewave fixture only: the same run declared `literature_uppeak0`; `.error_id` must be `core_snpm:circPhaseConventionInverted`, `.message` carries the numbers |');
    w('| `fl_Zfit`, `fl_R` | n x 2nCh (Hotelling) or n x 1 (circular-linear) | Freedman-Lane fit and residual; non-evaluable columns: `fl_Zfit` 0 and `fl_R` NaN. Empty for U^2. |');
    w('| `perm_type` | char | `free` (Freedman-Lane) or `labels` (U^2) |');
    w('| `stat`, `p_uncorrected`, `tfce` | 1 x nCh | observed maps, NaN outside `evaluable` (= `results.T.real_T`, `results.p.real`, `results.T.real_TFCE`) |');
    w('| `tfce_null_sorted` | n_draws x 1, **sorted descending** | replayed max-TFCE null (= `results.T.tMaxTFCE`) |');
    w('| `cluster_mass_null_sorted` | n_draws x 1, **sorted descending** | max cluster mass per permutation (0 when none), recorded inside `snpm_perm_correction`; every `Clusters(c).p` is computed against it |');
    w('| `p_correctedTFCE` | 1 x nCh | replayed TFCE-corrected p (= `results.p.correctedTFCE`) |');
    w('| `replay_checks`, `replay_detail`, `replay_ok` | cellstr, logical 1 x K, logical | each stage the export rebuilt independently from the input files and found equal to the pipeline: montage, import, count rule, covariate rules (point-biserial, usecov, collinear), gate masks (gate0, gate, complete, evaluable), excluded set, inversion check (and probe), Freedman-Lane context, engine statistic, observed maps, loop tag and size, TFCE null, mass null, cluster p, TFCE p, descriptive panel |');
    w('');
    w('## Not in these fixtures');
    w('');
    w('- **A global test for single-effect GLM presets.** The pipeline computes none (mixed2way has one per effect).');
    w('- **A global test for the circular tier.** None exists, by design; `circ_descriptive` is the no-inference substitute.');
    w('- The circular `per_subject` count mode, the `custom` convention, and Tier-1 linearised measures.');
    w('- Mixed-model analyses (a later phase).');
end
