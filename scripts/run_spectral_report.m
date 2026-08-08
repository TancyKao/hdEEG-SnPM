%% run_spectral_report.m  —  headless two-condition spectral report (no GUI; for NCI)
% Within-subject paired contrast (condition A vs B), per stage x band, on the
% 178-ch montage, TFCE + cluster correction, periodograms + topoplots, and a
% ready-to-open HTML report. Wraps export_report.m.
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_spectral_report.m')"

% ======================= CONFIG (edit me) =======================
ROOT     = fileparts(fileparts(mfilename('fullpath')));   % scripts/ -> repo root
if isempty(ROOT), ROOT = pwd; end
DATA     = fullfile(ROOT,'test_data','hdeeg_analysis_all_sub');   % folder of per-subject spectral .mat files
OUT      = fullfile(ROOT,'test_data','hdeeg_paired_report_2026-07');

opts = struct();
opts.permutations     = 5000;
opts.condA            = 'Condition A';
opts.condB            = 'Condition B';
opts.fmin             = 0.1;     opts.fmax = 40;
opts.add_nrem         = true;                 % add NREM = mean(N2,N3)
opts.seed             = 20260722;             % one seed for the whole sweep
opts.exclude_subjects = {};                   % EDIT ME: drop whole subjects, e.g. {'sub-XX'}
% EDIT ME: bad channels, NaN-masked for that subject only (pairwise deletion)
% across all bands, both conditions, all stages, and the global-PSD channel
% average. Format: {subject_id, channel_label; ...}, e.g.
%   opts.exclude_chan_subject = {'sub-XX','E197'; 'sub-YY','E184'};
% An entry that matches no record is a hard error, so leave it {} until you have
% identified your own artifact channels.
opts.exclude_chan_subject = {};
% opts.periodogram_only = true;               % fast path: only regenerate periodograms
% ================================================================

cd(ROOT); addpath(genpath(ROOT));
if ~exist(OUT,'dir'), mkdir(OUT); end
% No template copy: export_report reads templates/sleep_eeg_report.html from the
% toolbox. OUT must contain results only.

EXPECT_N = 27;   % matched subjects for this dataset; guards a silent subject loss

export_report(DATA, OUT, opts);

% n is only reported inside REPORT.js -- assert it rather than trusting the log.
js = fileread(fullfile(OUT,'REPORT.js'));
n  = str2double(regexp(js,'n_subj:\s*(\d+)','tokens','once'));
assert(n==EXPECT_N, 'REPORT.n_subj is %d, expected %d - subjects were lost.', n, EXPECT_N);
nc = str2double(regexp(js,'"low-delta":(\d+)','tokens','once'));
fprintf('Asserted REPORT.n_subj = %d (first per-cell n = %d)\n', n, nc);
fprintf('Done: spectral report -> %s\n', fullfile(OUT,'sleep_eeg_report_filled.html'));
