%% run_spectral_report.m  —  headless two-condition spectral report (no GUI; for NCI)
% Within-subject paired contrast (condition A vs B), per stage x band, on the
% 178-ch montage, TFCE + cluster correction, periodograms + topoplots, and a
% ready-to-open HTML report. Wraps export_report.m.
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_spectral_report.m')"

% ======================= CONFIG (edit me) =======================
ROOT     = '/path/to/hdEEG-SnPM';
DATA     = '/path/to/your-study/derivatives/eeg/spectral_power';   % folder of per-subject spectral .mat files
OUT      = fullfile(ROOT,'03_results','statistical_outputs','hdEEG_overnightPSA'); % must contain sleep_eeg_report.html (master in templates/sleep_eeg_report.html)

opts = struct();
opts.permutations     = 1000;
opts.condA            = 'Condition A';
opts.condB            = 'Condition B';
opts.fmin             = 0.1;     opts.fmax = 40;
opts.add_nrem         = true;                 % add NREM = mean(N2,N3)
opts.exclude_subjects = {'sub-XX'};           % drop bad-data subjects (e.g. 20 Hz artifact)
% opts.periodogram_only = true;               % fast path: only regenerate periodograms
% ================================================================

addpath(genpath(ROOT));
export_report(DATA, OUT, opts);
fprintf('Done: spectral report -> %s\n', fullfile(OUT,'sleep_eeg_report_filled.html'));
