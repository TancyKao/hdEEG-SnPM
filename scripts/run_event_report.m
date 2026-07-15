%% run_event_report.m  —  headless event-level (spindle/SO) report (no GUI; for NCI)
% Two-group contrast per parameter (density/amplitude/duration) x event band
% (slow spindle / fast spindle / slow wave) on the 256-ch montage, TFCE +
% cluster correction, group-mean + signed t-map topoplots, and an HTML report.
% Wraps export_event_report.m.
%
% INPUT layout in DATA (produced by db_to_group_table.py — see README_scripts.md):
%   subjects.csv                          (Subject + group + covariates)
%   eventStat_density_spindle_9-12.csv    (rows=Subject, cols=E1..E256[+Cz])
%   eventStat_density_spindle_12-15.csv
%   eventStat_density_sw_0.5-1.25.csv
%   eventStat_amplitude_*.csv, eventStat_duration_*.csv
%
% Run headless / on NCI:
%   matlab -batch "run('scripts/run_event_report.m')"

% ======================= CONFIG (edit me) =======================
ROOT = '/path/to/hdEEG-SnPM';
DATA = fullfile(ROOT,'test_data','synthetic_events');     % folder with subjects.csv + eventStat_*.csv
OUT  = fullfile(ROOT,'03_results','statistical_outputs','hdEEG_events');

opts = struct();
opts.group_col    = 'group';
opts.permutations = 1000;
opts.html         = fullfile(ROOT,'sleep_eeg_report_filled.html');   % template to inject into
% ================================================================

addpath(genpath(ROOT));
if ~exist(OUT,'dir'), mkdir(OUT); end
export_event_report(DATA, OUT, opts);
fprintf('Done: event report -> %s\n', fullfile(OUT,'event_report.html'));
