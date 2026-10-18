function func_genSnpmTable(S, uncorrsigch, correctTFCEsigch, SnPMsigch, chanlocs, outname, select_mode)

% select_mode ('label'|'positional'): in label mode the chanlocs align 1:1 with
% the statistic columns, so labels are assigned directly. The legacy positional
% 256-column case drops the 'Cz' label (UWM montage). Defaults to positional so
% any older caller keeps its behavior.
if nargin < 7 || isempty(select_mode), select_mode = 'positional'; end

% Check if S is a results_struct (from core_snpm_analysis) or direct SnPM output
if isfield(S, 'T') && isfield(S, 'p')
    % S is a results_struct from core_snpm_analysis
    T_real_T = array2table(round(S.T.real_T,4));
    T_p_real = array2table(round(S.p.real,4));
    T_p_correctedTFCE = array2table(round(S.p.correctedTFCE,4));
else
    % S is direct SnPM output (original format)
    T_real_T = array2table(round(S.real_T,4));
    T_p_real = array2table(round(S.p_real,4));
    T_p_correctedTFCE = array2table(round(S.p_correctedTFCE,4));
end

Chanlabels = string({chanlocs.labels});
if strcmp(select_mode, 'positional') && width(T_real_T) == 256
    Chanlabels(Chanlabels == "Cz") = []; % legacy UWM 256-column layout
end

% Every per-channel sheet is headed by the channel labels in the analysed
% order. The old guard on the p sheet compared length({T_p_real}) - always 1,
% the length of a one-cell array - with the channel count, so that sheet was
% never labelled and came out as Var1..VarN. A width mismatch is now reported
% instead of silently leaving generic headers.
T_real_T          = label_columns(T_real_T,          Chanlabels, 'statVal');
T_p_real          = label_columns(T_p_real,          Chanlabels, 'pVal_or_pCritVal');
T_p_correctedTFCE = label_columns(T_p_correctedTFCE, Chanlabels, 'pValTFCE');

output_xls = [outname, '.xlsx'];
writetable(T_real_T, output_xls,'Sheet','statVal')
writetable(T_p_real, output_xls,'Sheet','pVal_or_pCritVal')
writetable(T_p_correctedTFCE,output_xls,'Sheet','pValTFCE')

% Correlation only: effective N per channel. TWO quantities that are easy to
% confuse, so they are labelled row-wise in one sheet ('Metric' column):
%   n_analysed   the subjects the statistic was actually computed on. Under the
%                complete-column rule this is the same number on every retained
%                channel and 0 on every excluded one - a constant plus an
%                exclusion flag, not a per-channel quantity.
%   n_available  the matched pairs the channel HAD before the exclusion rule was
%                applied. This is the diagnostic that varies: it separates a
%                channel that was one subject short from one that was half empty,
%                which the post-mask column cannot show.
% Guarded so GLM/LMM callers (no per_channel_n field) are unaffected; n_available
% is optional, so an older results_struct (or core_snpm_circ, which sets only
% per_channel_n) still writes the single labelled row.
if isfield(S, 'per_channel_n') && ~isempty(S.per_channel_n)
    effN_rows   = {double(S.per_channel_n(:)')};
    effN_metric = {'n_analysed'};
    if isfield(S, 'n_available') && ~isempty(S.n_available) && ...
            numel(S.n_available) == numel(S.per_channel_n)
        effN_rows{end+1}   = double(S.n_available(:)');
        effN_metric{end+1} = 'n_available';
    end
    T_effN = label_columns(array2table(vertcat(effN_rows{:})), Chanlabels, 'effectiveN');
    T_effN = addvars(T_effN, string(effN_metric(:)), 'Before', 1, ...
        'NewVariableNames', 'Metric');
    writetable(T_effN, output_xls, 'Sheet', 'effectiveN')
end


% Channels dropped before the permutation loop because missing data made them
% unevaluable in some permutations (GLM tier, snpm_glm_fl_context). Always
% written when the field exists so a reader of the workbook can see that the
% analysed channel set is smaller than the montage - and see exactly which.
% One row per (analysis, channel), with that channel's own reason. Since the
% zero-variance rule (2026-10-14) a channel can leave for missing data, for
% having no variation left after the nuisance model, or because its observed
% statistic is undefined, and each GLM loop (omnibus, each post-hoc pair, each
% mixed2way effect and simple effect) lists its own; the Analysis column says
% which loop a row belongs to.
if isfield(S, 'excluded_channels') && isstruct(S.excluded_channels)
    T_excl = excluded_sheet(S);
    if isempty(T_excl)
        T_excl = table("None", "All channels evaluable (no missing data)", ...
            'VariableNames', {'ChannelLabel','Reason'});
    end
    writetable(T_excl, output_xls, 'Sheet', 'excludedChannels')
end

% Pairwise post-hoc comparisons (GLM tier): one row per pair, RUN OR NOT. A
% pair is tested on the rows of its two levels only, so it can rest on fewer
% subjects than the effect, and a pair that cannot be estimated is not run at
% all; neither fact is visible from the maps, so both are tabulated here.
T_ph = posthoc_sheet(S);
if ~isempty(T_ph)
    writetable(T_ph, output_xls, 'Sheet', 'posthoc')
end

if ~isempty(uncorrsigch)
    writetable(array2table(uncorrsigch),output_xls,'Sheet','uncorrectsigch')
elseif isempty(uncorrsigch)
    writetable(table({'No significant chans'}),output_xls,'Sheet','uncorrectsigch')
end

if ~isempty(correctTFCEsigch)
    writetable(array2table(correctTFCEsigch),output_xls,'Sheet','correctTFCEsigch')
elseif isempty(correctTFCEsigch)
    writetable(table({'No significant chans'}),output_xls,'Sheet','correctTFCEsigch')
end

if ~isempty(SnPMsigch)
    writetable(array2table(SnPMsigch),output_xls,'Sheet','correctSnPMsigch')
elseif isempty(SnPMsigch)
    writetable(table({'No significant chans'}),output_xls,'Sheet','correctSnPMsigch')
end

end

function T = label_columns(T, labels, sheet)
% Head a channels-in-columns table with the channel labels. Warns (does not
% throw) on a width mismatch so the rest of the workbook is still written.
if width(T) == numel(labels)
    T.Properties.VariableNames = cellstr(labels);
else
    warning('func_genSnpmTable:labelMismatch', ...
        'Sheet %s has %d columns but %d channel labels; left with generic headers.', ...
        sheet, width(T), numel(labels));
end
end

function T = excluded_sheet(S)
% ChannelIndex | ChannelLabel | Reason | Analysis, over every loop of the run.
    rows = {};
    if isfield(S, 'effects') && numel(S.effects) > 1
        for k = 1:numel(S.effects)
            e = S.effects(k);
            rows = add_excl_rows(rows, effect_name(e), e);
            rows = add_posthoc_excl_rows(rows, effect_name(e), e);
            if isfield(e, 'simple_effects')
                for s = 1:numel(e.simple_effects)
                    se = e.simple_effects(s);
                    rows = add_excl_rows(rows, effect_name(se), se);
                    rows = add_posthoc_excl_rows(rows, effect_name(se), se);
                end
            end
        end
    else
        nm = '';
        if isfield(S.excluded_channels, 'context'), nm = char(S.excluded_channels.context); end
        if isfield(S, 'glm') && isfield(S.glm, 'effect_label'), nm = char(S.glm.effect_label); end
        rows = add_excl_rows(rows, nm, S);
        rows = add_posthoc_excl_rows(rows, nm, S);
    end
    T = [];
    if isempty(rows), return; end
    T = cell2table(rows, 'VariableNames', {'ChannelIndex', 'ChannelLabel', 'Reason', 'Analysis'});
end

function rows = add_excl_rows(rows, name, e)
    if ~isfield(e, 'excluded_channels') || ~isstruct(e.excluded_channels), return; end
    ex = e.excluded_channels;
    if ~isfield(ex, 'n') || ex.n == 0, return; end
    if isfield(ex, 'reasons') && numel(ex.reasons) == numel(ex.index)
        r = ex.reasons;
    else
        r = repmat({char(ex.reason)}, 1, numel(ex.index));
    end
    for i = 1:numel(ex.index)
        rows(end+1, :) = {ex.index(i), char(ex.labels{i}), char(r{i}), name}; %#ok<AGROW>
    end
end

function rows = add_posthoc_excl_rows(rows, name, e)
    if ~isfield(e, 'posthoc'), return; end
    for i = 1:numel(e.posthoc)
        rows = add_excl_rows(rows, sprintf('%s: post-hoc %s', name, char(e.posthoc(i).label)), ...
            e.posthoc(i));
    end
end

function T = posthoc_sheet(S)
% Effect | Pair | Status | RunOn | Rows | Subjects | SubjectsLeftOut | TFCEsigChannels | ClusterSigChannels
    rows = {};
    if isfield(S, 'effects') && numel(S.effects) > 1
        for k = 1:numel(S.effects)
            e = S.effects(k);
            rows = add_effect_rows(rows, effect_name(e), e);
            if isfield(e, 'simple_effects')
                for s = 1:numel(e.simple_effects)
                    rows = add_effect_rows(rows, effect_name(e.simple_effects(s)), e.simple_effects(s));
                end
            end
        end
    else
        nm = 'main effect';
        if isfield(S, 'glm') && isfield(S.glm, 'effect_label'), nm = char(S.glm.effect_label); end
        rows = add_effect_rows(rows, nm, S);
    end
    T = [];
    if isempty(rows), return; end
    T = cell2table(rows, 'VariableNames', {'Effect', 'Pair', 'Status', 'RunOn', 'Rows', ...
        'Subjects', 'SubjectsLeftOut', 'TFCEsigChannels', 'ClusterSigChannels'});
end

function nm = effect_name(e)
    nm = '';
    if isfield(e, 'effect_label'), nm = char(e.effect_label); end
    if isempty(nm) && isfield(e, 'key'), nm = char(e.key); end
end

function rows = add_effect_rows(rows, name, e)
    if isfield(e, 'posthoc')
        for i = 1:numel(e.posthoc)
            q = e.posthoc(i);
            rows(end+1, :) = {name, char(q.label), 'run', fget(q, 'run_on', ''), ...
                fget(q, 'n_rows', NaN), fget(q, 'n_subjects', NaN), ...
                strjoin(fget(q, 'dropped_subjects', {}), ', '), ...
                numel(q.correctTFCEsigch), numel(q.SnPMsigch)}; %#ok<AGROW>
        end
    end
    if isfield(e, 'posthoc_skipped')
        for i = 1:numel(e.posthoc_skipped)
            q = e.posthoc_skipped(i);
            rows(end+1, :) = {name, char(q.label), ['not run: ' char(q.reason)], '', NaN, NaN, '', NaN, NaN}; %#ok<AGROW>
        end
    end
    % mixed2way: interaction simple effects that were not run (a group with
    % fewer than 2 subjects). Listed here with the skipped pairs so the
    % workbook states every test the report does not contain.
    if isfield(e, 'simple_skipped')
        for i = 1:numel(e.simple_skipped)
            q = e.simple_skipped(i);
            rows(end+1, :) = {name, char(q.label), ['not run: ' char(q.reason)], '', NaN, NaN, '', NaN, NaN}; %#ok<AGROW>
        end
    end
end

function v = fget(s, f, dflt)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
