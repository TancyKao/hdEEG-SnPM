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

Chanlabels={};
for i = 1:size(chanlocs,2)
    Chanlabels{i} = chanlocs(i).labels;
end

Chanlabels = string(Chanlabels);

if strcmp(select_mode, 'positional') && length(T_real_T.Properties.VariableNames) == 256
    Chanlabels(strcmp(Chanlabels, 'Cz')) = []; % for UWM
else
    T_real_T.Properties.VariableNames = Chanlabels;
end
    

if  length({T_p_real})== length(Chanlabels)
    T_p_real.Properties.VariableNames = Chanlabels;
end
T_p_correctedTFCE.Properties.VariableNames = Chanlabels;

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
    T_effN = array2table(vertcat(effN_rows{:}));
    if width(T_effN) == length(Chanlabels)
        T_effN.Properties.VariableNames = Chanlabels;
    end
    T_effN = addvars(T_effN, string(effN_metric(:)), 'Before', 1, ...
        'NewVariableNames', 'Metric');
    writetable(T_effN, output_xls, 'Sheet', 'effectiveN')
end


% Channels dropped before the permutation loop because missing data made them
% unevaluable in some permutations (GLM tier, snpm_glm_fl_context). Always
% written when the field exists so a reader of the workbook can see that the
% analysed channel set is smaller than the montage - and see exactly which.
if isfield(S, 'excluded_channels') && isstruct(S.excluded_channels)
    ex = S.excluded_channels;
    if isfield(ex, 'n') && ex.n > 0
        T_excl = table(ex.index(:), string(ex.labels(:)), ...
            repmat(string(ex.reason), numel(ex.index), 1), ...
            'VariableNames', {'ChannelIndex','ChannelLabel','Reason'});
    else
        T_excl = table("None", "All channels evaluable (no missing data)", ...
            'VariableNames', {'ChannelLabel','Reason'});
    end
    writetable(T_excl, output_xls, 'Sheet', 'excludedChannels')
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
