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
