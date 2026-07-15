function out_csv = export_results_csv(results_text, outputSname, params)
%EXPORT_RESULTS_CSV  Save the GUI "Results" section to a CSV next to the HTML report.
%
%   out_csv = export_results_csv(results_text, outputSname, params)
%
%   results_text is the same cell array of summary lines that the toolbox shows
%   in the GUI Results panel (and returns as the second output of
%   core_snpm_analysis / core_snpm_glm / core_snpm_lmm). It is written to
%   [outputSname '_report.csv'] so the file name matches the HTML report
%   [outputSname '_report.html'].
%
%   Two CSV-only edits versus the on-screen panel (params drives the first):
%     * a "Data files" row (input file name(s)) is inserted near the top;
%     * the "--- Statistical Results ---" block is dropped (those counts live in
%       the HTML report / .xlsx).
%
%   The CSV has two columns:
%     field - text before the first ':' on a line (or the whole line if none)
%     value - text after the first ':' (empty when there is no ':')

    if nargin < 3, params = struct(); end

    out_csv = [outputSname '_report.csv'];

    lines = normalize_lines(results_text);
    lines = drop_statistical_results(lines);
    lines = insert_data_files(lines, params);

    field = cell(numel(lines), 1);
    value = cell(numel(lines), 1);
    for i = 1:numel(lines)
        ln = lines{i};
        ci = find(ln == ':', 1, 'first');
        if isempty(ci)
            field{i} = strtrim(ln);
            value{i} = '';
        else
            field{i} = strtrim(ln(1:ci-1));
            value{i} = strtrim(ln(ci+1:end));
        end
    end

    T = table(string(field), string(value), 'VariableNames', {'field', 'value'});
    writetable(T, out_csv);
end

function lines = insert_data_files(lines, params)
%INSERT_DATA_FILES  Add a "Data files: ..." row after the title line.
    df = data_files_str(params);
    if isempty(df), return; end
    row = {['Data files: ' df]};
    if isempty(lines)
        lines = row;
    else
        lines = [lines(1); row; lines(2:end)];   % after the title line
    end
end

function s = data_files_str(params)
%DATA_FILES_STR  Input file name(s) from params, whichever pipeline set them.
    names = {};
    if isfield(params, 'data_file')  && ~isempty(params.data_file),  names{end+1} = base_name(params.data_file);  end
    if isfield(params, 'data1_file') && ~isempty(params.data1_file), names{end+1} = base_name(params.data1_file); end
    if isfield(params, 'data2_file') && ~isempty(params.data2_file), names{end+1} = base_name(params.data2_file); end
    names = names(~cellfun(@isempty, names));
    s = strjoin(names, '; ');
end

function b = base_name(p)
    p = char(string(p));
    [~, n, e] = fileparts(p);
    b = [n e];
end

function lines = drop_statistical_results(lines)
%DROP_STATISTICAL_RESULTS  Remove the "Statistical Results" header + count rows.
    drop = false(numel(lines), 1);
    prefixes = { ...
        'Uncorrected significant channels', ...
        'TFCE corrected significant channels', ...
        'Significant clusters N', ...
        'SnPM significant channels'};
    for i = 1:numel(lines)
        t = strtrim(lines{i});
        if strcmp(t, '--- Statistical Results ---')
            drop(i) = true; continue;
        end
        for k = 1:numel(prefixes)
            if startsWith(t, prefixes{k}), drop(i) = true; break; end
        end
    end
    lines = lines(~drop);
    lines = collapse_blank_runs(lines);
end

function out = collapse_blank_runs(lines)
%COLLAPSE_BLANK_RUNS  Squeeze consecutive blank lines to a single blank.
    out = cell(0, 1);
    prevBlank = false;
    for i = 1:numel(lines)
        isBlank = isempty(strtrim(lines{i}));
        if isBlank && prevBlank, continue; end
        out{end+1, 1} = lines{i}; %#ok<AGROW>
        prevBlank = isBlank;
    end
end

function lines = normalize_lines(results_text)
%NORMALIZE_LINES  Flatten results_text into a column cellstr of char lines.
    if isstring(results_text)
        results_text = cellstr(results_text);
    elseif ischar(results_text)
        results_text = cellstr(results_text);  % splits a char matrix by row
    elseif ~iscell(results_text)
        results_text = {char(string(results_text))};
    end

    results_text = results_text(:);
    lines = cell(0, 1);
    for i = 1:numel(results_text)
        el = results_text{i};
        if isstring(el) || (ischar(el) && size(el,1) > 1)
            el = cellstr(el);
        end
        if iscell(el)
            for j = 1:numel(el)
                lines{end+1, 1} = to_char(el{j}); %#ok<AGROW>
            end
        else
            lines{end+1, 1} = to_char(el); %#ok<AGROW>
        end
    end
end

function s = to_char(x)
    if ischar(x)
        s = x;
    elseif isstring(x)
        s = char(x);
    else
        s = char(string(x));
    end
    s = s(:).';  % force a row
end
