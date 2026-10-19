function txt = snpm_lmm_dropped_text(d)
% One-line summary of results_struct.lmm.dropped_rows (snpm_lmm_drop_missing_meta),
% shared by the results text, the HTML model card and the Excel summary sheet.
    if isempty(d) || ~isstruct(d) || ~isfield(d, 'n') || d.n == 0
        txt = 'Rows dropped for a missing model variable: none';
        return
    end
    parts = cell(1, 0);
    for k = 1:numel(d.columns)
        if d.per_column(k) > 0
            parts{end+1} = sprintf('%s %d', d.columns{k}, d.per_column(k)); %#ok<AGROW>
        end
    end
    txt = sprintf('Rows dropped for a missing model variable: %d of %d (%s)', ...
        d.n, d.n_input, strjoin(parts, ', '));
    if ~isempty(d.subjects_lost)
        txt = sprintf('%s; subjects lost entirely: %s', txt, strjoin(d.subjects_lost, ', '));
    end
end
