function info = snpm_circ_tier1_guard(params)
%SNPM_CIRC_TIER1_GUARD  Refuse an unadjusted analysis of the Tier-1 UNSIGNED
%   circular measure.
%
%   info = snpm_circ_tier1_guard(params)
%
% snpm_circ_linearise writes a one-line provenance stamp into a reserved
% 'circ_provenance' column of the file it produces (grammar in that file's
% header). This guard decodes the stamp and, when it says the data are
% circ_unsigned_distance, insists that the chosen analysis actually adjusts for
% the event count.
%
% WHY ONLY THE UNSIGNED MEASURE. |deviation from the up-state peak| is folded
% and cannot go below zero, so estimation noise pushes it UP. A subject whose
% preferred phase came from few events has a noisier angle and therefore a
% systematically larger unsigned deviation; any group difference in event count
% manufactures a group difference in the measure. Unadjusted, that path runs a
% Type-I error rate of 0.280 against a nominal 0.05 under an event-count
% imbalance -- the worst-calibrated route in the whole design.
% The SIGNED measure is NOT guarded: noise is symmetric about the pooled mean
% direction, so it is already nominal unadjusted (0.050), and forcing the
% covariate on it costs real power (0.827 -> 0.567).
%
% The only presets with a covariate slot are 'ancova' and 'regression'; every
% other preset is rejected outright. Where a slot exists, params.covariate_cols
% must actually carry the event-count column.
%
% INPUT  params (the orchestrator's params struct; unrelated fields ignored)
%        Optional params.circ_count_col names the count covariate column
%        explicitly, for when it is not called something with "count" in it.
% OUTPUT info : struct with .checked (logical), .stamp (struct or []),
%        .files (cellstr actually inspected). Returns .checked=false when no
%        file carries a provenance stamp (the ordinary, non-Tier-1 case).
%
% Errors with identifier core_snpm:circUnsignedNeedsCountCovariate.
%
% See also SNPM_CIRC_LINEARISE, CORE_SNPM_ANALYSIS.

    info = struct('checked', false, 'stamp', [], 'files', {{}});

    files  = {};
    sheets = {};
    for f = {'data_file', 'data1_file', 'data2_file'}
        if isfield(params, f{1}) && ~isempty(params.(f{1}))
            files{end+1}  = params.(f{1}); %#ok<AGROW>
            shf = strrep(f{1}, '_file', '_sheet');
            if isfield(params, shf) && ~isempty(params.(shf))
                sheets{end+1} = params.(shf); %#ok<AGROW>
            else
                sheets{end+1} = 'CSV File'; %#ok<AGROW>
            end
        end
    end

    stamp = [];
    for k = 1:numel(files)
        s = read_stamp(files{k}, sheets{k});
        if ~isempty(s)
            stamp = s;
            info.files{end+1} = files{k};
        end
    end
    if isempty(stamp), return; end

    info.checked = true;
    info.stamp   = stamp;

    if ~isfield(stamp, 'requires_count_covariate') || ...
            ~strcmp(strtrim(stamp.requires_count_covariate), '1')
        return   % signed measure (or a stamp that does not demand adjustment)
    end

    measure = 'circ_unsigned_distance';
    if isfield(stamp, 'measure'), measure = stamp.measure; end
    cfile = 'the event-count column';
    if isfield(stamp, 'count_file') && ~isempty(stamp.count_file) && ...
            ~strcmp(stamp.count_file, 'none')
        cfile = stamp.count_file;
    end

    preset = '';
    if isfield(params, 'comparison'), preset = char(params.comparison); end

    if ~ismember(preset, {'ancova', 'regression'})
        error('core_snpm:circUnsignedNeedsCountCovariate', ...
            ['This file is the Tier-1 UNSIGNED circular measure (%s, from ' ...
             'snpm_circ_linearise) but ''%s'' has no covariate slot, so the ' ...
             'event count cannot be adjusted for. |deviation from the up-state ' ...
             'peak| is folded at zero: a subject with fewer events has a noisier ' ...
             'angle and therefore a systematically LARGER value, so any group ' ...
             'difference in event count manufactures a group difference in the ' ...
             'measure (Type-I rate 0.280 against a nominal 0.05 under an event- ' ...
             'count imbalance). Use one of the TWO legal presets, ''ancova'' or ' ...
             '''regression'', with the event count (%s) in params.covariate_cols. ' ...
             'Alternatively re-run snpm_circ_linearise with measure = ''signed'', ' ...
             'which is symmetric about the pooled mean direction and needs no ' ...
             'count adjustment.'], measure, preset, cfile);
    end

    covs = {};
    if isfield(params, 'covariate_cols') && ~isempty(params.covariate_cols)
        covs = cellstr(string(params.covariate_cols));
    end
    if ~has_count_covariate(covs, params, stamp)
        error('core_snpm:circUnsignedNeedsCountCovariate', ...
            ['This file is the Tier-1 UNSIGNED circular measure (%s, from ' ...
             'snpm_circ_linearise). Preset ''%s'' does have a covariate slot, but ' ...
             'params.covariate_cols {%s} does not contain the event-count column ' ...
             '(expected something matching ''%s'', or set params.circ_count_col ' ...
             'explicitly). Without it the measure is not calibrated: the folded ' ...
             'deviation rises with estimation noise, so a group difference in ' ...
             'event count alone produces a Type-I rate of 0.280 against a nominal ' ...
             '0.05. The two legal presets are ''ancova'' and ''regression'', both ' ...
             'WITH the count covariate.'], measure, preset, ...
             strjoin(covs, ', '), cfile);
    end
end

% ======================================================================
function tf = has_count_covariate(covs, params, stamp)
    tf = false;
    if isempty(covs), return; end
    if isfield(params, 'circ_count_col') && ~isempty(params.circ_count_col)
        tf = any(strcmpi(covs, char(params.circ_count_col)));
        return
    end
    if any(contains(lower(covs), 'count')), tf = true; return; end
    if isfield(stamp, 'count_file') && ~isempty(stamp.count_file)
        [~, base] = fileparts(stamp.count_file);
        if ~isempty(base) && any(strcmpi(covs, base)), tf = true; end
    end
end

% ======================================================================
function stamp = read_stamp(file, sheet)
% Read ONLY the reserved circ_provenance column, if the file has one.
    stamp = [];
    if ~ischar(file) && ~isstring(file), return; end
    file = char(file);
    if isempty(file) || ~exist(file, 'file'), return; end
    try
        if isempty(sheet) || strcmpi(sheet, 'CSV File')
            o = detectImportOptions(file);
        else
            o = detectImportOptions(file, 'Sheet', sheet);
        end
        if ~ismember('circ_provenance', o.VariableNames), return; end
        o = setvartype(o, {'circ_provenance'}, 'char');
        o.SelectedVariableNames = {'circ_provenance'};
        Tc = readtable(file, o);
    catch ME
        warning('core_snpm:circStampUnreadable', ...
            ['Could not inspect ''%s'' for a circular provenance stamp (%s); ' ...
             'the Tier-1 unsigned-measure guard was NOT applied to it.'], ...
            file, ME.message);
        return
    end
    if isempty(Tc) || height(Tc) == 0, return; end
    s = Tc.circ_provenance;
    if iscell(s), s = s{1}; end
    s = char(string(s));
    if isempty(s), return; end
    stamp = decode_stamp(s);
end

% ======================================================================
function stamp = decode_stamp(s)
% Reference decoder for the stamp grammar documented in snpm_circ_linearise:
% key=value pairs joined by '|', keys/values from [A-Za-z0-9_.:+-].
    stamp = struct();
    kv = strsplit(s, '|');
    for i = 1:numel(kv)
        t = strsplit(kv{i}, '=');
        if numel(t) < 2, continue; end
        key = matlab.lang.makeValidName(strtrim(t{1}));
        stamp.(key) = strtrim(strjoin(t(2:end), '='));
    end
end
