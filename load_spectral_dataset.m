function [Tbl, info] = load_spectral_dataset(root, opts)
% Assemble a per-channel band-power analysis table from BIDS-style hd-EEG
% spectral-power derivative folders (EEG_processor output):
%   <folder>/sub-XX_..._run-N_desc-<stage>_powerspect.mat
% Each .mat holds an EEG struct whose EEG.features (or legacy EEG.bands) carries
% precomputed per-channel band power for 'absolute power'/'normalized power'.
%
% PRIMARY model (opts.folders): one folder per cohort/condition; the design
% factor 'level' = the folder's label. The same subject appearing in two folders
% (matched by the sub-XX token) gives a within-subject contrast; distinct
% subjects per folder give a between-group contrast. FALLBACK (no opts.folders):
% a single root containing condition-* subfolders (level = the part after
% 'condition-'). Pick a band + type; the loader extracts that channel vector
% from every file and assembles a long table the GLM/t-test path consumes.
%
% INPUTS
%   root : (fallback only) a directory containing condition-* subfolders.
%          Ignored when opts.folders is given (pass '' then).
%   opts : struct
%       .folders    (cellstr) one folder per factor level (PRIMARY model)
%       .labels     (cellstr) level label per folder (default = folder basename)
%       .band       (char)  band label to extract, e.g. 'delta','sigma' (default 'delta')
%       .type       (char)  'absolute' | 'normalized' (default 'absolute')
%       .conditions (cellstr) fallback: condition-* subfolder names to include
%       .stages     (cellstr) desc stages to include, e.g. {'n2','n3'} (default = all)
%       .subjects   (cellstr) restrict to these subject ids (default = all)
%       .subjects_csv (char) optional subject-metadata CSV to left-join on Subject
%
% OUTPUTS
%   Tbl  : table, one row per file: Subject, level, stage, run, + channel
%          columns (chanlocs labels E1..E178)
%   info : struct with .channel_cols, .chanlabels, .band, .type, .n_files,
%          .levels, .meta_cols (for passing to core_snpm_analysis)
%
% Then e.g.:
%   [T,info] = load_spectral_dataset(root, struct('band','delta','type','absolute'));
%   writetable(T, 'analysis.csv');
%   params.comparison='rmanova'; params.meta_cols=info.meta_cols;
%   params.subject_col='Subject'; params.condition_col='stage'; ...

    if nargin < 2, opts = struct(); end
    band  = getdef(opts, 'band', 'delta');
    type  = getdef(opts, 'type', 'absolute');

    % Build the list of source folders, each tagged with a factor-level label.
    % Primary model: opts.folders = {pathA, pathB, ...} -- one folder per cohort/
    % condition (EEG_processor output). The factor 'level' = the folder's label
    % (opts.labels, default = folder basename). Fallback: a single root containing
    % condition-* subfolders (label = the part after 'condition-').
    if isfield(opts,'folders') && ~isempty(opts.folders)
        folders = cellstr(opts.folders);
        if isfield(opts,'labels') && ~isempty(opts.labels)
            labels = cellstr(opts.labels);
            assert(numel(labels)==numel(folders), 'opts.labels must match opts.folders in length');
        else
            labels = cellfun(@folder_basename, folders, 'UniformOutput', false);
        end
    else
        if isfield(opts,'conditions') && ~isempty(opts.conditions)
            sub = cellstr(opts.conditions);
        else
            d = dir(fullfile(root, 'condition-*')); sub = {d([d.isdir]).name};
        end
        assert(~isempty(sub), 'No condition-* folders found under %s (or pass opts.folders).', root);
        folders = fullfile(root, sub);
        labels  = regexprep(sub, '^condition-', '');   % 'condition-a' -> 'a'
    end

    rows = {};            % accumulate row structs
    chanlabels = {};
    nfiles = 0;

    for ci = 1:numel(folders)
        cdir = folders{ci};
        files = dir(fullfile(cdir, '*powerspect*.mat'));   % PSG _powerspect.mat AND KDT _powerspect_final.mat
        for fi = 1:numel(files)
            meta = parse_bids_name(files(fi).name);
            if isempty(meta), continue; end
            if isfield(opts,'stages') && ~isempty(opts.stages) && ~ismember(meta.stage, opts.stages), continue; end
            if isfield(opts,'subjects') && ~isempty(opts.subjects) && ~ismember(meta.subject, opts.subjects), continue; end

            S = load(fullfile(cdir, files(fi).name));
            if ~isfield(S,'EEG') || isempty(S.EEG)
                warning('load_spectral_dataset:noEEG', 'Skipping %s (no EEG struct).', files(fi).name);
                continue;
            end
            EEG = S.EEG;
            vec = extract_band(EEG, band, type);    % 1 x nCh
            if isempty(vec), continue; end

            if isempty(chanlabels)
                chanlabels = {EEG.chanlocs.labels};
            end

            r = struct();
            r.Subject = string(meta.subject);
            r.level   = string(labels{ci});      % factor level = folder label
            r.stage   = string(meta.stage);
            r.run     = string(meta.run);
            r.vec     = vec;
            rows{end+1} = r; %#ok<AGROW>
            nfiles = nfiles + 1;
        end
    end

    assert(nfiles > 0, 'No matching files (band=%s type=%s).', band, type);

    Subject = arrayfun(@(k) rows{k}.Subject, 1:nfiles)';
    level   = arrayfun(@(k) rows{k}.level,   1:nfiles)';
    stage   = arrayfun(@(k) rows{k}.stage,   1:nfiles)';
    P = zeros(nfiles, numel(chanlabels));
    for k = 1:nfiles, P(k, :) = rows{k}.vec; end

    % Ignore run (and any other within-cell duplication): collapse rows sharing
    % (Subject, level, stage) by averaging the channel vectors. No-op for PSG
    % (one file per subject/stage); for KDT it averages run-1/run-4.
    key = strcat(Subject, '|', level, '|', stage);
    [~, ia, g] = unique(key, 'stable');
    if numel(ia) < nfiles
        nU = numel(ia); Pa = zeros(nU, size(P,2));
        for u = 1:nU, Pa(u,:) = mean(P(g==u, :), 1, 'omitnan'); end
        Subject = Subject(ia); level = level(ia); stage = stage(ia); P = Pa;
        fprintf('Collapsed %d files -> %d rows (Subject x level x stage), averaging duplicate runs.\n', nfiles, nU);
    end

    meta_tbl = table(Subject, level, stage);
    extra_meta = {};
    if isfield(opts,'subjects_csv') && ~isempty(opts.subjects_csv)
        [meta_tbl, extra_meta] = join_subjects(meta_tbl, opts.subjects_csv);
    end
    chan_tbl = array2table(P, 'VariableNames', matlab.lang.makeValidName(chanlabels));
    Tbl = [meta_tbl, chan_tbl];

    info = struct();
    info.channel_cols = chan_tbl.Properties.VariableNames;
    info.chanlabels   = chanlabels;
    info.band = band; info.type = type;
    info.n_files = nfiles;
    info.meta_cols = [{'Subject','level','stage'}, extra_meta];
    info.levels = cellstr(unique(level))';
    fprintf('Loaded %d files | band=%s type=%s | %d channels | %d subjects | levels=%s | stages=%s\n', ...
        nfiles, band, type, numel(chanlabels), numel(unique(Subject)), ...
        strjoin(cellstr(unique(level)),','), strjoin(cellstr(unique(stage)),','));
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function b = folder_basename(p)
    p = char(p);
    if ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end) = []; end
    [~, n, e] = fileparts(p);
    b = [n e];
end

function meta = parse_bids_name(name)
% Key each file on sub-XXX (subject) and desc-XXX (stage/segment), ignoring the
% variable middle (task-psg / task-kdtcomb, run, condition) and the file suffix.
% PSG stage = n1/n2/n3/rem; KDT segment = eyesopen/eyesclosed (the desc token).
    meta = [];
    sub  = regexp(name, 'sub-([^_]+)', 'tokens', 'once');
    cond = regexp(name, 'condition-([^_]+)', 'tokens', 'once');
    run  = regexp(name, 'run-([^_]+)', 'tokens', 'once');
    stg  = regexp(name, 'desc-([^_.]+)', 'tokens', 'once');   % up to next '_' or '.', suffix-agnostic
    if isempty(sub) || isempty(stg), return; end
    meta.subject = sub{1};
    if ~isempty(cond), meta.cond = cond{1}; else, meta.cond = ''; end
    if ~isempty(run),  meta.run  = run{1};  else, meta.run  = ''; end
    meta.stage = stg{1};
end

function vec = extract_band(EEG, band, type)
% Pull the per-channel band-power vector for (label==band, type) from the band
% struct (EEG.bands or EEG.features). Type is matched tolerantly so 'absolute'
% matches 'absolute'/'absolute power' and 'normalized'/'normalised' match
% 'normalized power' (this dataset stores '<type> power' strings).
    vec = [];
    B = bandsource(EEG);
    if isempty(B), return; end
    labels = {B.label};
    types  = {B.type};
    idx = find(strcmpi(labels, band) & type_matches(types, type), 1);
    if isempty(idx)
        % fall back: match label only, warn about available types
        li = find(strcmpi(labels, band));
        if isempty(li)
            warning('load_spectral_dataset:noBand', 'Band ''%s'' not found. Available: %s', ...
                band, strjoin(unique(labels), ', '));
            return;
        end
        idx = li(1);
        warning('load_spectral_dataset:noType', 'Type ''%s'' not found for band ''%s''; using ''%s''.', ...
            type, band, B(idx).type);
    end
    vec = B(idx).data(:)';   % 1 x nCh
end

function [meta_tbl, extra] = join_subjects(meta_tbl, csv)
% Left-join subject-level metadata (group/covariates/predictor) from a CSV onto
% meta_tbl on a normalized Subject key (lowercase, strip leading 'sub-'). Row
% order of meta_tbl is preserved; unmatched rows get NaN / <missing> / ''.
    S = readtable(csv, 'TextType','string', 'VariableNamingRule','preserve');
    assert(any(strcmpi(S.Properties.VariableNames,'Subject')), ...
        'subjects_csv must have a Subject column: %s', csv);
    keyfun = @(v) regexprep(lower(strtrim(string(v))), '^sub-?', '');
    kL = keyfun(meta_tbl.Subject);
    kR = keyfun(S.Subject);
    extra = setdiff(S.Properties.VariableNames, {'Subject'}, 'stable');
    [tf, loc] = ismember(kL, kR);
    out = S(max(loc,1), extra);                 % gather (unmatched temporarily -> row 1)
    if any(~tf)
        for j = 1:numel(extra)
            c = out.(extra{j});
            if     isnumeric(c),     c(~tf) = NaN;
            elseif isstring(c),      c(~tf) = missing;
            elseif iscategorical(c), c(~tf) = missing;
            elseif iscell(c),        c(~tf) = {''};
            end
            out.(extra{j}) = c;
        end
        warning('load_spectral_dataset:subjUnmatched', ...
            '%d of %d rows had no subjects_csv match.', sum(~tf), numel(tf));
    end
    meta_tbl = [meta_tbl, out];
end

function B = bandsource(EEG)
% Per-channel band-power struct array (label/type/freqrange/data). Prefer
% EEG.bands; fall back to EEG.features (some the example study exports use that name).
    B = [];
    if isfield(EEG,'bands') && ~isempty(EEG.bands), B = EEG.bands;
    elseif isfield(EEG,'features') && ~isempty(EEG.features), B = EEG.features; end
end

function tf = type_matches(types, want)
% Compare ignoring a trailing ' power' and unifying normalised/normalized.
    norm = @(s) strrep(regexprep(lower(strtrim(s)), '\s*power$', ''), 'normalised', 'normalized');
    w = norm(want);
    tf = cellfun(@(t) strcmp(norm(t), w), types);
end
