function [fig, G] = plot_global_spectrum(root, opts)
% Global power-spectrum comparison: mean PSD across all channels, averaged
% within each group (condition by default), plotted overlaid on a log power
% axis vs frequency. Reproduces the "Spectral density vs Frequency(Hz)"
% comparison figure for the the example study hd-EEG spectral dataset.
%
% For each file it reads EEG.data (nCh x nFreq) and EEG.freqs, averages PSD
% over channels -> one global spectrum per (subject, group). Subject spectra
% are then averaged within group and plotted (mean line; optional +/- SEM band).
%
% INPUTS
%   root : dataset directory containing condition-* folders
%   opts : struct
%       .groupby    'condition' (default) | 'stage'  -> what defines the lines
%       .stages     (cellstr) restrict stages (default all)
%       .conditions (cellstr) restrict conditions (default all)
%       .fmax       upper frequency for the x-axis (default 45 Hz)
%       .logy       true (default) -> log10 y-axis
%       .sem        true (default) -> shaded +/- SEM band
%       .savepath   if set, save PNG here
%       .title      figure title
%
% OUTPUTS
%   fig : figure handle
%   G   : struct array per group with .name, .freqs, .mean, .sem, .nsubj

    if nargin < 2, opts = struct(); end
    groupby = getdef(opts, 'groupby', 'condition');
    fmax    = getdef(opts, 'fmax', 45);
    logy    = getdef(opts, 'logy', true);
    showsem = getdef(opts, 'sem', true);

    d = dir(fullfile(root, 'condition-*'));
    conds = {d([d.isdir]).name};
    if isfield(opts,'conditions') && ~isempty(opts.conditions), conds = opts.conditions; end

    % accumulate: group -> subject -> mean spectrum
    acc = containers.Map('KeyType','char','ValueType','any');
    freqs = [];

    for ci = 1:numel(conds)
        files = dir(fullfile(root, conds{ci}, '*_powerspect.mat'));
        for fi = 1:numel(files)
            meta = parse_bids_name(files(fi).name);
            if isempty(meta), continue; end
            if isfield(opts,'stages') && ~isempty(opts.stages) && ~ismember(meta.stage, opts.stages), continue; end

            S = load(fullfile(root, conds{ci}, files(fi).name)); EEG = S.EEG;
            if isempty(freqs), freqs = double(EEG.freqs(:))'; end
            gspec = mean(double(EEG.data), 1, 'omitnan');   % 1 x nFreq, global mean PSD

            switch groupby
                case 'condition', gkey = meta.cond;
                case 'stage',     gkey = meta.stage;
                otherwise,        gkey = meta.cond;
            end
            skey = meta.subject;
            if ~isKey(acc, gkey), acc(gkey) = containers.Map('KeyType','char','ValueType','any'); end
            sub_map = acc(gkey);
            if isKey(sub_map, skey)   % average multiple files (e.g. stages) for that subject
                sub_map(skey) = [sub_map(skey); gspec];
            else
                sub_map(skey) = gspec;
            end
            acc(gkey) = sub_map;
        end
    end

    assert(~isempty(freqs), 'No files found under %s', root);
    fmask = freqs <= fmax;
    f = freqs(fmask);

    fig = figure('Color','w','Position',[60 60 720 420]); hold on;
    colors = lines(max(acc.Count,1));
    gnames = keys(acc); G = struct('name',{},'freqs',{},'mean',{},'sem',{},'nsubj',{});
    h = [];
    for gi = 1:numel(gnames)
        sub_map = acc(gnames{gi});
        subj = keys(sub_map);
        M = zeros(numel(subj), numel(freqs));
        for si = 1:numel(subj)
            sp = sub_map(subj{si});
            M(si, :) = mean(sp, 1, 'omitnan');   % collapse multiple files per subject
        end
        mu = mean(M, 1, 'omitnan'); se = std(M, 0, 1, 'omitnan') ./ sqrt(size(M,1));
        mu = mu(fmask); se = se(fmask);
        if showsem
            fill([f fliplr(f)], [mu+se fliplr(mu-se)], colors(gi,:), ...
                'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
        end
        h(gi) = plot(f, mu, 'Color', colors(gi,:), 'LineWidth', 2); %#ok<AGROW>
        G(gi) = struct('name',gnames{gi},'freqs',f,'mean',mu,'sem',se,'nsubj',numel(subj));
    end
    if logy, set(gca,'YScale','log'); end
    xlabel('Frequency (Hz)'); ylabel('Spectral density');
    legend(h, gnames, 'Interpreter','none', 'Location','northeast');
    if isfield(opts,'title') && ~isempty(opts.title), title(opts.title, 'Interpreter','none'); end
    box off; xlim([min(f) fmax]);

    if isfield(opts,'savepath') && ~isempty(opts.savepath)
        if ~exist(opts.savepath,'dir'), mkdir(opts.savepath); end
        fn = fullfile(opts.savepath, 'global_spectrum.png');
        print(fig, fn, '-dpng', '-r150');
        fprintf('Saved %s\n', fn);
    end
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d), if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end; end

function meta = parse_bids_name(name)
    meta = [];
    sub  = regexp(name, 'sub-([^_]+)', 'tokens', 'once');
    cond = regexp(name, 'condition-([^_]+)', 'tokens', 'once');
    stg  = regexp(name, 'desc-([^_]+)_powerspect', 'tokens', 'once');
    if isempty(sub) || isempty(stg), return; end
    meta.subject = sub{1};
    if ~isempty(cond), meta.cond = cond{1}; else, meta.cond = ''; end
    meta.stage = stg{1};
end
