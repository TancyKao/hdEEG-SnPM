function generateAnalysisReport(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch)
%GENERATEANALYSISREPORT  One self-contained HTML report per analysis (design v2).
%
% Writes <outputSname>_report.html for ONE run of any tier that reaches it:
% the legacy t-test / correlation flow (core_snpm_analysis), the GLM presets
% and the three-effect mixed2way report (core_snpm_glm), and the circular tier
% (core_snpm_circ). One layout, one stylesheet, one behaviour, whatever the
% statistic.
%
% LAYOUT is the approved design mock templates/snpm_report_v2_master.html:
% masthead -> sticky per-effect summary strip -> model card -> amber notice
% block (only when something needs saying) -> legend -> descriptive figures ->
% one collapsible .effect section per effect (effect map + Cluster / TFCE /
% Uncorrected toggle + post-hoc tabs) -> global test card -> footer.
% v2_css() and v2_js() below carry lines 5-166 and 167-185 of that file
% BYTE-FOR-BYTE; the mock stays the source of truth for the styling and a test
% diffs the two. Never hand-edit the CSS here -- edit the template and re-port.
%
% STRUCTURE of this file: normalise FIRST (build_model turns whatever the tier
% put in results_struct into one R model: R.effects{k} with a statistic vector,
% p-vectors, clusters, significant-channel sets, PNG names and an optional
% global), render SECOND (write_document reads R only). Adding a tier means
% teaching build_model about it, not touching the emitter.
%
% INPUTS  results_struct, params, base_filename, outputSname,
%         uncorrsigch       channels with uncorrected p <= alpha
%         correctTFCEsigch  channels with TFCE-corrected p <= alpha
%         SnPMsigch         channels in clusters with cluster p <= alpha
%
% ALPHA is the run's own (results_struct.alpha, published by every core, .05 in
% every shipped preset). It drives the cluster numbering here AND was handed to
% TopoplotSignificant_single, so the numbers in the table and on the PNG cannot
% disagree; the prose quotes it rather than a second hard-coded literal.
%
% Every number comes from results_struct; nothing is recomputed here. Cluster
% numbering is delegated to dependencies/snpm_order_clusters.m, the same helper
% TopoplotSignificant_single uses to number the discs it draws on the Cluster
% PNG, so a "Cluster 2" in the table is the "2" on the map.

    % Images are embedded as base64 data URIs by default, so the single .html
    % file survives being moved, emailed or imported into Notion. Set
    % params.report_embed_images = false to go back to relative <img src>
    % filenames (smaller file, PNGs must travel with it). The PNGs are written
    % next to outputSname, so that folder is where the encoder looks.
    embed_dir = fileparts(outputSname);
    embed_on  = true;
    if isstruct(params) && isfield(params, 'report_embed_images') ...
            && ~isempty(params.report_embed_images)
        embed_on = logical(params.report_embed_images(1));
    end
    embed_ctx(embed_dir, embed_on);

    report_filename = [outputSname '_report.html'];
    fid = fopen(report_filename, 'w');
    if fid == -1
        warning('generateAnalysisReport:open', 'Could not open %s for writing.', report_filename);
        return;
    end
    closer = onCleanup(@() fclose(fid));

    R = build_model(results_struct, params, base_filename, outputSname, ...
        uncorrsigch, correctTFCEsigch, SnPMsigch);
    write_document(fid, R);
end

% =========================================================================
% MODEL: results_struct (any tier) -> one R struct the emitter can render
% =========================================================================
function R = build_model(rs, params, base, outputSname, unc, tfce, clu)
    R = struct();
    R.rs    = rs;                       % only for channel labels
    % The RUN'S significance level (every core publishes results_struct.alpha).
    % The cluster table below numbers the clusters with it, and the same value
    % went to TopoplotSignificant_single, so a "3" in the table is the "3" on
    % the PNG. Falling back to .05 only for an older .mat re-rendered here.
    R.alpha = double(gf(rs, 'alpha', 0.05));
    R.nperm = 1000;
    if isfield(params,'permutations') && ~isempty(params.permutations)
        R.nperm = double(params.permutations);
    end
    R.comparison = '';
    if isfield(params,'comparison') && ischar(params.comparison)
        R.comparison = params.comparison;
    end
    % How clusters were scored (GLM / circular / mixed-model tiers only); ''
    % for the legacy t-tests, which use cluster extent, and for results saved
    % before the setting existed.
    R.cluster_stat_html = cluster_stat_html(rs, R.comparison);
    [statName, statSym] = stat_info(rs, params);
    R.statName = statName;
    R.kind = report_kind(rs);

    [~, csvbase] = fileparts(outputSname);
    R.csv = [csvbase '_report.csv'];

    R.title    = analysis_title(rs, params);
    R.methods  = methods_line(rs, params, statName, R.nperm, R.alpha);
    R.model    = model_card_rows(rs, params, base);
    R.notices  = build_notices(rs);
    R.effects  = build_effects(rs, params, base, statName, statSym, unc, tfce, clu);
    R.descriptive = build_descriptive(rs, params, base, statName);

    % Global card. The circular tier deliberately has none (no valid arithmetic
    % channel mean of an angle); a preset that simply does not compute one says
    % so rather than showing an empty card.
    R.global_note = '';
    anyg = false;
    for k = 1:numel(R.effects)
        if ~isempty(R.effects{k}.global), anyg = true; end
    end
    R.has_global = anyg;
    if ~anyg
        if strcmp(R.kind, 'circ')
            R.global_note = ['<b>No whole-head test.</b> An angle has no valid arithmetic ' ...
                'channel mean, so a channel-averaged omnibus would not be the quantity it ' ...
                'appears to be. The circular summary above is descriptive only.'];
        else
            R.global_note = ['<b>No whole-head test for this design.</b> This preset does not ' ...
                'compute a channel-averaged omnibus; read the per-channel maps below.'];
        end
    end
end

function k = report_kind(rs)
    if isfield(rs,'effects') && numel(rs.effects) > 1
        k = 'multi';
    elseif isfield(rs,'circ')
        k = 'circ';
    elseif isfield(rs,'lmm')
        k = 'lmm';
    elseif isfield(rs,'glm') && isfield(rs.glm,'preset')
        k = 'glm';
    else
        k = 'legacy';
    end
end

% ---- one effect model ---------------------------------------------------
function E = build_effects(rs, params, base, statName, statSym, unc, tfce, clu)
    E = {};
    if isfield(rs,'effects') && numel(rs.effects) > 1
        keymap = struct('group','g','condition','c','interaction','i');
        nse = 0;
        for k = 1:numel(rs.effects)
            eff = rs.effects(k);
            key = sprintf('e%d', k);
            if isfield(eff,'key') && isfield(keymap, lower(eff.key)), key = keymap.(lower(eff.key)); end
            E{end+1} = effect_from_glm_effect(rs, eff, key); %#ok<AGROW>
            se = [];
            if isfield(eff,'simple_effects'), se = eff.simple_effects; end
            for s = 1:numel(se)
                nse = nse + 1;
                E{end+1} = effect_from_glm_effect(rs, se(s), sprintf('s%d', nse)); %#ok<AGROW>
            end
        end
        return;
    end

    e = new_effect();
    e.key      = 'e';
    e.title    = effect_heading(rs, params, statName);
    e.sym      = statSym;
    e.signed   = ~any(strcmp(statName, {'F','circ'}));
    e.subtitle = sprintf('%s &middot; %s', stat_word(statSym, e.signed), perm_phrase(rs, params));
    e.T        = statvec(rs);
    e.punc     = pvec(rs, 'real');
    e.ptfce    = pvec(rs, 'correctedTFCE');
    e.Clusters = gf(rs, 'Clusters', []);
    e.sunc     = idxvec(unc); e.stfce = idxvec(tfce); e.sclu = idxvec(clu);
    e.lamp     = ~isempty(e.sclu) || ~isempty(e.stfce);
    e.map_base = base;
    [e.eff_png, e.eff_title] = single_effect_map(rs, params, base, statName);
    e.posthoc  = posthoc_list(gf(rs,'posthoc',[]));
    e.skipped  = skipped_list(gf(rs,'posthoc_skipped',[]));
    e.global   = global_from_legacy(rs, statSym);
    E = {e};
end

function e = effect_from_glm_effect(rs, eff, key)
    e = new_effect();
    e.key    = key;
    e.title  = eff_title(eff);
    sym = 'F'; if isfield(eff,'contrast_type') && strcmpi(eff.contrast_type,'t'), sym = 't'; end
    e.sym    = sym;
    e.signed = strcmp(sym, 't');
    e.subtitle = sprintf('%s &middot; %s permutation', stat_word(sym, e.signed), perm_word(eff));
    e.T      = statvec(eff);
    e.punc   = pvec(eff, 'real');
    e.ptfce  = pvec(eff, 'correctedTFCE');
    e.Clusters = gf(eff, 'Clusters', []);
    e.sunc  = idxvec(gf(eff,'uncorrsigch',[]));
    e.stfce = idxvec(gf(eff,'correctTFCEsigch',[]));
    e.sclu  = idxvec(gf(eff,'SnPMsigch',[]));
    e.lamp  = ~isempty(e.sclu) || ~isempty(e.stfce);
    e.map_base = gf(eff,'map_base','');
    if strcmp(sym,'F')
        e.eff_png = gf(eff,'fmap_png','');
        e.eff_title = 'Effect map (F)';
    else
        e.eff_png = gf(eff,'diff_png','');
        e.eff_title = 'Effect map (t)';
    end
    e.posthoc = posthoc_list(gf(eff,'posthoc',[]));
    e.skipped = skipped_list(gf(eff,'posthoc_skipped',[]));
    e.excl_html = excluded_effect_html(rs, eff);
    % mixed2way condition main effect: say what it averages over. With unequal
    % groups the equal-weight average is not the pooled one, and a reader who
    % expects the pooled effect would otherwise misread the map.
    if isfield(eff,'key') && strcmpi(char(eff.key), 'condition') && isfield(rs,'glm') ...
            && isfield(rs.glm,'preset') && strcmpi(char(rs.glm.preset), 'mixed2way')
        e.note_html = sprintf(['    <p class="methods" style="margin:14px 0 0">With unequal group ' ...
            'sizes the condition main effect is the average across groups with equal weight ' ...
            'per group, as in a standard mixed analysis of variance.</p>\n']);
    end
    % mixed2way interaction: simple effects that were NOT run (a group with
    % fewer than 2 subjects). Said here, next to the simple-effect sections,
    % so a missing section is not read as "nothing significant".
    ssk = gf(eff,'simple_skipped',[]);
    for i = 1:numel(ssk)
        e.note_html = [e.note_html sprintf(['    <p class="methods" style="margin:14px 0 0;' ...
            'color:#b8862a"><b>%s was not run:</b> %s.</p>\n'], ...
            esc(char(ssk(i).label)), esc(char(ssk(i).reason)))];
    end
    e.global = [];
    % A global the engine declared undefined (round-4 spec 3a: a perfectly
    % fitted whole-head average, 'perfect fit, statistic undefined') keeps its
    % card: NaN with the reason, so a missing number is never read as "no
    % whole-head test was run".
    greason = char(gf(eff,'global_reason',''));
    if isfield(eff,'global_stat') && ~isempty(eff.global_stat) ...
            && (~isnan(eff.global_stat) || ~isempty(greason))
        g = struct('label', eff_title(eff), 'sym', sym, 'stat', eff.global_stat, ...
            'p', gf(eff,'global_pval',NaN), 'nch', gf(eff,'global_n_channels',NaN), ...
            'df1', gf(eff,'global_df1',NaN), 'df2', gf(eff,'global_df2',NaN), 'flag', '', ...
            'reason', greason);
        e.global = g;
    end
end

function e = new_effect()
    e = struct('key','', 'title','', 'subtitle','', 'sym','t', 'signed',true, ...
        'T',[], 'punc',[], 'ptfce',[], 'Clusters',[], 'sunc',[], 'stfce',[], 'sclu',[], ...
        'lamp',false, 'map_base','', 'eff_png','', 'eff_title','Effect map', ...
        'posthoc',{{}}, 'skipped',{{}}, 'global',[], 'excl_html','', 'note_html','');
end

% Pairwise comparisons that were NOT run: {label, reason} per entry.
function K = skipped_list(sk)
    K = {};
    for i = 1:numel(sk)
        K{end+1} = struct('label', char(gf(sk(i),'label','')), ...
            'reason', char(gf(sk(i),'reason',''))); %#ok<AGROW>
    end
end

function P = posthoc_list(ph)
    P = {};
    for i = 1:numel(ph)
        q = struct();
        q.label = gf(ph(i),'label', sprintf('Contrast %d', i));
        q.T     = statvec(ph(i));
        q.punc  = pvec(ph(i),'real');
        q.ptfce = pvec(ph(i),'correctedTFCE');
        q.Clusters = gf(ph(i),'Clusters',[]);
        q.sunc  = idxvec(gf(ph(i),'uncorrsigch',[]));
        q.stfce = idxvec(gf(ph(i),'correctTFCEsigch',[]));
        q.sclu  = idxvec(gf(ph(i),'SnPMsigch',[]));
        q.map_base = gf(ph(i),'map_base','');
        % What the pair was run on: its own rows, which can be fewer subjects
        % than the effect's header states.
        q.run_on  = char(gf(ph(i),'run_on',''));
        q.dropped = gf(ph(i),'dropped_subjects',{});
        if ~iscell(q.dropped), q.dropped = cellstr(string(q.dropped)); end
        q.name  = tab_name(q.label, i);
        % This pair's own exclusions (zero-variance rule 4: a pair can exclude
        % a channel its omnibus tested).
        q.excl_html = '';
        ex = gf(ph(i),'excluded_channels',[]);
        if isstruct(ex) && isfield(ex,'n') && ex.n > 0
            q.excl_html = excl_reason_html(ex);
        end
        P{end+1} = q; %#ok<AGROW>
    end
end

function g = global_from_legacy(rs, statSym)
    g = [];
    hasS = isfield(rs,'global_stat') && ~isempty(rs.global_stat);
    hasP = isfield(rs,'global_pval') && ~isempty(rs.global_pval);
    if ~hasS && ~hasP, return; end
    g = struct('label','Whole head', 'sym', statSym, 'stat', NaN, 'p', NaN, ...
        'nch', NaN, 'df1', NaN, 'df2', NaN, 'flag', '', 'reason', '');
    if hasS, g.stat = rs.global_stat(1); end
    if hasP, g.p    = rs.global_pval(1); end
    if isfield(rs,'global_n_channels') && ~isempty(rs.global_n_channels)
        g.nch = rs.global_n_channels;
    end
    % A global the engine itself declared uninterpretable (correlation +
    % datatype 'normalize': the channel average of a per-subject z-scored row is
    % analytically zero, so the statistic is rounding noise) must never be
    % painted as a significant headline number.
    if isfield(rs,'global_not_interpretable') && rs.global_not_interpretable
        g.flag = 'not interpretable';
    end
end

% ---- descriptive block --------------------------------------------------
% kind: 'cells' (mixed2way grid + interaction), 'scatter' (correlation),
% 'means' (A/B/effect topos or per-level means), 'circ', 'none'.
function D = build_descriptive(rs, params, base, statName)
    D = struct('kind','none', 'title','Descriptive means', 'sub','', ...
        'figs',{{}}, 'rows',{{}}, 'cols',{{}}, 'cells',{{}}, 'note','');
    % Source space has no scalp layout, so the run deliberately writes no
    % images. Listing figure cards would fill the report with grey "No data"
    % slots for files nobody intended to produce; the notice block says so and
    % the significant-voxel CSV is the deliverable.
    if isfield(rs,'is_source') && rs.is_source, return; end
    d = gf(rs,'descriptive',[]);

    if isfield(rs,'effects') && numel(rs.effects) > 1 && ~isempty(d)
        D.kind = 'cells';
        D.sub  = 'Cell means per design cell and the interaction plot';
        D.rows = gf(d,'cell_rows',{});
        D.cols = gf(d,'cell_cols',{});
        pngs = gf(d,'cellmean_png',{}); labs = gf(d,'cellmean_labels',{});
        D.cells = {pngs, labs};
        ip = gf(d,'interaction_png','');
        if ~isempty(ip)
            nch = numel(gf(d,'interaction_channels',[]));
            D.figs = {fig_spec(ip, 'Interaction plot', ...
                sprintf(['Channel-averaged signal per cell, mean &plusmn; SE, one line per ' ...
                'group &mdash; averaged over <b>%d</b> channel(s): the interaction cluster, ' ...
                'else the TFCE-significant set, else all channels.'], nch), 'wide')};
        end
        return;
    end

    if has(lower(base_comparison(params)), 'correlation')
        D.kind = 'scatter';
        D.title = 'Descriptive';
        D.sub   = 'The association behind the map';
        sp = ''; nch = 0; chlabel = ''; is_rank = false;
        if ~isempty(d)
            sp = gf(d,'scatter_png','');
            sch = gf(d,'scatter_channels',[]);
            nch = numel(sch);
            chlabel = gf(d,'scatter_channel_set','');
            is_rank = isfield(d,'scatter_is_rank') && d.scatter_is_rank;
        end
        ylab = 'Data 2';
        if isfield(rs,'behavior_col') && ~isempty(rs.behavior_col), ylab = char(rs.behavior_col); end
        skipmsg = gf(d,'scatter_skipped_reason','');
        % The scatter is absent by design when an axis had no usable spread
        % (datatype 'normalize' flattens the channel mean). Say so IN THIS CARD:
        % without it the section shows only mean topographies under the heading
        % "The association behind the map" and a reader of the HTML alone cannot
        % tell a suppressed figure from one that failed to render. The engine's
        % own sentence is reused verbatim so console and report never diverge.
        % The 'Data Type caveat' notice does not cover this -- it only fires when
        % Data 2 is a single broadcast behavioural score, not for the
        % montage-vs-montage correlation, which degenerates under 'normalize'
        % just the same.
        if isempty(sp) && ~isempty(skipmsg)
            D.figs{end+1} = notecard_spec(sprintf('Scatter &mdash; EEG vs %s', esc(ylab)), ...
                esc(char(skipmsg)));
        end
        if ~isempty(sp)
            if isempty(chlabel), chlabel = 'the analysed channels'; end
            if is_rank
                % Spearman: the engine correlates RANKS, so the plotted axes are
                % ranks and the fitted line is the relation the map reports.
                % Saying so is the difference between a caption and a claim.
                cap = sprintf(['Per subject: the <b>rank</b> of the signal averaged over ' ...
                    '<b>%d</b> channel(s) &mdash; %s &mdash; against the <b>rank</b> of %s. ' ...
                    'Spearman correlates ranks, so both axes are ranks and the line is the ' ...
                    'fit to them. Descriptive; the inference is the map.'], ...
                    nch, esc(chlabel), esc(ylab));
                ttl = sprintf('Scatter (ranks) &mdash; EEG vs %s', esc(ylab));
            else
                cap = sprintf(['Per subject: the signal averaged over <b>%d</b> channel(s) ' ...
                    '&mdash; %s &mdash; against %s. Descriptive; the inference is the map.'], ...
                    nch, esc(chlabel), esc(ylab));
                ttl = sprintf('Scatter &mdash; EEG vs %s', esc(ylab));
            end
            D.figs{end+1} = fig_spec(sp, ttl, cap, 'wide');
        end
        % Mean topographies. plot_topoInd writes them for the correlation path
        % too (<data1>_topo.png, and <data2>_topo.png unless Data 2 is a
        % broadcast behavioural score); without these cards those PNGs sat in
        % the output folder unreferenced. The r-map is NOT repeated here -- it
        % is the effect map of the effect section below.
        [d1, d2] = ab_basenames(base);
        hideB = isfield(rs,'hide_condition_b') && rs.hide_condition_b;
        if hideB
            D.figs{end+1} = fig_spec([d1 '_topo.png'], 'Mean topography &mdash; EEG', ...
                'Channel means of the EEG measure across subjects.', '');
        else
            D.figs{end+1} = fig_spec([d1 '_topo.png'], 'Mean topography &mdash; Data 1', '', '');
            D.figs{end+1} = fig_spec([d2 '_topo.png'], 'Mean topography &mdash; Data 2', '', '');
        end
        return;
    end

    if isfield(rs,'circ_descriptive') && isstruct(rs.circ_descriptive)
        % Circular tier: the per-group table (rendered from rs.circ_descriptive)
        % plus the mean-DIRECTION topographies. For the two group tests
        % <base>_topo.png is the circ_dist difference of mean directions -- a
        % descriptive map, not the statistic -- so it belongs here, not in the
        % effect section.
        D.kind = 'circ';
        D.title = 'Circular summary';
        D.sub   = 'Descriptive only &mdash; no whole-head test';
        [d1, d2] = ab_basenames(base);
        hideB = isfield(rs,'hide_condition_b') && rs.hide_condition_b;
        D.figs{end+1} = fig_spec([d1 '_topo.png'], 'Mean direction &mdash; A (deg)', '', '');
        if ~hideB
            D.figs{end+1} = fig_spec([d2 '_topo.png'], 'Mean direction &mdash; B (deg)', '', '');
            D.figs{end+1} = fig_spec([base '_topo.png'], ...
                'Mean-direction difference (circ_dist, deg)', '', '');
        end
        return;
    end

    % Per-level means for an omnibus-F GLM preset.
    if isfield(rs,'glm') && isfield(rs.glm,'group_mean_png') && ~isempty(rs.glm.group_mean_png)
        D.kind = 'means';
        D.sub  = 'Per-level mean topographies';
        pngs = rs.glm.group_mean_png; labs = gf(rs.glm,'group_labels',{});
        for i = 1:numel(pngs)
            lab = sprintf('Level %d', i);
            if i <= numel(labs) && ~isempty(labs{i}), lab = labs{i}; end
            D.figs{end+1} = fig_spec(pngs{i}, sprintf('Mean &mdash; %s', esc(lab)), '', '');
        end
        return;
    end

    % Regression / single-map GLM preset: the slope map is the effect map, so
    % there is nothing descriptive to add.
    if isfield(rs,'glm') && isfield(rs.glm,'effect_map_only') && rs.glm.effect_map_only
        return;
    end

    % Legacy t-test layout: condition means (+ the effect map lives in the
    % effect section, so only the means are descriptive).
    if strcmp(statName, 't')
        D.kind = 'means';
        D.sub  = 'Condition / group mean topographies';
        [d1, d2] = ab_basenames(base);
        hideB = isfield(rs,'hide_condition_b') && rs.hide_condition_b;
        if hideB
            D.figs{end+1} = fig_spec([d1 '_topo.png'], 'Mean topography', '', '');
        else
            D.figs{end+1} = fig_spec([d1 '_topo.png'], 'Mean topography &mdash; A', '', '');
            D.figs{end+1} = fig_spec([d2 '_topo.png'], 'Mean topography &mdash; B', '', '');
        end
    end
end

function s = fig_spec(png, title, note, cls)
    s = struct('png', png, 'title', title, 'note', note, 'cls', cls);
end

% A figure-shaped card that carries prose INSTEAD of an image, for a figure the
% run deliberately did not produce. write_fig renders .emptynote in the slot's
% place, so the reason sits exactly where the missing figure would have been
% rather than as a grey "No data" placeholder that reads like a broken link.
function s = notecard_spec(title, text)
    s = struct('png', '', 'title', title, 'note', '', 'cls', 'wide', 'emptynote', text);
end

function c = base_comparison(params)
    c = '';
    if isfield(params,'comparison') && ischar(params.comparison), c = params.comparison; end
end

function [d1, d2] = ab_basenames(base)
    if has(base, ' VS ')
        parts = strsplit(base, ' VS '); d1 = parts{1}; d2 = parts{2};
    else
        d1 = [base '_Data1']; d2 = [base '_Data2'];
    end
end

function [png, ttl] = single_effect_map(rs, params, base, statName)
% The standalone statistic map shown at the top of an effect section.
%   legacy t / correlation : <base>_topo.png (plot_topoInd's third card)
%   GLM omnibus F          : glm.fmap_png
%   circular group tests   : none -- <base>_topo.png is the descriptive
%                            mean-DIRECTION difference, not the statistic
%   circular-linear        : <base>_topo.png IS the F map (fed as 'diffdata')
%   source space           : no scalp layout, so no map at all
    png = ''; ttl = 'Effect map';
    if isfield(rs,'is_source') && rs.is_source, return; end
    if isfield(rs,'glm') && isfield(rs.glm,'fmap_png') && ~isempty(rs.glm.fmap_png)
        png = rs.glm.fmap_png; ttl = 'Effect map (F)'; return;
    end
    comp = base_comparison(params);
    if has(comp, 'circ')
        if isfield(rs,'hide_condition_b') && rs.hide_condition_b
            png = [base '_topo.png']; ttl = 'Effect map (circular-linear F)';
        end
        return;
    end
    png = [base '_topo.png'];
    switch statName
        case 'F', ttl = 'Effect map (F)';
        case 'r', ttl = 'Correlation r-map';
        otherwise, ttl = 'Effect map (t)';
    end
end

% =========================================================================
% STATISTIC RESOLUTION + TITLES
% =========================================================================
function [name, sym] = stat_info(rs, params)
    name = 't';
    if isfield(rs,'glm') && isfield(rs.glm,'contrast_type') && ~isempty(rs.glm.contrast_type)
        if strcmpi(rs.glm.contrast_type, 'F'), name = 'F'; else, name = 't'; end
    elseif isfield(params,'comparison') && ischar(params.comparison)
        c = lower(params.comparison);
        if     has(c,'correlation'), name = 'r';
        elseif has(c,'circ'),        name = 'circ';
        elseif has(c,'anova') || has(c,'ancova') || has(c,'regression') || has(c,'mixed')
            name = 'F';
        end
    end
    switch name
        case 'F',    sym = 'F';
        case 'r',    sym = 'r';
        case 'circ'
            % The circular tier does NOT emit a generic "stat": two of its three
            % analyses report an F (Hotelling T^2 on the (cos,sin) embedding, and
            % the circular-linear model F), and only Watson's secondary test
            % reports U^2. Naming the column honestly matters because the reader
            % has to know which reference distribution the p came from.
            sym = 'F';
            if isfield(params,'comparison') && ischar(params.comparison) && ...
                    strcmp(params.comparison, 'circ_phase_group_u2')
                sym = 'U&sup2;';
            end
        otherwise,   sym = 't';
    end
end

function t = analysis_title(rs, params)
    comp = lower(base_comparison(params));
    switch comp
        case 'pairedt',      t = 'Paired t-test &mdash; condition A vs B';
        case 'unpairedt',    t = 'Independent-groups t-test &mdash; group A vs B';
        case 'onesamplet',   t = 'One-sample t-test &mdash; condition vs 0';
        case 'correlationp', t = 'Pearson correlation';
        case 'correlations', t = 'Spearman correlation';
        case 'anova1',       t = 'One-way ANOVA &mdash; between subjects';
        case 'ancova',       t = 'ANCOVA &mdash; group effect adjusted for covariates';
        case 'regression',   t = 'Regression &mdash; continuous predictor';
        case 'rmanova',      t = 'Repeated-measures ANOVA &mdash; within subject';
        case 'mixed2way',    t = 'Two-way mixed ANOVA &mdash; group &times; condition';
        case 'mixedmodel',   t = 'Linear mixed model';
        case 'circ_phase_group',    t = 'Preferred phase &mdash; Hotelling T&sup2; between groups';
        case 'circ_phase_group_u2', t = 'Preferred phase &mdash; Watson U&sup2; between groups';
        case 'circ_corranglinear',  t = 'Circular-linear association';
        otherwise
            t = 'SnPM analysis';
            if ~isempty(comp), t = ['SnPM analysis &mdash; ' esc(base_comparison(params))]; end
    end
    if any(strcmp(comp, {'correlationp','correlations'}))
        if isfield(rs,'behavior_col') && ~isempty(rs.behavior_col)
            t = [t ' &mdash; EEG vs ' esc(char(rs.behavior_col))];
        else
            t = [t ' &mdash; EEG vs Data 2'];
        end
    end
end

% Heading of the single effect section of a single-effect analysis.
function t = effect_heading(rs, params, statName)
    if isfield(rs,'glm') && isfield(rs.glm,'effect_label') && ~isempty(rs.glm.effect_label)
        % The circular tier parks its whole provenance line in effect_label
        % ("Hotelling T^2 ... | zero convention ... | count covariate: ...").
        % That belongs in the model card, not in a section heading and a summary
        % chip -- keep the part before the first pipe.
        lab = char(rs.glm.effect_label);
        bar = strfind(lab, '|');
        if ~isempty(bar), lab = strtrim(lab(1:bar(1)-1)); end
        t = esc(lab);
        return;
    end
    switch statName
        case 'r',    t = 'Correlation';
        case 'F',    t = 'Omnibus effect';
        case 'circ', t = 'Phase effect';
        otherwise
            comp = lower(base_comparison(params));
            if strcmp(comp,'onesamplet'), t = 'Condition vs zero';
            elseif strcmp(comp,'pairedt'), t = 'Condition difference (A &minus; B)';
            elseif strcmp(comp,'unpairedt'), t = 'Group difference (A &minus; B)';
            else, t = 'Effect';
            end
    end
end

function w = stat_word(sym, signed)
    if signed
        w = sprintf('signed %s', sym);
    else
        w = sprintf('omnibus %s', sym);
    end
end

function w = perm_word(eff)
    pt = ''; if isfield(eff,'perm_type'), pt = eff.perm_type; end
    switch lower(pt)
        case 'free',   w = 'between-subjects';
        case 'within', w = 'within-subject';
        otherwise
            if isempty(pt), w = 'residual'; else, w = pt; end
    end
end

function s = perm_phrase(rs, params)
    comp = lower(base_comparison(params));
    if isfield(rs,'glm') && isfield(rs.glm,'perm_type') && ~isempty(rs.glm.perm_type)
        s = [perm_word(rs.glm) ' permutation'];
        return;
    end
    % 'unpaired' CONTAINS 'paired' -- test it first or every independent-groups
    % run is captioned as a within-subject sign flip.
    if     has(comp,'unpaired'),    s = 'group-relabelling permutation';
    elseif has(comp,'paired') || has(comp,'onesample'), s = 'within-subject sign-flip permutation';
    elseif has(comp,'correlation') && has_fl_covariates(rs)
        s = 'Freedman&ndash;Lane permutation (approximate)';
    elseif has(comp,'correlation'), s = 'subject-relabelling permutation';
    elseif has(comp,'circ'),        s = 'Freedman&ndash;Lane permutation';
    else,                           s = 'permutation';
    end
end

function t = eff_title(eff)
    key = ''; if isfield(eff,'key'), key = eff.key; end
    switch lower(key)
        case 'group',       t = 'Group main effect';
        case 'condition',   t = 'Condition main effect';
        case 'interaction', t = 'Interaction (group &times; condition)';
        otherwise
            lab = '';
            if isfield(eff,'effect_label') && ~isempty(eff.effect_label), lab = char(eff.effect_label); end
            if isempty(lab), lab = char(key); end
            if strncmpi(key, 'condition@', 10)
                t = ['Simple effect: ' esc(lab)];
            else
                t = esc(lab);
            end
    end
end

% =========================================================================
% METHODS SENTENCE (masthead)
% =========================================================================
function s = methods_line(rs, params, statName, nperm, alpha)
    comp = '';
    if isfield(params,'comparison') && ischar(params.comparison), comp = lower(params.comparison); end
    % A MULTI-EFFECT REPORT HAS NO SINGLE STATISTIC. mixed2way carries three
    % effects with three error terms and three permutation schemes; the
    % top-level rs.glm.contrast_type describes only the interaction, so naming
    % it here ("omnibus F-test", or "t-test" when the interaction happens to be
    % a 1-df contrast) mis-describes the other two. Name the DESIGN and let each
    % effect section state its own statistic.
    multi = isfield(rs,'effects') && numel(rs.effects) > 1;
    switch statName
        case 'F',    statword = 'omnibus F-test';
        case 'r',    statword = 'correlation';
        case 'circ', statword = 'circular statistic';
        otherwise,   statword = 't-test';
    end
    if multi, statword = design_word(comp); end
    phrase = 'Permutation';
    if     has(comp,'unpair'),     phrase = 'Independent-groups permutation';
    elseif has(comp,'onesample'),  phrase = 'One-sample permutation';
    elseif has(comp,'pair'),       phrase = 'Paired permutation';
    elseif has(comp,'rmanova'),    phrase = 'Repeated-measures permutation';
    elseif has(comp,'mixed'),      phrase = 'Mixed-design permutation';
    elseif has(comp,'ancova'),     phrase = 'ANCOVA permutation';
    elseif has(comp,'anova'),      phrase = 'One-way permutation';
    elseif has(comp,'regression'), phrase = 'Regression permutation';
    elseif has(comp,'circ'),       phrase = 'Circular permutation';
    elseif has(comp,'correlation') && has_fl_covariates(rs)
        % Never "exact" with covariates: Freedman-Lane permutes residuals,
        % which are only approximately exchangeable.
        phrase = 'Freedman&ndash;Lane partial-correlation permutation (approximate)';
    end

    % TAIL. Only the legacy t-test / correlation paths have one. An F statistic
    % (GLM omnibus, RM-ANOVA, interaction, LMM factor effect) and every circular
    % statistic are non-negative omnibus quantities with no sidedness, so the
    % tail wording is omitted for them rather than guessed; the GLM's 1-row t
    % contrasts are evaluated two-sided (|t| against a max-|t| null in
    % snpm_perm_correction) and read 'two-tailed'.
    %
    % params.tail is VALIDATED against the three real tail keywords rather than
    % trusted, because core_snpm_glm used to park the contrast type ('t'/'F')
    % in it purely to give this line something to print -- 'F' is not 'both',
    % so every GLM report claimed "one-tailed". The overload is gone, but an
    % older results .mat re-rendered through this function can still carry it.
    tail = '';
    if ~multi && ~any(strcmp(statName, {'F', 'circ'}))
        tail = 'two-tailed';
        if isfield(params,'tail') && ~isempty(params.tail) && ischar(params.tail) ...
                && any(strcmpi(params.tail, {'left','right'}))
            tail = 'one-tailed';
        end
    end
    if isempty(tail), tail_clause = ''; else, tail_clause = [tail ', ']; end
    dt = '';
    if isfield(params,'datatype') && ~isempty(params.datatype) && ischar(params.datatype)
        dt = [' of the ' esc(params.datatype) ' data'];
    end

    if multi
        s = sprintf(['Per-channel <b>%s</b>%s; per-effect statistics listed below. ' ...
            '<b>%s permutations</b>, family-wise corrected with <b>TFCE</b> and ' ...
            '<b>cluster-based</b> methods at <b>&alpha; = %s</b>.'], ...
            statword, dt, fmt_int(nperm), fmt_alpha(alpha));
    else
        s = sprintf(['Per-channel <b>%s</b>%s. <b>%s permutations</b>, %sfamily-wise corrected ' ...
            'with <b>TFCE</b> and <b>cluster-based</b> methods at <b>&alpha; = %s</b>.'], ...
            statword, dt, fmt_int(nperm), tail_clause, fmt_alpha(alpha));
    end
    s = [upper(phrase(1)) phrase(2:end) ' test. ' s];

    nc = obs_clause(rs, params);
    if ~isempty(nc), s = [s ' &nbsp;<b>' nc '</b>.']; end
end

% The DESIGN a multi-effect report covers, for the masthead sentence.
function w = design_word(comp)
    switch comp
        case 'mixed2way', w = 'two-way mixed ANOVA';
        otherwise,        w = 'multi-effect design';
    end
end

% OBSERVATIONS, IN THE UNIT THE DESIGN IS ACTUALLY COUNTED IN.
%
% This used to print data1_size(1) whenever the two arms matched, so a balanced
% 16-vs-16 independent-groups test read "16" -- half the sample, and silent
% about the balance that matters most for a between-groups claim. It also
% printed ROWS for the repeated-measures and mixed designs, where a row is a
% subject x condition cell, not a subject.
%
% The rules, one per design family:
%   two independent arms  "16 vs 16 (32 subjects)"  -- always both arms, even
%                         when equal, plus the total
%   paired / one-sample   "24 pairs" / "24 subjects"
%   correlation           "27 subjects"
%   GLM presets           "48 rows, 24 subjects" when rows > subjects, else
%                         "24 subjects"; plus "(12 / 12 per group)" when the
%                         design has a grouping factor
%   LMM                   "620 rows, 24 subjects"
% Everything is read from published fields (rs.data_summary, rs.glm, rs.lmm) --
% nothing is inferred from the map or the channel count.
function nc = obs_clause(rs, params)
    nc = '';
    comp = lower(base_comparison(params));

    % ---- LMM: trial rows nested in subjects ------------------------------
    if isfield(rs,'lmm') && isstruct(rs.lmm) && isfield(rs.lmm,'n_trials') ...
            && ~isempty(rs.lmm.n_trials)
        nsub = gf(rs.lmm, 'n_subjects', []);
        if isempty(nsub) || ~isfinite(nsub)
            nc = sprintf('%s rows', fmt_int(rs.lmm.n_trials));
        else
            nc = sprintf('%s rows, %s subjects', fmt_int(rs.lmm.n_trials), fmt_int(nsub));
        end
        return;
    end

    % ---- GLM presets, including the mixed2way multi-effect report --------
    % Keyed on glm.preset: the circular tier also writes an rs.glm (effect_label
    % only) and must not be counted as a GLM design.
    if isfield(rs,'glm') && isstruct(rs.glm) && isfield(rs.glm,'preset') ...
            && isfield(rs.glm,'n_obs') && ~isempty(rs.glm.n_obs)
        nrow = double(rs.glm.n_obs);
        nsub = gf(rs.glm, 'n_subjects', []);
        if isempty(nsub) || ~isfinite(nsub), nsub = nrow; end
        if nsub < nrow
            nc = sprintf('%s rows, %s subjects', fmt_int(nrow), fmt_int(nsub));
        else
            nc = sprintf('%s subjects', fmt_int(nsub));
        end
        gs = gf(rs.glm, 'group_sizes', []);
        if numel(gs) > 1
            parts = cell(1, numel(gs));
            for i = 1:numel(gs), parts{i} = fmt_int(gs(i)); end
            nc = sprintf('%s (%s per group)', nc, strjoin(parts, ' / '));
        end
        return;
    end

    % ---- legacy t-test / correlation / circular: two data matrices -------
    if ~isfield(rs,'data_summary'), return; end
    ds = rs.data_summary;
    n1 = []; n2 = [];
    if isfield(ds,'data1_size') && ~isempty(ds.data1_size), n1 = ds.data1_size(1); end
    if isfield(ds,'data2_size') && ~isempty(ds.data2_size), n2 = ds.data2_size(1); end
    if isempty(n1), return; end

    two_arm = any(strcmp(comp, {'unpairedt','circ_phase_group','circ_phase_group_u2'}));
    if two_arm && ~isempty(n2) && n2 > 0
        nc = sprintf('%d vs %d (%d subjects)', n1, n2, n1 + n2);
        return;
    end
    if strcmp(comp, 'pairedt')
        nc = sprintf('%d pairs', n1);
        return;
    end
    if any(strcmp(comp, {'onesamplet','correlationp','correlations','circ_corranglinear'}))
        nc = sprintf('%d subjects', n1);
        return;
    end
    % Unknown comparison: state both matrices rather than guess a unit.
    if ~isempty(n2) && n2 > 0 && n2 ~= n1
        nc = sprintf('%d vs %d', n1, n2);
    else
        nc = sprintf('%d', n1);
    end
end

% =========================================================================
% MODEL CARD: the fitted formula + the specification rows
% =========================================================================
function M = model_card_rows(rs, params, base)
    comp = lower(base_comparison(params));
    dv = 'SIGNAL';
    if isfield(rs,'lmm') && isfield(rs.lmm,'dv') && ~isempty(rs.lmm.dv), dv = upper(char(rs.lmm.dv)); end
    v = @(x) ['<span class="v">' x '</span>'];
    r = @(x) ['<span class="r">' x '</span>'];
    gcol = pget(params,'group_col','group');
    ccol = pget(params,'condition_col','condition');
    pcol = pget(params,'predictor_col','predictor');
    covs = covariate_names(params);

    rows = {};
    switch comp
        case 'pairedt'
            f = [v(dv) ' ~ ' r('condition') ' + Error(Subject)'];
            rows{end+1} = {'Design', 'Paired / within subject', true};
        case 'unpairedt'
            f = [v(dv) ' ~ ' r('group')];
            rows{end+1} = {'Design', 'Two independent groups', true};
        case 'onesamplet'
            f = [v(dv) ' ~ ' r('1') '  (tested against 0)'];
            rows{end+1} = {'Design', 'One sample vs zero', true};
        case {'correlationp','correlations'}
            y = 'Data 2';
            if isfield(rs,'behavior_col') && ~isempty(rs.behavior_col), y = char(rs.behavior_col); end
            if strcmp(comp,'correlations'), est = 'Spearman'; else, est = 'Pearson'; end
            f = ['corr( ' v(dv) ' , ' r(esc(y)) ' )'];
            rows{end+1} = {'Design', [est ' correlation, per channel'], true};
            if ~isempty(covs)
                f = [f ' &nbsp;|&nbsp; ' esc(strjoin(covs, ', '))];
                rows{end+1} = {'Partialled out', strjoin(covs, ', '), false};
            end
        case 'anova1'
            f = [v(dv) ' ~ ' r(esc(gcol))];
            rows{end+1} = {'Design', 'One-way, between subjects', true};
        case 'ancova'
            f = [v(dv) ' ~ ' r(esc(gcol)) ' + ' esc(strjoin_or(covs,'covariates'))];
            rows{end+1} = {'Design', 'ANCOVA, between subjects', true};
        case 'regression'
            f = [v(dv) ' ~ ' r(esc(pcol))];
            if ~isempty(covs), f = [f ' + ' esc(strjoin(covs, ' + '))]; end
            rows{end+1} = {'Design', 'Linear regression', true};
        case 'rmanova'
            f = [v(dv) ' ~ ' r(esc(ccol)) ' + Error(Subject)'];
            rows{end+1} = {'Design', 'Repeated measures, within subject', true};
        case 'mixed2way'
            f = [v(dv) ' ~ ' r(esc(gcol)) ' &times; ' r(esc(ccol)) ' + Error(Subject)'];
            rows{end+1} = {'Design', 'Two-way mixed (between &times; within)', true};
        case {'circ_phase_group','circ_phase_group_u2'}
            f = ['( cos ' v('&theta;') ' , sin ' v('&theta;') ' ) ~ ' r('group')];
            rows{end+1} = {'Design', 'Preferred phase, two independent groups', true};
        case 'circ_corranglinear'
            f = [v('&theta;') ' ~ ' r('linear measure') '  (circular-linear)'];
            rows{end+1} = {'Design', 'Circular-linear association, paired by subject', true};
        otherwise
            f = [v(dv) ' ~ ' r('effect')];
    end

    if isfield(rs,'glm') && isfield(rs.glm,'effect_label') && ~isempty(rs.glm.effect_label) ...
            && ~strcmp(comp,'mixed2way')
        rows{end+1} = {'Effect tested', char(rs.glm.effect_label), false};
    end
    if any(strcmp(comp, {'anova1','ancova','mixed2way'})) && ~isempty(gcol)
        rows{end+1} = {'Between-subjects factor', gcol, true};
    end
    if any(strcmp(comp, {'rmanova','mixed2way'})) && ~isempty(ccol)
        rows{end+1} = {'Within-subject factor', ccol, true};
    end
    if strcmp(comp,'regression') && ~isempty(pcol)
        rows{end+1} = {'Predictor', pcol, true};
    end
    if ~isempty(covs) && ~any(strcmp(comp, {'correlationp','correlations'}))
        rows{end+1} = {'Covariates', strjoin(covs, ', '), false};
    end
    nc = obs_clause(rs, params);
    if ~isempty(nc), rows{end+1} = {'Observations', nc, false}; end
    if isfield(params,'datatype') && ~isempty(params.datatype)
        rows{end+1} = {'Data type', char(params.datatype), true};
    end
    rows{end+1} = {'Montage', montage_label(rs, params), true};
    nperm = 1000;
    if isfield(params,'permutations') && ~isempty(params.permutations), nperm = double(params.permutations); end
    rows{end+1} = {'Permutations &middot; &alpha;', ...
        sprintf('%s &middot; %s', fmt_int(nperm), fmt_alpha(gf(rs,'alpha',0.05))), false};
    rows{end+1} = {'Correction', 'TFCE + cluster', true};
    csh = cluster_stat_html(rs, comp);
    if ~isempty(csh), rows{end+1} = {'Cluster statistic', [upper(csh(1)) csh(2:end)], false}; end
    rows{end+1} = {'Dataset', strrep(base, '_', ' '), false};

    M = struct('formula', f, 'rows', {rows});
end

% A legacy correlation run that partialled covariates out (core_snpm_analysis
% sets rs.covariates.method only on that route).
function tf = has_fl_covariates(rs)
    tf = isfield(rs,'covariates') && isstruct(rs.covariates) ...
        && isfield(rs.covariates,'method') && ~isempty(rs.covariates.method);
end

function s = strjoin_or(c, dflt)
    if isempty(c), s = dflt; else, s = strjoin(c, ' + '); end
end

function v = pget(params, f, dflt)
    v = dflt;
    if isfield(params, f) && ~isempty(params.(f)), v = char(string(params.(f))); end
end

function c = covariate_names(params)
    c = {};
    if isfield(params,'covariate_cols') && ~isempty(params.covariate_cols)
        c = cellstr(string(params.covariate_cols));
    elseif isfield(params,'use_covariates') && params.use_covariates ...
            && isfield(params,'covariate_file') && ~isempty(params.covariate_file)
        c = {'covariate file'};
    end
    c = c(:)';
end

% Recording system as the reader knows it (the registry's display name, not the
% legacy GUI key '164 channels'), plus how many channels were actually analysed.
function s = montage_label(rs, params)
    s = '';
    if isfield(params,'channels') && ~isempty(params.channels)
        s = char(string(params.channels));
        try
            m = snpm_montage_registry(s);
            s = m.display;
        catch
            % unknown key: show what the caller asked for, verbatim
        end
    end
    n = 0;
    if isfield(rs,'chanlocs'), n = numel(rs.chanlocs); end
    if n > 0
        if isempty(s), s = sprintf('%d channels analysed', n);
        else,          s = sprintf('%s &mdash; %d channels analysed', s, n);
        end
    elseif isempty(s)
        s = 'not stated';
    end
end

% =========================================================================
% NOTICES (amber block under the model card)
% Anything that qualifies EVERY number in the report and would otherwise die
% in a headless console: excluded channels, an unverified Data 2 pairing, the
% normalize/broadcast constraint, an uninterpretable global, the circular
% rotation actually applied. Rendered only when the list is non-empty.
% The literal strings 'Excluded channels', 'Data Type caveat' and
% 'rotation +X.X deg' are asserted by the test suites -- do not reword them.
% =========================================================================
function N = build_notices(rs)
    N = {};

    % -- excluded channels (single set, or the union across effects) --------
    if isfield(rs,'effects') && numel(rs.effects) > 1
        [idx, labs, per_effect] = union_excluded(rs);
        if ~isempty(idx) && ~per_effect
            N{end+1} = notice('Excluded channels', sprintf([ ...
                '<b>%d</b> of %d channel(s) were not tested. The <b>same channel set is used ' ...
                'for every effect</b> in this report &mdash; group, condition, interaction, ' ...
                'the simple effects and all post-hoc comparisons, and the whole-head averages: ' ...
                'a channel is analysed only if it is present in every row of the table. ' ...
                'Excluded: <span style="font-family:var(--mono);font-size:11.5px">%s</span>'], ...
                numel(idx), effect_nchannels(rs), esc(strjoin(labs, ', '))));
        elseif ~isempty(idx)
            % Zero-variance exclusions are made per test (each effect, simple
            % effect and post-hoc pair has its own rows and nuisance model), so
            % the sets differ; the missing-data set and the whole-head averages
            % are still shared.
            N{end+1} = notice('Excluded channels', sprintf([ ...
                '<b>%d</b> of %d channel(s) were not tested in at least one test of this ' ...
                'report. Channels missing in any row of the table are excluded from every ' ...
                'effect and from the whole-head averages. Channels with no variation to ' ...
                'permute, or whose observed statistic is undefined (a perfect fit), are ' ...
                'excluded only from the test where that happens; each effect and post-hoc ' ...
                'below lists its own. Whole-head averages are over the channels with no ' ...
                'missing data. Excluded somewhere: <span style="font-family:var(--mono);' ...
                'font-size:11.5px">%s</span>'], ...
                numel(idx), effect_nchannels(rs), esc(strjoin(labs, ', '))));
        end
    elseif isfield(rs,'excluded_channels') && isstruct(rs.excluded_channels) ...
            && isfield(rs.excluded_channels,'n') && rs.excluded_channels.n > 0
        ex = rs.excluded_channels;
        labs = ex.labels; if ~iscell(labs), labs = cellstr(string(labs)); end
        nch = 0; if isfield(ex,'n_channels'), nch = ex.n_channels; end
        % The engine's own reason string, not a hard-coded "missing data": the
        % correlation tier also excludes non-finite cells (log10(0) = -Inf under
        % datatype 'logscale'), which has a different remedy, and the Excel
        % sheet already prints ex.reason.
        reason = 'missing data: channel not usable in every analysed subject';
        if isfield(ex,'reason') && ~isempty(ex.reason), reason = char(ex.reason); end
        if isfield(ex,'reasons') && numel(unique(ex.reasons)) > 1 || is_zero_var_only(ex)
            % Zero-variance exclusions (2026-10-14) carry their own reason per
            % channel, and their consequence differs from the missing-data one:
            % a channel with nothing to permute, or whose observed statistic is
            % undefined, is not tested -- it is not a row-set problem.
            N{end+1} = notice('Excluded channels', sprintf([ ...
                '<b>%d</b> of %d channel(s) were not tested and carry no statistic and no ' ...
                'p-value in the maps and tables. By reason: %s'], ...
                ex.n, nch, excl_reason_html(ex)));
        else
            N{end+1} = notice('Excluded channels', sprintf([ ...
                '<b>%d</b> of %d channel(s) were not tested. Reason: %s. They cannot be evaluated ' ...
                'on the same row set in every permutation, so they carry no statistic and no ' ...
                'p-value in either the observed map or the null. Channels: <span style="font-family:' ...
                'var(--mono);font-size:11.5px">%s</span>'], ...
                ex.n, nch, esc(reason), esc(strjoin(labs(:)', ', '))));
        end
    end

    % -- how Data 2 was matched (correlation) -------------------------------
    if isfield(rs,'data2_matched') && ~isempty(rs.data2_matched)
        switch char(rs.data2_matched)
            case 'positional'
                N{end+1} = notice('Data 2 matched by POSITION', [ ...
                    'Data 2''s column names are not recognised channel labels, so channel ' ...
                    '<i>i</i> of Data 1 was paired with Data 2''s <i>i</i>-th data column ' ...
                    '<b>on trust</b>. If the two files do not list their channels in the same ' ...
                    'order, every value in this report is mis-paired. Re-export Data 2 with ' ...
                    'the recording system''s labels (E1, Cz, ...) to have the pairing verified.']);
            case 'behavior'
                lab = 'a single behavioural measure';
                if isfield(rs,'behavior_col') && ~isempty(rs.behavior_col)
                    lab = sprintf('the behavioural measure <b>%s</b>', esc(char(rs.behavior_col)));
                end
                N{end+1} = notice('Data 2', sprintf([ ...
                    '%s &mdash; one value per subject, correlated against every channel of ' ...
                    'Data 1. There is no channel-to-channel pairing in this analysis.'], lab));
        end
    end

    % -- datatype 'normalize' + a broadcast behavioural score ---------------
    % Not a preference note: per-subject z-scoring across channels forces the
    % sd-weighted sum of the r-map to zero, so a uniformly signed association is
    % arithmetically excluded and the symmetric +/- map is the constraint, not a
    % finding.
    if isfield(rs,'normalize_broadcast_caveat') && ~isempty(rs.normalize_broadcast_caveat)
        N{end+1} = notice('Data Type caveat', esc(char(rs.normalize_broadcast_caveat)));
    end

    % -- partial correlation: how the covariates were controlled -------------
    if has_fl_covariates(rs)
        cv = rs.covariates;
        k = gf(cv, 'k_eff', NaN); df = gf(cv, 'error_df', NaN);
        N{end+1} = notice('Covariate control', sprintf([ ...
            'Partial correlation with <b>%d</b> covariate(s) (error df n &minus; 2 &minus; k = ' ...
            '<b>%d</b>). Both measures were residualized on the covariates and each ' ...
            'permutation re-projects the permuted EEG residuals off the covariates before ' ...
            'correlating (<b>Freedman&ndash;Lane</b>). This permutation test is ' ...
            '<b>approximate</b>, not exact.'], k, df));
        if isfield(cv,'conservative_note') && ~isempty(cv.conservative_note)
            N{end+1} = notice('Few error degrees of freedom', esc(char(cv.conservative_note)));
        end
        if isfield(cv,'spearman_caveat') && ~isempty(cv.spearman_caveat)
            N{end+1} = notice('Partial Spearman caveat', esc(char(cv.spearman_caveat)));
        end
    end

    % -- a global the engine declared uninterpretable ------------------------
    if isfield(rs,'global_not_interpretable') && rs.global_not_interpretable
        g = ['The channel-averaged signal is analytically zero under this transform, so the ' ...
            'global value below is floating-point noise rather than a whole-head result.'];
        if isfield(rs,'global_note') && ~isempty(rs.global_note), g = char(rs.global_note); end
        N{end+1} = notice('Global test not interpretable', esc(g));
    end

    % -- source space: no scalp images at all --------------------------------
    if isfield(rs,'is_source') && rs.is_source
        N{end+1} = notice('Source space', ['This run is in source space (cortical voxels), ' ...
            'which has no scalp layout, so no topography is drawn anywhere in this report ' ...
            '&mdash; that is the design, not a missing file. The per-voxel deliverable is the ' ...
            'sibling <span style="font-family:var(--mono);font-size:11.5px">_sigvoxels.csv</span>, ' ...
            'which carries every voxel''s statistic, p-values, coordinates and cluster ' ...
            'membership.']);
    end

    % -- circular: the rotation actually applied ----------------------------
    if isfield(rs,'circ') && isfield(rs.circ,'rotation_deg')
        N{end+1} = notice('Angle convention', sprintf([ ...
            'Angles were rotated by <b>%+.1f deg</b> (declared convention <b>%s</b>) so that ' ...
            '0 is the slow-oscillation up-state peak. A wrong declaration rotates every angle ' ...
            'invisibly. %s, not cluster extent.'], ...
            rs.circ.rotation_deg, esc(char(rs.circ.convention)), circ_cluster_note(rs)));
    end

    % -- circular-linear: subjects dropped for a missing measure value -------
    if isfield(rs,'circ') && isfield(rs.circ,'measure_nan_removed') && ...
            rs.circ.measure_nan_removed > 0
        ids = cellstr(string(rs.circ.dropped_subjects));
        N{end+1} = notice('Subjects dropped', sprintf([ ...
            '<b>%d</b> subject(s) have no finite value for the linear measure and were dropped ' ...
            'from the analysis at every channel before anything was computed: %s. Missing angles ' ...
            'are handled separately, by deleting that subject at that channel only.'], ...
            rs.circ.measure_nan_removed, esc(strjoin(ids(:)', ', '))));
    end
end

function n = notice(title, body)
    n = struct('title', title, 'body', body);
end

% Per-effect exclusion statement, rendered inside the effect section so it sits
% next to the maps it describes. Every effect of a multi-effect (mixed2way)
% report uses the one channel set, so this repeats the top notice's list in
% place. (Until 2026-10-14 the group effect could keep a channel the
% within-subject effects dropped, and an effect that excluded nothing printed
% "all channels were tested for this effect (the exclusions noted at the top
% apply to other effects)". With one shared set that case cannot arise -- an
% effect with nothing excluded implies a report with nothing excluded -- and
% the branch has been removed.)
function h = excluded_effect_html(rs, eff)
    h = '';
    ex = [];
    if isfield(eff,'excluded_channels') && isstruct(eff.excluded_channels)
        ex = eff.excluded_channels;
    end
    n = 0; labs = {};
    if ~isempty(ex) && isfield(ex,'n'), n = ex.n; end
    if n == 0, return; end
    if isfield(ex,'labels')
        labs = ex.labels; if ~iscell(labs), labs = cellstr(string(labs)); end
        labs = labs(:)';
    end
    nch = effect_nchannels(rs);
    if isfield(ex,'reasons') && numel(ex.reasons) == numel(ex.index)
        what = excl_reason_html(ex);
    else
        what = esc(strjoin(labs, ', '));
    end
    h = sprintf(['    <p class="methods" style="margin:14px 0 0;color:#b8862a"><b>Excluded ' ...
        'for this effect:</b> %d of %d channel(s) &mdash; %s. They carry no statistic and no ' ...
        'p-value in the maps and tables below.</p>\n'], n, nch, what);
end

% 'E1, E2 (reason A); E7 (reason B)' as HTML, channel labels in the mono font.
function h = excl_reason_html(ex)
    labs = ex.labels; if ~iscell(labs), labs = cellstr(string(labs)); end
    labs = labs(:)';
    if isfield(ex,'reasons') && numel(ex.reasons) == numel(labs)
        r = ex.reasons(:)';
    else
        r = repmat({char(ex.reason)}, 1, numel(labs));
    end
    u = unique(r, 'stable');
    parts = cell(1, numel(u));
    for k = 1:numel(u)
        parts{k} = sprintf(['<span style="font-family:var(--mono);font-size:11.5px">%s</span> ' ...
            '(%s)'], esc(strjoin(labs(strcmp(r, u{k})), ', ')), esc(u{k}));
    end
    h = strjoin(parts, '; ');
end

% True when every listed reason is a zero-variance one (none is missing data).
function tf = is_zero_var_only(ex)
    zv = {'no variation left after the nuisance model', 'perfect fit, statistic undefined', ...
          'identical values in both conditions', 'every value is zero', ...
          'same value in every subject', 'perfect separation, t undefined', ...
          'constant across subjects', 'fully explained by the covariates'};
    tf = false;
    if ~isfield(ex,'reasons') || isempty(ex.reasons), return; end
    tf = all(ismember(ex.reasons, zv));
end

function n = effect_nchannels(rs)
    n = 0;
    if isfield(rs,'chanlocs'), n = numel(rs.chanlocs); end
    if n == 0 && isfield(rs,'effects') && ~isempty(rs.effects)
        e = rs.effects(1);
        if isfield(e,'excluded_channels') && isfield(e.excluded_channels,'n_channels')
            n = e.excluded_channels.n_channels;
        end
    end
end

function [idx, labs, per_effect] = union_excluded(rs)
% per_effect: true when the loops do not all exclude the same channels (the
% zero-variance rule excludes per test), so the "one channel set for every
% effect" sentence would be false.
    idx = []; sets = {};
    if isfield(rs,'effects')
        for k = 1:numel(rs.effects)
            sets{end+1} = eff_excl_idx(rs.effects(k)); %#ok<AGROW>
            sets = [sets, posthoc_excl_sets(rs.effects(k))]; %#ok<AGROW>
            se = [];
            if isfield(rs.effects(k),'simple_effects'), se = rs.effects(k).simple_effects; end
            for s = 1:numel(se)
                sets{end+1} = eff_excl_idx(se(s)); %#ok<AGROW>
                sets = [sets, posthoc_excl_sets(se(s))]; %#ok<AGROW>
            end
        end
    end
    per_effect = false;
    for k = 1:numel(sets)
        idx = [idx, sets{k}]; %#ok<AGROW>
        if k > 1 && ~isequal(sort(sets{k}), sort(sets{1})), per_effect = true; end
    end
    idx = unique(idx(:)');
    labs = cell(1, numel(idx));
    for i = 1:numel(idx), labs{i} = chan_label(rs, idx(i)); end
end

function S = posthoc_excl_sets(eff)
    S = {};
    if ~isfield(eff,'posthoc'), return; end
    for i = 1:numel(eff.posthoc)
        S{end+1} = eff_excl_idx(eff.posthoc(i)); %#ok<AGROW>
    end
end

function v = eff_excl_idx(eff)
    v = [];
    if isfield(eff,'excluded_channels') && isstruct(eff.excluded_channels) && ...
            isfield(eff.excluded_channels,'index')
        v = reshape(eff.excluded_channels.index, 1, []);
    end
end

% =========================================================================
% RENDER: R -> HTML, in the order of templates/snpm_report_v2_master.html
% =========================================================================
function write_document(fid, R)
    fprintf(fid, '%s', '<!DOCTYPE html>');
    fprintf(fid, '\n%s', '<html lang="en"><head><meta charset="utf-8">');
    fprintf(fid, '\n%s', '<meta name="viewport" content="width=device-width, initial-scale=1">');
    fprintf(fid, '\n<title>SnPM report &mdash; %s</title>\n', R.title);
    fprintf(fid, '%s\n', v2_css());
    fprintf(fid, '%s\n', v2_js());
    fprintf(fid, '%s', '</head>');
    fprintf(fid, '\n<body>\n<div class="wrap">\n\n');

    write_masthead(fid, R);
    write_strip(fid, R);
    write_model(fid, R);
    write_notices(fid, R);
    write_legend(fid, R);
    write_descriptive(fid, R);
    for k = 1:numel(R.effects)
        write_effect(fid, R, R.effects{k});
    end
    write_global(fid, R);
    write_footer(fid);

    fprintf(fid, '\n</div>\n</body>\n</html>\n');
end

function write_masthead(fid, R)
    fprintf(fid, '<header class="masthead">\n  <div>\n');
    fprintf(fid, '    <span class="eyebrow">SnPM Topographic Statistics</span>\n');
    fprintf(fid, '    <h1>%s</h1>\n', R.title);
    fprintf(fid, '    <p class="methods">%s</p>\n', R.methods);
    fprintf(fid, '  </div>\n  <div class="actions">\n');
    fprintf(fid, ['    <a class="btn" href="%s" download><svg viewBox="0 0 24 24" fill="none" ' ...
        'stroke="currentColor" stroke-width="1.8"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 ' ...
        '1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>Results CSV</a>\n'], url(R.csv));
    fprintf(fid, ['    <button class="btn" onclick="window.print()"><svg viewBox="0 0 24 24" ' ...
        'fill="none" stroke="currentColor" stroke-width="1.8"><path d="M6 9V3h12v6M6 18H4a2 2 ' ...
        '0 0 1-2-2v-4a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v4a2 2 0 0 1-2 2h-2M6 14h12v7H6z"/></svg>' ...
        'Print / PDF</button>\n']);
    fprintf(fid, '  </div>\n</header>\n\n');
end

% Sticky per-effect summary strip: traffic light + counts + global headline.
function write_strip(fid, R)
    fprintf(fid, '<nav class="strip">\n');
    for k = 1:numel(R.effects)
        e = R.effects{k};
        fprintf(fid, '  <a class="chipfx" href="#fx-%s"><span class="lamp %s"></span>', ...
            e.key, lamp_cls(e.lamp));
        fprintf(fid, '<span class="nm">%s</span>', e.title);
        fprintf(fid, '<span class="cts">Clu %d &middot; TFCE %d &middot; Unc %d</span>', ...
            numel(e.sclu), numel(e.stfce), numel(e.sunc));
        if ~isempty(e.global)
            g = e.global;
            fprintf(fid, '<span class="gp%s">%s = %s &middot; %s</span>', ...
                psig_cls(g, R.alpha), g.sym, fmt_stat(g.stat), global_p_text(g, R.nperm));
        end
        fprintf(fid, '</a>\n');
    end
    fprintf(fid, '  <a class="chipfx" href="#global"><span class="nm">Global test</span></a>\n');
    fprintf(fid, '</nav>\n\n');
end

function write_model(fid, R)
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, '  <div class="head"><h2>Model</h2><span class="sub">The design fitted at every channel</span></div>\n');
    fprintf(fid, '  <div class="model-card">\n');
    fprintf(fid, '    <div class="formula"><span class="fk">Fitted per channel</span>\n');
    fprintf(fid, '      <code>%s</code>\n    </div>\n', R.model.formula);
    fprintf(fid, '    <div class="specs">\n');
    for i = 1:numel(R.model.rows)
        row = R.model.rows{i};
        cls = 'val'; if row{3}, cls = 'val tag'; end
        fprintf(fid, '      <div class="row"><span class="k">%s</span><span class="%s">%s</span></div>\n', ...
            row{1}, cls, esc_soft(row{2}));
    end
    fprintf(fid, '    </div>\n  </div>\n</section>\n\n');
end

% Amber notice block. Not in the mock (which shows a clean run), so it borrows
% the .legend card and states its colour inline; every string in it is a
% validity caveat that must survive a headless run.
function write_notices(fid, R)
    if isempty(R.notices), return; end
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, ['  <div class="head"><h2>Notes on this run</h2>' ...
        '<span class="sub">Read before interpreting the maps</span></div>\n']);
    fprintf(fid, ['  <div class="legend" style="border-color:#b8862a;background:#fdf8ee;' ...
        'display:block">\n']);
    for i = 1:numel(R.notices)
        n = R.notices{i};
        fprintf(fid, ['    <div class="grp" style="align-items:flex-start;margin:6px 0">' ...
            '<span class="k" style="color:#b8862a;white-space:nowrap">%s</span>' ...
            '<span>%s</span></div>\n'], n.title, n.body);
    end
    fprintf(fid, '  </div>\n</section>\n\n');
end

function write_legend(fid, R)
    signed = false; unsigned = false;
    for k = 1:numel(R.effects)
        if R.effects{k}.signed, signed = true; else, unsigned = true; end
    end
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, ['  <div class="head"><h2>Legend</h2><span class="sub">Statistics, colour ' ...
        'scales &amp; markers</span></div>\n']);
    fprintf(fid, '  <div class="legend">\n');
    if signed
        lab = 'Signed t';
        if strcmp(R.statName,'r'), lab = 'Correlation r'; end
        fprintf(fid, '%s\n', ['    <div class="grp"><span class="k">' lab '</span>' ...
            '<span class="tmapbar" style="background:linear-gradient(90deg,#2f6db0,#eef2f6 50%,#c2543a)"></span>' ...
            '<span style="color:var(--cool)">&minus;</span><span style="color:var(--warm)">+</span></div>']);
        fprintf(fid, '    <span class="sep"></span>\n');
    end
    if unsigned
        fprintf(fid, '%s\n', ['    <div class="grp"><span class="k">Omnibus statistic</span>' ...
            '<span class="tmapbar" style="background:linear-gradient(90deg,#f3f6f8,#f4b63a 45%,#c2543a 78%,#7a1f10)"></span>' ...
            '<span style="color:var(--muted)">low</span><span style="color:var(--warm)">high</span></div>']);
        fprintf(fid, '    <span class="sep"></span>\n');
    end
    fprintf(fid, '%s\n', ['    <div class="grp"><span class="k">Channel dots</span>' ...
        '<span style="display:inline-flex;align-items:center;gap:6px"><span class="dd" style="background:#1b2733"></span>uncorrected</span>' ...
        '<span style="display:inline-flex;align-items:center;gap:6px"><span class="dd" style="background:#fff;border:1px solid var(--ink-2)"></span>corrected</span></div>']);
    fprintf(fid, '    <span class="sep"></span>\n');
    fprintf(fid, ['    <div class="grp"><span class="k">Cluster number</span><span class="cnum">1</span>' ...
        '<span>at cluster centroid, matches the table</span></div>\n']);
    fprintf(fid, '  </div>\n</section>\n\n');
end

% =========================================================================
% DESCRIPTIVE BLOCK
% =========================================================================
function write_descriptive(fid, R)
    D = R.descriptive;
    if strcmp(D.kind,'none'), return; end
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, '  <div class="head"><h2>%s</h2><span class="sub">%s</span></div>\n', D.title, D.sub);
    switch D.kind
        case 'cells'
            fprintf(fid, '  <div class="desc-flex">\n');
            write_cellgrid(fid, D);
            for i = 1:numel(D.figs), write_fig(fid, D.figs{i}, 4); end
            fprintf(fid, '  </div>\n');
        case 'circ'
            fprintf(fid, '  <div class="desc-flex" style="grid-template-columns:1fr 1fr">\n');
            write_circ_table(fid, R);
            for i = 1:numel(D.figs), write_fig(fid, D.figs{i}, 4); end
            fprintf(fid, '  </div>\n');
        otherwise
            n = max(numel(D.figs), 1);
            fprintf(fid, '  <div class="desc-flex" style="grid-template-columns:repeat(%d,1fr)">\n', n);
            for i = 1:numel(D.figs), write_fig(fid, D.figs{i}, 4); end
            fprintf(fid, '  </div>\n');
    end
    fprintf(fid, '</section>\n\n');
end

% Group x condition cell means. Any grid size: rows/cols come from the engine.
function write_cellgrid(fid, D)
    rows = D.rows; cols = D.cols;
    pngs = {}; labs = {};
    if numel(D.cells) == 2, pngs = D.cells{1}; labs = D.cells{2}; end
    if isempty(rows) || isempty(cols) || isempty(pngs), return; end
    fprintf(fid, '    <div class="fig">\n');
    fprintf(fid, ['      <div class="cap"><span class="t">Cell mean topographies</span>' ...
        '<span class="f">&hellip;_cell_&lt;row&gt;_&lt;col&gt;.png</span></div>\n']);
    fprintf(fid, '      <div class="cellgrid" style="grid-template-columns:56px repeat(%d,1fr)">\n', numel(cols));
    fprintf(fid, '        <span></span>');
    for c = 1:numel(cols), fprintf(fid, '<span class="ch">%s</span>', esc(cols{c})); end
    fprintf(fid, '\n');
    for r = 1:numel(rows)
        fprintf(fid, '        <span class="rh">%s</span>\n', esc(rows{r}));
        for c = 1:numel(cols)
            want = sprintf('%s / %s', rows{r}, cols{c});
            png = '';
            for i = 1:numel(labs)
                if strcmp(labs{i}, want) && i <= numel(pngs), png = pngs{i}; break; end
            end
            fprintf(fid, '        %s\n', slot_html(png, ''));
        end
    end
    fprintf(fid, '      </div>\n    </div>\n');
end

function write_circ_table(fid, R)
    d = R.rs.circ_descriptive;
    if ~isfield(d,'group') || isempty(d.group), return; end
    fprintf(fid, '    <div class="tablewrap">\n');
    fprintf(fid, ['      <div class="thead-note"><span>Per-group circular descriptives</span>' ...
        '<span>no inference &mdash; the maps carry the test</span></div>\n']);
    fprintf(fid, ['      <table class="sum"><thead><tr><th>Group</th><th>n</th>' ...
        '<th>Mean dir (deg)</th><th>SD (deg)</th><th>R</th><th>Rayleigh z</th>' ...
        '<th>Rayleigh p</th></tr></thead><tbody>\n']);
    for k = 1:numel(d.group)
        g = d.group(k);
        fprintf(fid, ['        <tr><td class="band">%s</td><td>%d</td><td>%.1f</td>' ...
            '<td>%.1f</td><td>%.3f</td><td>%.2f</td><td>%s</td></tr>\n'], ...
            esc(g.label), g.n, g.mean_deg, g.sd_deg, g.R, g.rayleigh_z, fmt_p(g.rayleigh_p, 0));
    end
    fprintf(fid, '      </tbody></table>\n');
    if isfield(d,'note') && ~isempty(d.note)
        fprintf(fid, '      <div class="empty-note" style="padding:12px 16px;text-align:left">%s</div>\n', ...
            esc(char(d.note)));
    end
    fprintf(fid, '    </div>\n');
end

% =========================================================================
% ONE EFFECT SECTION
% =========================================================================
function write_effect(fid, R, e)
    opencls = ''; if e.lamp, opencls = ' open'; end
    fprintf(fid, '<section class="effect%s" id="fx-%s">\n', opencls, e.key);
    fprintf(fid, '  <div class="ehead" onclick="toggleFx(''fx-%s'')">\n', e.key);
    fprintf(fid, ['    <svg class="chev" viewBox="0 0 24 24" fill="none" stroke="currentColor" ' ...
        'stroke-width="2"><path d="m9 6 6 6-6 6"/></svg>\n']);
    fprintf(fid, '    <h2>%s</h2><span class="esub">%s</span>\n', e.title, e.subtitle);
    fprintf(fid, ['    <span class="lamp %s"></span><span class="ecount">Clu %d &middot; ' ...
        'TFCE %d &middot; Unc %d</span>\n'], lamp_cls(e.lamp), numel(e.sclu), numel(e.stfce), numel(e.sunc));
    fprintf(fid, '  </div>\n  <div class="ebody">\n');

    if ~isempty(e.excl_html), fprintf(fid, '%s', e.excl_html); end
    if ~isempty(e.note_html), fprintf(fid, '%s', e.note_html); end

    if ~isempty(e.eff_png)
        fprintf(fid, '    <div class="blk">\n      <div class="bh">Effect map</div>\n');
        fprintf(fid, '      <div class="effect-map">\n');
        write_fig(fid, fig_spec(e.eff_png, e.eff_title, effect_map_note(e), ''), 8);
        fprintf(fid, '      </div>\n    </div>\n');
    end

    fprintf(fid, '    <div class="blk">\n      <div class="bh">Significant channels</div>\n');
    write_seg(fid, e.key, numel(e.sclu), numel(e.stfce), numel(e.sunc));
    write_cluster_view(fid, R, e);
    write_tfce_view(fid, R, e);
    write_unc_view(fid, R, e);
    fprintf(fid, '    </div>\n');

    write_posthoc(fid, R, e);

    fprintf(fid, '  </div>\n</section>\n\n');
end

function s = effect_map_note(e)
    if e.signed
        s = map_note_signed(sprintf('%s &lt; 0', e.sym), sprintf('%s &gt; 0', e.sym));
    else
        s = ['<span class="dl"><span class="dd" style="background:#f4b63a"></span>low</span>' ...
             '<span class="dl"><span class="dd" style="background:#7a1f10"></span>high</span>'];
    end
end

function s = map_note_signed(lo, hi)
    s = ['<span class="dl"><span class="dd" style="background:#2f6db0"></span>' lo '</span>' ...
         '<span class="dl"><span class="dd" style="background:#c2543a"></span>' hi '</span>'];
end

function write_seg(fid, fx, nC, nT, nU)
    fprintf(fid, '      <div class="seg">\n');
    fprintf(fid, ['        <button id="btn_%s_clu" class="on" onclick="showView(''%s'',''clu'')">' ...
        'Cluster <span class="ct%s">%d</span></button>\n'], fx, fx, none_cls(nC), nC);
    fprintf(fid, ['        <button id="btn_%s_tfce" onclick="showView(''%s'',''tfce'')">' ...
        'TFCE <span class="ct%s">%d</span></button>\n'], fx, fx, none_cls(nT), nT);
    fprintf(fid, ['        <button id="btn_%s_unc" onclick="showView(''%s'',''unc'')">' ...
        'Uncorrected <span class="ct%s">%d</span></button>\n'], fx, fx, none_cls(nU), nU);
    fprintf(fid, '      </div>\n');
end

function write_cluster_view(fid, R, e)
    fprintf(fid, '      <div id="panel_%s_clu" class="vpanel" data-view-name="Cluster-corrected">\n', e.key);
    [ord, sgn] = snpm_order_clusters(e.Clusters, e.T, R.alpha);
    if isempty(ord)
        % Only point at "the effect map above" when there is one (a post-hoc
        % pair has no effect map of its own).
        mapnote = '';
        if isfield(e, 'eff_png') && ~isempty(e.eff_png)
            mapnote = ' The effect map above shows the uncorrected pattern.';
        end
        csnote = '';
        if ~isempty(R.cluster_stat_html), csnote = [' Cluster statistic: ' R.cluster_stat_html '.']; end
        fprintf(fid, '        <div class="tablewrap"><div class="empty-note">%s</div></div>\n', ...
            sprintf('<b>No significant clusters</b> at &alpha; = %s%s.%s%s', fmt_alpha(R.alpha), ...
            smallest_clause(min_cluster_p(e.Clusters), R.nperm, 'cluster p'), mapnote, csnote));
        fprintf(fid, '      </div>\n');
        return;
    end
    fprintf(fid, '        <div class="vbody">\n');
    if ~isempty(e.map_base)
        % A post-hoc pair passes its own caption and a signed legend prefix
        % (clu_title / map_legend); an effect uses the generic ones.
        write_fig(fid, fig_spec([e.map_base ' Cluster.png'], ...
            gf(e, 'clu_title', 'Cluster significance map'), ...
            [gf(e, 'map_legend', '') ...
             '<span class="dl"><span class="cnum" style="width:14px;height:14px;font-size:9px">1</span>' ...
             'cluster number at centroid</span><span class="dl"><span class="dd" style="background:#fff">' ...
             '</span>cluster-corrected</span>'], ''), 10);
    end
    write_cluster_table(fid, R, e, ord, sgn);
    fprintf(fid, '        </div>\n      </div>\n');
end

function write_cluster_table(fid, R, e, ord, sgn)
    total = 0;
    for j = 1:numel(ord), total = total + numel(cluster_chs(e.Clusters, ord(j))); end
    if e.signed, right = 'positive and negative clusters listed separately';
    else,        right = 'grouped by cluster, ascending cluster p'; end
    fprintf(fid, '          <div class="tablewrap">\n');
    fprintf(fid, ['            <div class="thead-note"><span><b>%d</b> channels in <b>%d</b> ' ...
        'significant cluster(s)</span><span>%s</span></div>\n'], total, numel(ord), right);
    fprintf(fid, ['            <table class="sum"><thead><tr><th>Channel</th><th>%s</th>' ...
        '<th>cluster p</th></tr></thead><tbody>\n'], e.sym);
    shown_pos = false; shown_neg = false;
    for j = 1:numel(ord)
        ci = ord(j);
        if e.signed && sgn(j) > 0 && ~shown_pos
            fprintf(fid, ['              <tr class="sign-row"><td colspan="3">Positive clusters ' ...
                '<span class="sgn">(%s &gt; 0)</span></td></tr>\n'], e.sym);
            shown_pos = true;
        elseif e.signed && sgn(j) < 0 && ~shown_neg
            fprintf(fid, ['              <tr class="sign-row"><td colspan="3">Negative clusters ' ...
                '<span class="sgn">(%s &lt; 0)</span></td></tr>\n'], e.sym);
            shown_neg = true;
        end
        ch = cluster_chs(e.Clusters, ci);
        cp = cluster_p(e.Clusters, ci);
        sv = arrayfun(@(c) abs(stat_at(e, c)), ch);
        [~, o] = sort(sv, 'descend'); ch = ch(o);
        fprintf(fid, ['              <tr class="grp-row"><td colspan="3"><span class="cnum">%d' ...
            '</span>Cluster %d<span class="cp">%s &middot; %d ch</span></td></tr>\n'], ...
            j, j, p_phrase(cp, R.nperm), numel(ch));
        for k = 1:numel(ch)
            fprintf(fid, ['              <tr><td class="band"><span class="dotkey cor"></span>%s' ...
                '</td><td class="stat-v">%s</td><td class="p-sig">%s</td></tr>\n'], ...
                esc(chan_label(R.rs, ch(k))), fmt_stat(stat_at(e, ch(k))), fmt_p(cp, R.nperm));
        end
    end
    fprintf(fid, '            </tbody></table>\n');
    if ~isempty(R.cluster_stat_html)
        fprintf(fid, '            <div class="thead-note"><span>Cluster statistic: %s</span></div>\n', ...
            R.cluster_stat_html);
    end
    fprintf(fid, '          </div>\n');
end

function write_tfce_view(fid, R, e)
    fprintf(fid, ['      <div id="panel_%s_tfce" class="vpanel" data-view-name="TFCE-corrected" ' ...
        'style="display:none">\n'], e.key);
    if isempty(e.stfce)
        fprintf(fid, '        <div class="tablewrap"><div class="empty-note">%s</div></div>\n', ...
            sprintf('<b>No TFCE-corrected significant channels</b> at &alpha; = %s%s.', ...
            fmt_alpha(R.alpha), smallest_clause(min_finite(e.ptfce), R.nperm, 'p')));
        fprintf(fid, '      </div>\n');
        return;
    end
    fprintf(fid, '        <div class="vbody">\n');
    if ~isempty(e.map_base)
        % Same optional per-pair fields as write_cluster_view, plus tfce_note
        % for the table's right-hand note.
        if isfield(e, 'map_legend')
            legend_html = e.map_legend;     % signed t legend replaces "uncorrected"
        else
            legend_html = '<span class="dl"><span class="dd" style="background:#1b2733"></span>uncorrected</span>';
        end
        write_fig(fid, fig_spec([e.map_base ' TFCE.png'], ...
            gf(e, 'tfce_title', 'TFCE significance map'), ...
            [legend_html ...
             '<span class="dl"><span class="dd" style="background:#fff"></span>TFCE-corrected</span>'], ''), 10);
    end
    write_flat_table(fid, R, e, e.stfce, e.ptfce, 'TFCE-corrected', ...
        gf(e, 'tfce_note', 'sorted by ascending p'), 'cor', Inf, 10);
    fprintf(fid, '        </div>\n      </div>\n');
end

function write_unc_view(fid, R, e)
    fprintf(fid, ['      <div id="panel_%s_unc" class="vpanel" data-view-name="Uncorrected" ' ...
        'style="display:none">\n'], e.key);
    if isempty(e.sunc)
        fprintf(fid, '        <div class="tablewrap"><div class="empty-note">%s</div></div>\n', ...
            sprintf('<b>No uncorrected significant channels</b> at &alpha; = %s%s.', ...
            fmt_alpha(R.alpha), smallest_clause(min_finite(e.punc), R.nperm, 'p')));
    else
        write_flat_table(fid, R, e, e.sunc, e.punc, 'uncorrected', ...
            'sorted by ascending p &middot; full list in the results CSV', 'unc', 20, 8);
    end
    fprintf(fid, '      </div>\n');
end

% Channel | statistic | p, ascending p. CAP caps the visible rows (extra rows
% become tr.xrow behind the "Show all N channels" button); Inf = never capped
% (cluster and TFCE tables are complete by design).
function write_flat_table(fid, R, e, chs, pv, label, right, dotcls, cap, ind)
    pad = repmat(' ', 1, ind);
    chs = chs(:).';
    ps = arrayfun(@(c) p_at(pv, c), chs);
    [ps, o] = sort(ps); chs = chs(o);
    fprintf(fid, '%s<div class="tablewrap">\n', pad);
    fprintf(fid, '%s  <div class="thead-note"><span><b>%d</b> %s significant channels</span><span>%s</span></div>\n', ...
        pad, numel(chs), label, right);
    fprintf(fid, '%s  <table class="sum"><thead><tr><th>Channel</th><th>%s</th><th>p</th></tr></thead><tbody>\n', ...
        pad, e.sym);
    for k = 1:numel(chs)
        xcls = ''; if k > cap, xcls = ' class="xrow"'; end
        fprintf(fid, ['%s    <tr%s><td class="band"><span class="dotkey %s"></span>%s</td>' ...
            '<td class="stat-v">%s</td><td class="p-sig">%s</td></tr>\n'], ...
            pad, xcls, dotcls, esc(chan_label(R.rs, chs(k))), fmt_stat(stat_at(e, chs(k))), ...
            fmt_p(ps(k), R.nperm));
    end
    fprintf(fid, '%s  </tbody></table>\n', pad);
    if numel(chs) > cap
        fprintf(fid, '%s  <button class="morebtn" onclick="expandTable(this)">Show all %d channels</button>\n', ...
            pad, numel(chs));
    end
    fprintf(fid, '%s</div>\n', pad);
end

% ---- post-hoc pairwise contrasts as tabs --------------------------------
function write_posthoc(fid, R, e)
    P = e.posthoc;
    K = e.skipped;
    if isempty(P) && isempty(K), return; end
    grp = ['ph' e.key];
    names = cell(1, numel(P));
    for i = 1:numel(P), names{i} = P{i}.name; end
    nlist = ['[''' strjoin(names, ''',''') ''']'];
    fprintf(fid, '    <div class="blk">\n');
    fprintf(fid, ['      <div class="bh">Post-hoc pairwise contrasts &mdash; signed t per ' ...
        'comparison</div>\n']);
    % A pair that was not run has no tab below. Say so here, in the report, so
    % its absence is not mistaken for "nothing significant".
    for i = 1:numel(K)
        fprintf(fid, ['      <p class="methods" style="margin:0 0 10px;color:#b8862a"><b>Post-hoc ' ...
            '%s was not run:</b> %s.</p>\n'], esc(K{i}.label), esc(K{i}.reason));
    end
    if isempty(P)
        fprintf(fid, '    </div>\n');
        return;
    end
    fprintf(fid, '      <div class="tabs">\n');
    for i = 1:numel(P)
        oncls = ''; if i == 1, oncls = ' class="on"'; end
        n = numel(P{i}.sclu);   % the Cluster view opens first
        fprintf(fid, ['        <button id="tbtn_%s_%s"%s onclick="showTab(''%s'',''%s'',%s)">' ...
            '%s <span class="ct%s">%d</span></button>\n'], ...
            grp, P{i}.name, oncls, grp, P{i}.name, nlist, esc(P{i}.label), none_cls(n), n);
    end
    fprintf(fid, '      </div>\n');
    for i = 1:numel(P)
        q = P{i};
        hide = ' style="display:none"'; if i == 1, hide = ''; end
        fprintf(fid, '      <div id="tab_%s_%s" class="tabpanel" data-tab-name="Post-hoc: %s"%s>\n', ...
            grp, q.name, esc(q.label), hide);
        % What THIS pair was run on. It is tested on the rows of its two levels
        % only, so the count can be smaller than the one in the report header.
        if ~isempty(q.run_on)
            left = '';
            if ~isempty(q.dropped)
                left = sprintf(' Left out because they have only one of the two levels: %s.', ...
                    esc(strjoin(q.dropped(:)', ', ')));
            end
            fprintf(fid, ['        <p class="methods" style="margin:0 0 10px">Tested on the rows ' ...
                'of these two levels only: %s.%s</p>\n'], esc(q.run_on), left);
        end
        if ~isempty(q.excl_html)
            fprintf(fid, ['        <p class="methods" style="margin:0 0 10px;color:#b8862a"><b>' ...
                'Excluded for this comparison:</b> %s. They carry no statistic and no p-value ' ...
                'for this pair.</p>\n'], q.excl_html);
        end
        % Same Cluster | TFCE | Uncorrected switch as the effect, keyed per
        % pair so showView() toggles this pair's panels only. The PNGs are the
        % '<map_base> Cluster.png' / ' TFCE.png' pair the core already wrote.
        qe = q; qe.key = [e.key '_' q.name]; qe.sym = 't'; qe.signed = true;
        qe.clu_title  = sprintf('%s &mdash; Cluster map', esc(q.label));
        qe.tfce_title = sprintf('%s &mdash; TFCE map', esc(q.label));
        qe.map_legend = map_note_signed('t &lt; 0', 't &gt; 0');
        qe.tfce_note  = 'pairwise contrast &mdash; read alongside the omnibus effect';
        write_seg(fid, qe.key, numel(q.sclu), numel(q.stfce), numel(q.sunc));
        write_cluster_view(fid, R, qe);
        write_tfce_view(fid, R, qe);
        write_unc_view(fid, R, qe);
        fprintf(fid, '      </div>\n');
    end
    fprintf(fid, '    </div>\n');
end

% =========================================================================
% GLOBAL TEST CARD + FOOTER
% =========================================================================
function write_global(fid, R)
    fprintf(fid, '<section class="section" id="global" style="scroll-margin-top:70px">\n');
    fprintf(fid, ['  <div class="head"><h2>Global test</h2><span class="sub">Whole-head ' ...
        'statistic on the channel-averaged signal</span></div>\n']);
    fprintf(fid, '  <div class="global-card">\n');
    if R.has_global
        first = true;
        for k = 1:numel(R.effects)
            g = R.effects{k}.global;
            if isempty(g), continue; end
            if ~first, fprintf(fid, '    <div class="gsep"></div>\n'); end
            first = false;
            fprintf(fid, ['    <div class="gstat"><span class="gk">%s</span>' ...
                '<span class="gv">%s = %s</span><span class="gp%s">%s</span></div>\n'], ...
                g.label, g.sym, fmt_stat(g.stat), psig_cls(g, R.alpha), global_p_text(g, R.nperm));
        end
        fprintf(fid, '    <div class="gnote">%s</div>\n', global_note(R));
    else
        fprintf(fid, '    <div class="gnote" style="margin-left:0">%s</div>\n', R.global_note);
    end
    fprintf(fid, '  </div>\n</section>\n\n');
end

% The p-value text of a global card, or the engine's reason when it declared
% the statistic undefined (then stat and p are NaN).
function s = global_p_text(g, nperm)
    if isfield(g,'reason') && ~isempty(g.reason)
        s = esc(char(g.reason));
    else
        s = p_phrase(g.p, nperm);
    end
end

% What the channel average was actually taken over. NOT cosmetic: the global is
% a whole-head claim only when it rests on the whole montage. In a mixed2way
% report every effect is averaged over the SAME channel set (the one set fixed
% before any test: finite in every row of the table), so the counts are equal
% and collapse to one sentence. The unequal-counts wording at the end is kept
% for any other multi-effect caller whose effects do have different sets; the
% current pipelines do not produce that case.
function s = global_note(R)
    nch = effect_nchannels(R.rs);
    labs = {}; cnts = []; flagged = false;
    for k = 1:numel(R.effects)
        g = R.effects{k}.global;
        if isempty(g), continue; end
        if ~isempty(g.flag), flagged = true; end
        if ~isnan(g.nch)
            labs{end+1} = R.effects{k}.title; %#ok<AGROW>
            cnts(end+1) = g.nch;              %#ok<AGROW>
        end
    end
    if flagged
        s = ['<b>Not interpretable.</b> The value(s) above are stated only because the run ' ...
            'computed them &mdash; see the note near the top of this report.'];
        return;
    end
    if isempty(cnts)
        s = ['One whole-head test on the <b>channel-averaged</b> signal &mdash; the omnibus ' ...
            'answer before any per-channel localisation.'];
        return;
    end
    if isscalar(unique(cnts))
        if nch > 0 && cnts(1) < nch
            s = sprintf(['Computed on the mean across the <b>%d of %d</b> channels complete in ' ...
                'every analysed unit &mdash; not the whole head. Averaging each unit over its ' ...
                'own available channels would let a scalp gradient masquerade as an effect; ' ...
                'the channel-wise maps are unaffected.'], cnts(1), nch);
        else
            s = sprintf(['Computed on the mean across all <b>%d</b> analysed channels &mdash; ' ...
                'a single whole-head test, complementary to the channel-wise maps.'], cnts(1));
        end
        return;
    end
    parts = cell(1, numel(cnts));
    for i = 1:numel(cnts), parts{i} = sprintf('%s %d', labs{i}, cnts(i)); end
    s = sprintf(['Each effect is averaged over the channels <b>it</b> can evaluate, so the ' ...
        'counts differ (of %d): %s. Averaging every unit over one common set is what keeps a ' ...
        'scalp gradient out of the whole-head number; the channel-wise maps are unaffected.'], ...
        nch, strjoin(parts, ' &middot; '));
end

function write_footer(fid)
    [~, embed_on] = embed_ctx();
    if embed_on
        note = 'Self-contained report &mdash; images are embedded; the PNG files beside it are copies.';
    else
        note = 'Self-contained report &mdash; keep the PNG images in the same folder as this HTML.';
    end
    fprintf(fid, '<footer class="foot"><span>%s</span><span>Generated %s</span></footer>\n', ...
        note, esc(char(datetime('now','Format','dd-MMM-uuuu HH:mm'))));
end

% =========================================================================
% FIGURE / SLOT PRIMITIVES  (the mock's .fig card + onerror placeholder)
% =========================================================================
function write_fig(fid, F, ind)
    pad = repmat(' ', 1, ind);
    note_only = isfield(F,'emptynote') && ~isempty(F.emptynote);
    fprintf(fid, '%s<div class="fig">\n', pad);
    if note_only
        fprintf(fid, '%s  <div class="cap"><span class="t">%s</span><span class="f">%s</span></div>\n', ...
            pad, F.title, '(not produced)');
        fprintf(fid, '%s  <div class="empty-note">%s</div>\n', pad, F.emptynote);
        fprintf(fid, '%s</div>\n', pad);
        return;
    end
    fprintf(fid, '%s  <div class="cap"><span class="t">%s</span><span class="f">%s</span></div>\n', ...
        pad, F.title, esc(F.png));
    fprintf(fid, '%s  %s\n', pad, slot_html(F.png, F.cls));
    if ~isempty(F.note)
        fprintf(fid, '%s  <div class="map-note">%s</div>\n', pad, F.note);
    end
    fprintf(fid, '%s</div>\n', pad);
end

function s = slot_html(png, cls)
    if ~isempty(cls), cls = [' ' cls]; end
    if isempty(png), png = '(not produced)'; end
    % src is either the base64 payload of the PNG (default: the report is one
    % portable file) or the relative filename (report_embed_images = false, or
    % the PNG is missing -- in which case the src 404s and the onerror handler
    % below reveals the "No data" placeholder, exactly as before). The .fn
    % caption always names the real file, embedded or not.
    src = embed_src(png);
    s = sprintf(['<figure class="slot%s"><img src="%s" alt="" onerror="this.classList.add(''err'')">' ...
        '<figcaption class="ph"><span class="nd">No data</span><span class="fn">%s</span>' ...
        '</figcaption></figure>'], cls, src, esc(png));
end

% Per-report image-embedding context, set once by generateAnalysisReport and
% read by embed_src. Held here rather than threaded through the ~10 write_fig
% call sites, which take only (fid, F, ind).
function [dirpath, on] = embed_ctx(setdir, seton)
    persistent D ON
    if nargin == 2
        D = setdir; ON = logical(seton);
    end
    if isempty(D),  D  = '';    end
    if isempty(ON), ON = false; end
    dirpath = D; on = ON;
end

% Filename -> value for the img src attribute.
function src = embed_src(png)
    src = url(png);
    [dirpath, on] = embed_ctx();
    if ~on || isempty(png), return; end
    f = png;
    if ~isfile(f), f = fullfile(dirpath, png); end
    if ~isfile(f), return; end          % missing -> relative src -> onerror
    b = read_bytes(f);
    if isempty(b), return; end
    b64 = reshape(char(matlab.net.base64encode(b(:)')), 1, []);
    src = ['data:image/png;base64,' b64];
end

function b = read_bytes(f)
    b = uint8([]);
    fidp = fopen(f, 'r');
    if fidp == -1, return; end
    c = onCleanup(@() fclose(fidp));
    b = fread(fidp, Inf, '*uint8');
end

% =========================================================================
% DATA HELPERS
% =========================================================================
% The cluster-statistic wording for the report, from the same source as the
% results text (snpm_cluster_statistic), e.g. 'sum of t&#178; (or F) over the
% cluster'. '' when the run did not record a setting: the legacy t-tests use
% cluster extent, and a results file saved before the setting existed has none.
function h = cluster_stat_html(rs, comp)
    h = '';
    if ~isstruct(rs) || ~isfield(rs, 'cluster_statistic') || isempty(rs.cluster_statistic)
        return;
    end
    what = '';
    if strncmpi(char(comp), 'circ', 4), what = 'the statistic'; end
    [~, line] = snpm_cluster_statistic(rs.cluster_statistic, what);
    h = regexprep(line, '^Cluster statistic: ', '');
    h = strrep(esc(h), char(178), '&#178;');
end

% The cluster sentence of the circular convention notice. A results file saved
% before the setting existed was scored with the mean, so it says so.
function s = circ_cluster_note(rs)
    if isfield(rs, 'cluster_statistic') && ~isempty(rs.cluster_statistic)
        cs = snpm_cluster_statistic(rs.cluster_statistic);
    else
        cs = 'mean';
    end
    if strcmp(cs, 'sum')
        s = 'Cluster statistic: <b>sum</b> of the statistic over the cluster';
    else
        s = ['Cluster statistic: <b>mean</b> (sum of the statistic over the cluster, ' ...
             'divided by the number of channels)'];
    end
end

function v = gf(s, f, dflt)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function v = statvec(s)
    v = [];
    if ~isstruct(s) || ~isfield(s,'T'), return; end
    T = s.T;
    if isstruct(T) && isfield(T,'real_T'), v = double(T.real_T(:).');
    elseif isnumeric(T),                   v = double(T(:).');
    end
end

function v = pvec(s, fieldname)
    v = [];
    if isstruct(s) && isfield(s,'p') && isstruct(s.p) && isfield(s.p, fieldname)
        v = double(s.p.(fieldname)(:).');
    end
end

function v = idxvec(x)
    if isempty(x), v = []; else, v = double(x(:).'); end
end

function v = stat_at(e, ch)
    v = NaN;
    if ch >= 1 && ch <= numel(e.T), v = e.T(ch); end
end

function p = p_at(pv, ch)
    p = NaN;
    if ~isempty(pv) && ch >= 1 && ch <= numel(pv), p = pv(ch); end
end

function lab = chan_label(rs, ch)
    lab = sprintf('Ch %d', ch);
    if isstruct(rs) && isfield(rs,'chanlocs') && numel(rs.chanlocs) >= ch && isfield(rs.chanlocs,'labels')
        l = rs.chanlocs(ch).labels;
        if ~isempty(l), lab = l; end
    end
end

function pv = cluster_p(C, i)
    pv = NaN;
    if isstruct(C) && isfield(C,'p')
        v = C(i).p;
        if ~isempty(v), pv = v(1); end
    end
end

function ch = cluster_chs(C, i)
    ch = [];
    if isstruct(C) && isfield(C,'channels')
        v = C(i).channels;
        if ~isempty(v), ch = double(v(:).'); end
    end
end

function p = min_cluster_p(C)
    p = NaN;
    if isempty(C) || ~isstruct(C) || ~isfield(C,'p'), return; end
    v = [];
    for i = 1:numel(C)
        q = C(i).p;
        if ~isempty(q), v(end+1) = q(1); end %#ok<AGROW>
    end
    p = min_finite(v);
end

function m = min_finite(v)
    m = NaN;
    v = v(isfinite(v));
    if ~isempty(v), m = min(v); end
end

% " (smallest cluster p = .11)" -- the honest headline when a panel is empty.
function s = smallest_clause(p, nperm, what)
    if isnan(p), s = ''; return; end
    s = sprintf(' (smallest %s = %s)', what, fmt_p(p, nperm));
end

% =========================================================================
% FORMATTING HELPERS
% =========================================================================
% "p = .004" / "p &lt; .001" -- the equals sign is wrong in front of a bound.
function s = p_phrase(p, nperm)
    v = fmt_p(p, nperm);
    if strncmp(v, '&lt;', 4)
        s = ['p &lt; ' strtrim(v(5:end))];
    else
        s = ['p = ' v];
    end
end

function s = fmt_p(p, nperm)
    if isempty(p) || isnan(p)
        s = '&mdash;'; return;
    end
    if p <= 0
        if nperm >= 1, s = ['&lt; ' fmt_pnum(1/nperm)]; else, s = '&lt; .001'; end
        return;
    end
    s = fmt_pnum(p);
end

function s = fmt_pnum(p)
    if p >= 1, s = '1.00'; return; end
    if p < 0.001, s = '&lt;.001'; return; end
    s = sprintf('%.3f', p);
    if strncmp(s, '0', 1), s = s(2:end); end
end

function s = fmt_stat(v)
    if isempty(v) || isnan(v), s = '&mdash;'; return; end
    if v < 0, sgn = '&minus;'; else, sgn = ''; end
    s = sprintf('%s%.2f', sgn, abs(v));
end

function s = fmt_int(n)
    s = sprintf('%d', round(double(n)));
end

% The run's significance level in the report's own style (no leading zero), so
% every "at alpha = ..." sentence quotes the value the maps were actually
% thresholded at rather than a hard-coded .05.
function s = fmt_alpha(a)
    if isempty(a) || ~isfinite(a), a = 0.05; end
    s = sprintf('%.3g', double(a));
    if strncmp(s, '0.', 2), s = s(2:end); end
end

function c = none_cls(n)
    if n > 0, c = ''; else, c = ' none'; end
end

function c = lamp_cls(tf)
    if tf, c = 'sig'; else, c = 'ns'; end
end

% Highlight class for a global (whole-head) p-value. Thresholded at the RUN's
% alpha, not a hard-coded .05 -- a run at alpha = .01 must not paint a p = .03
% global as significant while every map card says otherwise.
function c = psig_cls(g, alpha)
    if nargin < 2 || isempty(alpha) || ~isfinite(alpha), alpha = 0.05; end
    c = '';
    if isstruct(g) && isfield(g,'p') && ~isempty(g.p) && ~isnan(g.p) && g.p <= alpha ...
            && (~isfield(g,'flag') || isempty(g.flag))
        c = ' p-sig';
    end
end

function t = has(s, sub)
    t = ischar(s) && ~isempty(strfind(s, sub)); %#ok<STREMP>
end

% Sanitised, unique tab id fragment (letters/digits only -- it goes into an
% element id and a JS string list).
function n = tab_name(label, i)
    n = regexprep(char(string(label)), '[^A-Za-z0-9]', '');
    if numel(n) > 16, n = n(1:16); end
    % '_' before the index: without it 'abc1'+1 and 'abc'+11 give the same id.
    n = sprintf('%s_%d', n, i);
end

% HTML-escape. The double quote is escaped too: engine strings (a behavioural
% column name, a group level, a dataset basename) reach attribute contexts such
% as the effect-map alt text and the chip title, where a bare " closes the
% attribute and everything after it becomes markup.
function s = esc(s)
    if ~ischar(s) && ~isstring(s), s = char(string(s)); end
    s = char(s);
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    s = strrep(s, '"', '&quot;');
end

% Like esc, but leaves an already-formed HTML entity (&mdash;, &#183;) alone --
% used for the model-card values, which mix engine strings with entities this
% file wrote itself.
function s = esc_soft(s)
    s = char(string(s));
    s = regexprep(s, '&(?!([A-Za-z]+|#[0-9]+);)', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    s = strrep(s, '"', '&quot;');
end

% Filename -> URL for src/href. Report filenames legitimately contain spaces
% ("<a> VS <b> Cluster.png"), which must be percent-encoded in an attribute.
function s = url(s)
    s = char(string(s));
    s = strrep(s, '%', '%25');
    s = strrep(s, ' ', '%20');
    s = strrep(s, '#', '%23');
    s = strrep(s, '?', '%3F');
    % A quote in a filename would close the src="..."/href="..." attribute it is
    % written into; percent-encode it rather than letting it become markup.
    s = strrep(s, '"', '%22');
    s = strrep(s, '&', '&amp;');
end

% =========================================================================
% STATIC ASSETS
% Lines 5-166 (<style>...</style>) and 167-185 (<script>...</script>) of
% templates/snpm_report_v2_master.html, ported BYTE-FOR-BYTE. The template is
% the source of truth for the design; a test diffs these against it. If you
% need a style change, change the template and re-port -- do not edit here.
% =========================================================================
function s = v2_css()
    c = {};
    c{end+1} = '<style>';
    c{end+1} = ':root{--ink:#1b2733;--ink-2:#445261;--muted:#7c8a98;--line:#e3e9ef;--line-2:#eef2f6;--bg:#f5f7f9;--panel:#ffffff;--accent:#1f6f8b;--accent-soft:#e8f1f4;--sig:#1f8a5b;--sig-soft:#e7f4ec;--warm:#c2543a;--cool:#2f6db0;--mono:ui-monospace,"SF Mono",Menlo,Consolas,"Roboto Mono",monospace;--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;--serif:Charter,"Iowan Old Style",Georgia,"Times New Roman",serif}';
    c{end+1} = '*{box-sizing:border-box}html{-webkit-text-size-adjust:100%;scroll-behavior:smooth}';
    c{end+1} = 'body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased}';
    c{end+1} = 'a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}';
    c{end+1} = '.wrap{max-width:1100px;margin:0 auto;padding:40px 32px 96px}';
    c{end+1} = '.eyebrow{font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);font-weight:600}';
    c{end+1} = '.masthead{display:flex;justify-content:space-between;align-items:flex-start;gap:32px;flex-wrap:wrap;padding-bottom:26px;border-bottom:1px solid var(--line)}';
    c{end+1} = '.masthead h1{font-family:var(--serif);font-weight:600;font-size:29px;line-height:1.2;letter-spacing:-.01em;margin:9px 0 11px;max-width:760px}';
    c{end+1} = '.methods{color:var(--ink-2);font-size:13.5px;max-width:820px;margin:0}.methods b{color:var(--ink);font-weight:600}';
    c{end+1} = '.actions{display:flex;gap:8px;flex-shrink:0}';
    c{end+1} = '.btn{font:inherit;font-size:12.5px;font-weight:600;color:var(--ink-2);background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:8px 13px;cursor:pointer;display:inline-flex;align-items:center;gap:6px;transition:border-color .15s,color .15s;text-decoration:none}';
    c{end+1} = '.btn:hover{border-color:var(--accent);color:var(--accent);text-decoration:none}.btn svg{width:14px;height:14px}';
    c{end+1} = '';
    c{end+1} = '/* ---- sticky per-effect summary strip ---- */';
    c{end+1} = '.strip{position:sticky;top:0;z-index:20;background:rgba(245,247,249,.94);backdrop-filter:blur(6px);margin:0 -32px;padding:10px 32px;border-bottom:1px solid var(--line);display:flex;gap:10px;flex-wrap:wrap}';
    c{end+1} = '.chipfx{display:flex;align-items:center;gap:10px;background:var(--panel);border:1px solid var(--line);border-radius:9px;padding:8px 14px;font-size:12px;color:var(--ink-2);text-decoration:none;transition:border-color .15s;box-shadow:0 1px 2px rgba(20,40,60,.04)}';
    c{end+1} = '.chipfx:hover{border-color:var(--accent);text-decoration:none}';
    c{end+1} = '.chipfx .lamp{width:9px;height:9px;border-radius:50%;flex:none}';
    c{end+1} = '.lamp.sig{background:var(--sig)}.lamp.ns{background:#c3ccd4}';
    c{end+1} = '.chipfx .nm{font-weight:650;color:var(--ink)}';
    c{end+1} = '.chipfx .cts{font-family:var(--mono);font-size:11px;color:var(--muted);font-variant-numeric:tabular-nums}';
    c{end+1} = '.chipfx .gp{font-family:var(--mono);font-size:11px;font-variant-numeric:tabular-nums;color:var(--ink-2)}';
    c{end+1} = '.chipfx .gp.p-sig{color:var(--sig);font-weight:700}';
    c{end+1} = '';
    c{end+1} = '.section{margin-top:40px}';
    c{end+1} = '.section>.head{display:flex;align-items:baseline;gap:14px;margin-bottom:16px;flex-wrap:wrap}';
    c{end+1} = '.section>.head h2{font-family:var(--serif);font-size:19px;font-weight:600;letter-spacing:-.01em;margin:0}';
    c{end+1} = '.section>.head .sub{font-size:12.5px;color:var(--muted)}';
    c{end+1} = '';
    c{end+1} = '/* ---- model card ---- */';
    c{end+1} = '.model-card{background:var(--panel);border:1px solid var(--line);border-radius:12px;box-shadow:0 1px 2px rgba(20,40,60,.04);overflow:hidden}';
    c{end+1} = '.formula{padding:22px 26px;border-bottom:1px solid var(--line-2);background:linear-gradient(180deg,#fbfdfe,#fff)}';
    c{end+1} = '.formula .fk{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600;display:block;margin-bottom:9px}';
    c{end+1} = '.formula code{font-family:var(--mono);font-size:19px;color:var(--ink);font-weight:500}';
    c{end+1} = '.formula code .v{color:var(--accent)}.formula code .r{color:var(--warm)}';
    c{end+1} = '.specs{display:grid;grid-template-columns:repeat(2,1fr)}';
    c{end+1} = '@media(max-width:680px){.specs{grid-template-columns:1fr}}';
    c{end+1} = '.specs .row{display:flex;justify-content:space-between;gap:16px;padding:11px 26px;border-bottom:1px solid var(--line-2);font-size:13.5px}';
    c{end+1} = '.specs .row:nth-child(odd){border-right:1px solid var(--line-2)}';
    c{end+1} = '.specs .k{color:var(--ink-2)}';
    c{end+1} = '.specs .val{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink);font-weight:500;text-align:right}';
    c{end+1} = '.specs .val.tag{font-family:var(--sans);font-weight:600}';
    c{end+1} = '';
    c{end+1} = '/* ---- legend ---- */';
    c{end+1} = '.legend{display:flex;flex-wrap:wrap;gap:12px 26px;align-items:center;padding:14px 18px;background:var(--panel);border:1px solid var(--line);border-radius:10px;font-size:12.5px;color:var(--ink-2)}';
    c{end+1} = '.legend .grp{display:flex;align-items:center;gap:9px}.legend .grp .k{font-weight:600;color:var(--ink)}';
    c{end+1} = '.tmapbar{width:120px;height:13px;border-radius:3px;border:1px solid rgba(0,0,0,.12)}';
    c{end+1} = '.legend .sep{width:1px;height:26px;background:var(--line)}';
    c{end+1} = '.legend .dd{width:10px;height:10px;border-radius:50%;flex:none}';
    c{end+1} = '.cnum{display:inline-flex;align-items:center;justify-content:center;width:17px;height:17px;border-radius:50%;background:#fff;border:1.5px solid var(--accent);color:var(--accent);font-family:var(--mono);font-size:10px;font-weight:700;flex:none}';
    c{end+1} = '';
    c{end+1} = '/* ---- figures ---- */';
    c{end+1} = '.fig{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;box-shadow:0 1px 2px rgba(20,40,60,.04)}';
    c{end+1} = '.fig>.cap{display:flex;align-items:baseline;justify-content:space-between;gap:10px;padding:12px 16px;border-bottom:1px solid var(--line-2)}';
    c{end+1} = '.fig>.cap .t{font-size:13px;font-weight:650}.fig>.cap .f{font-family:var(--mono);font-size:10.5px;color:var(--muted)}';
    c{end+1} = '.slot{position:relative;margin:0;background:#f1f5f7;aspect-ratio:1/1}';
    c{end+1} = '.slot.wide{aspect-ratio:16/10}';
    c{end+1} = '.slot img{position:absolute;inset:0;width:100%;height:100%;object-fit:contain;background:#fff;display:block}';
    c{end+1} = '.slot img.err{display:none}.slot img.err + .ph{display:flex}';
    c{end+1} = '.slot .ph{position:absolute;inset:0;display:none;flex-direction:column;align-items:center;justify-content:center;gap:6px;color:#9fb0bd;text-align:center;padding:12px}';
    c{end+1} = '.slot .ph .nd{font-size:12px;font-weight:600;color:#8b9aa7}';
    c{end+1} = '.slot .ph .fn{font-family:var(--mono);font-size:9px;color:#a7b6c2;word-break:break-all;max-width:92%}';
    c{end+1} = '.map-note{padding:9px 16px;border-top:1px solid var(--line-2);font-size:11px;color:var(--muted);display:flex;align-items:center;gap:14px;flex-wrap:wrap}';
    c{end+1} = '.map-note .dl{display:inline-flex;align-items:center;gap:6px}';
    c{end+1} = '.map-note .dd{width:9px;height:9px;border-radius:50%;border:1px solid rgba(0,0,0,.25);flex:none}';
    c{end+1} = '';
    c{end+1} = '/* ---- descriptive means grid ---- */';
    c{end+1} = '.cellgrid{display:grid;grid-template-columns:56px repeat(3,1fr);gap:8px;align-items:center;padding:18px}';
    c{end+1} = '.cellgrid .ch{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--ink-2);font-weight:650;text-align:center}';
    c{end+1} = '.cellgrid .rh{font-size:12px;font-weight:650;color:var(--ink-2);text-align:right;padding-right:6px}';
    c{end+1} = '.cellgrid .slot{border-radius:6px;overflow:hidden}';
    c{end+1} = '.desc-flex{display:grid;grid-template-columns:1.4fr 1fr;gap:22px;align-items:start}';
    c{end+1} = '@media(max-width:900px){.desc-flex{grid-template-columns:1fr}}';
    c{end+1} = '';
    c{end+1} = '/* ---- effect sections (collapsible) ---- */';
    c{end+1} = '.effect{background:var(--panel);border:1px solid var(--line);border-radius:12px;box-shadow:0 1px 2px rgba(20,40,60,.04);margin-top:26px;overflow:hidden;scroll-margin-top:70px}';
    c{end+1} = '.effect>.ehead{display:flex;align-items:center;gap:14px;padding:16px 20px;cursor:pointer;user-select:none}';
    c{end+1} = '.effect>.ehead:hover{background:#fbfdfe}';
    c{end+1} = '.effect>.ehead .chev{width:20px;height:20px;flex:none;color:var(--muted);transition:transform .2s}';
    c{end+1} = '.effect.open>.ehead .chev{transform:rotate(90deg)}';
    c{end+1} = '.effect>.ehead h2{font-family:var(--serif);font-size:18px;font-weight:600;margin:0}';
    c{end+1} = '.effect>.ehead .esub{font-size:12px;color:var(--muted)}';
    c{end+1} = '.effect>.ehead .lamp{width:9px;height:9px;border-radius:50%;flex:none;margin-left:auto}';
    c{end+1} = '.effect>.ehead .ecount{font-family:var(--mono);font-size:11.5px;color:var(--ink-2);font-variant-numeric:tabular-nums;white-space:nowrap}';
    c{end+1} = '.effect>.ebody{display:none;padding:6px 20px 24px;border-top:1px solid var(--line-2)}';
    c{end+1} = '.effect.open>.ebody{display:block}';
    c{end+1} = '.blk{margin-top:20px}';
    c{end+1} = '.blk>.bh{font-size:12px;font-weight:650;color:var(--ink-2);margin-bottom:10px;letter-spacing:.02em}';
    c{end+1} = '.effect-map{max-width:380px}';
    c{end+1} = '';
    c{end+1} = '/* ---- segmented toggles ---- */';
    c{end+1} = '.seg{display:inline-flex;background:var(--line-2);border-radius:9px;padding:4px;gap:3px;flex-wrap:wrap}';
    c{end+1} = '.seg button{font:inherit;font-size:13px;font-weight:600;color:var(--ink-2);background:transparent;border:0;border-radius:6px;padding:8px 16px;cursor:pointer;transition:.15s;display:flex;align-items:center;gap:8px}';
    c{end+1} = '.seg button.on{background:var(--panel);color:var(--accent);box-shadow:0 1px 2px rgba(20,40,60,.1)}';
    c{end+1} = '.seg button .ct{font-family:var(--mono);font-size:11px;font-weight:600;color:#fff;background:var(--muted);border-radius:20px;padding:1px 8px;min-width:22px;text-align:center}';
    c{end+1} = '.seg button.on .ct{background:var(--accent)}.seg button .ct.none{background:var(--line);color:var(--muted)}';
    c{end+1} = '.vpanel{margin-top:18px}';
    c{end+1} = '.vbody{display:grid;grid-template-columns:340px 1fr;gap:20px;align-items:start}';
    c{end+1} = '@media(max-width:860px){.vbody{grid-template-columns:1fr}}';
    c{end+1} = '';
    c{end+1} = '/* ---- tables ---- */';
    c{end+1} = '.tablewrap{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}';
    c{end+1} = '.tablewrap .thead-note{padding:11px 16px;border-bottom:1px solid var(--line-2);font-size:12px;color:var(--muted);display:flex;justify-content:space-between;gap:12px;flex-wrap:wrap}';
    c{end+1} = '.tablewrap .thead-note b{color:var(--ink)}';
    c{end+1} = 'table.sum{border-collapse:collapse;width:100%;font-size:13.5px}';
    c{end+1} = 'table.sum th,table.sum td{padding:8px 16px;text-align:right;border-bottom:1px solid var(--line-2)}';
    c{end+1} = 'table.sum th:first-child,table.sum td:first-child{text-align:left}';
    c{end+1} = 'table.sum thead th{font-size:11px;letter-spacing:.03em;color:var(--ink-2);font-weight:650;background:#fafcfd}';
    c{end+1} = 'table.sum tbody td{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink-2)}';
    c{end+1} = 'table.sum tbody td.band{font-family:var(--sans);color:var(--ink);font-weight:600;display:flex;align-items:center;gap:8px}';
    c{end+1} = 'table.sum tbody tr:last-child td{border-bottom:0}';
    c{end+1} = 'table.sum .p-sig{color:var(--sig);font-weight:700}table.sum .stat-v{color:var(--ink);font-weight:600}';
    c{end+1} = '.dotkey{width:9px;height:9px;border-radius:50%;flex:none}';
    c{end+1} = '.dotkey.unc{background:#1b2733}.dotkey.cor{background:#fff;border:1px solid var(--ink-2)}';
    c{end+1} = 'tr.sign-row td{background:#f3f6f9;color:var(--ink-2);font-family:var(--sans);font-weight:650;font-size:11px;text-align:left;letter-spacing:.06em;text-transform:uppercase}';
    c{end+1} = 'tr.sign-row td .sgn{font-family:var(--mono);text-transform:none;letter-spacing:0}';
    c{end+1} = 'tr.grp-row td{background:var(--accent-soft);color:var(--accent);font-family:var(--sans);font-weight:650;font-size:12px;text-align:left}';
    c{end+1} = 'tr.grp-row td .cp{font-family:var(--mono);font-weight:600;color:var(--ink-2);margin-left:8px}';
    c{end+1} = 'tr.grp-row td .cnum{margin-right:8px;vertical-align:-3px}';
    c{end+1} = 'tr.xrow{display:none}.tablewrap.expanded tr.xrow{display:table-row}';
    c{end+1} = '.morebtn{display:block;width:100%;font:inherit;font-size:12px;font-weight:600;color:var(--accent);background:#fbfdfe;border:0;border-top:1px solid var(--line-2);padding:10px;cursor:pointer}';
    c{end+1} = '.morebtn:hover{background:var(--accent-soft)}';
    c{end+1} = '.empty-note{padding:26px 18px;text-align:center;color:var(--muted);font-size:13px}';
    c{end+1} = '.empty-note b{color:var(--ink-2)}';
    c{end+1} = '';
    c{end+1} = '/* ---- post-hoc tabs ---- */';
    c{end+1} = '.tabs{display:flex;gap:0;border-bottom:1px solid var(--line);flex-wrap:wrap}';
    c{end+1} = '.tabs button{font:inherit;font-size:13px;font-weight:600;color:var(--ink-2);background:transparent;border:0;border-bottom:2px solid transparent;padding:9px 16px;cursor:pointer;display:flex;align-items:center;gap:8px;margin-bottom:-1px}';
    c{end+1} = '.tabs button.on{color:var(--accent);border-bottom-color:var(--accent)}';
    c{end+1} = '.tabs button .ct{font-family:var(--mono);font-size:10.5px;font-weight:600;color:#fff;background:var(--muted);border-radius:20px;padding:1px 7px}';
    c{end+1} = '.tabs button.on .ct{background:var(--accent)}.tabs button .ct.none{background:var(--line);color:var(--muted)}';
    c{end+1} = '.tabpanel{margin-top:18px}';
    c{end+1} = '';
    c{end+1} = '/* ---- global test ---- */';
    c{end+1} = '.global-card{display:flex;align-items:stretch;background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;flex-wrap:wrap;box-shadow:0 1px 2px rgba(20,40,60,.04)}';
    c{end+1} = '.gstat{padding:20px 26px;display:flex;flex-direction:column;gap:6px;min-width:140px}';
    c{end+1} = '.gstat .gk{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600}';
    c{end+1} = '.gstat .gv{font-family:var(--mono);font-size:24px;font-weight:600;color:var(--ink);font-variant-numeric:tabular-nums;line-height:1.1}';
    c{end+1} = '.gstat .gp{font-family:var(--mono);font-size:13px;color:var(--muted)}';
    c{end+1} = '.gstat .gp.p-sig{color:var(--sig);font-weight:700}';
    c{end+1} = '.gsep{width:1px;background:var(--line)}';
    c{end+1} = '.gnote{margin-left:auto;align-self:center;padding:18px 26px;font-size:12.5px;color:var(--ink-2);max-width:340px}.gnote b{color:var(--ink)}';
    c{end+1} = '.foot{margin-top:56px;padding-top:20px;border-top:1px solid var(--line);font-size:11.5px;color:var(--muted);display:flex;justify-content:space-between;gap:16px;flex-wrap:wrap}';
    c{end+1} = '@media(max-width:820px){.wrap{padding:28px 18px 72px}.masthead h1{font-size:24px}.strip{margin:0 -18px;padding:10px 18px}}';
    c{end+1} = '';
    c{end+1} = '/* ---- print ---- */';
    c{end+1} = '@media print{';
    c{end+1} = '  body{background:#fff}.wrap{max-width:none;padding:0}';
    c{end+1} = '  .actions,.seg,.morebtn,.tabs{display:none !important}';
    c{end+1} = '  .strip{position:static;background:#fff;margin:0;padding:10px 0}';
    c{end+1} = '  .effect>.ebody{display:block !important}';
    c{end+1} = '  .effect>.ehead .chev{display:none}';
    c{end+1} = '  .vpanel{display:block !important;margin-top:16px;page-break-inside:avoid}';
    c{end+1} = '  .vpanel::before{content:attr(data-view-name);display:block;font-weight:700;font-size:12px;color:var(--accent);margin-bottom:8px}';
    c{end+1} = '  .tabpanel{display:block !important;margin-top:16px;page-break-inside:avoid}';
    c{end+1} = '  .tabpanel::before{content:attr(data-tab-name);display:block;font-weight:700;font-size:12px;color:var(--accent);margin-bottom:8px}';
    c{end+1} = '  tr.xrow{display:none !important}';
    c{end+1} = '  .fig,.tablewrap,.model-card,.global-card,.legend,.effect{box-shadow:none;break-inside:avoid}';
    c{end+1} = '  .section{break-inside:avoid}';
    c{end+1} = '}';
    c{end+1} = '</style>';
    s = strjoin(c, newline);
end

function s = v2_js()
    c = {};
    c{end+1} = '<script>';
    c{end+1} = 'function showView(fx,w){';
    c{end+1} = '  var v=["clu","tfce","unc"];';
    c{end+1} = '  for(var i=0;i<v.length;i++){';
    c{end+1} = '    var p=document.getElementById("panel_"+fx+"_"+v[i]); if(p)p.style.display=(v[i]===w)?"block":"none";';
    c{end+1} = '    var b=document.getElementById("btn_"+fx+"_"+v[i]);   if(b)b.className=(v[i]===w)?"on":"";';
    c{end+1} = '  }';
    c{end+1} = '}';
    c{end+1} = 'function showTab(grp,w,names){';
    c{end+1} = '  for(var i=0;i<names.length;i++){';
    c{end+1} = '    var p=document.getElementById("tab_"+grp+"_"+names[i]); if(p)p.style.display=(names[i]===w)?"block":"none";';
    c{end+1} = '    var b=document.getElementById("tbtn_"+grp+"_"+names[i]); if(b)b.className=(names[i]===w)?"on":"";';
    c{end+1} = '  }';
    c{end+1} = '}';
    c{end+1} = 'function toggleFx(id){document.getElementById(id).classList.toggle("open");}';
    c{end+1} = 'function expandTable(btn){';
    c{end+1} = '  var tw=btn.closest(".tablewrap"); tw.classList.add("expanded"); btn.style.display="none";';
    c{end+1} = '}';
    c{end+1} = '</script>';
    s = strjoin(c, newline);
end
