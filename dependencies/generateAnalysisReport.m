function generateAnalysisReport(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch)
%GENERATEANALYSISREPORT  One self-contained, journal-styled HTML report.
% Writes <outputSname>_report.html for a SINGLE statistical comparison
% (no faceting): a calm, card-based page with the mean topographies and the
% statistic map, a three-way Uncorrected / TFCE / Cluster significance toggle
% (each view lists its significant channels with statistic + p-value), an
% optional post-hoc table for omnibus-F presets, and a short "Reading the
% statistics" explainer. Stat-aware (t / F / correlation r / circular).
%
% INPUTS  results_struct, params, base_filename, outputSname,
%         uncorrsigch       channels with uncorrected p <= .05
%         correctTFCEsigch  channels with TFCE-corrected p <= .05
%         SnPMsigch         channels in clusters with cluster p <= .05
%
% All values come from results_struct; nothing is recomputed.
%   p.real            uncorrected per-channel p
%   p.correctedTFCE   TFCE-corrected per-channel p
%   Clusters(k).channels / .p   cluster membership + cluster-level p
%   T.real_T          statistic per channel
%   chanlocs(i).labels  channel name (e.g. "E129")
%   glm.contrast_type ('t'/'F') / .effect_label ; posthoc(.label,...)

    report_filename = [outputSname '_report.html'];
    fid = fopen(report_filename, 'w');
    if fid == -1
        warning('generateAnalysisReport:open', 'Could not open %s for writing.', report_filename);
        return;
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    nperm = 1000;
    if isfield(params,'permutations') && ~isempty(params.permutations)
        nperm = double(params.permutations);
    end
    [statName, statSym] = stat_info(results_struct, params);

    % Multi-effect report (two-way mixed ANOVA): one section per effect (group,
    % condition, interaction), a descriptive block and a per-effect global
    % table. Single-effect presets fall through to the original layout below.
    if isfield(results_struct, 'effects') && numel(results_struct.effects) > 1
        write_head(fid);
        fprintf(fid, '<body>\n<div class="wrap">\n');
        write_masthead(fid, results_struct, params, base_filename, statName, nperm);
        write_legend_multi(fid);
        write_descriptive(fid, results_struct);
        for k = 1:numel(results_struct.effects)
            write_effect_one(fid, results_struct, results_struct.effects(k), k, nperm);
        end
        write_global_multi(fid, results_struct, nperm);
        write_footer(fid);
        fprintf(fid, '</div>\n</body>\n</html>\n');
        return;
    end

    write_head(fid);
    fprintf(fid, '<body>\n<div class="wrap">\n');

    write_masthead(fid, results_struct, params, base_filename, statName, nperm);
    write_legend(fid, statName);
    write_topographies(fid, results_struct, base_filename, statName);
    write_significance(fid, results_struct, base_filename, statSym, nperm, ...
        uncorrsigch, correctTFCEsigch, SnPMsigch);
    write_posthoc(fid, results_struct, nperm);
    write_global(fid, results_struct, statName, nperm);
    write_footer(fid);

    fprintf(fid, '</div>\n</body>\n</html>\n');
end

% =========================================================================
% STATISTIC RESOLUTION
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
        case 'circ', sym = 'stat';
        otherwise,   sym = 't';
    end
end

% =========================================================================
% HEAD : inline CSS + the (only) bit of JS, the view toggle
% =========================================================================
function write_head(fid)
    css = [ ...
'<!DOCTYPE html>' newline '<html lang="en"><head><meta charset="utf-8">' newline ...
'<meta name="viewport" content="width=device-width, initial-scale=1">' newline ...
'<title>SnPM Analysis Report</title>' newline '<style>' newline ...
':root{--ink:#1b2733;--ink-2:#445261;--muted:#7c8a98;--line:#e3e9ef;--line-2:#eef2f6;--bg:#f5f7f9;--panel:#ffffff;--accent:#1f6f8b;--accent-soft:#e8f1f4;--sig:#1f8a5b;--sig-soft:#e7f4ec;--warm:#c2543a;--cool:#2f6db0;--amber:#b8862a;--mono:ui-monospace,"SF Mono",Menlo,Consolas,"Roboto Mono",monospace;--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;--serif:Charter,"Iowan Old Style",Georgia,"Times New Roman",serif}' newline ...
'*{box-sizing:border-box}html{-webkit-text-size-adjust:100%}' newline ...
'body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility}' newline ...
'.wrap{max-width:1100px;margin:0 auto;padding:40px 32px 96px}' newline ...
'.eyebrow{font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);font-weight:600}' newline ...
'.masthead{display:flex;justify-content:space-between;align-items:flex-start;gap:32px;flex-wrap:wrap}' newline ...
'.masthead h1{font-family:var(--serif);font-weight:600;font-size:29px;line-height:1.2;letter-spacing:-.01em;margin:9px 0 11px;max-width:760px;color:var(--ink)}' newline ...
'.methods{color:var(--ink-2);font-size:13.5px;max-width:800px;margin:0}.methods b{color:var(--ink);font-weight:600}' newline ...
'.actions{display:flex;gap:8px;flex-shrink:0}' newline ...
'.btn{font:inherit;font-size:12.5px;font-weight:600;color:var(--ink-2);background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:8px 13px;cursor:pointer;display:inline-flex;align-items:center;gap:6px;transition:border-color .15s,color .15s}' newline ...
'.btn:hover{border-color:var(--accent);color:var(--accent)}.btn svg{width:14px;height:14px}' newline ...
'.legend{display:flex;flex-wrap:wrap;gap:10px 26px;align-items:center;margin:22px 0 4px;padding:14px 18px;background:var(--panel);border:1px solid var(--line);border-radius:10px;font-size:12.5px;color:var(--ink-2)}' newline ...
'.legend .grp{display:flex;align-items:center;gap:9px}.legend .grp .k{font-weight:600;color:var(--ink)}' newline ...
'.chip{display:inline-flex;align-items:center;gap:6px}' newline ...
'.swatch{width:13px;height:13px;border-radius:3px;border:1px solid rgba(0,0,0,.08);flex:none}' newline ...
'.tmapbar{width:120px;height:13px;border-radius:3px;border:1px solid rgba(0,0,0,.1)}' newline ...
'.legend .sep{width:1px;height:26px;background:var(--line)}' newline ...
'.section{margin-top:46px}' newline ...
'.section>.head{display:flex;align-items:baseline;gap:14px;margin-bottom:16px;flex-wrap:wrap}' newline ...
'.section>.head h2{font-size:18px;font-weight:650;letter-spacing:-.01em;margin:0}' newline ...
'.section>.head .sub{font-size:12.5px;color:var(--muted)}' newline ...
'.topo-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:20px}' newline ...
'@media(max-width:820px){.topo-grid{grid-template-columns:1fr}}' newline ...
'.fig{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}' newline ...
'.fig>.cap{display:flex;align-items:baseline;justify-content:space-between;gap:10px;padding:12px 16px;border-bottom:1px solid var(--line-2)}' newline ...
'.fig>.cap .t{font-size:13px;font-weight:650}.fig>.cap .f{font-family:var(--mono);font-size:10.5px;color:var(--muted)}' newline ...
'.slot{position:relative;margin:0;background:#f1f5f7;aspect-ratio:1/1}' newline ...
'.slot img{position:absolute;inset:0;width:100%;height:100%;object-fit:contain;background:#fff;display:block}' newline ...
'.slot img.err{display:none}' newline ...
'.slot .ph{position:absolute;inset:0;display:flex;flex-direction:column;align-items:center;justify-content:center;gap:8px;color:#9fb0bd;text-align:center;padding:14px}' newline ...
'.slot .ph svg{width:30px;height:30px;opacity:.6}' newline ...
'.slot .ph .fn{font-family:var(--mono);font-size:10px;color:#a7b6c2;word-break:break-all;max-width:90%}' newline ...
'.seg{display:inline-flex;background:var(--line-2);border-radius:9px;padding:4px;gap:3px;flex-wrap:wrap}' newline ...
'.seg button{font:inherit;font-size:13px;font-weight:600;color:var(--ink-2);background:transparent;border:0;border-radius:6px;padding:9px 18px;cursor:pointer;transition:.15s;display:flex;align-items:center;gap:8px}' newline ...
'.seg button.on{background:var(--panel);color:var(--accent);box-shadow:0 1px 2px rgba(20,40,60,.1)}' newline ...
'.seg button .ct{font-family:var(--mono);font-size:11px;font-weight:600;color:#fff;background:var(--muted);border-radius:20px;padding:1px 8px;min-width:22px;text-align:center}' newline ...
'.seg button.on .ct{background:var(--accent)}.seg button .ct.none{background:var(--line);color:var(--muted)}' newline ...
'.vpanel{margin-top:24px;animation:fade .25s ease}' newline ...
'@keyframes fade{from{opacity:0;transform:translateY(4px)}to{opacity:1;transform:none}}' newline ...
'.vbody{display:grid;grid-template-columns:380px 1fr;gap:22px;align-items:start}' newline ...
'@media(max-width:860px){.vbody{grid-template-columns:1fr}}' newline ...
'.tablewrap{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden}' newline ...
'.tablewrap .thead-note{padding:11px 16px;border-bottom:1px solid var(--line-2);font-size:12px;color:var(--muted);display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap}' newline ...
'.tablewrap .thead-note b{color:var(--ink)}' newline ...
'table.sum{border-collapse:collapse;width:100%;font-size:13.5px}' newline ...
'table.sum th,table.sum td{padding:9px 16px;text-align:right;border-bottom:1px solid var(--line-2)}' newline ...
'table.sum th:first-child,table.sum td:first-child{text-align:left}' newline ...
'table.sum thead th{font-size:11px;letter-spacing:.03em;color:var(--ink-2);font-weight:650;white-space:nowrap;background:#fafcfd}' newline ...
'table.sum tbody td{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink-2)}' newline ...
'table.sum tbody td.band{font-family:var(--sans);color:var(--ink);font-weight:600;white-space:nowrap}' newline ...
'table.sum tbody tr.sig td.band{box-shadow:inset 3px 0 0 var(--sig)}' newline ...
'table.sum tbody tr:last-child td{border-bottom:0}' newline ...
'table.sum .p-sig{color:var(--sig);font-weight:700}table.sum .stat-v{color:var(--ink);font-weight:600}' newline ...
'tr.grp-row td{background:var(--accent-soft);color:var(--accent);font-family:var(--sans);font-weight:650;font-size:12px;text-align:left;letter-spacing:.01em}' newline ...
'tr.grp-row td .cp{font-family:var(--mono);font-weight:600;color:var(--ink-2);margin-left:8px}' newline ...
'.empty-note{padding:34px 18px;text-align:center;color:var(--muted);font-size:13.5px}' newline ...
'.empty-note svg{width:26px;height:26px;opacity:.5;display:block;margin:0 auto 8px}' newline ...
'.global-card{display:flex;align-items:stretch;background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;flex-wrap:wrap}' newline ...
'.gstat{padding:22px 30px;display:flex;flex-direction:column;gap:7px;min-width:160px}' newline ...
'.gstat .gk{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600}' newline ...
'.gstat .gv{font-family:var(--mono);font-size:30px;font-weight:600;color:var(--ink);font-variant-numeric:tabular-nums;line-height:1}' newline ...
'.gstat .gv.p-sig{color:var(--sig)}' newline ...
'.gsep{width:1px;background:var(--line)}' newline ...
'.gnote{margin-left:auto;align-self:center;padding:20px 30px;font-size:12.5px;color:var(--ink-2);max-width:340px;line-height:1.5}.gnote b{color:var(--ink)}' newline ...
'.map-note{padding:9px 16px;border-top:1px solid var(--line-2);font-size:11px;color:var(--muted);display:flex;align-items:center;gap:14px;flex-wrap:wrap}' newline ...
'.map-note .dl{display:inline-flex;align-items:center;gap:6px}' newline ...
'.map-note .dd{width:9px;height:9px;border-radius:50%;border:1px solid rgba(0,0,0,.25);flex:none}' newline ...
'.foot{margin-top:60px;padding-top:20px;border-top:1px solid var(--line);font-size:11.5px;color:var(--muted);display:flex;justify-content:space-between;gap:16px;flex-wrap:wrap}' newline ...
'@media(max-width:820px){.wrap{padding:28px 18px 72px}.masthead h1{font-size:24px}}' newline ...
'@media print{:root{--bg:#fff}body{font-size:11px}.wrap{max-width:none;padding:0}.actions,.seg{display:none!important}.vpanel{display:block!important;animation:none;break-inside:avoid;margin-top:22px}.vpanel::before{content:attr(data-view-name);display:block;font-family:var(--serif);font-size:16px;font-weight:600;margin:0 0 12px;color:#000}.fig,.tablewrap,.legend,.global-card{break-inside:avoid;box-shadow:none}}' newline ...
'</style>' newline ...
'<script>function showView(w){var v=["unc","tfce","clu"];for(var i=0;i<v.length;i++){var k=v[i];document.getElementById("panel_"+k).style.display=(k===w)?"block":"none";document.getElementById("btn_"+k).className=(k===w)?"on":"";}}' newline ...
'function showViewG(g,w){var v=["unc","tfce","clu"];for(var i=0;i<v.length;i++){var k=v[i];var p=document.getElementById("panel_"+g+"_"+k);if(p)p.style.display=(k===w)?"block":"none";var b=document.getElementById("btn_"+g+"_"+k);if(b)b.className=(k===w)?"on":"";}}</script>' newline ...
'</head>' newline];
    fprintf(fid, '%s', css);
end

% =========================================================================
% MASTHEAD
% =========================================================================
function write_masthead(fid, rs, params, base, statName, nperm)
    title = strrep(base, '_', ' ');
    fprintf(fid, '<header class="masthead">\n  <div>\n');
    fprintf(fid, '    <span class="eyebrow">SnPM Topographic Statistics</span>\n');
    fprintf(fid, '    <h1>%s</h1>\n', esc(title));
    fprintf(fid, '    <p class="methods">%s</p>\n', methods_line(rs, params, statName, nperm));
    fprintf(fid, '  </div>\n');
    fprintf(fid, ['  <div class="actions"><button class="btn" onclick="window.print()">' ...
        '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8">' ...
        '<path d="M6 9V3h12v6M6 18H4a2 2 0 0 1-2-2v-4a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v4a2 2 0 0 1-2 2h-2M6 14h12v7H6z"/>' ...
        '</svg>Print / PDF</button></div>\n']);
    fprintf(fid, '</header>\n');
end

function s = methods_line(rs, params, statName, nperm)
    % comparison phrase + statistical word
    comp = '';
    if isfield(params,'comparison') && ischar(params.comparison), comp = lower(params.comparison); end
    switch statName
        case 'F',    statword = 'omnibus F-test';
        case 'r',    statword = 'correlation';
        case 'circ', statword = 'circular statistic';
        otherwise,   statword = 't-test';
    end
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
    end

    % tail
    tail = 'two-tailed';
    if isfield(params,'tail') && ~isempty(params.tail) && ischar(params.tail)
        if ~strcmpi(params.tail,'both'), tail = 'one-tailed'; end
    end
    % data type
    dt = '';
    if isfield(params,'datatype') && ~isempty(params.datatype) && ischar(params.datatype)
        dt = [' of the ' esc(params.datatype) ' data'];
    end

    s = sprintf('%s <b>%s</b>%s, tested channel by channel. <b>%s permutations</b>, %s, family-wise corrected with <b>TFCE</b> and <b>cluster-based</b> methods at <b>&alpha; = .05</b>.', ...
        phrase, statword, dt, fmt_int(nperm), tail);

    nc = n_clause(rs);
    if ~isempty(nc), s = [s ' &nbsp;<b>' nc '</b>.']; end

    if isfield(rs,'glm') && isfield(rs.glm,'effect_label') && ~isempty(rs.glm.effect_label)
        s = [s ' &nbsp;<span style="color:var(--muted)">Effect: ' esc(rs.glm.effect_label) '</span>'];
    end
end

function nc = n_clause(rs)
    nc = '';
    if ~isfield(rs,'data_summary'), return; end
    ds = rs.data_summary;
    n1 = []; n2 = [];
    if isfield(ds,'data1_size') && ~isempty(ds.data1_size), n1 = ds.data1_size(1); end
    if isfield(ds,'data2_size') && ~isempty(ds.data2_size), n2 = ds.data2_size(1); end
    if ~isempty(n1) && ~isempty(n2) && n2 ~= n1
        nc = sprintf('n = %d vs %d', n1, n2);
    elseif ~isempty(n1)
        nc = sprintf('n = %d', n1);
    end
end

% =========================================================================
% LEGEND
% =========================================================================
function write_legend(fid, statName)
    fprintf(fid, '<div class="legend">\n');
    if strcmp(statName,'F')
        bar = 'linear-gradient(90deg,#f3f6f8,#0e4d63)';
        fprintf(fid, '  <div class="grp"><span class="k">Statistic</span><span>F-statistic (omnibus)</span></div>\n');
        fprintf(fid, '  <span class="sep"></span>\n');
        fprintf(fid, '  <div class="grp"><span class="k">F-map</span><span class="chip"><span class="tmapbar" style="background:%s"></span></span><span style="color:var(--muted)">low</span><span style="color:var(--muted)">&middot;</span><span style="color:var(--warm)">high F</span></div>\n', bar);
    elseif strcmp(statName,'r')
        bar = 'linear-gradient(90deg,#2f6db0,#cfe0ef 46%,#f2efe6 50%,#edc6b6 54%,#c2543a)';
        fprintf(fid, '  <div class="grp"><span class="k">Statistic</span><span>correlation (r)</span></div>\n');
        fprintf(fid, '  <span class="sep"></span>\n');
        fprintf(fid, '  <div class="grp"><span class="k">r-map</span><span class="chip"><span class="tmapbar" style="background:%s"></span></span><span style="color:var(--cool)">negative</span><span style="color:var(--muted)">&middot;</span><span style="color:var(--warm)">positive</span></div>\n', bar);
    else
        bar = 'linear-gradient(90deg,#2f6db0,#cfe0ef 46%,#f2efe6 50%,#edc6b6 54%,#c2543a)';
        fprintf(fid, '  <div class="grp"><span class="k">Statistic</span><span>t-statistic &middot; A vs B</span></div>\n');
        fprintf(fid, '  <span class="sep"></span>\n');
        fprintf(fid, '  <div class="grp"><span class="k">T-map</span><span class="chip"><span class="tmapbar" style="background:%s"></span></span><span style="color:var(--cool)">B&gt;A</span><span style="color:var(--muted)">&middot;</span><span style="color:var(--warm)">warm = A&gt;B</span></div>\n', bar);
    end
    fprintf(fid, '  <span class="sep"></span>\n');
    fprintf(fid, '  <div class="grp"><span class="k">Significant</span><span class="chip"><span class="swatch" style="background:var(--sig-soft);border-color:var(--sig)"></span>corrected channel</span></div>\n');
    fprintf(fid, '</div>\n');
end

% =========================================================================
% TOPOGRAPHIES
% =========================================================================
function write_topographies(fid, rs, base, statName)
    % Omnibus-F presets: one mean topography per group + the standalone F-map
    % (the per-group means and F-map are saved by core_snpm_glm).
    if strcmp(statName, 'F') && has_group_topos(rs)
        write_topographies_F(fid, rs);
        return;
    end
    % Regression (continuous predictor): no A/B means exist, so show only the
    % single signed slope map (core_snpm_glm sets glm.effect_map_only).
    if isfield(rs,'glm') && isfield(rs.glm,'effect_map_only') && rs.glm.effect_map_only
        effTitle = 'Effect (t-map)';
        if isfield(rs.glm,'effect_label') && ~isempty(rs.glm.effect_label)
            effTitle = [rs.glm.effect_label ' (t-map)'];
        end
        fprintf(fid, '<section class="section">\n  <div class="head"><h2>Topography</h2><span class="sub">Per-channel effect (slope) map</span></div>\n');
        fprintf(fid, '  <div class="topo-grid">\n');
        write_fig_card(fid, [base '_topo.png'], esc(effTitle));
        fprintf(fid, '  </div>\n</section>\n');
        return;
    end
    % t / r presets: unchanged two-means + effect-map layout.
    if has(base, ' VS ')
        parts = strsplit(base, ' VS '); d1 = parts{1}; d2 = parts{2};
    else
        d1 = [base '_Data1']; d2 = [base '_Data2'];
    end
    switch statName
        case 'F',    statTitle = 'F-map (omnibus)';
        case 'r',    statTitle = 'r-map';
        case 'circ', statTitle = 'Statistic map';
        otherwise,   statTitle = 't-map (A vs B)';
    end
    % One-sample vs 0: condition B is a zero map -> show only the condition
    % mean + the effect (t-vs-0) map.
    hideB = isfield(rs,'hide_condition_b') && rs.hide_condition_b;
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Topographies</h2><span class="sub">Group / condition means and the statistic map</span></div>\n');
    fprintf(fid, '  <div class="topo-grid">\n');
    if hideB
        write_fig_card(fid, [d1 '_topo.png'], 'Condition mean');
        write_fig_card(fid, [base '_topo.png'], 'Effect (t vs 0)');
    else
        write_fig_card(fid, [d1 '_topo.png'], 'Mean topography &mdash; A');
        write_fig_card(fid, [d2 '_topo.png'], 'Mean topography &mdash; B');
        write_fig_card(fid, [base '_topo.png'], statTitle);
    end
    fprintf(fid, '  </div>\n</section>\n');
end

function tf = has_group_topos(rs)
    tf = isfield(rs,'glm') && isfield(rs.glm,'group_mean_png') && ~isempty(rs.glm.group_mean_png);
end

% F (omnibus): one card per group mean, then the F-map card.
function write_topographies_F(fid, rs)
    labels = {}; pngs = {}; fmap = '';
    if isfield(rs.glm,'group_labels'),   labels = rs.glm.group_labels;   end
    if isfield(rs.glm,'group_mean_png'), pngs   = rs.glm.group_mean_png; end
    if isfield(rs.glm,'fmap_png'),       fmap   = rs.glm.fmap_png;       end

    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Topographies</h2><span class="sub">Per-group means and the omnibus F-map</span></div>\n');
    fprintf(fid, '  <div class="topo-grid">\n');
    for i = 1:numel(pngs)
        lab = sprintf('Group %d', i);
        if i <= numel(labels) && ~isempty(labels{i}), lab = labels{i}; end
        write_fig_card(fid, pngs{i}, sprintf('Mean topography &mdash; %s', esc(lab)));
    end
    if ~isempty(fmap)
        write_fig_card(fid, fmap, 'F-map (omnibus)');
    end
    fprintf(fid, '  </div>\n</section>\n');
end

function write_fig_card(fid, fname, title)
    fprintf(fid, '    <div class="fig"><div class="cap"><span class="t">%s</span><span class="f">%s</span></div>', title, esc(fname));
    write_slot(fid, fname);
    fprintf(fid, '</div>\n');
end

function write_slot(fid, fname)
    fprintf(fid, ['<figure class="slot"><img src="%s" alt="%s" loading="lazy" ' ...
        'onload="var p=this.parentNode.querySelector(''.ph'');if(p)p.style.display=''none''" ' ...
        'onerror="this.classList.add(''err'')">' ...
        '<figcaption class="ph"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">' ...
        '<rect x="3" y="3" width="18" height="18" rx="2"/><circle cx="8.5" cy="8.5" r="1.6"/><path d="M21 15l-5-5L5 21"/></svg>' ...
        '<span class="fn">%s</span></figcaption></figure>'], esc(fname), esc(fname), esc(fname));
end

% =========================================================================
% SIGNIFICANCE : three-way toggle (Uncorrected / TFCE / Cluster)
% =========================================================================
function write_significance(fid, rs, base, statSym, nperm, unc, tfce, clu)
    nU = numel(unc); nT = numel(tfce); nC = numel(clu);
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, '  <div class="head"><h2>Significant channels</h2><span class="sub">Choose the multiple-comparison correction</span></div>\n');
    fprintf(fid, '  <div class="seg" id="viewSeg">\n');
    fprintf(fid, '    <button id="btn_clu" class="on" onclick="showView(''clu'')">Cluster <span class="ct%s">%d</span></button>\n', none_cls(nC), nC);
    fprintf(fid, '    <button id="btn_tfce" onclick="showView(''tfce'')">TFCE <span class="ct%s">%d</span></button>\n', none_cls(nT), nT);
    fprintf(fid, '    <button id="btn_unc" onclick="showView(''unc'')">Uncorrected <span class="ct%s">%d</span></button>\n', none_cls(nU), nU);
    fprintf(fid, '  </div>\n');

    % --- Cluster panel (shown by default) ---
    fprintf(fid, '  <div id="panel_clu" class="vpanel" data-view-name="Cluster-corrected">\n    <div class="vbody">\n');
    write_map_card(fid, [base ' Cluster.png'], 'cluster-corrected');
    write_cluster_table(fid, rs, clu, statSym, nperm);
    fprintf(fid, '    </div>\n  </div>\n');

    % --- TFCE panel ---
    fprintf(fid, '  <div id="panel_tfce" class="vpanel" data-view-name="TFCE-corrected" style="display:none">\n    <div class="vbody">\n');
    write_map_card(fid, [base ' TFCE.png'], 'TFCE-corrected');
    write_flat_table(fid, rs, tfce, get_pvec(rs,'correctedTFCE'), statSym, 'TFCE-corrected', nperm);
    fprintf(fid, '    </div>\n  </div>\n');

    % --- Uncorrected panel : no map (uncorrected electrodes already show as
    %     black dots on the TFCE / Cluster maps); full-width channel list ---
    fprintf(fid, '  <div id="panel_unc" class="vpanel" data-view-name="Uncorrected" style="display:none">\n');
    write_flat_table(fid, rs, unc, get_pvec(rs,'real'), statSym, 'uncorrected', nperm);
    fprintf(fid, '  </div>\n');

    fprintf(fid, '</section>\n');
end

function write_map_card(fid, fname, kind)
    fprintf(fid, '      <div class="fig"><div class="cap"><span class="t">%s significance map</span><span class="f">%s</span></div>', cap_kind(kind), esc(fname));
    write_slot(fid, fname);
    fprintf(fid, ['<div class="map-note"><span class="dl"><span class="dd" style="background:#1b2733"></span>uncorrected</span>' ...
        '<span class="dl"><span class="dd" style="background:#fff"></span>%s</span></div>'], kind);
    fprintf(fid, '</div>\n');
end

function t = cap_kind(kind)
    if strncmpi(kind,'TFCE',4), t = 'TFCE'; else, t = 'Cluster'; end
end

% Post-hoc significance-map card. Same as write_map_card but the legend shows
% the effect DIRECTION instead of uncorrected/corrected (the dot meaning is
% already printed on the PNG as "B: uncorrectedSig / W: correctedSig"). The map
% is a signed t = mean_a - mean_b with a diverging jet colormap, so warm (red) =
% A > B and cool (blue) = A < B, for label "A vs B".
function write_posthoc_map_card(fid, fname, kind, lab)
    parts = strsplit(lab, ' vs ');
    if numel(parts) == 2, A = strtrim(parts{1}); B = strtrim(parts{2}); else, A = ''; B = ''; end
    fprintf(fid, '      <div class="fig"><div class="cap"><span class="t">%s significance map</span><span class="f">%s</span></div>', cap_kind(kind), esc(fname));
    write_slot(fid, fname);
    if ~isempty(A) && ~isempty(B)
        fprintf(fid, ['<div class="map-note"><span class="dl"><span class="dd" style="background:#c0392b"></span>%s &gt; %s</span>' ...
            '<span class="dl"><span class="dd" style="background:#2c6fbb"></span>%s &lt; %s</span></div>'], ...
            esc(A), esc(B), esc(A), esc(B));
    end
    fprintf(fid, '</div>\n');
end

% ---- flat table (uncorrected / TFCE) : Channel | statistic | p, ascending p ----
function write_flat_table(fid, rs, chs, pvec, statSym, label, nperm)
    chs = chs(:).';
    if isempty(chs)
        fprintf(fid, '      <div class="tablewrap">'); write_empty(fid, label); fprintf(fid, '</div>\n');
        return;
    end
    ps = arrayfun(@(c) get_p(pvec, c), chs);
    [ps, ord] = sort(ps);
    chs = chs(ord);
    fprintf(fid, '      <div class="tablewrap">\n');
    fprintf(fid, '        <div class="thead-note"><span><b>%d</b> %s significant channels</span><span>sorted by ascending p</span></div>\n', numel(chs), label);
    fprintf(fid, '        <table class="sum"><thead><tr><th>Channel</th><th>%s</th><th>p</th></tr></thead><tbody>\n', statSym);
    for k = 1:numel(chs)
        c = chs(k);
        pcls = ''; if ps(k) <= 0.05, pcls = 'p-sig'; end
        fprintf(fid, '          <tr class="sig"><td class="band">%s</td><td class="stat-v">%s</td><td class="%s">%s</td></tr>\n', ...
            esc(chan_label(rs,c)), fmt_stat(stat_val(rs,c)), pcls, fmt_p(ps(k), nperm));
    end
    fprintf(fid, '        </tbody></table>\n      </div>\n');
end

% ---- cluster table : channels grouped under their cluster, cluster-level p ----
function write_cluster_table(fid, rs, clu, statSym, nperm)
    have = isfield(rs,'Clusters') && ~isempty(rs.Clusters);
    idx = []; cps = [];
    if have
        C = rs.Clusters;
        nCl = numel(C);
        cps = nan(1,nCl); keep = false(1,nCl);
        for i = 1:nCl
            pv = cluster_p(C, i);
            ch = cluster_chs(C, i);
            if ~isnan(pv) && pv <= 0.05 && ~isempty(ch)
                cps(i) = pv; keep(i) = true;
            end
        end
        idx = find(keep);
        [~, ord] = sort(cps(idx));
        idx = idx(ord);
    end

    if isempty(idx)
        % Fallback: no per-cluster structure but channels supplied -> flat list
        if ~isempty(clu)
            write_flat_table(fid, rs, clu, get_pvec(rs,'real'), statSym, 'cluster-corrected', nperm);
        else
            fprintf(fid, '      <div class="tablewrap">'); write_empty(fid, 'cluster-corrected'); fprintf(fid, '</div>\n');
        end
        return;
    end

    total = 0;
    for j = 1:numel(idx), total = total + numel(cluster_chs(rs.Clusters, idx(j))); end
    fprintf(fid, '      <div class="tablewrap">\n');
    fprintf(fid, '        <div class="thead-note"><span><b>%d</b> channels in <b>%d</b> significant cluster(s)</span><span>grouped by cluster, ascending cluster p</span></div>\n', total, numel(idx));
    fprintf(fid, '        <table class="sum"><thead><tr><th>Channel</th><th>%s</th><th>cluster p</th></tr></thead><tbody>\n', statSym);
    for j = 1:numel(idx)
        ci = idx(j);
        ch = cluster_chs(rs.Clusters, ci);
        cp = cluster_p(rs.Clusters, ci);
        % order channels within the cluster by descending |statistic|
        sv = arrayfun(@(c) abs(stat_val(rs,c)), ch);
        [~, o] = sort(sv, 'descend');
        ch = ch(o);
        fprintf(fid, '          <tr class="grp-row"><td colspan="3">Cluster %d<span class="cp">p = %s &middot; %d ch</span></td></tr>\n', ...
            j, fmt_p(cp, nperm), numel(ch));
        for k = 1:numel(ch)
            c = ch(k);
            fprintf(fid, '          <tr class="sig"><td class="band">%s</td><td class="stat-v">%s</td><td class="p-sig">%s</td></tr>\n', ...
                esc(chan_label(rs,c)), fmt_stat(stat_val(rs,c)), fmt_p(cp, nperm));
        end
    end
    fprintf(fid, '        </tbody></table>\n      </div>\n');
end

function write_empty(fid, label)
    fprintf(fid, ['<div class="empty-note"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6">' ...
        '<circle cx="12" cy="12" r="9"/><path d="M8 12h8"/></svg>No %s significant channels at &alpha; = .05.</div>'], label);
end

% =========================================================================
% POST-HOC (omnibus-F presets)
% =========================================================================
function write_posthoc(fid, rs, nperm)
    if ~isfield(rs,'posthoc') || isempty(rs.posthoc), return; end
    ph = rs.posthoc;
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, '  <div class="head"><h2>Post-hoc pairwise contrasts</h2><span class="sub">Significant channels per comparison &mdash; choose the correction</span></div>\n');
    for i = 1:numel(ph)
        write_posthoc_one(fid, rs, ph(i), i, nperm);
    end
    fprintf(fid, '</section>\n');
end

% One pairwise contrast: a Cluster / TFCE / Uncorrected toggle whose panels
% reuse the omnibus channel-list helpers. Pairwise contrasts are 2-group t
% tests, so the statistic shown is t.
function write_posthoc_one(fid, rs, phe, idx, nperm, prefix)
    if nargin < 6 || isempty(prefix)
        gid = sprintf('ph%d', idx);
    else
        gid = sprintf('%s_ph%d', prefix, idx);
    end
    lab = '';
    if isfield(phe,'label') && ~isempty(phe.label), lab = phe.label; end

    unc  = getfield_or(phe, 'uncorrsigch', []);
    tfce = getfield_or(phe, 'correctTFCEsigch', []);
    clu  = getfield_or(phe, 'SnPMsigch', []);
    nU = numel(unc); nT = numel(tfce); nC = numel(clu);

    % Mini results struct so the existing table helpers work unchanged: they
    % read .T.real_T, .p.real/.correctedTFCE, .Clusters, .chanlocs.
    phrs = struct();
    if isfield(phe,'T'),        phrs.T = phe.T;               end
    if isfield(phe,'p'),        phrs.p = phe.p;               end
    if isfield(phe,'Clusters'), phrs.Clusters = phe.Clusters; end
    if isfield(rs,'chanlocs'),  phrs.chanlocs = rs.chanlocs;  end

    statSym = 't';
    % Base name of the per-contrast significance maps: "<base> Cluster.png" /
    % "<base> TFCE.png", produced by TopoplotSignificant_single (black dots =
    % uncorrected-sig, white dots = correction-sig). Empty -> table-only.
    mbase = getfield_or(phe, 'map_base', '');

    fprintf(fid, '  <div style="margin-top:30px">\n');
    fprintf(fid, ['    <div style="display:flex;align-items:baseline;gap:12px;margin-bottom:12px">' ...
        '<span style="font-family:var(--serif);font-size:16px;font-weight:600;color:var(--ink)">%s</span>' ...
        '<span style="font-size:12px;color:var(--muted)">pairwise t-contrast</span></div>\n'], esc(lab));

    fprintf(fid, '    <div class="seg" id="seg_%s">\n', gid);
    fprintf(fid, '      <button id="btn_%s_clu" class="on" onclick="showViewG(''%s'',''clu'')">Cluster <span class="ct%s">%d</span></button>\n', gid, gid, none_cls(nC), nC);
    fprintf(fid, '      <button id="btn_%s_tfce" onclick="showViewG(''%s'',''tfce'')">TFCE <span class="ct%s">%d</span></button>\n', gid, gid, none_cls(nT), nT);
    fprintf(fid, '      <button id="btn_%s_unc" onclick="showViewG(''%s'',''unc'')">Uncorrected <span class="ct%s">%d</span></button>\n', gid, gid, none_cls(nU), nU);
    fprintf(fid, '    </div>\n');

    % Cluster panel (default shown): significance map beside the cluster table.
    fprintf(fid, '    <div id="panel_%s_clu" class="vpanel" data-view-name="%s &mdash; Cluster-corrected">\n      <div class="vbody">\n', gid, esc(lab));
    if ~isempty(mbase), write_posthoc_map_card(fid, [mbase ' Cluster.png'], 'cluster-corrected', lab); end
    write_cluster_table(fid, phrs, clu, statSym, nperm);
    fprintf(fid, '      </div>\n    </div>\n');
    % TFCE panel
    fprintf(fid, '    <div id="panel_%s_tfce" class="vpanel" data-view-name="%s &mdash; TFCE-corrected" style="display:none">\n      <div class="vbody">\n', gid, esc(lab));
    if ~isempty(mbase), write_posthoc_map_card(fid, [mbase ' TFCE.png'], 'TFCE-corrected', lab); end
    write_flat_table(fid, phrs, tfce, get_pvec(phrs,'correctedTFCE'), statSym, 'TFCE-corrected', nperm);
    fprintf(fid, '      </div>\n    </div>\n');
    % Uncorrected panel: no map (its channels already show as black dots on the
    % Cluster / TFCE maps), full-width channel list.
    fprintf(fid, '    <div id="panel_%s_unc" class="vpanel" data-view-name="%s &mdash; Uncorrected" style="display:none">\n', gid, esc(lab));
    write_flat_table(fid, phrs, unc, get_pvec(phrs,'real'), statSym, 'uncorrected', nperm);
    fprintf(fid, '    </div>\n');

    fprintf(fid, '  </div>\n');
end

function v = getfield_or(s, f, dflt)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

% =========================================================================
% GLOBAL TEST : statistic + p on the channel-averaged signal (all good chans)
% =========================================================================
function write_global(fid, rs, statName, nperm)
    hasS = isfield(rs,'global_stat') && ~isempty(rs.global_stat);
    hasP = isfield(rs,'global_pval') && ~isempty(rs.global_pval);
    if ~hasS && ~hasP, return; end

    switch statName
        case 'F',    slab = 'F-statistic';
        case 'r',    slab = 'correlation r';
        case 'circ', slab = 'statistic';
        otherwise,   slab = 't-statistic';
    end

    if hasS, sval = fmt_stat(rs.global_stat(1)); else, sval = '&mdash;'; end
    if hasP
        pv = rs.global_pval(1);
        pval = fmt_p(pv, nperm);
        if pv <= 0.05, pcls = ' p-sig'; else, pcls = ''; end
    else
        pval = '&mdash;'; pcls = '';
    end

    % number of good channels averaged over
    nch = [];
    if isfield(rs,'T') && isfield(rs.T,'real_T') && ~isempty(rs.T.real_T)
        nch = numel(rs.T.real_T);
    elseif isfield(rs,'chanlocs')
        nch = numel(rs.chanlocs);
    end
    if ~isempty(nch)
        note = sprintf('Statistic and permutation p-value computed on the mean across all <b>%d</b> good channels &mdash; a single whole-head test, complementary to the channel-wise maps above.', nch);
    else
        note = 'Statistic and permutation p-value computed on the mean across all good channels &mdash; a single whole-head test, complementary to the channel-wise maps above.';
    end

    fprintf(fid, '<section class="section">\n');
    fprintf(fid, '  <div class="head"><h2>Global test</h2><span class="sub">Channel-averaged statistic</span></div>\n');
    fprintf(fid, '  <div class="global-card">\n');
    fprintf(fid, '    <div class="gstat"><div class="gk">%s</div><div class="gv">%s</div></div>\n', slab, sval);
    fprintf(fid, '    <div class="gsep"></div>\n');
    fprintf(fid, '    <div class="gstat"><div class="gk">p-value</div><div class="gv%s">%s</div></div>\n', pcls, pval);
    fprintf(fid, '    <div class="gnote">%s</div>\n', note);
    fprintf(fid, '  </div>\n</section>\n');
end

% =========================================================================
% FOOTER
% =========================================================================
function write_footer(fid)
    fprintf(fid, '<footer class="foot"><span>Self-contained report &mdash; keep the PNG images in the same folder as this HTML.</span><span>Generated %s</span></footer>\n', esc(datestr(now)));
end

% =========================================================================
% DATA HELPERS
% =========================================================================
function v = get_pvec(rs, fieldname)
    v = [];
    if isfield(rs,'p') && isfield(rs.p, fieldname), v = rs.p.(fieldname); end
end

function p = get_p(pvec, ch)
    p = NaN;
    if ~isempty(pvec) && ch >= 1 && ch <= numel(pvec), p = pvec(ch); end
end

function v = stat_val(rs, ch)
    v = NaN;
    if isfield(rs,'T') && isfield(rs.T,'real_T')
        t = rs.T.real_T;
        if ch >= 1 && ch <= numel(t), v = t(ch); end
    end
end

function lab = chan_label(rs, ch)
    lab = sprintf('Ch %d', ch);
    if isfield(rs,'chanlocs') && numel(rs.chanlocs) >= ch && isfield(rs.chanlocs,'labels')
        l = rs.chanlocs(ch).labels;
        if ~isempty(l), lab = l; end
    end
end

function pv = cluster_p(C, i)
    pv = NaN;
    if isfield(C,'p')
        v = C(i).p;
        if ~isempty(v), pv = v(1); end
    end
end

function ch = cluster_chs(C, i)
    ch = [];
    if isfield(C,'channels')
        v = C(i).channels;
        if ~isempty(v), ch = double(v(:).'); end
    end
end

% =========================================================================
% FORMATTING HELPERS
% =========================================================================
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

function c = none_cls(n)
    if n > 0, c = ''; else, c = ' none'; end
end

function t = has(s, sub)
    t = ischar(s) && ~isempty(strfind(s, sub)); %#ok<STREMP>
end

function s = esc(s)
    if ~ischar(s) && ~isstring(s), s = char(string(s)); end
    s = char(s);
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
end

% =========================================================================
% MULTI-EFFECT REPORT (two-way mixed ANOVA): group / condition / interaction
% =========================================================================
function write_legend_multi(fid)
    fprintf(fid, '<div class="legend">\n');
    fprintf(fid, '  <div class="grp"><span class="k">Report</span><span>Two-way mixed ANOVA &mdash; group &amp; condition main effects and their interaction</span></div>\n');
    fprintf(fid, '  <span class="sep"></span>\n');
    fprintf(fid, '  <div class="grp"><span class="k">Maps</span><span>each effect has its own statistic map and colour scale</span></div>\n');
    fprintf(fid, '  <span class="sep"></span>\n');
    fprintf(fid, ['  <div class="grp"><span class="k">Dots</span>' ...
        '<span class="chip"><span style="width:11px;height:11px;border-radius:50%%;background:#1b2733;display:inline-block"></span>uncorrected</span>' ...
        '<span class="chip"><span style="width:11px;height:11px;border-radius:50%%;background:#fff;border:1px solid #bbb;display:inline-block"></span>corrected</span></div>\n']);
    fprintf(fid, '</div>\n');
end

function write_descriptive(fid, rs)
    if ~isfield(rs, 'descriptive'), return; end
    d = rs.descriptive;
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Descriptive means</h2><span class="sub">Group &times; condition cell means and the interaction plot</span></div>\n');
    if isfield(d, 'cellmean_png') && ~isempty(d.cellmean_png)
        fprintf(fid, '  <div class="topo-grid">\n');
        for i = 1:numel(d.cellmean_png)
            lab = '';
            if isfield(d, 'cellmean_labels') && i <= numel(d.cellmean_labels), lab = d.cellmean_labels{i}; end
            write_fig_card(fid, d.cellmean_png{i}, esc(lab));
        end
        fprintf(fid, '  </div>\n');
    end
    if isfield(d, 'interaction_png') && ~isempty(d.interaction_png)
        nch = 0; if isfield(d, 'interaction_channels'), nch = numel(d.interaction_channels); end
        fprintf(fid, '  <div style="margin-top:20px;max-width:640px">\n');
        write_fig_card(fid, d.interaction_png, 'Group &times; condition interaction');
        if nch > 0
            fprintf(fid, '    <p class="methods" style="margin-top:8px">Mean &plusmn; SE power averaged over <b>%d</b> channel(s) &mdash; the interaction cluster, else TFCE-significant, else all channels.</p>\n', nch);
        end
        fprintf(fid, '  </div>\n');
    end
    fprintf(fid, '</section>\n');
end

function write_effect_one(fid, rs, eff, k, nperm)
    gid = sprintf('eff%d', k);
    [~, statSym] = eff_stat_info(eff);
    fprintf(fid, '<section class="section">\n');
    fprintf(fid, '  <div class="head"><h2>%s</h2><span class="sub">%s-statistic &middot; %s permutation</span></div>\n', ...
        esc(eff_title(eff)), statSym, esc(perm_word(eff)));

    write_effect_topographies(fid, eff);

    effrs = mini_rs(rs, eff);
    write_significance_grouped(fid, effrs, eff.map_base, gid, statSym, nperm, ...
        eff.uncorrsigch, eff.correctTFCEsigch, eff.SnPMsigch);

    if isfield(eff, 'posthoc') && ~isempty(eff.posthoc)
        fprintf(fid, '  <div style="margin-top:30px"><span style="font-family:var(--serif);font-size:16px;font-weight:600">Pairwise post-hoc</span></div>\n');
        for i = 1:numel(eff.posthoc)
            write_posthoc_one(fid, rs, eff.posthoc(i), i, nperm, gid);
        end
    end

    if isfield(eff, 'simple_effects') && ~isempty(eff.simple_effects)
        fprintf(fid, '  <div style="margin-top:34px"><div class="head"><h2 style="font-size:16px">Simple effects</h2><span class="sub">Condition within each group</span></div></div>\n');
        for s = 1:numel(eff.simple_effects)
            se = eff.simple_effects(s);
            sgid = sprintf('%s_se%d', gid, s);
            [~, sSym] = eff_stat_info(se);
            selab = eff_title(se);
            fprintf(fid, '  <div style="margin-top:22px"><span style="font-family:var(--serif);font-size:15px;font-weight:600">%s</span></div>\n', esc(selab));
            sers = mini_rs(rs, se);
            write_significance_grouped(fid, sers, se.map_base, sgid, sSym, nperm, ...
                se.uncorrsigch, se.correctTFCEsigch, se.SnPMsigch);
            if isfield(se, 'posthoc') && ~isempty(se.posthoc)
                for i = 1:numel(se.posthoc)
                    write_posthoc_one(fid, rs, se.posthoc(i), i, nperm, sgid);
                end
            end
        end
    end
    fprintf(fid, '</section>\n');
end

function write_effect_topographies(fid, eff)
    fprintf(fid, '  <div class="topo-grid">\n');
    if strcmpi(eff.contrast_type, 'F')
        pngs = {}; labels = {};
        if isfield(eff, 'group_mean_png'), pngs = eff.group_mean_png; end
        if isfield(eff, 'group_labels'),   labels = eff.group_labels; end
        for i = 1:numel(pngs)
            lab = sprintf('Level %d', i);
            if i <= numel(labels) && ~isempty(labels{i}), lab = labels{i}; end
            write_fig_card(fid, pngs{i}, sprintf('Mean &mdash; %s', esc(lab)));
        end
        if isfield(eff, 'fmap_png') && ~isempty(eff.fmap_png)
            write_fig_card(fid, eff.fmap_png, 'F-map');
        end
    else
        if isfield(eff, 'data1_png') && ~isempty(eff.data1_png)
            labels = {}; if isfield(eff, 'group_labels'), labels = eff.group_labels; end
            la = 'A'; lb = 'B';
            if numel(labels) >= 2, la = labels{1}; lb = labels{2}; end
            write_fig_card(fid, eff.data1_png, sprintf('Mean &mdash; %s', esc(la)));
            write_fig_card(fid, eff.data2_png, sprintf('Mean &mdash; %s', esc(lb)));
            if isfield(eff, 'diff_png') && ~isempty(eff.diff_png)
                write_fig_card(fid, eff.diff_png, 't-map');
            end
        elseif isfield(eff, 'diff_png') && ~isempty(eff.diff_png)
            write_fig_card(fid, eff.diff_png, sprintf('%s (t-map)', esc(eff.effect_label)));
        end
    end
    fprintf(fid, '  </div>\n');
end

% Grouped three-way significance toggle (reuses showViewG + panel_<gid>_<w>).
function write_significance_grouped(fid, rs, map_base, gid, statSym, nperm, unc, tfce, clu)
    nU = numel(unc); nT = numel(tfce); nC = numel(clu);
    fprintf(fid, '  <div class="seg" id="seg_%s">\n', gid);
    fprintf(fid, '    <button id="btn_%s_clu" class="on" onclick="showViewG(''%s'',''clu'')">Cluster <span class="ct%s">%d</span></button>\n', gid, gid, none_cls(nC), nC);
    fprintf(fid, '    <button id="btn_%s_tfce" onclick="showViewG(''%s'',''tfce'')">TFCE <span class="ct%s">%d</span></button>\n', gid, gid, none_cls(nT), nT);
    fprintf(fid, '    <button id="btn_%s_unc" onclick="showViewG(''%s'',''unc'')">Uncorrected <span class="ct%s">%d</span></button>\n', gid, gid, none_cls(nU), nU);
    fprintf(fid, '  </div>\n');

    fprintf(fid, '  <div id="panel_%s_clu" class="vpanel" data-view-name="Cluster-corrected">\n    <div class="vbody">\n', gid);
    if ~isempty(map_base), write_map_card(fid, [map_base ' Cluster.png'], 'cluster-corrected'); end
    write_cluster_table(fid, rs, clu, statSym, nperm);
    fprintf(fid, '    </div>\n  </div>\n');

    fprintf(fid, '  <div id="panel_%s_tfce" class="vpanel" data-view-name="TFCE-corrected" style="display:none">\n    <div class="vbody">\n', gid);
    if ~isempty(map_base), write_map_card(fid, [map_base ' TFCE.png'], 'TFCE-corrected'); end
    write_flat_table(fid, rs, tfce, get_pvec(rs,'correctedTFCE'), statSym, 'TFCE-corrected', nperm);
    fprintf(fid, '    </div>\n  </div>\n');

    fprintf(fid, '  <div id="panel_%s_unc" class="vpanel" data-view-name="Uncorrected" style="display:none">\n', gid);
    write_flat_table(fid, rs, unc, get_pvec(rs,'real'), statSym, 'uncorrected', nperm);
    fprintf(fid, '  </div>\n');
end

% Per-effect global (whole-head, channel-averaged) parametric test table.
function write_global_multi(fid, rs, nperm)
    if ~isfield(rs, 'effects') || isempty(rs.effects), return; end
    nch = [];
    if isfield(rs, 'chanlocs'), nch = numel(rs.chanlocs); end
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Global test</h2><span class="sub">Whole-head mixed ANOVA on the channel-averaged signal</span></div>\n');
    fprintf(fid, '  <div class="tablewrap">\n');
    if ~isempty(nch)
        fprintf(fid, '    <div class="thead-note"><span>Parametric test on the mean across all <b>%d</b> good channels</span><span>complementary to the maps</span></div>\n', nch);
    end
    fprintf(fid, '    <table class="sum"><thead><tr><th>Effect</th><th>stat</th><th>value</th><th>df</th><th>p</th></tr></thead><tbody>\n');
    for k = 1:numel(rs.effects)
        eff = rs.effects(k);
        [~, sSym] = eff_stat_info(eff);
        val = '&mdash;';
        if isfield(eff,'global_stat') && ~isempty(eff.global_stat) && ~isnan(eff.global_stat)
            val = fmt_stat(eff.global_stat);
        end
        dfs = '&mdash;';
        if isfield(eff,'global_df1') && ~isnan(eff.global_df1) && isfield(eff,'global_df2') && ~isnan(eff.global_df2)
            if strcmpi(sSym,'F'), dfs = sprintf('%g, %g', eff.global_df1, eff.global_df2);
            else,                 dfs = sprintf('%g', eff.global_df2); end
        end
        pcls = ''; ps = '&mdash;';
        if isfield(eff,'global_pval') && ~isempty(eff.global_pval) && ~isnan(eff.global_pval)
            ps = fmt_p(eff.global_pval, nperm);
            if eff.global_pval <= 0.05, pcls = 'p-sig'; end
        end
        fprintf(fid, '      <tr class="sig"><td class="band">%s</td><td class="stat-v">%s</td><td class="stat-v">%s</td><td>%s</td><td class="%s">%s</td></tr>\n', ...
            esc(eff_title(eff)), sSym, val, dfs, pcls, ps);
    end
    fprintf(fid, '    </tbody></table>\n  </div>\n</section>\n');
end

function [name, sym] = eff_stat_info(eff)
    name = 'F';
    if isfield(eff,'contrast_type') && strcmpi(eff.contrast_type,'t'), name = 't'; end
    if strcmpi(name,'F'), sym = 'F'; else, sym = 't'; end
end

function t = eff_title(eff)
    key = ''; if isfield(eff,'key'), key = eff.key; end
    switch lower(key)
        case 'group',       t = 'Group main effect';
        case 'condition',   t = 'Condition main effect';
        case 'interaction', t = 'Interaction (group by condition)';
        otherwise
            if isfield(eff,'effect_label') && ~isempty(eff.effect_label), t = eff.effect_label; else, t = key; end
    end
end

function w = perm_word(eff)
    pt = ''; if isfield(eff,'perm_type'), pt = eff.perm_type; end
    switch lower(pt)
        case 'free',   w = 'between-subjects';
        case 'within', w = 'within-subject';
        otherwise,     w = pt;
    end
end

function m = mini_rs(rs, eff)
    m = struct();
    if isfield(eff,'T'),        m.T = eff.T;               end
    if isfield(eff,'p'),        m.p = eff.p;               end
    if isfield(eff,'Clusters'), m.Clusters = eff.Clusters; end
    if isfield(rs,'chanlocs'),  m.chanlocs = rs.chanlocs;  end
end
