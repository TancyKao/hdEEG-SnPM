function generateAnalysisReport_lmm(results_struct, params, base_filename, outputSname, uncorrsigch, correctTFCEsigch, SnPMsigch)
%GENERATEANALYSISREPORT_LMM  One self-contained HTML report for the LMM path.
% Dedicated report for core_snpm_lmm (params.comparison == 'mixedmodel'); the
% shared generateAnalysisReport is for the GLM/t-test structs and is left
% untouched. Stat-aware via results_struct.lmm.effect_type: 'continuous' ->
% signed t (diverging legend), 'factor' -> omnibus F (hot legend).
%
% Sections: masthead, model card (the fitted formula + spec), legend, effect
% map, a three-way Cluster/TFCE/Uncorrected significance toggle, a descriptive
% panel (scatter for continuous, group-means for factor), a whole-head global
% test, footer. Layout/CSS ported from the approved design mockups.
%
% Reads:
%   results_struct.lmm.{formula,dv,fixed,random,effect,effect_type,perm_scheme,
%       tail,n_trials,n_subjects,group_levels,effect_map_png,descriptive_png,
%       descriptive_channels,descriptive_kind,global_stat,global_pval}
%   results_struct.{T.real_T, p.real, p.correctedTFCE, Clusters, chanlocs}
%   base_filename -> "<base> Cluster.png" / "<base> TFCE.png" significance maps

    report_filename = [outputSname '_report.html'];
    fid = fopen(report_filename, 'w');
    if fid == -1
        warning('generateAnalysisReport_lmm:open', 'Could not open %s for writing.', report_filename);
        return;
    end
    closer = onCleanup(@() fclose(fid)); %#ok<NASGU>

    nperm = 1000;
    if isfield(params, 'permutations') && ~isempty(params.permutations)
        nperm = double(params.permutations);
    end
    lmm = struct();
    if isfield(results_struct, 'lmm'), lmm = results_struct.lmm; end
    is_t    = ~(isfield(lmm, 'effect_type') && strcmpi(lmm.effect_type, 'factor'));
    statSym = 't'; if ~is_t, statSym = 'F'; end

    write_head(fid);
    fprintf(fid, '<body>\n<div class="wrap">\n');
    write_masthead(fid, lmm, nperm, is_t);
    write_model_card(fid, lmm, nperm);
    write_legend(fid, is_t);
    write_effect_map(fid, lmm, is_t);
    write_significance(fid, results_struct, base_filename, statSym, nperm, uncorrsigch, correctTFCEsigch, SnPMsigch);
    write_descriptive(fid, lmm);
    write_global(fid, lmm, statSym, nperm);
    write_footer(fid);
    fprintf(fid, '</div>\n</body>\n</html>\n');
end

% =========================================================================
% HEAD (inline CSS + toggle JS) — ported verbatim from the design mockups
% =========================================================================
function write_head(fid)
    css = strjoin({ ...
'<!DOCTYPE html>', '<html lang="en"><head><meta charset="utf-8">', ...
'<meta name="viewport" content="width=device-width, initial-scale=1">', ...
'<title>LMM Topographic Report</title>', '<style>', ...
':root{--ink:#1b2733;--ink-2:#445261;--muted:#7c8a98;--line:#e3e9ef;--line-2:#eef2f6;--bg:#f5f7f9;--panel:#ffffff;--accent:#1f6f8b;--accent-soft:#e8f1f4;--sig:#1f8a5b;--sig-soft:#e7f4ec;--warm:#c2543a;--cool:#2f6db0;--mono:ui-monospace,"SF Mono",Menlo,Consolas,"Roboto Mono",monospace;--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;--serif:Charter,"Iowan Old Style",Georgia,"Times New Roman",serif}', ...
'*{box-sizing:border-box}html{-webkit-text-size-adjust:100%}', ...
'body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.55;-webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility}', ...
'.wrap{max-width:1100px;margin:0 auto;padding:40px 32px 96px}', ...
'.eyebrow{font-size:11px;letter-spacing:.14em;text-transform:uppercase;color:var(--muted);font-weight:600}', ...
'.masthead{display:flex;justify-content:space-between;align-items:flex-start;gap:32px;flex-wrap:wrap;padding-bottom:26px;border-bottom:1px solid var(--line)}', ...
'.masthead h1{font-family:var(--serif);font-weight:600;font-size:29px;line-height:1.2;letter-spacing:-.01em;margin:9px 0 11px;max-width:760px}', ...
'.methods{color:var(--ink-2);font-size:13.5px;max-width:820px;margin:0}.methods b{color:var(--ink);font-weight:600}', ...
'.actions{display:flex;gap:8px;flex-shrink:0}', ...
'.btn{font:inherit;font-size:12.5px;font-weight:600;color:var(--ink-2);background:var(--panel);border:1px solid var(--line);border-radius:7px;padding:8px 13px;cursor:pointer;display:inline-flex;align-items:center;gap:6px;transition:border-color .15s,color .15s}', ...
'.btn:hover{border-color:var(--accent);color:var(--accent)}.btn svg{width:14px;height:14px}', ...
'.section{margin-top:40px}', ...
'.section>.head{display:flex;align-items:baseline;gap:14px;margin-bottom:16px;flex-wrap:wrap}', ...
'.section>.head h2{font-family:var(--serif);font-size:19px;font-weight:600;letter-spacing:-.01em;margin:0}', ...
'.section>.head .sub{font-size:12.5px;color:var(--muted)}', ...
'.model-card{background:var(--panel);border:1px solid var(--line);border-radius:12px;box-shadow:0 1px 2px rgba(20,40,60,.04);overflow:hidden}', ...
'.formula{padding:22px 26px;border-bottom:1px solid var(--line-2);background:linear-gradient(180deg,#fbfdfe,#fff)}', ...
'.formula .fk{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600;display:block;margin-bottom:9px}', ...
'.formula code{font-family:var(--mono);font-size:19px;color:var(--ink);font-weight:500;letter-spacing:-.01em}', ...
'.formula code .v{color:var(--accent)}.formula code .r{color:var(--warm)}', ...
'.specs{display:grid;grid-template-columns:repeat(2,1fr);gap:0}', ...
'@media(max-width:680px){.specs{grid-template-columns:1fr}}', ...
'.specs .row{display:flex;justify-content:space-between;gap:16px;padding:11px 26px;border-bottom:1px solid var(--line-2);font-size:13.5px}', ...
'.specs .row:nth-child(odd){border-right:1px solid var(--line-2)}', ...
'@media(max-width:680px){.specs .row:nth-child(odd){border-right:0}}', ...
'.specs .k{color:var(--ink-2)}', ...
'.specs .val{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink);font-weight:500;text-align:right}', ...
'.specs .val.tag{font-family:var(--sans);font-weight:600}', ...
'.specs .val .hint{color:var(--muted);font-weight:400;font-size:12px}', ...
'.legend{display:flex;flex-wrap:wrap;gap:12px 26px;align-items:center;padding:14px 18px;background:var(--panel);border:1px solid var(--line);border-radius:10px;font-size:12.5px;color:var(--ink-2)}', ...
'.legend .grp{display:flex;align-items:center;gap:9px}.legend .grp .k{font-weight:600;color:var(--ink)}', ...
'.chip{display:inline-flex;align-items:center;gap:6px}', ...
'.tmapbar{width:130px;height:13px;border-radius:3px;border:1px solid rgba(0,0,0,.12)}', ...
'.legend .sep{width:1px;height:26px;background:var(--line)}', ...
'.legend .dd{width:10px;height:10px;border-radius:50%;flex:none}', ...
'.fig{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;box-shadow:0 1px 2px rgba(20,40,60,.04)}', ...
'.fig>.cap{display:flex;align-items:baseline;justify-content:space-between;gap:10px;padding:12px 16px;border-bottom:1px solid var(--line-2)}', ...
'.fig>.cap .t{font-size:13px;font-weight:650}.fig>.cap .f{font-family:var(--mono);font-size:10.5px;color:var(--muted)}', ...
'.slot{position:relative;margin:0;background:#f1f5f7;aspect-ratio:1/1}', ...
'.slot.wide{aspect-ratio:16/10}', ...
'.slot img{position:absolute;inset:0;width:100%;height:100%;object-fit:contain;background:#fff;display:block}', ...
'.slot img.err{display:none}.slot img.err + .ph{display:flex}', ...
'.slot .ph{position:absolute;inset:0;display:none;flex-direction:column;align-items:center;justify-content:center;gap:8px;color:#9fb0bd;text-align:center;padding:14px}', ...
'.slot .ph svg{width:30px;height:30px;opacity:.6}', ...
'.slot .ph .fn{font-family:var(--mono);font-size:10px;color:#a7b6c2;word-break:break-all;max-width:90%}', ...
'.map-note{padding:9px 16px;border-top:1px solid var(--line-2);font-size:11px;color:var(--muted);display:flex;align-items:center;gap:16px;flex-wrap:wrap}', ...
'.map-note .dl{display:inline-flex;align-items:center;gap:6px}', ...
'.map-note .dd{width:9px;height:9px;border-radius:50%;border:1px solid rgba(0,0,0,.25);flex:none}', ...
'.effect-wrap{max-width:440px}', ...
'.seg{display:inline-flex;background:var(--line-2);border-radius:9px;padding:4px;gap:3px;flex-wrap:wrap}', ...
'.seg button{font:inherit;font-size:13px;font-weight:600;color:var(--ink-2);background:transparent;border:0;border-radius:6px;padding:9px 18px;cursor:pointer;transition:.15s;display:flex;align-items:center;gap:8px}', ...
'.seg button.on{background:var(--panel);color:var(--accent);box-shadow:0 1px 2px rgba(20,40,60,.1)}', ...
'.seg button .ct{font-family:var(--mono);font-size:11px;font-weight:600;color:#fff;background:var(--muted);border-radius:20px;padding:1px 8px;min-width:22px;text-align:center}', ...
'.seg button.on .ct{background:var(--accent)}.seg button .ct.none{background:var(--line);color:var(--muted)}', ...
'.vpanel{margin-top:22px;animation:fade .25s ease}', ...
'@keyframes fade{from{opacity:0;transform:translateY(4px)}to{opacity:1;transform:none}}', ...
'.vbody{display:grid;grid-template-columns:360px 1fr;gap:22px;align-items:start}', ...
'@media(max-width:860px){.vbody{grid-template-columns:1fr}}', ...
'.tablewrap{background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;box-shadow:0 1px 2px rgba(20,40,60,.04)}', ...
'.tablewrap .thead-note{padding:11px 16px;border-bottom:1px solid var(--line-2);font-size:12px;color:var(--muted);display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap}', ...
'.tablewrap .thead-note b{color:var(--ink)}', ...
'table.sum{border-collapse:collapse;width:100%;font-size:13.5px}', ...
'table.sum th,table.sum td{padding:9px 16px;text-align:right;border-bottom:1px solid var(--line-2)}', ...
'table.sum th:first-child,table.sum td:first-child{text-align:left}', ...
'table.sum thead th{font-size:11px;letter-spacing:.03em;color:var(--ink-2);font-weight:650;white-space:nowrap;background:#fafcfd}', ...
'table.sum tbody td{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink-2)}', ...
'table.sum tbody td.band{font-family:var(--sans);color:var(--ink);font-weight:600;white-space:nowrap;display:flex;align-items:center;gap:8px}', ...
'table.sum tbody tr:last-child td{border-bottom:0}', ...
'table.sum .p-sig{color:var(--sig);font-weight:700}table.sum .stat-v{color:var(--ink);font-weight:600}', ...
'.dotkey{width:9px;height:9px;border-radius:50%;flex:none}', ...
'.dotkey.unc{background:#1b2733}.dotkey.cor{background:#fff;border:1px solid var(--ink-2)}', ...
'tr.grp-row td{background:var(--accent-soft);color:var(--accent);font-family:var(--sans);font-weight:650;font-size:12px;text-align:left;letter-spacing:.01em}', ...
'tr.grp-row td .cp{font-family:var(--mono);font-weight:600;color:var(--ink-2);margin-left:8px}', ...
'.empty-note{padding:30px 18px;text-align:center;color:var(--muted);font-size:13.5px}', ...
'.desc-body{display:grid;grid-template-columns:1.2fr 1fr;gap:22px;align-items:start}', ...
'@media(max-width:860px){.desc-body{grid-template-columns:1fr}}', ...
'.desc-note{font-size:13px;color:var(--ink-2);line-height:1.6}', ...
'.desc-note .lead{font-family:var(--serif);font-size:15px;color:var(--ink);font-weight:600;display:block;margin-bottom:8px}', ...
'.desc-note ul{margin:14px 0 0;padding-left:18px}.desc-note li{margin:4px 0}', ...
'.desc-note code{font-family:var(--mono);color:var(--accent)}', ...
'.global-card{display:flex;align-items:stretch;background:var(--panel);border:1px solid var(--line);border-radius:12px;overflow:hidden;flex-wrap:wrap;box-shadow:0 1px 2px rgba(20,40,60,.04)}', ...
'.gstat{padding:22px 30px;display:flex;flex-direction:column;gap:7px;min-width:150px}', ...
'.gstat .gk{font-size:11px;letter-spacing:.06em;text-transform:uppercase;color:var(--muted);font-weight:600}', ...
'.gstat .gv{font-family:var(--mono);font-size:30px;font-weight:600;color:var(--ink);font-variant-numeric:tabular-nums;line-height:1}', ...
'.gstat .gv.p-sig{color:var(--sig)}.gsep{width:1px;background:var(--line)}', ...
'.gnote{margin-left:auto;align-self:center;padding:20px 30px;font-size:12.5px;color:var(--ink-2);max-width:360px;line-height:1.5}.gnote b{color:var(--ink)}', ...
'.foot{margin-top:56px;padding-top:20px;border-top:1px solid var(--line);font-size:11.5px;color:var(--muted);display:flex;justify-content:space-between;gap:16px;flex-wrap:wrap}', ...
'@media(max-width:820px){.wrap{padding:28px 18px 72px}.masthead h1{font-size:24px}}', ...
'@media print{body{background:#fff}.wrap{max-width:none;padding:0}.actions,.seg{display:none !important}.vpanel{display:block !important;margin-top:18px;animation:none;page-break-inside:avoid}.vpanel::before{content:attr(data-view-name);display:block;font-weight:700;font-size:12px;color:var(--accent);margin-bottom:8px;letter-spacing:.02em}.fig,.tablewrap,.model-card,.global-card,.legend{box-shadow:none;break-inside:avoid}.section{break-inside:avoid}}', ...
'</style>', ...
'<script>function showView(w){var v=["clu","tfce","unc"];for(var i=0;i<v.length;i++){var k=v[i];var p=document.getElementById("panel_"+k);if(p)p.style.display=(k===w)?"block":"none";var b=document.getElementById("btn_"+k);if(b)b.className=(k===w)?"on":"";}}</script>', ...
'</head>'}, newline);
    fprintf(fid, '%s\n', css);
end

% =========================================================================
% MASTHEAD
% =========================================================================
function write_masthead(fid, lmm, nperm, is_t)
    eff = lmm_field(lmm, 'effect', 'effect');
    dv  = lmm_field(lmm, 'dv', 'DV');
    tail = lmm_field(lmm, 'tail', 'both');
    if strcmpi(tail, 'both'), tailw = 'two-sided'; else, tailw = 'one-tailed'; end
    if is_t, kind = 'continuous &rarr; signed t'; else, kind = 'categorical &rarr; omnibus F'; end
    title = sprintf('Effect of %s on %s &mdash; scalp topography', esc(eff), esc(dv));
    fprintf(fid, '<header class="masthead">\n  <div>\n');
    fprintf(fid, '    <span class="eyebrow">LMM Topographic Statistics</span>\n');
    fprintf(fid, '    <h1>%s</h1>\n', title);
    fprintf(fid, ['    <p class="methods">Per-channel <b>linear mixed model</b>, tested channel by channel. ' ...
        '<b>%s permutations</b>, %s, family-wise corrected with <b>TFCE</b> and <b>cluster-based</b> methods at ' ...
        '<b>&alpha; = .05</b>. &nbsp;<span style="color:var(--muted)">Effect of interest: %s (%s)</span></p>\n'], ...
        fmt_int(nperm), tailw, esc(eff), kind);
    fprintf(fid, '  </div>\n');
    fprintf(fid, ['  <div class="actions"><button class="btn" onclick="window.print()">' ...
        '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8">' ...
        '<path d="M6 9V3h12v6M6 18H4a2 2 0 0 1-2-2v-4a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v4a2 2 0 0 1-2 2h-2M6 14h12v7H6z"/>' ...
        '</svg>Print / PDF</button></div>\n']);
    fprintf(fid, '</header>\n');
end

% =========================================================================
% MODEL CARD
% =========================================================================
function write_model_card(fid, lmm, nperm)
    formula = lmm_field(lmm, 'formula', '');
    dv      = lmm_field(lmm, 'dv', '');
    eff     = lmm_field(lmm, 'effect', '');
    is_t    = ~(isfield(lmm, 'effect_type') && strcmpi(lmm.effect_type, 'factor'));
    if is_t, hint = 'continuous &rarr; signed t'; else, hint = 'categorical &rarr; omnibus F'; end
    tail = lmm_field(lmm, 'tail', 'both');
    if strcmpi(tail, 'both'), tailw = 'two-sided'; else, tailw = tail; end
    perm = lmm_field(lmm, 'perm_scheme', '');
    rand_eff = lmm_field(lmm, 'random', '');
    grp_levels = '';
    if isfield(lmm, 'group_levels') && ~isempty(lmm.group_levels)
        gl = lmm.group_levels;
        if iscell(gl), grp_levels = strjoin(cellfun(@(x) char(string(x)), gl, 'UniformOutput', false), ', '); end
    end
    ntr = lmm_num(lmm, 'n_trials'); nsub = lmm_num(lmm, 'n_subjects');

    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Model</h2><span class="sub">The linear mixed model fitted at every channel</span></div>\n');
    fprintf(fid, '  <div class="model-card">\n');
    fprintf(fid, '    <div class="formula"><span class="fk">Fitted per channel</span>\n');
    fprintf(fid, '      <code>%s</code>\n', formula_html(formula, dv, eff));
    fprintf(fid, '    </div>\n    <div class="specs">\n');
    spec_row(fid, 'Dependent variable', esc(dv), true);
    spec_row(fid, 'Effect of interest', sprintf('%s <span class="hint">%s</span>', esc(eff), hint), true);
    spec_row(fid, 'Permutation scheme', esc(perm), true);
    spec_row(fid, 'Tail', esc(tailw), true);
    spec_row(fid, 'Trials &middot; Subjects', sprintf('%s &middot; %s', num_or_dash(ntr), num_or_dash(nsub)), false);
    if ~isempty(grp_levels), spec_row(fid, 'Groups', esc(grp_levels), true); end
    if ~isempty(rand_eff),   spec_row(fid, 'Random effects', esc(rand_eff), false); end
    spec_row(fid, 'Permutations', fmt_int(nperm), false);
    fprintf(fid, '    </div>\n  </div>\n</section>\n');
end

function spec_row(fid, k, val, isTag)
    if isTag, cls = 'val tag'; else, cls = 'val'; end
    fprintf(fid, '      <div class="row"><span class="k">%s</span><span class="%s">%s</span></div>\n', k, cls, val);
end

function s = formula_html(formula, dv, eff)
    s = esc(formula);
    parts = strsplit(s, ' ~ ');
    if numel(parts) >= 2 && ~isempty(dv)
        lhs = ['<span class="v">' parts{1} '</span>'];
        rhs = strjoin(parts(2:end), ' ~ ');
        if ~isempty(eff), rhs = strrep(rhs, esc(eff), ['<span class="r">' esc(eff) '</span>']); end
        s = [lhs ' ~ ' rhs];
    end
end

% =========================================================================
% LEGEND
% =========================================================================
function write_legend(fid, is_t)
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Legend</h2><span class="sub">Statistic, colour scale &amp; significant-channel markers</span></div>\n');
    fprintf(fid, '  <div class="legend">\n');
    if is_t
        bar = 'linear-gradient(90deg,#2f6db0,#eef2f6 50%,#c2543a)';
        fprintf(fid, '    <div class="grp"><span class="k">LMM effect (t)</span><span class="chip"><span class="tmapbar" style="background:%s"></span></span><span style="color:var(--cool)">&minus; negative</span><span style="color:var(--muted)">&middot;</span><span style="color:var(--warm)">positive +</span></div>\n', bar);
    else
        bar = 'linear-gradient(90deg,#f3f6f8,#f4b63a 45%,#c2543a 78%,#7a1f10)';
        fprintf(fid, '    <div class="grp"><span class="k">omnibus F</span><span class="chip"><span class="tmapbar" style="background:%s"></span></span><span style="color:var(--muted)">low</span><span style="color:var(--muted)">&middot;</span><span style="color:var(--warm)">high F</span></div>\n', bar);
    end
    fprintf(fid, '    <span class="sep"></span>\n');
    fprintf(fid, '    <div class="grp"><span class="k">Significant channel</span><span class="chip"><span class="dd" style="background:#1b2733"></span>uncorrected</span><span class="chip"><span class="dd" style="background:#fff;border:1px solid var(--ink-2)"></span>corrected</span></div>\n');
    fprintf(fid, '  </div>\n</section>\n');
end

% =========================================================================
% EFFECT MAP
% =========================================================================
function write_effect_map(fid, lmm, is_t)
    png = lmm_field(lmm, 'effect_map_png', '');
    if is_t
        sub = 'Signed t of the effect across the scalp'; cap = 'Effect map (t)';
        note = ['<span class="dl"><span class="dd" style="background:#2f6db0"></span>t &lt; 0</span>' ...
                '<span class="dl"><span class="dd" style="background:#c2543a"></span>t &gt; 0</span>'];
    else
        sub = 'Omnibus F of the effect across the scalp'; cap = 'Effect map (F)';
        note = ['<span class="dl"><span class="dd" style="background:#f4b63a"></span>low F</span>' ...
                '<span class="dl"><span class="dd" style="background:#7a1f10"></span>high F</span>'];
    end
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Effect map</h2><span class="sub">%s</span></div>\n', sub);
    fprintf(fid, '  <div class="effect-wrap">\n    <div class="fig">\n');
    fprintf(fid, '      <div class="cap"><span class="t">%s</span><span class="f">%s</span></div>\n', cap, esc(png));
    write_slot(fid, png, false);
    fprintf(fid, ['      <div class="map-note">%s<span class="dl"><span class="dd" style="background:#1b2733"></span>uncorrected</span>' ...
        '<span class="dl"><span class="dd" style="background:#fff"></span>corrected</span></div>\n'], note);
    fprintf(fid, '    </div>\n  </div>\n</section>\n');
end

function write_slot(fid, fname, wide)
    if wide, cls = 'slot wide'; else, cls = 'slot'; end
    fprintf(fid, ['      <figure class="%s"><img src="%s" alt="%s" onerror="this.classList.add(''err'')">' ...
        '<figcaption class="ph"><svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">' ...
        '<rect x="3" y="3" width="18" height="18" rx="2"/><circle cx="8.5" cy="8.5" r="1.5"/><path d="m21 15-5-5L5 21"/></svg>' ...
        '<span class="fn">%s</span></figcaption></figure>\n'], cls, esc(fname), esc(fname), esc(fname));
end

% =========================================================================
% SIGNIFICANCE — three-way toggle
% =========================================================================
function write_significance(fid, rs, base, statSym, nperm, unc, tfce, clu)
    nU = numel(unc); nT = numel(tfce); nC = numel(clu);
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Significant channels</h2><span class="sub">Choose the multiple-comparison correction</span></div>\n');
    fprintf(fid, '  <div class="seg" id="viewSeg">\n');
    fprintf(fid, '    <button id="btn_clu" class="on" onclick="showView(''clu'')">Cluster <span class="ct%s">%d</span></button>\n', none_cls(nC), nC);
    fprintf(fid, '    <button id="btn_tfce" onclick="showView(''tfce'')">TFCE <span class="ct%s">%d</span></button>\n', none_cls(nT), nT);
    fprintf(fid, '    <button id="btn_unc" onclick="showView(''unc'')">Uncorrected <span class="ct%s">%d</span></button>\n', none_cls(nU), nU);
    fprintf(fid, '  </div>\n');

    % Cluster
    fprintf(fid, '  <div id="panel_clu" class="vpanel" data-view-name="Cluster-corrected">\n    <div class="vbody">\n');
    write_map_card(fid, [base ' Cluster.png'], 'Cluster', 'cluster-corrected');
    write_cluster_table(fid, rs, clu, statSym, nperm);
    fprintf(fid, '    </div>\n  </div>\n');
    % TFCE
    fprintf(fid, '  <div id="panel_tfce" class="vpanel" data-view-name="TFCE-corrected" style="display:none">\n    <div class="vbody">\n');
    write_map_card(fid, [base ' TFCE.png'], 'TFCE', 'TFCE-corrected');
    write_flat_table(fid, rs, tfce, get_pvec(rs,'correctedTFCE'), statSym, 'TFCE-corrected', nperm, 'cor');
    fprintf(fid, '    </div>\n  </div>\n');
    % Uncorrected (no map)
    fprintf(fid, '  <div id="panel_unc" class="vpanel" data-view-name="Uncorrected" style="display:none">\n');
    write_flat_table(fid, rs, unc, get_pvec(rs,'real'), statSym, 'uncorrected', nperm, 'unc');
    fprintf(fid, '  </div>\n</section>\n');
end

function write_map_card(fid, fname, title, kind)
    fprintf(fid, '      <div class="fig"><div class="cap"><span class="t">%s significance map</span><span class="f">%s</span></div>\n', title, esc(fname));
    write_slot(fid, fname, false);
    fprintf(fid, ['      <div class="map-note"><span class="dl"><span class="dd" style="background:#1b2733"></span>uncorrected</span>' ...
        '<span class="dl"><span class="dd" style="background:#fff"></span>%s</span></div>\n'], kind);
    fprintf(fid, '      </div>\n');
end

function write_flat_table(fid, rs, chs, pvec, statSym, label, nperm, dotcls)
    chs = chs(:).';
    if isempty(chs)
        fprintf(fid, '      <div class="tablewrap"><div class="empty-note">No %s significant channels at &alpha; = .05.</div></div>\n', label);
        return;
    end
    ps = arrayfun(@(c) get_p(pvec, c), chs);
    [ps, ord] = sort(ps); chs = chs(ord);
    fprintf(fid, '      <div class="tablewrap">\n        <div class="thead-note"><span><b>%d</b> %s significant channels</span><span>sorted by ascending p</span></div>\n', numel(chs), label);
    fprintf(fid, '        <table class="sum"><thead><tr><th>Channel</th><th>%s</th><th>p</th></tr></thead><tbody>\n', statSym);
    for k = 1:numel(chs)
        c = chs(k); pcls = ''; if ps(k) <= 0.05, pcls = 'p-sig'; end
        fprintf(fid, '          <tr><td class="band"><span class="dotkey %s"></span>%s</td><td class="stat-v">%s</td><td class="%s">%s</td></tr>\n', ...
            dotcls, esc(chan_label(rs, c)), fmt_stat(stat_val(rs, c)), pcls, fmt_p(ps(k), nperm));
    end
    fprintf(fid, '        </tbody></table>\n      </div>\n');
end

function write_cluster_table(fid, rs, clu, statSym, nperm)
    have = isfield(rs, 'Clusters') && ~isempty(rs.Clusters);
    idx = []; cps = [];
    if have
        C = rs.Clusters; nCl = numel(C); cps = nan(1, nCl); keep = false(1, nCl);
        for i = 1:nCl
            pv = cluster_p(C, i); ch = cluster_chs(C, i);
            if ~isnan(pv) && pv <= 0.05 && ~isempty(ch), cps(i) = pv; keep(i) = true; end
        end
        idx = find(keep); [~, ord] = sort(cps(idx)); idx = idx(ord);
    end
    if isempty(idx)
        if ~isempty(clu)
            write_flat_table(fid, rs, clu, get_pvec(rs,'real'), statSym, 'cluster-corrected', nperm, 'cor');
        else
            fprintf(fid, '      <div class="tablewrap"><div class="empty-note">No cluster-corrected significant channels at &alpha; = .05.</div></div>\n');
        end
        return;
    end
    total = 0; for j = 1:numel(idx), total = total + numel(cluster_chs(rs.Clusters, idx(j))); end
    fprintf(fid, '      <div class="tablewrap">\n        <div class="thead-note"><span><b>%d</b> channels in <b>%d</b> significant cluster(s)</span><span>grouped by cluster, ascending cluster p</span></div>\n', total, numel(idx));
    fprintf(fid, '        <table class="sum"><thead><tr><th>Channel</th><th>%s</th><th>cluster p</th></tr></thead><tbody>\n', statSym);
    for j = 1:numel(idx)
        ci = idx(j); ch = cluster_chs(rs.Clusters, ci); cp = cluster_p(rs.Clusters, ci);
        sv = arrayfun(@(c) abs(stat_val(rs, c)), ch); [~, o] = sort(sv, 'descend'); ch = ch(o);
        fprintf(fid, '          <tr class="grp-row"><td colspan="3">Cluster %d<span class="cp">p = %s &middot; %d ch</span></td></tr>\n', j, fmt_p(cp, nperm), numel(ch));
        for k = 1:numel(ch)
            c = ch(k);
            fprintf(fid, '          <tr><td class="band"><span class="dotkey cor"></span>%s</td><td class="stat-v">%s</td><td class="p-sig">%s</td></tr>\n', ...
                esc(chan_label(rs, c)), fmt_stat(stat_val(rs, c)), fmt_p(cp, nperm));
        end
    end
    fprintf(fid, '        </tbody></table>\n      </div>\n');
end

% =========================================================================
% DESCRIPTIVE
% =========================================================================
function write_descriptive(fid, lmm)
    png = lmm_field(lmm, 'descriptive_png', '');
    if isempty(png), return; end
    is_t = ~(isfield(lmm, 'effect_type') && strcmpi(lmm.effect_type, 'factor'));
    dv  = lmm_field(lmm, 'dv', 'DV');
    eff = lmm_field(lmm, 'effect', 'effect');
    nch = lmm_num(lmm, 'descriptive_channels');
    if is_t
        cap  = sprintf('%s vs cluster-averaged %s', esc(dv), esc(eff));
        lead = 'Continuous effect &rarr; scatter';
        body = sprintf(['Each point is a trial. <b>%s</b> plotted against %s averaged over the significant ' ...
            'channels, with a fitted mixed-model regression line.'], esc(dv), esc(eff));
        b3 = 'Line is the fixed-effect fit; subject intercepts vary';
        b2 = sprintf('Positive slope &rarr; higher %s predicts higher %s', esc(eff), esc(dv));
    else
        cap  = sprintf('%s by %s &mdash; mean &plusmn; SE', esc(dv), esc(eff));
        lead = 'Categorical effect &rarr; DV by group';
        body = sprintf('<b>%s</b> averaged over the significant channels, summarised per <b>%s</b> as mean &plusmn; standard error.', esc(dv), esc(eff));
        b2 = 'One bar per level of the factor';
        b3 = 'Group means differ; the omnibus F tests any difference';
    end
    if ~isnan(nch), b1 = sprintf('Averaged over <code>%d</code> significant channels', nch);
    else,          b1 = 'Averaged over the significant channels'; end

    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Descriptive</h2><span class="sub">Pattern behind the effect</span></div>\n');
    fprintf(fid, '  <div class="model-card">\n');
    fprintf(fid, '    <div class="cap" style="display:flex;align-items:baseline;justify-content:space-between;gap:10px;padding:12px 20px;border-bottom:1px solid var(--line-2)"><span class="t" style="font-size:13px;font-weight:650">%s</span><span class="f" style="font-family:var(--mono);font-size:10.5px;color:var(--muted)">%s</span></div>\n', cap, esc(png));
    fprintf(fid, '    <div class="desc-body" style="padding:20px">\n');
    fprintf(fid, '      <div class="fig" style="box-shadow:none;border-color:var(--line-2)">\n');
    write_slot(fid, png, true);
    fprintf(fid, '      </div>\n');
    fprintf(fid, '      <div class="desc-note"><span class="lead">%s</span>%s<ul><li>%s</li><li>%s</li><li>%s</li></ul></div>\n', lead, body, b1, b2, b3);
    fprintf(fid, '    </div>\n  </div>\n</section>\n');
end

% =========================================================================
% GLOBAL TEST
% =========================================================================
function write_global(fid, lmm, statSym, nperm)
    if ~isfield(lmm, 'global_stat') || isempty(lmm.global_stat) || isnan(lmm.global_stat), return; end
    sval = sprintf('%s = %s', statSym, fmt_statnum(lmm.global_stat));
    pv = NaN; if isfield(lmm, 'global_pval') && ~isempty(lmm.global_pval), pv = lmm.global_pval; end
    if ~isnan(pv), pval = fmt_p(pv, nperm); else, pval = '&mdash;'; end
    if ~isnan(pv) && pv <= 0.05, pcls = ' p-sig'; else, pcls = ''; end
    eff = lmm_field(lmm, 'effect', 'effect');
    fprintf(fid, '<section class="section">\n  <div class="head"><h2>Global test</h2><span class="sub">Whole-head model on the channel-averaged signal</span></div>\n');
    fprintf(fid, '  <div class="global-card">\n');
    fprintf(fid, '    <div class="gstat"><span class="gk">Statistic</span><span class="gv">%s</span></div>\n', sval);
    fprintf(fid, '    <div class="gsep"></div>\n');
    fprintf(fid, '    <div class="gstat"><span class="gk">p-value</span><span class="gv%s">%s</span></div>\n', pcls, pval);
    fprintf(fid, '    <div class="gnote">A single whole-head LMM fitted on the <b>channel-averaged</b> signal, testing the %s effect across the entire montage.</div>\n', esc(eff));
    fprintf(fid, '  </div>\n</section>\n');
end

% =========================================================================
% FOOTER
% =========================================================================
function write_footer(fid)
    fprintf(fid, '<footer class="foot"><span>Self-contained report &mdash; keep the PNG images in the same folder as this HTML.</span><span>Generated %s</span></footer>\n', esc(datestr(now)));
end

% =========================================================================
% DATA + FORMAT HELPERS
% =========================================================================
function v = lmm_field(lmm, f, dflt)
    if isfield(lmm, f) && ~isempty(lmm.(f)), v = lmm.(f); else, v = dflt; end
    if ~ischar(v) && ~isstring(v), v = char(string(v)); else, v = char(v); end
end

function v = lmm_num(lmm, f)
    v = NaN; if isfield(lmm, f) && ~isempty(lmm.(f)) && isnumeric(lmm.(f)), v = double(lmm.(f)); end
end

function s = num_or_dash(v)
    if isnan(v), s = '&mdash;'; else, s = sprintf('%d', round(v)); end
end

function v = get_pvec(rs, fieldname)
    v = []; if isfield(rs,'p') && isfield(rs.p, fieldname), v = rs.p.(fieldname); end
end
function p = get_p(pvec, ch)
    p = NaN; if ~isempty(pvec) && ch >= 1 && ch <= numel(pvec), p = pvec(ch); end
end
function v = stat_val(rs, ch)
    v = NaN;
    if isfield(rs,'T') && isfield(rs.T,'real_T')
        t = rs.T.real_T; if ch >= 1 && ch <= numel(t), v = t(ch); end
    end
end
function lab = chan_label(rs, ch)
    lab = sprintf('Ch %d', ch);
    if isfield(rs,'chanlocs') && numel(rs.chanlocs) >= ch && isfield(rs.chanlocs,'labels')
        l = rs.chanlocs(ch).labels; if ~isempty(l), lab = l; end
    end
end
function pv = cluster_p(C, i)
    pv = NaN; if isfield(C,'p'), v = C(i).p; if ~isempty(v), pv = v(1); end, end
end
function ch = cluster_chs(C, i)
    ch = []; if isfield(C,'channels'), v = C(i).channels; if ~isempty(v), ch = double(v(:).'); end, end
end
function s = fmt_p(p, nperm)
    if isempty(p) || isnan(p), s = '&mdash;'; return; end
    if p <= 0
        if nperm >= 1, s = ['&lt; ' fmt_pnum(1/nperm)]; else, s = '&lt;.001'; end
        return;
    end
    s = fmt_pnum(p);
end
function s = fmt_pnum(p)
    if p >= 1, s = '1.00'; return; end
    if p < 0.001, s = '&lt;.001'; return; end
    s = sprintf('%.3f', p); if strncmp(s, '0', 1), s = s(2:end); end
end
function s = fmt_stat(v)
    if isempty(v) || isnan(v), s = '&mdash;'; return; end
    if v < 0, sgn = '&minus;'; else, sgn = ''; end
    s = sprintf('%s%.2f', sgn, abs(v));
end
function s = fmt_statnum(v)
    if v < 0, sgn = '&minus;'; else, sgn = ''; end
    s = sprintf('%s%.2f', sgn, abs(v));
end
function s = fmt_int(n)
    s = sprintf('%d', round(double(n)));
end
function c = none_cls(n)
    if n > 0, c = ''; else, c = ' none'; end
end
function s = esc(s)
    if ~ischar(s) && ~isstring(s), s = char(string(s)); end
    s = char(s);
    s = strrep(s, '&', '&amp;'); s = strrep(s, '<', '&lt;'); s = strrep(s, '>', '&gt;');
end
