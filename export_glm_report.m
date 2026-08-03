function export_glm_report(DATA, OUT, opts)
%EXPORT_GLM_REPORT  Faceted spectral GLM report (ANOVA / ANCOVA / regression).
%
%   export_glm_report(DATA, OUT, opts)
%
%   Sibling of export_report.m (the 2-condition paired-t path). For each
%   Power x Stage x Band cell it runs a GLM preset, draws the topoplots, and
%   fills one of two templates by the resulting statistic:
%     * F  (>=3-group ANOVA/ANCOVA, interaction)  -> anova_report_template.html
%           omnibus F-map + k group means + post-hoc pairwise t-maps.
%     * t  (regression, 2-group ANOVA/ANCOVA)      -> t_report_template.html
%           signed t-map (+ k group means for 2-group; none for regression).
%   The preset's contrast type (from snpm_glm_design) chooses t vs F, so ANCOVA
%   reuses the ANOVA (3+ group) or t (2-group) report with no template of its own.
%
%   DATA : folder from gen_synthetic_spectral_glm.m
%            subjects.csv  (Subject + group/group2/predictor/covariate columns)
%            <stage>_<power>_<band>.csv  (Subject + E1..E256)
%            gmPSD_<stage>_<power>.csv   (Subject + f_<hz>)   [periodogram]
%   OUT  : output folder (filled HTML + all PNGs)
%   opts : preset ('anova1'|'ancova'|'regression'), group_col, predictor_col,
%          covariate_cols, powers/stages/bands, permutations (200), alpha, title, html.
%
%   See also EXPORT_REPORT, SNPM_GLM_DESIGN, SNPM_GLM_STAT, SNPM_PERM_CORRECTION.

    if nargin < 3, opts = struct(); end
    preset    = getdef(opts, 'preset', 'anova1');
    group_col = getdef(opts, 'group_col', 'group');
    pred_col  = getdef(opts, 'predictor_col', 'cognition');
    cov_cols  = getdef(opts, 'covariate_cols', {});
    powers    = getdef(opts, 'powers', {'absolute','normalised'});
    stages    = getdef(opts, 'stages', {'n2','n3'});
    bands     = getdef(opts, 'bands',  {'delta','sigma','beta'});
    nperm     = getdef(opts, 'permutations', 200);
    alpha     = getdef(opts, 'alpha', 0.05);
    E = 0.5; H = 2;
    if ~exist(OUT, 'dir'), mkdir(OUT); end

    isGroup = ismember(preset, {'anova1','ancova'});
    isReg   = strcmp(preset, 'regression');
    subjects = readtable(fullfile(DATA,'subjects.csv'), 'TextType','string');
    if isGroup, glist = categories(categorical(subjects.(group_col))); else, glist = {}; end

    VAL = struct(); CH = struct(); PH = struct();
    isF_any = false; statCmap = 'jet';

    for pwi = 1:numel(powers), pw = powers{pwi};
      for si = 1:numel(stages), stg = stages{si};
        Grid = struct('band',{},'means',{},'stat',{},'sigU',{},'sigT',{},'clim',{});
        for bi = 1:numel(bands), bd = bands{bi};
            cellfile = fullfile(DATA, sprintf('%s_%s_%s.csv', stg, pw, bd));
            if ~isfile(cellfile), continue; end
            Tc = readtable(cellfile, 'VariableNamingRule','preserve');
            chan_cols = Tc.Properties.VariableNames(2:end);
            raw = Tc{:, 2:end};
            sid = string(Tc{:,1});

            if strcmpi(pw,'absolute'), Y = log10(max(raw, eps)); else, Y = raw; end   % guard log of <=0
            chsel = snpm_setup_channels('egi', {Y}, chan_cols);
            Y = chsel.data{1}; chanlocs = chsel.chanlocs; neighbors = chsel.neighbors;
            labels = {chanlocs.labels};

            % design from the preset (decides t vs F)
            [meta, dopts, grp] = build_meta(subjects, sid, preset, group_col, pred_col, cov_cols);
            D = snpm_glm_design(preset, meta, dopts);
            isF = strcmpi(D.contrast_type, 'F'); isF_any = isF_any || isF;
            statCmap = ternary(isF, 'hot', 'jet');
            statClim = ternary(isF, [], [-3 3]);

            % omnibus statistic + correction
            flctx = snpm_glm_fl_context(Y, D);
            if any(~flctx.evaluable)
                warning('export_glm_report:channelsExcluded', ...
                    ['%s/%s/%s: %d of %d channels excluded for missing data (not evaluable ' ...
                     'in every permutation): %s'], stg, pw, bd, sum(~flctx.evaluable), ...
                    numel(flctx.evaluable), strjoin(labels(~flctx.evaluable), ', '));
            end
            [stat, pm] = snpm_glm_stat(Y, D.X, D.C);
            [T, p, Cl] = snpm_perm_correction(stat, pm, @() permstat(flctx, D.X, D.C), ...
                neighbors, E, H, alpha, nperm, D.contrast_type, flctx.evaluable);
            iU = find(p.real <= alpha); iT = find(p.correctedTFCE <= alpha);
            iC = [Cl([Cl.p] <= alpha).channels]; minp = min([Cl.p]);
            VAL.(pw).(stg).(bd) = [numel(iU) numel(iT) numel(iC) minp];
            if numel(iU)+numel(iT)+numel(iC) > 0
                CH.(pw).(stg).(bd) = struct('u',{labels(iU)},'tfce',{labels(iT)},'cluster',{labels(iC)});
            end

            % group-mean topoplots (group presets only)
            gmeans = zeros(0, numel(labels));
            if isGroup
                ng = numel(glist); gmeans = nan(ng, numel(labels));
                for gi = 1:ng, gmeans(gi,:) = mean(Y(grp==string(glist{gi}),:),1,'omitnan'); end
                av = gmeans(~isnan(gmeans)); av = av(:);
                if isempty(av), clim = [0 1];
                else, clim = [quantile(av,0.02) quantile(av,0.98)];
                    if ~(clim(1)<clim(2)), clim = [min(av) max(av)]; end
                    if ~(clim(1)<clim(2)), clim = [0 1]; end
                end
                for gi = 1:ng
                    save_topo(gmeans(gi,:), chanlocs, clim, [], [], 'parula', ...
                        fullfile(OUT, sprintf('%s_%s_%s_g_%s_mean.png', stg, pw, bd, glist{gi})));
                end
            else
                clim = [];
            end

            % main stat map
            save_topo(T.real_T, chanlocs, statClim, iU, iT, statCmap, ...
                fullfile(OUT, sprintf('%s_%s_%s_stat.png', stg, pw, bd)));

            % post-hoc pairwise t-maps (F presets only)
            phlist = struct('label',{},'token',{},'tfce',{},'clu',{},'p',{},'ch',{});
            if isF
                for k = 1:numel(D.posthoc)
                    Cph = pad_contrast(D.posthoc(k).C, size(D.X,2));
                    [ts, tp] = snpm_glm_stat(Y, D.X, Cph);
                    [Tk, pk, Ck] = snpm_perm_correction(ts, tp, @() permstat(flctx, D.X, Cph), ...
                        neighbors, E, H, alpha, nperm, 't', flctx.evaluable);
                    iUp = find(pk.real <= alpha); iTp = find(pk.correctedTFCE <= alpha);
                    iCp = [Ck([Ck.p] <= alpha).channels];
                    tok = ['ph_' matlab.lang.makeValidName(D.posthoc(k).label)];
                    save_topo(Tk.real_T, chanlocs, [-3 3], iUp, iTp, 'jet', ...
                        fullfile(OUT, sprintf('%s_%s_%s_%s.png', stg, pw, bd, tok)));
                    ch = struct('u',{labels(iUp)},'tfce',{labels(iTp)},'cluster',{labels(iCp)});
                    phlist(end+1) = struct('label',D.posthoc(k).label,'token',tok, ...
                        'tfce',numel(iTp),'clu',numel(iCp),'p',min([Ck.p]),'ch',ch); %#ok<AGROW>
                end
                PH.(pw).(stg).(bd) = phlist;
            end

            Grid(end+1) = struct('band',bd,'means',gmeans,'stat',T.real_T, ...
                'sigU',iU,'sigT',iT,'clim',clim); %#ok<AGROW>
        end
        meanLabels = ternary(isGroup, glist(:)', {});
        make_grid_glm(Grid, chanlocs, stg, pw, meanLabels, statCmap, ...
            fullfile(OUT, sprintf('GRID_%s_%s.png', stg, pw)));
        if isGroup
            periodogram_glm(DATA, stg, pw, subjects, group_col, glist, alpha, ...
                fullfile(OUT, sprintf('global_periodogram_%s_%s.png', stg, pw)));
        end
      end
    end

    % ---- analysis metadata for the REPORT ----
    info = analysis_info(preset, glist, group_col, pred_col, cov_cols, subjects, isF_any, nperm, alpha, opts);
    tdir = fullfile(fileparts(mfilename('fullpath')), 'templates');   % templates live in repo-root/templates
    htmlin = getdef(opts, 'html', fullfile(tdir, ternary(isF_any,'anova_report_template.html','t_report_template.html')));
    js = build_report_js(info, powers, stages, bands, VAL, CH, PH);
    write_utf8(fullfile(OUT,'REPORT.js'), js);
    inject_report(htmlin, fullfile(OUT,'anova_report.html'), js);
    fprintf('Wrote %s  (stat=%s, template=%s)\n', fullfile(OUT,'anova_report.html'), ...
        ternary(isF_any,'F','t'), htmlin);
end

% =====================================================================
function [meta, dopts, grp] = build_meta(subjects, sid, preset, gcol, pcol, ccols)
    [tf, loc] = ismember(sid, string(subjects.Subject));
    assert(all(tf), 'export_glm_report:subjMismatch', 'Cell subjects not all in subjects.csv');
    S = subjects(loc, :);
    cols = {}; names = {}; dopts = struct(); grp = strings(0);
    if ismember(preset, {'anova1','ancova'})
        grp = string(S.(gcol)); cols{end+1}=categorical(grp); names{end+1}=gcol; dopts.group_col=gcol;
    end
    if strcmp(preset, 'regression')
        cols{end+1}=double(S.(pcol)); names{end+1}=pcol; dopts.predictor_col=pcol;
    end
    for i = 1:numel(ccols)
        v = S.(ccols{i});
        if isnumeric(v), cols{end+1}=double(v); else, cols{end+1}=categorical(string(v)); end %#ok<AGROW>
        names{end+1}=ccols{i}; %#ok<AGROW>
    end
    if ~isempty(ccols), dopts.covariate_cols = ccols; end
    meta = table(cols{:}, 'VariableNames', names);
end

function info = analysis_info(preset, glist, gcol, pcol, ccols, subjects, isF, nperm, alpha, opts)
    info = struct();
    info.nperm = nperm; info.alpha = alpha; info.nsubj = height(subjects);
    info.glist = glist; info.has_periodogram = ~isempty(glist);
    if ~isempty(glist)
        info.gn = arrayfun(@(i) sum(string(subjects.(gcol))==string(glist{i})), 1:numel(glist));
        info.mean_panels = arrayfun(@(i) struct('label',[glist{i} ' mean'],'token',['g_' glist{i} '_mean']), ...
            1:numel(glist), 'UniformOutput', false);
        info.condA = glist{1}; info.condB = glist{end};
    else
        info.gn = []; info.mean_panels = {}; info.condA = pcol; info.condB = '';
    end
    covtxt = ''; if ~isempty(ccols), covtxt = [' + ' strjoin(ccols, ' + ')]; end
    switch preset
        case 'anova1'
            if isF
                info.stat='F'; info.effect_label=sprintf('Group (%d levels)', numel(glist));
                info.legend='F statistic';
                info.methods=meth('One-way ANOVA of EEG power across groups', nperm, alpha);
            else
                info.stat='t'; info.effect_label=sprintf('%s vs %s', glist{2}, glist{1});
                info.legend=sprintf('%s > %s', glist{2}, glist{1});
                info.methods=meth(sprintf('Two-group comparison of EEG power (%s vs %s)', glist{2}, glist{1}), nperm, alpha);
            end
            info.model_info='power ~ group';
        case 'ancova'
            if isF
                info.stat='F'; info.effect_label=sprintf('Group (%d levels), adjusted', numel(glist));
                info.legend='F statistic';
            else
                info.stat='t'; info.effect_label=sprintf('%s vs %s, adjusted', glist{2}, glist{1});
                info.legend=sprintf('%s > %s', glist{2}, glist{1});
            end
            info.model_info=['power ~ group' covtxt];
            info.methods=meth(['ANCOVA of EEG power across groups, adjusting for ' strjoin(ccols,', ')], nperm, alpha);
        case 'regression'
            info.stat='t'; info.effect_label=sprintf('slope of %s', pcol);
            info.legend=sprintf('positive association with %s', pcol);
            info.model_info=['power ~ ' pcol covtxt];
            info.methods=meth(['Linear regression of EEG power on ' pcol], nperm, alpha);
        case 'mixed2way'
            info.stat=ternary(isF,'F','t');
            info.effect_label='Two-way mixed ANOVA (group × condition)';
            info.legend='F statistic';
            info.model_info='power ~ group * condition';
            info.methods=meth('Two-way mixed ANOVA testing whether the condition effect differs between groups', nperm, alpha);
        otherwise
            info.stat=ternary(isF,'F','t'); info.effect_label=preset; info.legend=''; info.model_info='';
            info.methods=meth('GLM of EEG power', nperm, alpha);
    end
    info.title = getdef(opts, 'title', ['Spectral power: ' info.effect_label]);
end
function s = meth(lead, nperm, alpha)
    s = sprintf(['%s, per stage and band. Permutation testing (%d permutations, Freedman-Lane), ' ...
        'corrected with <b>TFCE</b> and <b>cluster-based</b> methods at <b>α = %g</b>.'], lead, nperm, alpha);
end

function out = getdef(s, f, d), if isfield(s,f) && ~isempty(s.(f)), out=s.(f); else, out=d; end, end
function out = ternary(c, a, b), if c, out=a; else, out=b; end, end

% ---- GLM Freedman-Lane helpers (context builder shared with core_snpm_glm.m) ----
function [stat, pm] = permstat(flctx, X, C)
    Yp = snpm_glm_permute(flctx); [stat, pm] = snpm_glm_stat(Yp, X, C);
end
function Cpad = pad_contrast(C, p)
    Cpad = zeros(size(C,1), p); Cpad(:, 1:size(C,2)) = C;
end

% ---- topo plotting helpers (copied from export_report.m; de-dup later) ----
function save_topo(vec,chanlocs,clim,blackIdx,whiteIdx,cmap,fname)
    f=figure('Color','w','Position',[50 50 380 320],'Visible','off');
    try, topo_cell(vec,chanlocs,clim,blackIdx,whiteIdx,cmap); colorbar; catch, end
    try, print(f,fname,'-dpng','-r130'); catch, end; close(f);
end
function topo_cell(vec,chanlocs,clim,blackIdx,whiteIdx,cmap)
    mr=max([chanlocs.radius]);
    a={vec,chanlocs,'style','map','electrodes','off','whitebk','on','plotrad',mr,'headrad',0.65};
    if ~isempty(clim)&&clim(1)<clim(2), a=[a,{'maplimits',clim}]; end
    topoplot(a{:}); colormap(gca,cmap); hold on;
    [Y,X]=topo_xy(chanlocs);
    if ~isempty(blackIdx), plot(Y(blackIdx),X(blackIdx),'.','Color','k','MarkerSize',11); end
    if ~isempty(whiteIdx), plot(Y(whiteIdx),X(whiteIdx),'.','Color','w','MarkerSize',11); end
    set(gca,'XLim',[-0.62 0.62],'YLim',[-0.62 0.62]);
end
function [Y,X]=topo_xy(ch)
    Th=pi/180*[ch.theta]; Rd=[ch.radius]; [x,y]=pol2cart(Th,Rd);
    sq=0.5/max(Rd); X=x*sq; Y=y*sq;
end

function make_grid_glm(G,chanlocs,stg,pw,meanLabels,statCmap,fname)
    nB=numel(G); if nB==0, return; end
    ng=numel(meanLabels); ncol=ng+1;
    f=figure('Color','w','Position',[20 20 220*ncol max(220,170*nB)],'Visible','off');
    tl=tiledlayout(nB,ncol,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s  %s', upper(stg), pw),'Interpreter','none');
    for r=1:nB
        for gi=1:ng
            nexttile; try, topo_cell(G(r).means(gi,:),chanlocs,G(r).clim,[],[],'parula'); catch, end
            if r==1, title([meanLabels{gi} ' mean'],'Interpreter','none'); end
            if gi==1, ylabel(G(r).band,'Interpreter','none'); end
        end
        nexttile; try, topo_cell(G(r).stat,chanlocs,ternary(strcmp(statCmap,'jet'),[-3 3],[]),G(r).sigU,G(r).sigT,statCmap); catch, end
        if r==1, title(ternary(strcmp(statCmap,'hot'),'F-map','stat map')); end
        if ng==0, ylabel(G(r).band,'Interpreter','none'); end
    end
    try, print(f,fname,'-dpng','-r110'); catch, end; close(f);
end

function periodogram_glm(DATA, stg, pw, subjects, group_col, glist, alpha, fname)
    pf = fullfile(DATA, sprintf('gmPSD_%s_%s.csv', stg, pw));
    if ~isfile(pf), return; end
    Tp = readtable(pf, 'VariableNamingRule','preserve');
    fcols = Tp.Properties.VariableNames(2:end);
    freqs = cellfun(@(s) str2double(erase(s,'f_')), fcols);
    P = Tp{:, 2:end};
    [~,loc] = ismember(string(Tp{:,1}), string(subjects.Subject)); grp = string(subjects.(group_col)(loc));
    ng = numel(glist); nf = numel(freqs);
    pbin = ones(1, nf);
    for j = 1:nf, try, pbin(j) = anova1(P(:,j), cellstr(grp), 'off'); catch, end, end
    f=figure('Color','w','Position',[40 40 560 460],'Visible','off'); hold on;
    cmap = lines(ng); yl = [min(P(P>0)) max(P(:))];
    for j = find(pbin <= alpha)
        patch([freqs(j)-0.25 freqs(j)+0.25 freqs(j)+0.25 freqs(j)-0.25], ...
              [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], 'EdgeColor','none');
    end
    hs = gobjects(1,ng);
    for gi = 1:ng
        m = mean(P(grp==string(glist{gi}),:),1,'omitnan');
        hs(gi) = semilogy(freqs, m, '-', 'Color', cmap(gi,:), 'LineWidth', 1.5);
    end
    set(gca,'YScale','log'); xlim([min(freqs) max(freqs)]); xlabel('Frequency (Hz)'); ylabel('Power');
    legend(hs, glist, 'Location','northeast','Interpreter','none'); box on;
    title(sprintf('Global periodogram - %s (%s)', upper(stg), pw),'Interpreter','none');
    try, print(f,fname,'-dpng','-r110'); catch, end; close(f);
end

% =====================================================================
function js = build_report_js(info, powers, stages, bands, VAL, CH, PH)
    nl = sprintf('\n'); q = @(s) ['"' s '"'];
    mp = strjoin(cellfun(@(m) sprintf('{label:"%s",token:"%s"}', m.label, m.token), info.mean_panels, 'UniformOutput', false), ', ');
    pwArr = strjoin(cellfun(@(k) sprintf('{key:%s,name:%s}', q(k), q(capit(k))), powers, 'UniformOutput', false), ', ');
    stArr = strjoin(cellfun(@(k) sprintf('{key:%s,name:%s}', q(k), q(upper(k))), stages, 'UniformOutput', false), ', ');
    bdArr = strjoin(cellfun(@(k) band_obj(k), bands, 'UniformOutput', false), ', ');

    js = ['const REPORT = {' nl];
    js = [js '  title:   ' q(info.title) ',' nl];
    js = [js '  methods: ' q(info.methods) ',' nl];
    js = [js '  n_subj:  ' num2str(info.nsubj) ',' nl];
    if ~isempty(info.glist)
        grpArr = strjoin(arrayfun(@(i) sprintf('{key:"%s",name:"%s",n:%d,desc:""}', info.glist{i}, info.glist{i}, info.gn(i)), ...
            1:numel(info.glist), 'UniformOutput', false), ', ');
        js = [js '  groups: [ ' grpArr ' ],' nl];
    end
    js = [js '  analysis: { stat:' q(info.stat) ', effect_label:' q(info.effect_label) ', ' ...
             'legend:' q(info.legend) ', model_info:' q(info.model_info) ' },' nl];
    js = [js '  mean_panels: [ ' mp ' ],' nl];
    js = [js '  stat_token: "stat",' nl];
    js = [js '  has_periodogram: ' ternary(info.has_periodogram,'true','false') ',' nl];
    js = [js '  conditions: { A:' q(info.condA) ', B:' q(info.condB) ' },' nl];
    js = [js '  powers: [ ' pwArr ' ],' nl];
    js = [js '  stages: [ ' stArr ' ],' nl];
    js = [js '  bands:  [ ' bdArr ' ],' nl];
    js = [js '  values: ' nest3(VAL, powers, stages, bands, @val_fmt) ',' nl];
    js = [js '  channels: ' nest3(CH, powers, stages, bands, @ch_fmt)];
    if strcmpi(info.stat,'F')
        js = [js ',' nl '  posthoc: ' nest3(PH, powers, stages, bands, @ph_fmt) nl];
    else
        js = [js nl];
    end
    js = [js '};'];
end
function s = nest3(A, powers, stages, bands, fmt)
    nl = sprintf('\n'); parts = {};
    for pi = 1:numel(powers), pw = powers{pi};
        sp = {};
        for si = 1:numel(stages), st = stages{si};
            sb = {};
            for bi = 1:numel(bands), bd = bands{bi};
                if isfield(A, pw) && isfield(A.(pw), st) && isfield(A.(pw).(st), bd)
                    v = fmt(A.(pw).(st).(bd));
                    if ~isempty(v), sb{end+1} = [bd ':' v]; end %#ok<AGROW>
                end
            end
            if ~isempty(sb), sp{end+1} = [st ':{' strjoin(sb,', ') '}']; end %#ok<AGROW>
        end
        parts{end+1} = ['  ' pw ': {' strjoin(sp,', ') '}']; %#ok<AGROW>
    end
    s = ['{' nl strjoin(parts, [',' nl]) nl '  }'];
end
function s = val_fmt(v), s = sprintf('[%d,%d,%d,%.3f]', v(1),v(2),v(3),v(4)); end
function s = ch_fmt(c), s = sprintf('{u:%s, tfce:%s, cluster:%s}', jarr(c.u), jarr(c.tfce), jarr(c.cluster)); end
function s = ph_fmt(ph)
    if isempty(ph), s=''; return; end
    items = cell(1, numel(ph));
    for i = 1:numel(ph)
        e = ph(i);
        chs = sprintf('{u:%s,tfce:%s,cluster:%s}', jarr(e.ch.u), jarr(e.ch.tfce), jarr(e.ch.cluster));
        items{i} = sprintf('{label:"%s", token:"%s", tfce:%d, clu:%d, p:%.3f, ch:%s}', ...
            e.label, e.token, e.tfce, e.clu, e.p, chs);
    end
    s = ['[' strjoin(items, ', ') ']'];
end
function s = jarr(c)
    if isempty(c), s='[]'; return; end
    s = ['[' strjoin(cellfun(@(x) ['"' char(x) '"'], c, 'UniformOutput', false), ',') ']'];
end
function s = band_obj(k), [nm,gk,rg]=band_meta(k); s=sprintf('{key:"%s",name:"%s",gk:"%s",range:"%s"}', k,nm,gk,rg); end
function [nm, gk, rg] = band_meta(k)
    switch k
        case 'low-delta', nm='Low-delta'; gk='δ'; rg='0.5-1 Hz';
        case 'delta',     nm='Delta';     gk='δ'; rg='1-4 Hz';
        case 'theta',     nm='Theta';     gk='θ'; rg='4-8 Hz';
        case 'alpha',     nm='Alpha';     gk='α'; rg='8-12 Hz';
        case 'sigma',     nm='Sigma';     gk='σ'; rg='12-16 Hz';
        case 'beta',      nm='Beta';      gk='β'; rg='16-30 Hz';
        case 'gamma',     nm='Gamma';     gk='γ'; rg='30-45 Hz';
        otherwise,        nm=capit(k);    gk='';  rg='';
    end
end
function s = capit(k), s = k; if ~isempty(k), s(1)=upper(s(1)); end, end

function inject_report(htmlin, htmlout, js)
    h = fileread(htmlin);
    p1 = strfind(h, 'const REPORT = {'); assert(~isempty(p1), 'template REPORT marker missing');
    pend = strfind(h, 'END OF DATA');    assert(~isempty(pend), 'template END OF DATA marker missing');
    sl = strfind(h, '/*'); sl = sl(sl < pend(1)); pbanner = sl(end);
    newh = [h(1:p1(1)-1) js sprintf('\n') h(pbanner:end)];
    write_utf8(htmlout, newh);
end
function write_utf8(fname, txt)
    fid = fopen(fname, 'w', 'n', 'UTF-8'); fwrite(fid, unicode2native(txt, 'UTF-8')); fclose(fid);
end
