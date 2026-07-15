function export_report(DATA, OUT, opts)
% Populate the design-assistant report (sleep_eeg_report.html) with REAL data.
% Renders every image in the exact filename scheme the HTML expects, computes
% the REPORT data object (per-cell [u,tfce,clu,p] + real significant-channel
% labels), and writes both REPORT.js and a ready-to-open filled HTML copy.
%
% Two power variants:
%   absolute   - absolute band power (paired t on log10)
%   normalised - RELATIVE power = band / sum(7 bands) per channel (the file's
%                stored "normalized" values are degenerate all-1.0; relative
%                power is the meaningful stand-in until that source is fixed)
%
% Within-subject paired contrast condition-a vs condition-b (same subjects),
% per stage x band, on the 178-ch montage (NeighborMatrix_178), TFCE + cluster
% correction. Mirrors the per-cell batch runner but with the report's filename scheme,
% a per-frequency p-value periodogram, and REPORT emission.
%
% USAGE
%   DATA='.../hdeeg_analysis_all_sub';
%   OUT ='.../AnalyzeTools/SnPM_2025/test_data';   % where sleep_eeg_report.html lives
%   export_report(DATA, OUT, struct('permutations',1000));

    if nargin < 3, opts = struct(); end
    nperm  = getdef(opts,'permutations',1000);
    alpha  = getdef(opts,'alpha',0.05);
    powers = getdef(opts,'powers',{'absolute','normalised'});
    condA  = getdef(opts,'condA','');     % display names (default to the level labels below)
    condB  = getdef(opts,'condB','');
    fmax   = getdef(opts,'fmax',40);
    fmin   = getdef(opts,'fmin',0.1);
    htmlin = getdef(opts,'html', fullfile(OUT,'sleep_eeg_report.html'));
    comparison = getdef(opts,'comparison','pairedT');   % 2-level: pairedT | onesampleT | unpairedT
    paired     = ismember(lower(comparison), {'pairedt','onesamplet'});
    folders    = getdef(opts,'folders',{});             % multi-folder model (each folder = one level)
    labels     = getdef(opts,'labels',{});
    E=0.5; H=2;
    if ~exist(OUT,'dir'), mkdir(OUT); end

    Lc=load('EEG178chanlocs'); chanlocs=Lc.EEG178chanlocs;
    Ln=load('NeighborMatrix_178'); neighbors=Ln.neighbors;
    chanlabels={chanlocs.labels};

    fprintf('Scanning dataset...\n');
    [R, bandlabels, bandrange, freqs] = scan_dataset(DATA, folders, labels);
    levels = unique({R.cond}, 'stable');
    assert(numel(levels) >= 2, 'export_report needs >=2 levels (conditions/folders); found %d.', numel(levels));
    labelA = levels{1}; labelB = levels{2};
    if numel(levels) > 2
        warning('export_report:twoLevel', 'Dashboard uses the first 2 levels (%s vs %s) of %d found.', labelA, labelB, numel(levels));
    end
    if isempty(condA), condA = labelA; end
    if isempty(condB), condB = labelB; end
    if getdef(opts,'add_nrem',true), R = add_nrem(R); end   % NREM = mean(N2,N3)
    ex = getdef(opts,'exclude_subjects',{});                % drop bad-data subjects (e.g. sub-XX 20Hz artifact)
    if ~isempty(ex)
        n0=numel(R); R = R(~ismember({R.subject}, ex));
        fprintf('Excluded subjects {%s}: removed %d records\n', strjoin(ex,','), n0-numel(R));
    end
    order={'n1','n2','n3','nrem','rem'};
    stages = getdef(opts,'stages', order(ismember(order, unique({R.stage}))));
    nB = numel(bandlabels);
    fprintf('  %d records, %d bands, %d stages\n', numel(R), nB, numel(stages));

    % Fast path: regenerate only the periodograms (no permutation stats / topos / REPORT)
    if getdef(opts,'periodogram_only',false)
        for pwi=1:numel(powers)
            for si=1:numel(stages)
                periodogram_fig(R, stages{si}, powers{pwi}, freqs, fmin, fmax, bandrange, ...
                    fullfile(OUT,sprintf('global_periodogram_%s_%s.png',stages{si},powers{pwi})), condA, condB);
            end
        end
        fprintf('periodogram_only: regenerated %d periodograms in %s\n', numel(powers)*numel(stages), OUT);
        results_struct=struct(); results_text={}; %#ok<NASGU>
        return
    end

    % REPORT accumulators: values.(pw).(stage).(band)=[u tfce clu p]; channels likewise
    VAL = struct(); CH = struct();

    for pwi=1:numel(powers)
        pw = powers{pwi};
        fprintf('== power: %s ==\n', pw);
        for si=1:numel(stages)
            stg = stages{si};
            periodogram_fig(R, stg, pw, freqs, fmin, fmax, bandrange, fullfile(OUT,sprintf('global_periodogram_%s_%s.png',stg,pw)), condA, condB);
            G = struct('band',{},'meanA',{},'meanB',{},'tmap',{},'sigU',{},'sigT',{},'clim',{});
            for bi=1:nB
                [A,B] = assemble_AB(R, stg, bi, pw, labelA, labelB, paired);
                nA=size(A,1); nBn=size(B,1); ns=nA;
                vU=0; vT=0; vC=0; vP=1; labU={}; labT={}; labC={};
                if (paired && nA>=3) || (~paired && nA>=3 && nBn>=3)
                    if strcmpi(pw,'absolute'), As=log10(A); Bs=log10(B); else, As=A; Bs=B; end
                    if paired, df=nA-1; else, df=nA+nBn-2; end
                    thr=tinv(1-alpha/2,df);
                    [T,p]=snpm_single_threshold_with_TFCE(As,Bs,neighbors,E,H,alpha,comparison,'both',nperm);
                    Cl=snpm_cluster_analysis(As,Bs,thr,neighbors,alpha,comparison,'both',nperm);
                    iU=find(p.real<=alpha); iT=find(p.correctedTFCE<=alpha);
                    iC=[Cl(find([Cl.p]<=alpha)).channels]; %#ok<FNDSB>
                    vU=numel(iU); vT=numel(iT); vC=numel(iC); vP=min([Cl.p]);
                    labU=chanlabels(iU); labT=chanlabels(iT); labC=chanlabels(iC);
                    meanA=mean(A,1,'omitnan'); meanB=mean(B,1,'omitnan');
                    allv=[meanA meanB]; allv=allv(~isnan(allv));
                    clim=[prctile(allv,2) prctile(allv,98)]; if ~(clim(1)<clim(2)), clim=[min(allv) max(allv)]; end
                    save_topo(meanA,chanlocs,clim,[],[],'parula',fullfile(OUT,sprintf('%s_%s_%s_condA_mean.png',stg,pw,bandlabels{bi})));
                    save_topo(meanB,chanlocs,clim,[],[],'parula',fullfile(OUT,sprintf('%s_%s_%s_condB_mean.png',stg,pw,bandlabels{bi})));
                    save_topo(T.real_T,chanlocs,[-3 3],iU,iT,'jet',fullfile(OUT,sprintf('%s_%s_%s_Tmap.png',stg,pw,bandlabels{bi})));
                    k=numel(G)+1; G(k)=struct('band',bandlabels{bi},'meanA',meanA,'meanB',meanB,'tmap',T.real_T,'sigU',iU,'sigT',iT,'clim',clim);
                    fprintf('  %-4s %-10s n=%d TFCE=%d clu=%d\n', stg, bandlabels{bi}, ns, vT, vC);
                end
                VAL.(pwfield(pw)).(stg).(bf(bandlabels{bi})) = [vU vT vC vP];
                if vU+vT+vC>0
                    CH.(pwfield(pw)).(stg).(bf(bandlabels{bi})) = struct('u',{labU},'tfce',{labT},'cluster',{labC});
                end
            end
            make_grid(G, chanlocs, stg, pw, fullfile(OUT,sprintf('GRID_%s_%s.png',stg,pw)), condA, condB);
        end
    end

    % subject count for display: matched (paired) vs total (unpaired)
    if paired
        nsubj = numel(intersect(unique({R(strcmpi({R.cond},labelA)).subject}), ...
                                unique({R(strcmpi({R.cond},labelB)).subject})));
    else
        nsubj = numel(unique({R.subject}));
    end
    js = build_report_js(stages, bandlabels, bandrange, powers, VAL, CH, nsubj, condA, condB, nperm);
    fid=fopen(fullfile(OUT,'REPORT.js'),'w'); fwrite(fid, js); fclose(fid);
    fprintf('Wrote %s\n', fullfile(OUT,'REPORT.js'));

    % filled HTML copy with REPORT replaced
    try
        inject_report(htmlin, fullfile(OUT,'sleep_eeg_report_filled.html'), js);
        fprintf('Wrote %s\n', fullfile(OUT,'sleep_eeg_report_filled.html'));
    catch ME
        warning('Could not inject into HTML (%s); paste REPORT.js manually.', ME.message);
    end
    fprintf('Done.\n');
end

% ============================================================ data
function [R,bandlabels,bandrange,freqs] = scan_dataset(DATA, folders, labels)
    % Folder list (each folder = one level, cond = its label) or, as a fallback,
    % a single root containing condition-* subfolders (cond = part after 'condition-').
    if nargin>=2 && ~isempty(folders)
        dirs = cellstr(folders);
        if nargin>=3 && ~isempty(labels), labs = cellstr(labels);
        else, labs = cellfun(@dir_label, dirs, 'UniformOutput', false); end
    else
        d=dir(fullfile(DATA,'condition-*')); subn={d([d.isdir]).name};
        dirs = fullfile(DATA, subn);
        labs = regexprep(subn, '^condition-', '');
    end
    R=struct('subject',{},'cond',{},'stage',{},'abs',{},'rel',{},'gmPSD',{});
    bandlabels={}; bandrange={}; freqs=[];
    for ci=1:numel(dirs)
        ff=dir(fullfile(dirs{ci},'*powerspect*.mat'));
        for fi=1:numel(ff)
            m=parse_name(ff(fi).name); if isempty(m), continue; end
            S=load(fullfile(dirs{ci},ff(fi).name));
            if ~isfield(S,'EEG'), continue; end
            BND = bandsource(S.EEG);             % EEG.bands or EEG.features (same schema)
            if isempty(BND), continue; end
            EEG=S.EEG; lab={BND.label}; isabs=is_abs_type({BND.type}); ulab=unique(lab,'stable');
            if isempty(bandlabels)
                bandlabels=ulab;
                for b=1:numel(ulab)
                    ia=find(strcmpi(lab,ulab{b}) & isabs,1);
                    bandrange{b}=BND(ia).freqrange; %#ok<AGROW>
                end
            end
            if isempty(freqs), freqs=double(EEG.freqs(:))'; end
            nCh=EEG.nbchan; ab=nan(numel(ulab),nCh);
            for b=1:numel(ulab)
                ia=find(strcmpi(lab,ulab{b}) & isabs,1);
                if ~isempty(ia), ab(b,:)=BND(ia).data(:)'; end
            end
            rel = ab ./ sum(ab,1,'omitnan');     % relative power per channel
            k=numel(R)+1;
            R(k).subject=m.subject; R(k).cond=labs{ci}; R(k).stage=m.stage;   % cond = folder/condition label
            R(k).abs=ab; R(k).rel=rel; R(k).gmPSD=mean(double(EEG.data),1,'omitnan');
        end
    end
end

function R=add_nrem(R)
% Synthesise an 'nrem' stage = mean(N2,N3) per subject x condition (band power
% and global PSD averaged; relative power recomputed from the combined absolute).
    isN=ismember(lower({R.stage}),{'n2','n3'}); idxN=find(isN);
    if isempty(idxN), return; end
    pairs=strcat({R(idxN).subject},'|',{R(idxN).cond}); up=unique(pairs);
    for i=1:numel(up)
        sel=idxN(strcmp(pairs,up{i}));
        ab=mean(cat(3,R(sel).abs),3,'omitnan');
        gm=mean(cat(1,R(sel).gmPSD),1,'omitnan');
        k=numel(R)+1;
        R(k).subject=R(sel(1)).subject; R(k).cond=R(sel(1)).cond; R(k).stage='nrem';
        R(k).abs=ab; R(k).rel=ab./sum(ab,1,'omitnan'); R(k).gmPSD=gm;
    end
end

function [A,B]=assemble_AB(R,stg,bi,pw,labelA,labelB,paired)
% Two levels for a stage x band x power cell. paired=true matches subjects across
% the two levels (within); paired=false takes them independently (between groups).
    Rs=R(strcmpi({R.stage},stg)); isAbs=strcmpi(pw,'absolute'); nCh=size(Rs(1).abs,2);
    if paired
        sa={Rs(strcmpi({Rs.cond},labelA)).subject}; sb={Rs(strcmpi({Rs.cond},labelB)).subject};
        common=intersect(sa,sb,'stable'); n=numel(common);
        A=nan(n,nCh); B=nan(n,nCh);
        for s=1:n
            ka=find(strcmpi({Rs.subject},common{s})&strcmpi({Rs.cond},labelA),1);
            kb=find(strcmpi({Rs.subject},common{s})&strcmpi({Rs.cond},labelB),1);
            if isAbs, A(s,:)=Rs(ka).abs(bi,:); B(s,:)=Rs(kb).abs(bi,:);
            else,     A(s,:)=Rs(ka).rel(bi,:); B(s,:)=Rs(kb).rel(bi,:); end
        end
    else
        RA=Rs(strcmpi({Rs.cond},labelA)); RB=Rs(strcmpi({Rs.cond},labelB));
        A=nan(numel(RA),nCh); B=nan(numel(RB),nCh);
        for s=1:numel(RA), if isAbs, A(s,:)=RA(s).abs(bi,:); else, A(s,:)=RA(s).rel(bi,:); end, end
        for s=1:numel(RB), if isAbs, B(s,:)=RB(s).abs(bi,:); else, B(s,:)=RB(s).rel(bi,:); end, end
    end
end

function b=dir_label(p)
    p=char(p); if ~isempty(p) && (p(end)=='/'||p(end)=='\'), p(end)=[]; end
    [~,n,e]=fileparts(p); b=[n e];
end

% ============================================================ figures
function periodogram_fig(R,stg,pw,freqs,fmin,fmax,bandrange,fname,condA,condB)
    Rs=R(strcmpi({R.stage},stg));
    sa={Rs(strcmpi({Rs.cond},'a')).subject}; sb={Rs(strcmpi({Rs.cond},'b')).subject};
    common=intersect(sa,sb); ns=numel(common); nF=numel(freqs);
    A=nan(ns,nF); B=nan(ns,nF);
    for s=1:ns
        ka=find(strcmpi({Rs.subject},common{s})&strcmpi({Rs.cond},'a'),1);
        kb=find(strcmpi({Rs.subject},common{s})&strcmpi({Rs.cond},'b'),1);
        a=Rs(ka).gmPSD; b=Rs(kb).gmPSD;
        if strcmpi(pw,'normalised'), a=a/sum(a,'omitnan'); b=b/sum(b,'omitnan'); end
        A(s,:)=a; B(s,:)=b;
    end
    fm=freqs>=fmin & freqs<=fmax; f=freqs(fm);
    mA=mean(A(:,fm),1,'omitnan'); mB=mean(B(:,fm),1,'omitnan');
    isAbs=strcmpi(pw,'absolute');

    % --- paired t-test in 0.5 Hz bins -> significant bin edges (p<=0.05, uncorrected) ---
    bw=0.5; binedges=floor(fmin/bw)*bw : bw : fmax; sig=[];     % rows of [lo hi] significant bins
    for bi=1:numel(binedges)-1
        lo=binedges(bi); hi=binedges(bi+1); sel=freqs>=lo & freqs<hi;
        if ~any(sel), continue; end
        aB=mean(A(:,sel),2,'omitnan'); bB=mean(B(:,sel),2,'omitnan');   % mean PSD in bin per subject
        if isAbs, d=log10(aB)-log10(bB); else, d=aB-bB; end
        d=d(~isnan(d)&isfinite(d));
        if numel(d)>2
            [~,p]=ttest(d);
            if p<=0.05, sig=[sig; max(lo,fmin) min(hi,fmax)]; end %#ok<AGROW>
        end
    end

    % smooth the displayed curves (log-space moving average, ~1 Hz window);
    % the significance test above uses the raw 0.5 Hz bins, not the smoothed curve
    sm=@(y) 10.^movmean(log10(max(y,realmin)),9,'omitnan');
    mAp=sm(mA); mBp=sm(mB);

    fig=figure('Color','w','Position',[40 40 560 430],'Visible','off');
    semilogy(f,mAp,'-','Color',[0.18 0.43 0.69],'LineWidth',1.8); hold on;
    semilogy(f,mBp,'-','Color',[0.76 0.33 0.23],'LineWidth',1.8);
    yl=ylim;
    % gray shading behind the curves for significant 0.5 Hz bins
    for k=1:size(sig,1)
        patch([sig(k,1) sig(k,2) sig(k,2) sig(k,1)],[yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.85 0.85],'FaceAlpha',0.35,'EdgeColor','none','HandleVisibility','off');
    end
    uistack(findobj(gca,'Type','line'),'top'); ylim(yl);
    ylabel('Spectral density'); xlabel('Frequency (Hz)');
    legend({condA,condB},'Box','off','Location','northeast');
    title(sprintf('Global periodogram - %s (%s)',upper(stg),pw),'Interpreter','none');
    xlim([fmin fmax]); box off;
    try, print(fig,fname,'-dpng','-r150'); catch, end; close(fig);
end

function save_topo(vec,chanlocs,clim,blackIdx,whiteIdx,cmap,fname)
    f=figure('Color','w','Position',[50 50 380 320],'Visible','off');
    try, topo_cell(vec,chanlocs,clim,blackIdx,whiteIdx,cmap); colorbar; catch, end
    try, print(f,fname,'-dpng','-r130'); catch, end; close(f);
end
function topo_cell(vec,chanlocs,clim,blackIdx,whiteIdx,cmap)
    % map + manually-overlaid significance dots: black = uncorrected,
    % white = corrected (drawn on top). Manual overlay because topoplot's
    % emarker2 supports only one marker set.
    mr=max([chanlocs.radius]);   % plot all 178 channels (radius up to ~0.67)
    % headrad 0.65 draws a proper head outline (visible nose + ears) just inside the channel rim
    a={vec,chanlocs,'style','map','electrodes','off','whitebk','on','plotrad',mr,'headrad',0.65};
    if ~isempty(clim)&&clim(1)<clim(2), a=[a,{'maplimits',clim}]; end
    topoplot(a{:}); colormap(gca,cmap); hold on;
    [Y,X]=topo_xy(chanlocs);
    if ~isempty(blackIdx), plot(Y(blackIdx),X(blackIdx),'.','Color','k','MarkerSize',11); end
    if ~isempty(whiteIdx), plot(Y(whiteIdx),X(whiteIdx),'.','Color','w','MarkerSize',11); end
    set(gca,'XLim',[-0.62 0.62],'YLim',[-0.62 0.62]);  % margin so nose/ears aren't clipped
end
function [Y,X]=topo_xy(ch)
    % match topoplot's projection with plotrad = max channel radius (enlarged head):
    % squeeze so the outermost channel sits at rmax=0.5
    Th=pi/180*[ch.theta]; Rd=[ch.radius]; [x,y]=pol2cart(Th,Rd);
    sq=0.5/max(Rd); X=x*sq; Y=y*sq;
end
function make_grid(G,chanlocs,stg,pw,fname,condA,condB)
    nB=numel(G); if nB==0, return; end
    f=figure('Color','w','Position',[20 20 720 max(220,170*nB)],'Visible','off');
    tl=tiledlayout(nB,3,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s  %s  (%s mean | %s mean | T-map)',upper(stg),pw,condA,condB),'Interpreter','none');
    for r=1:nB
        nexttile; topo_cell(G(r).meanA,chanlocs,G(r).clim,[],[],'parula');
        if r==1, title([condA ' mean']); end; ylabel(G(r).band,'Interpreter','none');
        nexttile; topo_cell(G(r).meanB,chanlocs,G(r).clim,[],[],'parula'); if r==1, title([condB ' mean']); end
        nexttile; topo_cell(G(r).tmap,chanlocs,[-3 3],G(r).sigU,G(r).sigT,'jet'); if r==1, title('T-map'); end
    end
    try, print(f,fname,'-dpng','-r130'); catch, end; close(f);
end

% ============================================================ REPORT JS
function js = build_report_js(stages,bandlabels,bandrange,powers,VAL,CH,nsubj,condA,condB,nperm)
    % Greek symbols as JS \uXXXX escapes (built from code points -> no non-ASCII in source)
    gkcode=containers.Map({'low-delta','delta','theta','alpha','sigma','beta','gamma'}, ...
        {948,948,952,945,963,946,947});
    L=sprintf('const REPORT = {\n');
    L=[L sprintf('  title: "Spectral power in sleep: %s vs %s",\n',condA,condB)];
    L=[L sprintf('  methods: "Within-subject paired contrast of EEG power across sleep stages and %d frequency bands. Permutation testing (%d permutations, two-tailed), corrected with <b>TFCE</b> and <b>cluster-based</b> methods at <b>\\u03B1 = .05</b>.",\n',numel(bandlabels),nperm)];
    L=[L sprintf('  n_subj: %d,\n',nsubj)];
    L=[L sprintf('  conditions: { A:"%s", B:"%s" },\n',condA,condB)];
    pp=cellfun(@(p)sprintf('{key:"%s",name:"%s"}',p,powname(p)),powers,'uni',0);
    L=[L sprintf('  powers: [ %s ],\n', strjoin(pp,', '))];
    ss=cellfun(@(s)sprintf('{key:"%s",name:"%s"}',s,upper(s)),stages,'uni',0);
    L=[L sprintf('  stages: [ %s ],\n', strjoin(ss,', '))];
    bb={};
    for b=1:numel(bandlabels)
        nm=[upper(bandlabels{b}(1)) bandlabels{b}(2:end)];
        rg=bandrange{b}; code=948; if isKey(gkcode,bandlabels{b}), code=gkcode(bandlabels{b}); end
        g=sprintf('\\u%04X',code);
        bb{end+1}=sprintf('{key:"%s",name:"%s",gk:"%s",range:"%g\\u2013%g Hz"}',bandlabels{b},nm,g,rg(1),rg(2)); %#ok<AGROW>
    end
    L=[L sprintf('  bands: [ %s ],\n', strjoin(bb,', '))];
    % values
    L=[L sprintf('  values: {\n')];
    for pi=1:numel(powers)
        pw=powers{pi}; L=[L sprintf('   %s: {\n',pw)];
        for si=1:numel(stages)
            s=stages{si}; cells={};
            for b=1:numel(bandlabels)
                v=VAL.(pwfield(pw)).(s).(bf(bandlabels{b}));
                cells{end+1}=sprintf('"%s":[%d,%d,%d,%.3f]',bandlabels{b},v(1),v(2),v(3),v(4)); %#ok<AGROW>
            end
            term=','; if si==numel(stages), term=''; end
            L=[L sprintf('    %s:{ %s }%s\n',s,strjoin(cells,', '),term)];
        end
        term=','; if pi==numel(powers), term=''; end
        L=[L sprintf('   }%s\n',term)];
    end
    L=[L sprintf('  },\n')];
    % channels
    L=[L sprintf('  channels: {\n')];
    for pi=1:numel(powers)
        pw=powers{pi}; L=[L sprintf('   %s: {',pw)];
        entries={};
        if isfield(CH,pwfield(pw))
            sfs=fieldnames(CH.(pwfield(pw)));
            for si=1:numel(sfs)
                s=sfs{si}; bandsf=fieldnames(CH.(pwfield(pw)).(s)); be={};
                for bi=1:numel(bandsf)
                    c=CH.(pwfield(pw)).(s).(bandsf{bi});
                    be{end+1}=sprintf('"%s":{u:%s,tfce:%s,cluster:%s}', unbf(bandsf{bi}), jarr(c.u), jarr(c.tfce), jarr(c.cluster)); %#ok<AGROW>
                end
                entries{end+1}=sprintf('%s:{ %s }', s, strjoin(be,', ')); %#ok<AGROW>
            end
        end
        term=','; if pi==numel(powers), term=''; end
        L=[L sprintf(' %s }%s\n', strjoin(entries,', '), term)];
    end
    L=[L sprintf('  }\n};\n')];
    js=L;
end

function inject_report(htmlin, htmlout, js)
    h=fileread(htmlin);
    p1=strfind(h,'const REPORT = {'); assert(~isempty(p1),'REPORT start not found');
    pend=strfind(h,'END OF DATA'); assert(~isempty(pend),'END marker not found');
    % find the '/*' that begins the END-OF-DATA banner (last '/*' before pend)
    slashes=strfind(h,'/*'); slashes=slashes(slashes<pend(1)); pbanner=slashes(end);
    newh=[h(1:p1(1)-1) js sprintf('\n') h(pbanner:end)];
    fid=fopen(htmlout,'w'); fwrite(fid,newh); fclose(fid);
end

% ============================================================ helpers
function v=getdef(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end; end
function B = bandsource(EEG)
% Per-channel band-power struct array (label/type/freqrange/data). Prefer
% EEG.bands; fall back to EEG.features (some the example study exports use that field name).
    B=[];
    if isfield(EEG,'bands') && ~isempty(EEG.bands), B=EEG.bands;
    elseif isfield(EEG,'features') && ~isempty(EEG.features), B=EEG.features; end
end
function tf = is_abs_type(types)
% Match absolute-power entries whether typed 'absolute' or 'absolute power'.
    tf = cellfun(@(t) strncmpi(t,'absolute',8), types);
end
function m=parse_name(name)
    m=[]; s=regexp(name,'sub-([^_]+)','tokens','once'); c=regexp(name,'condition-([^_]+)','tokens','once');
    g=regexp(name,'desc-([^_]+)_powerspect','tokens','once');
    if isempty(s)||isempty(g), return; end; m.subject=s{1}; if ~isempty(c),m.cond=c{1};else,m.cond='';end; m.stage=g{1};
end
function s=pwfield(pw), s=pw; end
function s=powname(pw), if strcmpi(pw,'absolute'), s='Absolute'; else, s='Relative'; end; end
function f=bf(b), f=matlab.lang.makeValidName(b); end           % band -> struct field
function b=unbf(f) %#ok<DEFNU>  % struct field -> band label (reverse common cases)
    map=containers.Map({'lowDelta','low_delta','delta','theta','alpha','sigma','beta','gamma'}, ...
        {'low-delta','low-delta','delta','theta','alpha','sigma','beta','gamma'});
    if isKey(map,f), b=map(f); else, b=strrep(f,'_','-'); end
end
function s=jarr(c)
    if isempty(c), s='[]'; return; end
    q=cellfun(@(x)['"' x '"'],c,'uni',0); s=['[' strjoin(q,',') ']'];
end
