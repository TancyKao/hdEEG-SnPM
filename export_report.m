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
% correction. Adds the report's filename scheme, a per-frequency p-value
% periodogram, and REPORT emission on top of the plain per-cell analysis.
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
    % Template resolves from the toolbox itself (like export_glm_report.m), NOT
    % from OUT. Copying the blank template into an output directory left a
    % 41 KB page carrying the same title and plausible DEMO numbers sitting
    % next to the real report, indistinguishable at a glance. Nothing that is
    % not a result should ever land in an output directory.
    htmlin = getdef(opts,'html', default_template());
    comparison = getdef(opts,'comparison','pairedT');   % 2-level: pairedT | onesampleT | unpairedT
    paired     = ismember(lower(comparison), {'pairedt','onesamplet'});
    folders    = getdef(opts,'folders',{});             % multi-folder model (each folder = one level)
    labels     = getdef(opts,'labels',{});
    exchan     = getdef(opts,'exclude_chan_subject',{}); % {subject, channel} pairs -> NaN mask
    seed       = getdef(opts,'seed',20260722);
    E=0.5; H=2;
    if ~exist(OUT,'dir'), mkdir(OUT); end
    if ~isempty(exchan)
        assert(iscell(exchan) && size(exchan,2)==2, ...
            'opts.exclude_chan_subject must be an N x 2 cell array {subject, channel}.');
    end

    Lc=load('EEG178chanlocs'); chanlocs=Lc.EEG178chanlocs;
    Ln=load('NeighborMatrix_178'); neighbors=Ln.neighbors;
    chanlabels={chanlocs.labels};

    fprintf('Scanning dataset...\n');
    [R, bandlabels, bandrange, freqs, maskinfo] = scan_dataset(DATA, folders, labels, exchan);
    % A mask entry that matched nothing (typo'd subject or channel label) must be
    % fatal: the methods line would otherwise claim a mask that was never applied.
    for e=1:numel(maskinfo)
        assert(maskinfo(e).nrec > 0, 'export_report:maskNoMatch', ...
            ['exclude_chan_subject entry {''%s'',''%s''} matched no record. ' ...
             'Check the subject id and the channel label against the data files.'], ...
            maskinfo(e).subject, maskinfo(e).channel);
    end
    levels = unique({R.cond}, 'stable');
    assert(numel(levels) >= 2, 'export_report needs >=2 levels (conditions/folders); found %d.', numel(levels));
    labelA = levels{1}; labelB = levels{2};
    if numel(levels) > 2
        warning('export_report:twoLevel', 'Dashboard uses the first 2 levels (%s vs %s) of %d found.', labelA, labelB, numel(levels));
    end
    if isempty(condA), condA = labelA; end
    if isempty(condB), condB = labelB; end
    if getdef(opts,'add_nrem',true), R = add_nrem(R); end   % NREM = mean(N2,N3)
    ex = getdef(opts,'exclude_subjects',{});                % drop bad-data subjects entirely
    if ~isempty(ex)
        % R.subject is the token AFTER 'sub-' (parse_name), so a literal
        % ismember against {'sub-12'} never matches. Normalise BOTH sides, and
        % fail loudly if a requested exclusion removed nothing -- a silent no-op
        % here means the analysis quietly kept a subject it was told to drop.
        ex = cellstr(ex);
        exn = cellfun(@norm_subj, ex, 'UniformOutput', false);
        rsub = cellfun(@norm_subj, {R.subject}, 'UniformOutput', false);
        for e=1:numel(exn)
            assert(any(strcmp(rsub, exn{e})), 'export_report:excludeNoMatch', ...
                'exclude_subjects entry ''%s'' matched no record.', ex{e});
        end
        keep = ~ismember(rsub, exn);
        n0=numel(R); R = R(keep);
        fprintf('Excluded subjects {%s}: removed %d of %d records\n', strjoin(ex,','), n0-numel(R), n0);
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
                    fullfile(OUT,sprintf('global_periodogram_%s_%s.png',stages{si},powers{pwi})), condA, condB, labelA, labelB);
            end
        end
        fprintf('periodogram_only: regenerated %d periodograms in %s\n', numel(powers)*numel(stages), OUT);
        results_struct=struct(); results_text={}; %#ok<NASGU>
        return
    end

    % REPORT accumulators: values.(pw).(stage).(band)=[u tfce clu p]; channels,
    % per-cluster detail (CLU) and per-cell n (NN) keyed the same way.
    VAL = struct(); CH = struct(); CLU = struct(); NN = struct();

    % One seed for the WHOLE sweep. Re-seeding per cell would hand all 70 cells
    % an identical sign-flip matrix.
    rng(seed,'twister');
    fprintf('RNG seed %d (twister), %d permutations per cell\n', seed, nperm);

    for pwi=1:numel(powers)
        pw = powers{pwi};
        fprintf('== power: %s ==\n', pw);
        for si=1:numel(stages)
            stg = stages{si};
            periodogram_fig(R, stg, pw, freqs, fmin, fmax, bandrange, fullfile(OUT,sprintf('global_periodogram_%s_%s.png',stg,pw)), condA, condB, labelA, labelB);
            G = struct('band',{},'meanA',{},'meanB',{},'tmap',{},'sigU',{},'sigT',{},'clim',{},'clu',{});
            for bi=1:nB
                [A,B] = assemble_AB(R, stg, bi, pw, labelA, labelB, paired);
                nA=size(A,1); nBn=size(B,1); ns=nA;
                vU=0; vT=0; vC=0; vP=1; labU={}; labT={}; labC={}; C=empty_clu();
                if (paired && nA>=3) || (~paired && nA>=3 && nBn>=3)
                    if strcmpi(pw,'absolute'), As=log10(A); Bs=log10(B); else, As=A; Bs=B; end
                    if paired, df=nA-1; else, df=nA+nBn-2; end
                    thr=tinv(1-alpha/2,df);
                    [T,p]=snpm_single_threshold_with_TFCE(As,Bs,neighbors,E,H,alpha,comparison,'both',nperm);
                    Cl=snpm_cluster_analysis(As,Bs,thr,neighbors,alpha,comparison,'both',nperm);
                    iU=find(p.real<=alpha); iT=find(p.correctedTFCE<=alpha);
                    C  = rank_clusters(Cl, alpha, chanlabels, condA, condB, pw, T.real_T, A, B, paired, neighbors);
                    iC = [C.idx];
                    vU=numel(iU); vT=numel(iT); vC=numel(iC); vP=min_cluster_p(Cl);
                    labU=chanlabels(iU); labT=chanlabels(iT); labC=chanlabels(iC);
                    meanA=mean(A,1,'omitnan'); meanB=mean(B,1,'omitnan');
                    allv=[meanA meanB]; allv=allv(~isnan(allv));
                    clim=[prctile(allv,2) prctile(allv,98)]; if ~(clim(1)<clim(2)), clim=[min(allv) max(allv)]; end
                    save_topo(meanA,chanlocs,clim,[],[],'parula',fullfile(OUT,sprintf('%s_%s_%s_condA_mean.png',stg,pw,bandlabels{bi})));
                    save_topo(meanB,chanlocs,clim,[],[],'parula',fullfile(OUT,sprintf('%s_%s_%s_condB_mean.png',stg,pw,bandlabels{bi})));
                    save_topo(T.real_T,chanlocs,[-3 3],iU,iT,'jet',fullfile(OUT,sprintf('%s_%s_%s_Tmap.png',stg,pw,bandlabels{bi})),C);
                    % Companion map with the CLUSTER layer instead of the TFCE
                    % layer: black = uncorrected, white = cluster-significant,
                    % plus the numbered box. Blanking .ring suppresses the
                    % outline (topo_cell guards on ~isempty(rg)) while the
                    % number is still drawn, so the white dots stay readable.
                    Cn = C; if ~isempty(Cn), [Cn.ring] = deal([]); end
                    save_topo(T.real_T,chanlocs,[-3 3],iU,iC,'jet', ...
                        fullfile(OUT,sprintf('%s_%s_%s_Tmap_cluster.png',stg,pw,bandlabels{bi})),Cn);
                    k=numel(G)+1; G(k)=struct('band',bandlabels{bi},'meanA',meanA,'meanB',meanB,'tmap',T.real_T,'sigU',iU,'sigT',iT,'clim',clim,'clu',C);
                    fprintf('  %-4s %-10s n=%d TFCE=%d clu=%d (%d cluster%s, min p=%.4f)\n', ...
                        stg, bandlabels{bi}, ns, vT, vC, numel(C), plural(numel(C)), vP);
                else
                    fprintf('  %-4s %-10s n=%d SKIPPED (too few subjects)\n', stg, bandlabels{bi}, ns);
                end
                VAL.(pwfield(pw)).(stg).(bf(bandlabels{bi})) = [vU vT vC vP];
                NN.(pwfield(pw)).(stg).(bf(bandlabels{bi}))  = ns;
                if vU+vT+vC>0
                    CH.(pwfield(pw)).(stg).(bf(bandlabels{bi})) = struct('u',{labU},'tfce',{labT},'cluster',{labC});
                end
                if ~isempty(C)
                    CLU.(pwfield(pw)).(stg).(bf(bandlabels{bi})) = C;
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
    js = build_report_js(stages, bandlabels, bandrange, powers, VAL, CH, CLU, NN, ...
                         nsubj, condA, condB, nperm, seed, maskinfo);
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
function [R,bandlabels,bandrange,freqs,maskinfo] = scan_dataset(DATA, folders, labels, exchan)
    % Folder list (each folder = one level, cond = its label) or, as a fallback,
    % a single root containing condition-* subfolders (cond = part after 'condition-').
    %
    % exchan: N x 2 cell {subject, channel-label}. Bad-channel handling is
    % PAIRWISE DELETION via a NaN mask (no interpolation, no channel/subject
    % removal). The mask is applied to ALL bands, BOTH conditions and ALL stages
    % (blanket), and to the channel-averaged global PSD:
    %   * all bands, because relative power shares a per-channel denominator, so
    %     one contaminated band corrupts every band's share;
    %   * blanket across stages, because add_nrem averages N2/N3 with 'omitnan',
    %     so a stage-specific mask would silently turn NREM into an N3-only value
    %     at that channel while every other channel is an N2/N3 mean;
    %   * gmPSD, because a channel orders of magnitude hot distorts the reported
    %     global periodograms.
    % The mask is applied to `ab` BEFORE `rel = ab./sum(ab)`; masking afterwards
    % would compute relative power on a 6-band denominator (silent scale error).
    if nargin<4, exchan={}; end
    maskinfo = struct('subject',{},'channel',{},'nrec',{});
    for e=1:size(exchan,1)
        maskinfo(e).subject = char(exchan{e,1});
        maskinfo(e).channel = char(exchan{e,2});
        maskinfo(e).nrec    = 0;
    end
    if isempty(exchan), exsub = {}; else, exsub = cellfun(@norm_subj, exchan(:,1), 'UniformOutput', false); end
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
            % ---- bad-channel NaN mask: BEFORE the relative-power denominator ----
            D = double(EEG.data);
            hits = find(strcmp(exsub, norm_subj(m.subject)));
            for hh = 1:numel(hits)
                e  = hits(hh);
                ic = chan_index(EEG, exchan{e,2}, nCh);
                if isempty(ic)
                    warning('export_report:maskChan','Channel %s not found in %s; not masked.', ...
                        exchan{e,2}, ff(fi).name);
                    continue
                end
                ab(:,ic) = NaN;                 % all bands, this subject x channel
                D(ic,:)  = NaN;                 % keep it out of the gmPSD channel average
                maskinfo(e).nrec = maskinfo(e).nrec + 1;
            end
            rel = ab ./ sum(ab,1,'omitnan');     % relative power per channel
            k=numel(R)+1;
            R(k).subject=m.subject; R(k).cond=labs{ci}; R(k).stage=m.stage;   % cond = folder/condition label
            R(k).abs=ab; R(k).rel=rel; R(k).gmPSD=mean(D,1,'omitnan');
        end
    end
    for e=1:numel(maskinfo)
        fprintf('Bad-channel mask: subject %s channel %s -> NaN in %d records (all bands, both conditions, all stages)\n', ...
            maskinfo(e).subject, maskinfo(e).channel, maskinfo(e).nrec);
    end
end

function s = norm_subj(x)
% 'sub-12' / 'sub12' / '12' -> '12'
    s = lower(strtrim(char(x)));
    s = regexprep(s,'^sub[-_]?','');
end

function ic = chan_index(EEG, label, nCh)
% Index of a channel label within this file's montage. Prefers the file's own
% chanlocs; falls back to the E<n> numbering when chanlocs are absent.
    ic = [];
    label = char(label);
    if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs) && isfield(EEG.chanlocs,'labels')
        ic = find(strcmpi({EEG.chanlocs.labels}, label));
    end
    if isempty(ic)
        t = regexp(label,'^[Ee](\d+)$','tokens','once');
        if ~isempty(t)
            n = str2double(t{1});
            if n>=1 && n<=nCh, ic = n; end
        end
    end
    ic = ic(:)';
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

% ============================================================ clusters
function C = empty_clu()
% 0x0 struct array with the per-cluster schema (so [C.idx] is [] and
% numel(C)==0 when nothing survives).
    C = struct('num',{},'p',{},'nch',{},'dir',{},'dirtxt',{}, ...
               'mdiff',{},'mdifft',{},'munit',{},'tpeak',{},'tpeakch',{}, ...
               'idx',{},'ring',{},'labels',{});
end

function b = cluster_boundary(idx, neighbors)
% Members with at least one neighbour OUTSIDE the cluster. Ringing every member
% of a 123-of-178-channel cluster fills the head with white and destroys the
% black-uncorrected / white-TFCE dot distinction the legend describes; an
% outline conveys the same extent and stays readable in the small grid tiles.
% Degenerates to the whole cluster for small or thin clusters, which is correct.
    if isempty(neighbors), b = idx; return; end
    inC = false(1, size(neighbors,1)); inC(idx) = true;
    keep = false(1, numel(idx));
    for i = 1:numel(idx)
        nb = neighbors(idx(i),:); nb = nb(~isnan(nb) & nb>0);
        keep(i) = isempty(nb) || any(~inC(nb));
    end
    b = idx(keep);
    % b is empty only when the cluster covers the whole montage: there is no
    % boundary to draw, so nothing is ringed and the number alone marks it.
end

function C = rank_clusters(Cl, alpha, chanlabels, condA, condB, pw, tmap, A, B, paired, neighbors)
% Significant clusters, one struct element each, numbered WITHIN this cell
% (power x stage x band) -- never renumbered across cells.
%   * empty placeholder elements (a threshold with no suprathreshold channel;
%     snpm_cluster_analysis emits one with .channels=[] and .p=1) are dropped
%     BEFORE numbering;
%   * ranked by ascending corrected p, then descending cluster size, then
%     ascending minimum channel index. p is a deterministic monotone step
%     function of size here, so size is the finer resolution of the same
%     evidence; the third key only guarantees a deterministic order.
%   * sign(threshold) > 0 means the positive tail, i.e. condition A > condition B.
    C = empty_clu();
    if isempty(Cl), return; end
    keep = find(arrayfun(@(c) ~isempty(c.channels) && ~isempty(c.p) && c.p<=alpha, Cl(:)'));
    if isempty(keep), return; end

    % Mean difference on TWO scales, because for absolute power the test is run
    % on log10 while the readable units are the raw ones:
    %   dch  - original units (raw power, or share of total power)
    %   dcht - the scale the t-test actually used (log10 for absolute power)
    isAbs = strcmpi(pw,'absolute');
    if isAbs, At=log10(A); Bt=log10(B); munit='power units';
    else,     At=A;        Bt=B;        munit='share of total power'; end
    if paired
        dch  = mean(A-B,1,'omitnan');
        dcht = mean(At-Bt,1,'omitnan');
    else
        dch  = mean(A,1,'omitnan')  - mean(B,1,'omitnan');
        dcht = mean(At,1,'omitnan') - mean(Bt,1,'omitnan');
    end

    pv = arrayfun(@(c) c.p, Cl(keep));
    kv = arrayfun(@(c) numel(c.channels), Cl(keep));
    mv = arrayfun(@(c) min(c.channels), Cl(keep));
    [~,ord] = sortrows([pv(:) -kv(:) mv(:)], [1 2 3]);
    keep = keep(ord);

    for i = 1:numel(keep)
        c   = Cl(keep(i));
        idx = sort(c.channels(:))';
        [~,j] = max(abs(tmap(idx)));
        C(i).num     = i;
        C(i).p       = c.p;
        C(i).nch     = numel(idx);
        C(i).dir     = sign(c.threshold);
        C(i).dirtxt  = dir_text(sign(c.threshold), condA, condB, pw);
        C(i).mdiff   = mean(dch(idx),'omitnan');
        C(i).mdifft  = mean(dcht(idx),'omitnan');
        C(i).munit   = munit;
        C(i).tpeak   = tmap(idx(j));
        C(i).tpeakch = chanlabels{idx(j)};
        C(i).idx     = idx;
        C(i).ring    = cluster_boundary(idx, neighbors);
        C(i).labels  = chanlabels(idx);
    end
end

function s = dir_text(d, condA, condB, pw)
% Relative power is a SHARE of total power: a relative increase is compatible
% with an absolute decrease, so the wording must not claim "greater power".
    if d >= 0, hi=condA; lo=condB; else, hi=condB; lo=condA; end
    if strcmpi(pw,'absolute')
        s = sprintf('%s > %s', hi, lo);
    else
        s = sprintf('%s > %s (greater share of total power)', hi, lo);
    end
end

function v = min_cluster_p(Cl)
% Smallest corrected p over REAL clusters (placeholder elements with no
% channels excluded); 1 when there is no real cluster.
    v = 1;
    if isempty(Cl), return; end
    real = find(arrayfun(@(c) ~isempty(c.channels) && ~isempty(c.p), Cl(:)'));
    if isempty(real), return; end
    pv = arrayfun(@(k) Cl(k).p, real);
    v = min(pv);
end

function s = plural(n), if n==1, s=''; else, s='s'; end, end

function b=dir_label(p)
    p=char(p); if ~isempty(p) && (p(end)=='/'||p(end)=='\'), p(end)=[]; end
    [~,n,e]=fileparts(p); b=[n e];
end

% ============================================================ figures
function periodogram_fig(R,stg,pw,freqs,fmin,fmax,bandrange,fname,condA,condB,labelA,labelB)
    % labelA/labelB are the LEVEL keys in R.cond (folder/condition labels);
    % condA/condB are the display names. They are not interchangeable.
    if nargin<11 || isempty(labelA), labelA='a'; end
    if nargin<12 || isempty(labelB), labelB='b'; end
    Rs=R(strcmpi({R.stage},stg));
    sa={Rs(strcmpi({Rs.cond},labelA)).subject}; sb={Rs(strcmpi({Rs.cond},labelB)).subject};
    common=intersect(sa,sb); ns=numel(common); nF=numel(freqs);
    A=nan(ns,nF); B=nan(ns,nF);
    for s=1:ns
        ka=find(strcmpi({Rs.subject},common{s})&strcmpi({Rs.cond},labelA),1);
        kb=find(strcmpi({Rs.subject},common{s})&strcmpi({Rs.cond},labelB),1);
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

function save_topo(vec,chanlocs,clim,blackIdx,whiteIdx,cmap,fname,clu)
    if nargin<8, clu=[]; end
    f=figure('Color','w','Position',[50 50 380 320],'Visible','off');
    try, topo_cell(vec,chanlocs,clim,blackIdx,whiteIdx,cmap,clu); colorbar; catch, end
    try, print(f,fname,'-dpng','-r130'); catch, end; close(f);
end
function topo_cell(vec,chanlocs,clim,blackIdx,whiteIdx,cmap,clu)
    % map + manually-overlaid significance dots: black = uncorrected,
    % white = corrected (drawn on top). Manual overlay because topoplot's
    % emarker2 supports only one marker set.
    % clu (optional): struct array of significant clusters -> the cluster
    % OUTLINE (.ring, boundary members) is ringed. Outline rather than every
    % member so a 123-of-178-channel cluster does not bury the significance
    % dots. Cluster NUMBERS are drawn outside the head with a leader line, and
    % ONLY when there are 2+ clusters (see label_clusters).
    if nargin<7, clu=[]; end
    mr=max([chanlocs.radius]);   % plot all 178 channels (radius up to ~0.67)
    % headrad 0.65 draws a proper head outline (visible nose + ears) just inside the channel rim
    a={vec,chanlocs,'style','map','electrodes','off','whitebk','on','plotrad',mr,'headrad',0.65};
    if ~isempty(clim)&&clim(1)<clim(2), a=[a,{'maplimits',clim}]; end
    topoplot(a{:}); colormap(gca,cmap); hold on;
    [Y,X]=topo_xy(chanlocs);

    % Label geometry first: the leader lines must go UNDER the significance
    % markers, so they are drawn before the dots.
    [lab, rlab] = label_clusters(clu, Y, X);
    for k=1:numel(lab)
        plot([lab(k).ax lab(k).tx],[lab(k).ay lab(k).ty], ...
            '-','Color',[0.45 0.45 0.45],'LineWidth',0.75);
    end

    if ~isempty(blackIdx), plot(Y(blackIdx),X(blackIdx),'.','Color','k','MarkerSize',11); end
    if ~isempty(whiteIdx), plot(Y(whiteIdx),X(whiteIdx),'.','Color','w','MarkerSize',11); end
    for c=1:numel(clu)
        ix=clu(c).idx; if isempty(ix), continue; end
        % .ring is authoritative when present (empty = whole-montage cluster,
        % nothing to outline); only a struct without the field falls back.
        if isfield(clu,'ring'), rg=clu(c).ring; else, rg=ix; end
        % outline the cluster (open marker, distinct from the filled sig dots)
        if ~isempty(rg)
            plot(Y(rg),X(rg),'o','MarkerEdgeColor','w','MarkerFaceColor','none', ...
                'MarkerSize',8,'LineWidth',1.1);
        end
    end
    % numbers last, on the white margin: gray reads there without a backing box
    for k=1:numel(lab)
        text(lab(k).ax,lab(k).ay,sprintf('%d',lab(k).num), ...
            'HorizontalAlignment',lab(k).halign,'VerticalAlignment',lab(k).valign, ...
            'FontSize',9,'FontWeight','bold','Color',[0.35 0.35 0.35]);
    end

    if isempty(lab)
        set(gca,'XLim',[-0.62 0.62],'YLim',[-0.62 0.62]);  % margin so nose/ears aren't clipped
    else
        m = rlab + 0.10;                                   % room for the label glyph itself
        set(gca,'XLim',[-m m],'YLim',[-m m]);
    end
end

function [lab, rlab] = label_clusters(clu, Y, X)
% Anchor points for the cluster numbers, placed OUTSIDE the head outline.
%
% Gray text only reads reliably on white. Under a significant cluster the map
% is saturated (t ~ +3 dark red or -3 dark blue), so a gray label inside the
% head would need a backing box -- the very thing we are removing. Outside the
% head the background is plain white.
%
% Drawn ONLY when 2+ clusters are present: with a single cluster the number is
% redundant (the white dots ARE the cluster and the adjacent table names it).
%
% Note the plot order: horizontal axis is Y, vertical axis is X.
    lab = struct('num',{},'ax',{},'ay',{},'tx',{},'ty',{},'halign',{},'valign',{});
    rlab = 0;
    if isempty(clu), return; end
    have = find(arrayfun(@(c) ~isempty(c.idx), clu(:)'));
    if numel(have) < 2, return; end

    % Label ring radius, derived EMPIRICALLY from the projected channel
    % positions. Do NOT derive it from headrad/plotrad: topoplot's internal
    % scaling relative to topo_xy's squeeze is not worth deriving analytically
    % and would break silently if either changed. 1.42 clears the ear bumps,
    % which stick out further than the head circle (1.30 clipped them).
    rlab = 1.42 * max(hypot(Y,X));

    % direction of each cluster from the head centre
    th = zeros(1,numel(have));
    for i=1:numel(have)
        ix = clu(have(i)).idx;
        th(i) = atan2(mean(X(ix),'omitnan'), mean(Y(ix),'omitnan'));
    end

    % Collision handling: two clusters at a similar bearing would put their
    % labels on top of each other, and a leader line pointing at the wrong
    % cluster is actively misleading. Sort by angle and enforce a minimum
    % separation, nudging later labels round the ring.
    [th, ord] = sort(th); have = have(ord);
    minsep = deg2rad(22);
    for i=2:numel(th)
        if th(i) - th(i-1) < minsep, th(i) = th(i-1) + minsep; end
    end
    % wrap-around: last vs first, sharing the ring
    if numel(th)>1 && (th(1) + 2*pi - th(end)) < minsep
        push = (minsep - (th(1) + 2*pi - th(end))) / 2;
        th(end) = th(end) - push; th(1) = th(1) + push;
    end

    for i=1:numel(have)
        c  = clu(have(i));
        ax = rlab*cos(th(i)); ay = rlab*sin(th(i));
        % leader line ends at the cluster member NEAREST THE LABEL, not the
        % centroid: a line to the centroid can cut across a neighbouring cluster
        ix = c.idx;
        [~,j] = min(hypot(Y(ix)-ax, X(ix)-ay)); j = ix(j);
        k = numel(lab)+1;
        lab(k).num = c.num;
        lab(k).ax = ax;  lab(k).ay = ay;
        lab(k).tx = Y(j); lab(k).ty = X(j);
        % text grows outward, away from the head, rather than back over it
        if     ax >  0.15*rlab, lab(k).halign = 'left';
        elseif ax < -0.15*rlab, lab(k).halign = 'right';
        else,                   lab(k).halign = 'center';
        end
        if     ay >  0.15*rlab, lab(k).valign = 'bottom';
        elseif ay < -0.15*rlab, lab(k).valign = 'top';
        else,                   lab(k).valign = 'middle';
        end
    end
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
        nexttile; topo_cell(G(r).tmap,chanlocs,[-3 3],G(r).sigU,G(r).sigT,'jet',G(r).clu); if r==1, title('T-map'); end
    end
    try, print(f,fname,'-dpng','-r130'); catch, end; close(f);
end

% ============================================================ REPORT JS
function js = build_report_js(stages,bandlabels,bandrange,powers,VAL,CH,CLU,NN,nsubj,condA,condB,nperm,seed,maskinfo)
    % Greek symbols as JS \uXXXX escapes (built from code points -> no non-ASCII in source)
    gkcode=containers.Map({'low-delta','delta','theta','alpha','sigma','beta','gamma'}, ...
        {948,948,952,945,963,946,947});
    L=sprintf('const REPORT = {\n');
    L=[L sprintf('  title: %s,\n', jstr(sprintf('Spectral power in sleep: %s vs %s',condA,condB)))];
    meth = sprintf(['Within-subject paired contrast of EEG power across sleep stages and %d frequency bands. ' ...
        'Permutation testing (%d permutations, two-tailed, seed %d), corrected with <b>TFCE</b> and ' ...
        '<b>cluster-based</b> methods at <b>&alpha; = .05</b>. Absolute power is tested on ' ...
        'log<sub>10</sub> power; relative power is tested on the untransformed share of total power. ' ...
        'Cluster-corrected results are reported per cluster: each numbered cluster is a separate ' ...
        'corrected test.'], numel(bandlabels), nperm, seed);
    if ~isempty(maskinfo)
        mtxt=cell(1,numel(maskinfo));
        for e=1:numel(maskinfo), mtxt{e}=sprintf('%s/%s',maskinfo(e).subject,maskinfo(e).channel); end
        meth=[meth sprintf([' Artifactual channels excluded for the affected subject only ' ...
            '(set to missing in every band, both conditions and all stages; pairwise deletion): %s.'], ...
            strjoin(mtxt,', '))];
    end
    L=[L sprintf('  methods: %s,\n', jstr(meth))];
    L=[L sprintf('  n_subj: %d,\n',nsubj)];
    L=[L sprintf('  conditions: { A:%s, B:%s },\n',jstr(condA),jstr(condB))];
    % provenance: seed, permutation count and the bad-channel mask actually applied
    L=[L sprintf('  seed: %d,\n',seed)];
    L=[L sprintf('  permutations: %d,\n',nperm)];
    mi={};
    for e=1:numel(maskinfo)
        mi{end+1}=sprintf('{subject:%s,channel:%s,records:%d}', ...
            jstr(maskinfo(e).subject), jstr(maskinfo(e).channel), maskinfo(e).nrec); %#ok<AGROW>
    end
    L=[L sprintf('  excluded_channels: [ %s ],\n', strjoin(mi,', '))];
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
    L=[L sprintf('  },\n')];

    % per-cluster detail: same nesting as `channels`. Cells with no significant
    % cluster are simply absent (as in `channels`), so an older report template
    % that ignores this block still renders.
    L=[L sprintf('  clusters: {\n')];
    for pi=1:numel(powers)
        pw=powers{pi}; L=[L sprintf('   %s: {',pw)];
        entries={};
        if isfield(CLU,pwfield(pw))
            sfs=fieldnames(CLU.(pwfield(pw)));
            for si=1:numel(sfs)
                s=sfs{si}; bandsf=fieldnames(CLU.(pwfield(pw)).(s)); be={};
                for bi=1:numel(bandsf)
                    be{end+1}=sprintf('"%s":%s', unbf(bandsf{bi}), ...
                        jclu(CLU.(pwfield(pw)).(s).(bandsf{bi}))); %#ok<AGROW>
                end
                entries{end+1}=sprintf('%s:{ %s }', s, strjoin(be,', ')); %#ok<AGROW>
            end
        end
        term=','; if pi==numel(powers), term=''; end
        L=[L sprintf(' %s }%s\n', strjoin(entries,', '), term)];
    end
    L=[L sprintf('  },\n')];

    % per-cell n (assemble_AB intersects subject lists per stage, so a subject
    % missing one stage silently reduces n for that cell only)
    L=[L sprintf('  n_cell: {\n')];
    for pi=1:numel(powers)
        pw=powers{pi}; L=[L sprintf('   %s: {\n',pw)];
        for si=1:numel(stages)
            s=stages{si}; cells={};
            for b=1:numel(bandlabels)
                cells{end+1}=sprintf('"%s":%d',bandlabels{b},NN.(pwfield(pw)).(s).(bf(bandlabels{b}))); %#ok<AGROW>
            end
            term=','; if si==numel(stages), term=''; end
            L=[L sprintf('    %s:{ %s }%s\n',s,strjoin(cells,', '),term)];
        end
        term=','; if pi==numel(powers), term=''; end
        L=[L sprintf('   }%s\n',term)];
    end
    L=[L sprintf('  },\n')];

    % DESCRIPTIVE ONLY - Benjamini-Hochberg FDR q across the per-cell minimum
    % cluster p-values (one per power x stage x band cell). Not the primary
    % inference: each cell's p is already family-wise corrected across channels.
    [FQ, nq] = bh_over_cells(VAL, powers, stages, bandlabels);
    L=[L sprintf('  fdr_note: "DESCRIPTIVE ONLY: Benjamini-Hochberg q across the %d per-cell minimum cluster p-values. Primary inference is the per-cell FWE-corrected p.",\n', nq)];
    L=[L sprintf('  fdr_q: {\n')];
    for pi=1:numel(powers)
        pw=powers{pi}; L=[L sprintf('   %s: {\n',pw)];
        for si=1:numel(stages)
            s=stages{si}; cells={};
            for b=1:numel(bandlabels)
                cells{end+1}=sprintf('"%s":%.3f',bandlabels{b},FQ.(pwfield(pw)).(s).(bf(bandlabels{b}))); %#ok<AGROW>
            end
            term=','; if si==numel(stages), term=''; end
            L=[L sprintf('    %s:{ %s }%s\n',s,strjoin(cells,', '),term)];
        end
        term=','; if pi==numel(powers), term=''; end
        L=[L sprintf('   }%s\n',term)];
    end
    L=[L sprintf('  }\n};\n')];
    js=L;
end

function [FQ, n] = bh_over_cells(VAL, powers, stages, bandlabels)
% Benjamini-Hochberg step-up adjusted p (q) over every power x stage x band
% cell's minimum cluster p. Descriptive only.
    pv=[]; key={};
    for pi=1:numel(powers)
        for si=1:numel(stages)
            for b=1:numel(bandlabels)
                v=VAL.(pwfield(powers{pi})).(stages{si}).(bf(bandlabels{b}));
                pv(end+1)=v(4); %#ok<AGROW>
                key(end+1,:)={pwfield(powers{pi}), stages{si}, bf(bandlabels{b})}; %#ok<AGROW>
            end
        end
    end
    n=numel(pv);
    [sp,ord]=sort(pv(:));
    q=sp .* n ./ (1:n)';
    for i=n-1:-1:1, q(i)=min(q(i),q(i+1)); end
    q=min(q,1);
    qq=zeros(n,1); qq(ord)=q;
    FQ=struct();
    for i=1:n
        FQ.(key{i,1}).(key{i,2}).(key{i,3})=qq(i);
    end
end

function p = default_template()
% templates/sleep_eeg_report.html next to this file. Falls back to the MATLAB
% path so a relocated/deployed copy still resolves.
    p = fullfile(fileparts(mfilename('fullpath')), 'templates', 'sleep_eeg_report.html');
    if ~exist(p,'file')
        w = which('sleep_eeg_report.html');
        if ~isempty(w), p = w; end
    end
end

function inject_report(htmlin, htmlout, js)
    assert(exist(htmlin,'file')==2, 'export_report:noTemplate', ...
        'Report template not found: %s', htmlin);
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
% EEG.bands; fall back to EEG.features (some exports use that field name instead).
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
    q=cellfun(@jstr,c,'uni',0); s=['[' strjoin(q,',') ']'];
end
function s=jstr(x)
% Quoted JS string literal with backslash and double-quote escaped. Condition
% names are user-supplied and are embedded verbatim in title/conditions/dirtxt.
    x=char(x);
    x=strrep(x,'\','\\');
    x=strrep(x,'"','\"');
    x=regexprep(x,'[\r\n\t]',' ');
    s=['"' x '"'];
end
function s=jclu(C)
% JSON array of per-cluster records for one cell.
    if isempty(C), s='[]'; return; end
    e=cell(1,numel(C));
    for i=1:numel(C)
        e{i}=sprintf(['{n:%d,p:%.4f,k:%d,dir:%d,dirtxt:%s,mdiff:%.6g,munit:%s,' ...
                      'mdifft:%.6g,tpeak:%.3f,tpeakch:%s,ch:%s}'], ...
            C(i).num, C(i).p, C(i).nch, C(i).dir, jstr(C(i).dirtxt), ...
            C(i).mdiff, jstr(C(i).munit), C(i).mdifft, ...
            C(i).tpeak, jstr(C(i).tpeakch), jarr(C(i).labels));
    end
    s=['[' strjoin(e,',') ']'];
end
