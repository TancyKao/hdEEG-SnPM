function export_event_report(DATA, OUT, opts)
% STOPGAP event-analysis report (repurposes the spectral template):
%   tabs  = PARAMETER  (Density / Amplitude / Duration)   <- "By parameter"
%   rows  = event type (slow spindle 9-12, fast spindle 12-15, slow wave SO)
%   power toggle REMOVED (no Absolute/Normalised); periodogram section hidden.
% Two-group contrast (group_col) -> signed t-map (warm = group A > B) + per-group
% mean topos. Reuses core_snpm_analysis (anova1) for the SnPM stats.
%
% USAGE: export_event_report(DATA, OUT, struct('group_col','group','permutations',1000,'html',template))

    if nargin<3, opts=struct(); end
    grpcol = getd(opts,'group_col','group');
    nperm  = getd(opts,'permutations',1000);
    htmlin = getd(opts,'html', fullfile(OUT,'sleep_eeg_report_filled.html'));
    if ~exist(OUT,'dir'), mkdir(OUT); end
    scratch=fullfile(OUT,'_scratch'); if ~exist(scratch,'dir'), mkdir(scratch); end

    Lc=load('egi256_chanlocsCluster'); cl=Lc.EEG_256chansClusters; sub=cl(ismember({cl.type},'EEG'));
    chanlocs=sub(1:256);                                 % drop Cz to match 256-channel vectors
    mlabels=arrayfun(@(x)sprintf('E%d',x),1:256,'UniformOutput',false);

    meta=readtable(fullfile(DATA,'subjects.csv'),'TextType','string','VariableNamingRule','preserve');
    g=categorical(meta.(grpcol)); levels=categories(g); assert(numel(levels)>=2,'need >=2 groups');
    A=char(levels{1}); B=char(levels{2});                % warm = A>B

    params ={'density','amplitude','duration'};          % TABS
    pnames ={'Density','Amplitude','Duration'};
    events ={'spindle_9-12','spindle_12-15','sw_0.5-1.25'};  % ROWS
    enames ={'Slow spindle','Fast spindle','Slow wave'};
    erange ={'9-12 Hz','12-15 Hz','0.5-1.25 Hz'};
    PW='val';                                            % single dummy power token (toggle hidden)
    metacols=[{'Subject'} setdiff(meta.Properties.VariableNames,{'Subject'},'stable')];

    VAL=struct(); CH=struct(); VAL.(PW)=struct(); CH.(PW)=struct();
    for pi=1:numel(params)
        par=params{pi}; G=struct('band',{},'meanA',{},'meanB',{},'tmap',{},'sigU',{},'sigT',{},'clim',{});
        for ei=1:numel(events)
            ev=events{ei};
            M=readtable(fullfile(DATA,sprintf('eventStat_%s_%s.csv',par,ev)),'TextType','string','VariableNamingRule','preserve');
            M=M(:,['Subject', mlabels]);
            T=innerjoin(meta,M,'Keys','Subject');
            csv=fullfile(scratch,sprintf('%s_%s.csv',par,ev)); writetable(T,csv);
            gp=struct('data_file',csv,'data_sheet','CSV File','output_path',scratch, ...
                'channels','164 channels','datatype','absolute','comparison','anova1', ...
                'group_col',grpcol,'permutations',nperm,'meta_cols',{metacols});
            rs=core_snpm_analysis(gp);
            tmap=rs.T.real_T;
            X=T{:,mlabels}; gi=categorical(T.(grpcol));
            meanA=mean(X(gi==A,:),1,'omitnan'); meanB=mean(X(gi==B,:),1,'omitnan');
            if nansum_(tmap.*(meanA-meanB))<0, tmap=-tmap; end   % orient warm = A>B
            iU=rs.uncorrsigch; iT=rs.correctTFCEsigch; iC=rs.SnPMsigch;
            clim=robustclim([meanA meanB]);
            save_topo(meanA,chanlocs,clim,[],[],'parula',fullfile(OUT,sprintf('%s_%s_%s_condA_mean.png',par,PW,ev)));
            save_topo(meanB,chanlocs,clim,[],[],'parula',fullfile(OUT,sprintf('%s_%s_%s_condB_mean.png',par,PW,ev)));
            save_topo(tmap,chanlocs,[-3 3],iU,iT,'jet',fullfile(OUT,sprintf('%s_%s_%s_Tmap.png',par,PW,ev)));
            G(ei)=struct('band',enames{ei},'meanA',meanA,'meanB',meanB,'tmap',tmap,'sigU',iU,'sigT',iT,'clim',clim);
            VAL.(PW).(par).(bf(ev))=[numel(iU) numel(iT) numel(iC) safe_minp(rs)];
            if numel(iU)+numel(iT)+numel(iC)>0
                CH.(PW).(par).(bf(ev))=struct('u',{mlabels(iU)},'tfce',{mlabels(iT)},'cluster',{mlabels(iC)});
            end
            fprintf('  %-9s %-13s: TFCE=%d cluster=%d\n', par, ev, numel(iT), numel(iC));
        end
        make_grid(G,chanlocs,pnames{pi},fullfile(OUT,sprintf('GRID_%s_%s.png',par,PW)),A,B);
    end

    js=build_js(params,pnames,events,enames,erange,PW,VAL,CH,height(meta),A,B,grpcol);
    fid=fopen(fullfile(OUT,'REPORT.js'),'w'); fwrite(fid,js); fclose(fid);
    inject_event(htmlin, fullfile(OUT,'event_report.html'), js);
    fprintf('Done. Event report: %s\n', fullfile(OUT,'event_report.html'));
end

% ---- helpers ----
function v=getd(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end; end
function f=bf(b), f=matlab.lang.makeValidName(b); end
function y=nansum_(x), y=sum(x(~isnan(x))); end
function c=robustclim(v), v=v(~isnan(v)); c=[prctile(v,2) prctile(v,98)]; if ~(c(1)<c(2)), c=[min(v) max(v)]; end; end
function p=safe_minp(rs), if isfield(rs,'Clusters')&&~isempty(rs.Clusters), p=min([rs.Clusters.p]); else, p=1; end; end

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
    Th=pi/180*[chanlocs.theta]; Rd=[chanlocs.radius]; [x,y]=pol2cart(Th,Rd); sq=0.5/max(Rd); X=x*sq; Y=y*sq;
    if ~isempty(blackIdx), plot(Y(blackIdx),X(blackIdx),'.k','MarkerSize',11); end
    if ~isempty(whiteIdx), plot(Y(whiteIdx),X(whiteIdx),'.w','MarkerSize',11); end
    set(gca,'XLim',[-0.62 0.62],'YLim',[-0.62 0.62]);
end
function make_grid(G,chanlocs,parname,fname,A,B)
    nB=numel(G); if nB==0, return; end
    f=figure('Color','w','Position',[20 20 720 max(220,170*nB)],'Visible','off');
    tl=tiledlayout(nB,3,'TileSpacing','compact','Padding','compact');
    title(tl,sprintf('%s  (%s mean | %s mean | t-map)',parname,A,B),'Interpreter','none');
    for r=1:nB
        nexttile; topo_cell(G(r).meanA,chanlocs,G(r).clim,[],[],'parula'); if r==1, title([A ' mean']); end; ylabel(G(r).band);
        nexttile; topo_cell(G(r).meanB,chanlocs,G(r).clim,[],[],'parula'); if r==1, title([B ' mean']); end
        nexttile; topo_cell(G(r).tmap,chanlocs,[-3 3],G(r).sigU,G(r).sigT,'jet'); if r==1, title('t-map'); end
    end
    try, print(f,fname,'-dpng','-r130'); catch, end; close(f);
end

function js=build_js(params,pnames,events,enames,erange,PW,VAL,CH,nsubj,A,B,grpcol)
    L=sprintf('const REPORT = {\n');
    L=[L sprintf('  title: "Sleep events: %s vs %s (NREM2+3)",\n',A,B)];
    L=[L sprintf('  methods: "Per-channel spindle/slow-wave parameters by group (%s). Two-group SnPM, 1000 permutations, TFCE + cluster, \\u03B1=.05. Tabs = parameter; rows = event type.",\n',grpcol)];
    L=[L sprintf('  n_subj: %d,\n',nsubj)];
    L=[L sprintf('  conditions: { A:"%s", B:"%s" },\n',A,B)];
    L=[L sprintf('  powers: [ {key:"%s",name:"value"} ],\n', PW)];      % single, hidden toggle
    pp=arrayfun(@(i)sprintf('{key:"%s",name:"%s"}',params{i},pnames{i}),1:numel(params),'uni',0);
    L=[L sprintf('  stages: [ %s ],\n', strjoin(pp,', '))];             % TABS = parameters
    bb=arrayfun(@(i)sprintf('{key:"%s",name:"%s",gk:"",range:"%s"}',events{i},enames{i},erange{i}),1:numel(events),'uni',0);
    L=[L sprintf('  bands: [ %s ],\n', strjoin(bb,', '))];              % ROWS = event types
    % values: values[PW][param][event]
    L=[L sprintf('  values: {\n   %s: {\n', PW)];
    for pi=1:numel(params)
        cells=arrayfun(@(i)sprintf('"%s":[%d,%d,%d,%.3f]',events{i}, VAL.(PW).(params{pi}).(bf(events{i}))), 1:numel(events),'uni',0);
        term=','; if pi==numel(params), term=''; end
        L=[L sprintf('    %s:{ %s }%s\n', params{pi}, strjoin(cells,', '), term)];
    end
    L=[L sprintf('   }\n  },\n')];
    % channels: channels[PW][param][event]
    L=[L sprintf('  channels: {\n   %s: {\n', PW)];
    for pi=1:numel(params)
        be={};
        if isfield(CH.(PW),params{pi})
            for i=1:numel(events)
                if isfield(CH.(PW).(params{pi}),bf(events{i}))
                    c=CH.(PW).(params{pi}).(bf(events{i}));
                    be{end+1}=sprintf('"%s":{u:%s,tfce:%s,cluster:%s}',events{i},jarr(c.u),jarr(c.tfce),jarr(c.cluster)); %#ok<AGROW>
                end
            end
        end
        term=','; if pi==numel(params), term=''; end
        L=[L sprintf('    %s:{ %s }%s\n', params{pi}, strjoin(be,', '), term)];
    end
    L=[L sprintf('   }\n  }\n};\n')];
    js=L;
end
function s=jarr(c), if isempty(c), s='[]'; return; end; q=cellfun(@(x)['"' char(x) '"'],c,'uni',0); s=['[' strjoin(q,',') ']']; end

function inject_event(htmlin, htmlout, js)
    h=fileread(htmlin);
    p1=strfind(h,'const REPORT = {'); pend=strfind(h,'END OF DATA');
    sl=strfind(h,'/*'); sl=sl(sl<pend(1)); pb=sl(end);
    h=[h(1:p1(1)-1) js sprintf('\n') h(pb:end)];
    % rename the per-stage section heading -> "By parameter"
    h=strrep(h,'>By sleep stage<','>By parameter<');
    h=strrep(h,'Stage \ Band','Parameter \ Event');
    h=strrep(h,'Band grid','Event grid'); h=strrep(h,'per-band detail','per-event detail');
    % hide periodogram section + the power toggle (no Absolute/Normalised)
    hide=['<style>#powerSeg{display:none}</style>' ...
          '<script>(function(){var g=document.getElementById("perioGrid");if(g){var s=g.closest("section");if(s)s.style.display="none";}})();</script>'];
    if contains(h,'</body>'), h=strrep(h,'</body>',[hide '</body>']); else, h=[h hide]; end
    fid=fopen(htmlout,'w'); fwrite(fid,h); fclose(fid);
end
