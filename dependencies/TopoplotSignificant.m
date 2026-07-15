function TopoplotSignificant(topodata,comparison, alpha, roma, chanlocs,insidegoodch,data1_file, data2_file, Exportfig, titletext)


if contains(data1_file, '.')
    dat1 = extractBefore(data1_file,'_EEG.');
    dat2 = extractBefore(data2_file,'_EEG.');
else
    dat1 = data1_file;
    dat2 = data2_file;
end



numbands = size(topodata,2);
topofig = figure('position',[20 50 600 150*numbands]);
subplot('position',[.1 .96 .8 .02])

% if strcmp(dat1, 'THC_CBD')
%     dat1 = 'THC/CBD';
% end

figlabel = [topodata(1).stage,' ',titletext];

colordotlabel1 = ['B: ucSig']; 
colordotlabel2 = ['W: corSig']; 

UncMarkerColor = [0.9100, 0.7610, 0.0900];

%text(.5,.5,figlabel,'FontSize',12,'fontweight','bold','horizontalalignment','center');
text(0.8,-.5,colordotlabel1,'FontSize',10,'horizontalalignment','right');
text(0.8,-1.10,colordotlabel2,'FontSize',10,'horizontalalignment','right');

if (strcmp(comparison,'pairedT') ==1) | (strcmp(comparison,'unpairedT') ==1)
    colorlabel1 = ['warm: ', dat1, ' > ' dat2]; 
else
    colorlabel1 = ['warm: positive corr']; 
end



colormap(jet);
axis off


for kfreq = 1 : numbands
    
    pltdata = topodata(kfreq).real_T;
    sigch_unP       = find(topodata(kfreq).p_real <= alpha);  % uncorrected p<0.05
    sigch_TFCEP = find(topodata(kfreq).p_correctedTFCE <= alpha); % snpm suprathreshold p<0.05
    sigch_clusterP = topodata(kfreq).p_cluster_chans; % snpm cluster p value
    
    
    figure(topofig);
    subplot('position',[.15 (1-.07)*(numbands-kfreq)/numbands .7 (1-.07)/numbands]);

    if ~exist(roma)
        topoplot(pltdata,chanlocs,'headrad',0.57,'style','map','electrodes','on',...
                                  'maplimits','minmax','whitebk','on');
    else
        topoplot(pltdata,chanlocs,'headrad',0.57,'style','map','electrodes','on',...
                               'colormap', roma,'maplimits','minmax','whitebk','on');
    end

    caxis manual
    if strcmp(comparison,'pairedT') | strcmp(comparison,'unpairedT')
        caxis([-3 3]);
    elseif strcmp(comparison,'correlationP') | strcmp(comparison,'correlationS')
        caxis([-1 1]);
    end
    

    set(gca,'Xlim',[-.55 .55],'Ylim',[-.59 .59]);
    electrodes.x=get(findobj(gca,'Marker','.'),'XData'); 
    electrodes.y=get(findobj(gca,'Marker','.'),'YData'); 
    electrodes.z=get(findobj(gca,'Marker','.'),'ZData');
    delete(findobj(gca,'Marker','.'));

    % if you use inside good, be careful to pull correct indicies 
    if ~isequal(length(electrodes.x),length(insidegoodch))
        disp('error with plotting sig channels');
    else

        if sum(sigch_unP)>0
           [~,chi] = intersect(insidegoodch,sigch_unP);
           hold on;
           if ~exist(roma)
               hh = scatter(electrodes.x(chi),electrodes.y(chi),electrodes.z(chi),...
                    'filled','SizeData',25,'Cdata',[0 0 0],'MarkerEdgeColor',[0.5 0.5 0.5],'linewidth',.5);
           else
               hh = scatter(electrodes.x(chi),electrodes.y(chi),electrodes.z(chi),...
                    'filled','SizeData',25,'Cdata',UncMarkerColor,'MarkerEdgeColor',UncMarkerColor,'linewidth',.5);

           end
        
        end

        if contains(titletext, 'TFCE')==1
            if sum(sigch_TFCEP)>0
               [~,chi2] = intersect(insidegoodch,sigch_TFCEP);
               hold on;
               hh = scatter(electrodes.x(chi2),electrodes.y(chi2),electrodes.z(chi2),...
                    'filled','SizeData',25,'Cdata',[1 1 1],'MarkerEdgeColor',[1 1 1],'linewidth',.5);
            
            end
        elseif contains(titletext, 'Cluster') ==1
            if sum(sigch_clusterP)>0
               [~,chi2] = intersect(insidegoodch,sigch_clusterP);
               hold on;
               hh = scatter(electrodes.x(chi2),electrodes.y(chi2),electrodes.z(chi2),...
                    'filled','SizeData',22,'Cdata',[1 1 1],'MarkerEdgeColor',[1 1 1],'linewidth',.5);
            
            end
        end

    end

 
    text(-.8, 0, char(topodata(kfreq).freqband), 'FontSize', 18, 'rotation', 90, 'HorizontalAlignment', 'center', 'Interpreter', 'none');

    hh=colorbar;
    
    t=get(hh,'Limits');
    set(hh,'Ticks',linspace(t(1),t(2),3))
    set(hh,'position',[.68 (1-.05)*(numbands-kfreq)/numbands+0.01 .03 (1-.3)/(numbands)]);
    %freqband_display = strrep(topodata(kfreq).freqband, '_', '_');
    %text(-.8,0,char(freqband_display),'FontSize',18,'rotation',90,'HorizontalAlignment','center');

    figure(topofig);
    subplot('position',[.15 (1-.07)*(numbands-kfreq)/numbands .7 (1-.07)/numbands]);
    
    clear hh        

    set(gcf,'InvertHardCopy','off','color','w');

    
    sgtitle({figlabel,colorlabel1}, 'fontsize',14);
    
end

Exportpng = [Exportfig,'_',titletext,'_',num2str(alpha),'.png'];
set(topofig,'color','w','paperpositionmode','auto');
%print(topofig,'-dpng','-r300',[filepath,topodata(1).stage,'_',titletext,'_',num2str(alpha),'.png']);
print(topofig,'-dpng','-r300',Exportpng);
