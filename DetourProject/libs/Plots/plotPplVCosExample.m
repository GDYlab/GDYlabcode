function plotPplVCosExample(ax,ses,rplens,it,idir,pcorr,cmap,colors,colors2)
% this function plot an example for the population vector similarity
% for a given detoured track between detour and pre-detour session
% it will highlight the region of detour segment in detour session 
% and other middle segments in pre-detour session
% inputs:    ax, axis to make the plot
%            ses, ses give the track limit across sessions
%            rplens, segment length of detoured tracks
%            it, idir, track and direction to plot
%            pcorr, structure contain population vector similarities
%            cmap,colors,colors2, colormap and colors for the plot

tralinewidth = 5;
% we assume track 12 has the same stationary segments (start and end 50 cm) length
% we assume track 34 has the same stationary segments (start and end 50 cm) length
mapt13tot24 = [2,2,4,4]; 

% find detour session and pre-detour session for this track
detses = Det_FindDetTSes(it,ses);
predetses = detses - 1;
% x axis will be detour and y axis will be predetour
corrtemp = pcorr(idir,detses,predetses).crcef;
corrtemp = corrtemp';

meshnow = pcorr(idir,detses,predetses).plfmesh;
plotx = meshnow(1:size(corrtemp,2));
ploty = meshnow(1:size(corrtemp,1));
pcolor(plotx,ploty,corrtemp)
shading flat; axis tight;
colormap(cmap)
caxis([0 1])

% plot detses tracks
for dit = 1:length(ses(detses).tra_p)
    tralimnow = ses(detses).tralim(dit,:);
    trnum = ses(detses).tra_p(dit);
    if ismember(trnum,[1,2,3,4])
        plot(ax,tralimnow,[-1,-1],'-','linewidth',tralinewidth,'color',colors{trnum})
    else
        detseg = [tralimnow(1)+rplens{dit}(1),tralimnow(2)-rplens{dit}(2)];
        seg1 = [tralimnow(1),tralimnow(1)+rplens{dit}(1)];
        segend = [tralimnow(2)-rplens{dit}(2),tralimnow(2)];
        plot(ax,detseg,[-1,-1],'-','linewidth',tralinewidth,'color',colors{trnum})
        plot(ax,seg1,[-1,-1],'-','linewidth',tralinewidth,'color',colors{dit})
        plot(ax,segend,[-1,-1],'-','linewidth',tralinewidth,'color',colors{dit})
        plot(ax,detseg([1 1]),ploty([1 end]),':','color',[0.6 0.6 0.6],'linewidth',0.5)
        plot(ax,detseg([2 2]),ploty([1 end]),':','color',[0.6 0.6 0.6],'linewidth',0.5)
    end
end
xlabel(ax,'Detour linear position (cm)')

% plot predet ses tracks
rectcolors = {colors2{3},colors2{2},colors2{4},colors2{1}};
for pit = 1:length(ses(predetses).tra_p)
    tralimnow = ses(predetses).tralim(pit,:);
    trnum = ses(predetses).tra_p(pit);
    if ismember(trnum,[1,2,3,4])
        plot(ax,[-1 -1],tralimnow,'-','linewidth',tralinewidth,'color',colors{trnum})
        scseg = [tralimnow(1)+rplens{mapt13tot24(pit)}(1),tralimnow(2)-rplens{mapt13tot24(pit)}(2)];
        plot(ax,plotx([1 end]),scseg([1 1]),':','color',[0.6 0.6 0.6],'linewidth',0.5)
        plot(ax,plotx([1 end]),scseg([2 2]),':','color',[0.6 0.6 0.6],'linewidth',0.5)
        rectangle(ax,'Position',[detseg(1) scseg(1) diff(detseg) diff(scseg)],...
            'EdgeColor',rectcolors{pit},'linewidth',2);
    end
end
ylabel(ax,'Pre-detour linear position (cm)')

cb = addcolorbar(ax,-0.02,0.02);
cb.Ticks = [0 1];
cb.Label.String = 'PV cosine similarity';
cb.Label.Position(1) = cb.Label.Position(1) - 0.2;

hd(1) = plot(ax,nan,nan,'s','markeredgecolor',colors2{1},'linewidth',2);
hd(2) = plot(ax,nan,nan,'s','markeredgecolor',colors2{4},'linewidth',2);
hd(3) = plot(ax,nan,nan,'s','markeredgecolor',colors2{2},'linewidth',2);
hd(4) = plot(ax,nan,nan,'s','markeredgecolor',colors2{3},'linewidth',2);
legend(hd,{'DetvsMobile','DetvsT3','DetvsOpposite','DetvsT1'},'Location','northwestoutside')

end