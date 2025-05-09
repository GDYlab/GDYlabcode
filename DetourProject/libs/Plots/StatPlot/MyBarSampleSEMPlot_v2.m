function MyBarSampleSEMPlot_v2(d,ax,varargin)
%MyBarSampleSEMPlot compare group of observations specified in d and plot in axis ax
% for each group, it plots sample dots, bar of mean, and sem 

% input   d: data to compare, it must be a cell array with length equal or larger than 1,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot
%
%         plotx: x value for each group
%         barwidth : bar width
% example:  
% %          figure
% %          a = rand(100,1);b = rand(20,20); c = rand(1,15);
% %          pab = ranksum(a(:),b(:)); pac = ranksum(a,c); pbc = ranksum(b(:),c);
% %          MyBoxPlot({a,b,c},gca,1)
% %          sigstar({[1,2],[1,3],[2,3]},[pab,pac,pbc])

%% check inputs
args.barwidth = 0.25;
args.plotx = 1:length(d); 
args.colors = Color_FrenchDisp_modify;
args.markersize = 2;
args.linewidth = 1.5;
args.baralpha = 0.6;
args.barplot = 1;


args = parseArgs(varargin, args);

barwidth = args.barwidth;
plotx = args.plotx;
colors = args.colors;
markersize = args.markersize;
linewidth = args.linewidth;


if ~iscell(colors)
    tmpcolors = colors;
    clear colors
    colors = cell(1,2);
    colors{1} = tmpcolors;
    colors{2} = tmpcolors;
end
% % colors = Color_Budapest;

%% do the plot
hold(ax,'on')

for i = 1:length(d)
    colori = rem(i,length(colors));
    colori(colori == 0) = length(colors);
    ptemp = d{i};
    ptemp = ptemp(:);
    ptemp(isinf(ptemp)) = nan;
    ptemp(isnan(ptemp)) = [];
    
    if isempty(ptemp)
       continue 
    end
    if length(ptemp) <= 1
        continue
    end
    ptemp = berow(ptemp);
    % plot mean
    if args.barplot
    bar(ax,plotx(i),nanmean(ptemp),barwidth,'facecolor',colors{colori},'edgecolor',colors{colori},'facealpha',args.baralpha)
    end
    % plot data sample
    xnow = plotx(i) + (rand(1,length(ptemp))-0.5)/5;
    plot(ax,xnow,ptemp,'o','markerfacecolor',0.8*colors{colori},'markeredgecolor',0.8*colors{colori},'markersize',markersize)
    % add SEM  
    seml = nanmean(ptemp) - nansem(ptemp);
    semh = nanmean(ptemp) + nansem(ptemp);
    plot(ax,[plotx(i) plotx(i)],[seml semh],'-','color',0.4*colors{colori},'linewidth',linewidth)
    plot(ax,[plotx(i)-barwidth*0.35 plotx(i)+barwidth*0.35],[seml seml],'-','color',0.4*colors{colori},'linewidth',linewidth)
    plot(ax,[plotx(i)-barwidth*0.35 plotx(i)+barwidth*0.35],[semh semh],'-','color',0.4*colors{colori},'linewidth',linewidth)


end

xticks(ax,plotx)
% % xlim(ax,[0 length(d)+1])
% % set(gcf,'CurrentAxes',ax)
end

