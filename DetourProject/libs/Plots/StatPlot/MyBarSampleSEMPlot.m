function MyBarSampleSEMPlot(d,ax,colors,plotx,barwidth)
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

if nargin < 5
   barwidth = 0.25;
end

if nargin < 4
   plotx = 1:length(d); 
end

if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
end


if nargin < 3
colors = Color_FrenchDisp;
end

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
    bar(ax,plotx(i),nanmean(ptemp),barwidth,'facecolor',colors{colori},'edgecolor',colors{colori})
    % add SEM  
    seml = nanmean(ptemp) - nansem(ptemp);
    semh = nanmean(ptemp) + nansem(ptemp);
    plot(ax,[plotx(i) plotx(i)],[seml semh],'-','color',0.6*colors{colori},'linewidth',1.5)
    plot(ax,[plotx(i)-barwidth*0.3 plotx(i)+barwidth*0.3],[seml seml],'-','color',0.6*colors{colori},'linewidth',1.5)
    plot(ax,[plotx(i)-barwidth*0.3 plotx(i)+barwidth*0.3],[semh semh],'-','color',0.6*colors{colori},'linewidth',1.5)
    % plot data sample
    xnow = plotx(i) + (rand(1,length(ptemp))-0.5)/5;
    plot(ax,xnow,ptemp,'o','color',0.8*colors{colori},'markersize',2,'linewidth',2)

end

xticks(ax,plotx)
% % xlim(ax,[0 length(d)+1])
% setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

