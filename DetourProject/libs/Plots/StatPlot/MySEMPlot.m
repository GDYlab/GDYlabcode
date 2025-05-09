function MySEMPlot(d,ax,volinpdf,baseline,colors,plotx,barlen,varargin)
%MySEMPlot compare group of observations specified in d and plot in axis ax

% input   d: data to compare, it must be a cell array with length equal or larger than 1,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot
%
%         volinpdf: if plot the volin pdf or not
% example:  
% %          figure
% %          a = rand(100,1);b = rand(20,20); c = rand(1,15);
% %          [~,pab] = ttest2(a(:),b(:)); [~,pac] = ttest2(a,c); [~,pbc] = ttest2(b(:),c);
% %          MySEMPlot({a,b,c},gca,1)
% %          sigstar({[1,2],[1,3],[2,3]},[pab,pac,pbc])

%% check inputs
if nargin < 7
   barlen = 0.1;
end

if nargin < 6
   plotx = 1:length(d); 
end

if nargin < 4
   baseline = 0; 
end

if nargin < 3
   volinpdf = 0; 
end


if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
end

% % if length(d) < 2
% %    error('input data must have more than 2 groups') 
% % end


if nargin < 5
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

if volinpdf
    barlen = barlen/2;
end

args.markersize = 3; % arrow size
args.linewidth = 1.5; % linewidth
args.connectmean = 0; % connect mean or not
args = parseArgs(varargin, args);
%% do the plot
hold(ax,'on')
allmean = nan(1,length(d));
for i = 1:length(d)
    colori = rem(i,length(colors));
    colori(colori == 0) = length(colors);
    ptemp = d{i};
    ptemp = ptemp(:);
    ptemp(isnan(ptemp)) = [];
    
    if isempty(ptemp)
       continue 
    end
    if length(ptemp) <= 4
        continue
    end

    if i == baseline
        plot(ax,[min(plotx) max(plotx)],[mean(ptemp) mean(ptemp)],'--','color',[0.6 0.6 0.6])
    end
    
    % do the violin plot 
    if volinpdf
        rangetemp = [nanmin(ptemp), nanmax(ptemp)];
        
        if diff(rangetemp) == 0
            rangetemp(2) = rangetemp(2)+1;
            rangetemp(1) = rangetemp(1)-1;
        end
        [hv,hbin] = hist(ptemp,rangetemp(1):diff(rangetemp)/20:rangetemp(2));
        hv = hv./max(hv)*0.4;
        xx = [berow(hv),fliplr(berow(-hv))];
        xx = plotx(i)+xx;
        ytemp = berow(hbin);
        yy = [ytemp,fliplr(ytemp)];
        fill(xx, yy, 1-0.55*(1-colors{colori}),'LineStyle','none','FaceAlpha',0.6);
    end
    
    
    plot(ax,[plotx(i) plotx(i)],[mean(ptemp) mean(ptemp)],'s','markerfacecolor',colors{colori},...
        'markeredgecolor',colors{colori},'linewidth',args.linewidth,'markersize',args.markersize)
    plot(ax,[plotx(i) plotx(i)],mean(ptemp)+[-sem(ptemp) sem(ptemp)],'-','color',colors{colori},'linewidth',args.linewidth)
    plot(ax,plotx(i)+[-barlen barlen]/2,[mean(ptemp)-sem(ptemp) mean(ptemp)-sem(ptemp)],...
        '-','color',colors{colori},'linewidth',args.linewidth)
    plot(ax,plotx(i)+[-barlen barlen]/2,[mean(ptemp)+sem(ptemp) mean(ptemp)+sem(ptemp)],...
        '-','color',colors{colori},'linewidth',args.linewidth)
    allmean(i) = mean(ptemp);
end
if args.connectmean
   plot(ax,plotx,allmean,'-k','linewidth',args.linewidth)
end
xticks(ax,plotx)
% % xlim(ax,[0 length(d)+1])
setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

