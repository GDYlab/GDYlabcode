function MyBoxPlot(d,ax,volinpdf,baseline,colors,plotx,boxlen,rangebin,rawd,vlscale)
%MyBoxPlot compare group of observations specified in d and plot in axis ax

% input   d: data to compare, it must be a cell array with length equal or larger than 1,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot
%
%         volinpdf: if plot the volin pdf or not
%         rawd : if plot raw data points
% example:  
% %          figure
% %          a = rand(100,1);b = rand(20,20); c = rand(1,15);
% %          pab = ranksum(a(:),b(:)); pac = ranksum(a,c); pbc = ranksum(b(:),c);
% %          MyBoxPlot({a,b,c},gca,1)
% %          sigstar({[1,2],[1,3],[2,3]},[pab,pac,pbc])

%% check inputs
if nargin < 10
   vlscale = 0.4;
end

if nargin < 9
   rawd = 0;
end

if nargin < 8
   rangebin = [];
end

if nargin < 7
   boxlen = 0.15;
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
colors = Color_FrenchDisp_modify;
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
    boxlen = boxlen/2;
end
%% do the plot
hold(ax,'on')
allmin = [];
allmax = [];
for i = 1:length(d)
    colori = rem(i,length(colors));
    colori(colori == 0) = length(colors);
    ptemp = d{i};
    ptemp = ptemp(:);
    ptemp(isinf(ptemp)) = nan;
    ptemp(isnan(ptemp)) = [];
    ptemp = real(ptemp);
    
    if isempty(ptemp)
       continue 
    end
    if length(ptemp) <= 2
        continue
    end
    ptemp = berow(ptemp);
    [ptemp,~] = sort(ptemp,'ascend');
    dlen = length(ptemp);
    
    boxedge = [ptemp(max(1,round(dlen/4))),ptemp(round(3*dlen/4))];
    wisrange = [ptemp(max(1,round(dlen*0.05))),ptemp(round(dlen*0.95))];
    pdfrange = [ptemp(max(1,round(dlen*0.01))),ptemp(round(dlen*0.99))];
    allmin = min([allmin,pdfrange(1)]);
    allmax = max([allmax,pdfrange(2)]);
    mediap = median(ptemp);
    
    if i == baseline
        plot(ax,[min(plotx) max(plotx)],[mediap mediap],'--','color',[0.6 0.6 0.6],'linewidth',1.5)
    end
    
    if rawd
       xnow = plotx(i) + (rand(1,length(ptemp))-0.5)/5;
       plot(ax,xnow,ptemp,'o','color',1-0.8*(1-colors{colori}),'markersize',2)
    end
    
    % do the violin plot 
    if volinpdf
        rangetemp = [nanmin(ptemp), nanmax(ptemp)];
        
        if diff(rangetemp) == 0
            rangetemp(2) = rangetemp(2)+1;
            rangetemp(1) = rangetemp(1)-1;
        end
        if isempty(rangebin)
            [hv,hbin] = hist(ptemp,rangetemp(1):diff(rangetemp)/20:rangetemp(2));
        else
            [hv,hbin] = hist(ptemp,rangebin);
        end
        hv = hv./max(hv)*vlscale;
        xx = [berow(hv),fliplr(berow(-hv))];
        xx = plotx(i)+xx;
        ytemp = berow(hbin);
        yy = [ytemp,fliplr(ytemp)];
        fill(ax,xx, yy, 1-0.55*(1-colors{colori}),'LineStyle','none','FaceAlpha',0.6);
    end
    
    
    rectangle(ax,'Position',[plotx(i)-boxlen,boxedge(1),boxlen*2,diff(boxedge)],'facecolor',colors{colori},'edgecolor',colors{colori})
    plot(ax,[plotx(i) plotx(i)],[boxedge(1) wisrange(1)],'-','color',colors{colori},'linewidth',1.5)
    plot(ax,[plotx(i) plotx(i)],[boxedge(2) wisrange(2)],'-','color',colors{colori},'linewidth',1.5)
    plot(ax,[plotx(i)-boxlen plotx(i)+boxlen],[mediap mediap],'-','color',[1 1 1],'linewidth',2);


    
end
try
ylim(ax,[allmin allmax])
end
xticks(ax,plotx)
% % xlim(ax,[0 length(d)+1])
setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

