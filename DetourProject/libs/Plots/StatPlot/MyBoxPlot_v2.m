function MyBoxPlot_v2(d,ax,varargin)
% volinpdf,baseline,colors,plotx,boxlen,rangebin,rawd,vlscale)
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

args.vlscale = 0.4; 
args.rawd = 0; 
args.volinpdf = 1; % color of the arrow
args.baseline = 0;
args.sigplot = 1;
args.rangebin = [];
args.boxlen = 0.15;
args.meanlinewidth = 2;
args.plotx = 1:length(d);
args.colors = Color_FrenchDisp_modify;

args = parseArgs(varargin, args);

vlscale = args.vlscale;
rawd = args.rawd;
volinpdf = args.volinpdf;
baseline = args.baseline;
sigplot = args.sigplot;
rangebin = args.rangebin;
boxlen = args.boxlen;
plotx = args.plotx;
colors = args.colors;

if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
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
    if isempty(ptemp)
        continue
    end
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
        plot(ax,[min(plotx) max(plotx)],[mediap mediap],'--','color',[0.6 0.6 0.6],'linewidth',1)
    end
    
    if rawd
       xnow = plotx(i) + (rand(1,length(ptemp))-0.5)/5;
%        plot(ax,xnow,ptemp,'o','color',1-0.8*(1-colors{colori}),'markersize',3)
       plot(ax,xnow,ptemp,'o','color',0.5*colors{colori},'markersize',3)
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
    plot(ax,[plotx(i)-boxlen plotx(i)+boxlen],[mediap mediap],'-','color',[1 1 1],'linewidth',args.meanlinewidth);


    
end
try
    ylim(ax,[allmin allmax])
end

if sigplot
    % find groups with large sample numbers
    actgroup = 1:length(d);
    tmpd = d;
    todel = [];
    for ig = 1:length(d)
        dnow = d{ig};
        validd = dnow(~isnan(dnow));
        if length(validd) < 5
            todel = cat(1,todel,ig);
        end
    end
    actgroup(todel) = [];
    tmpd(todel) = [];
    
    try
        [g,p] = AllGroupRankSum(tmpd,actgroup);
        sigstar(g,p)
        %myoverlapsigstar_v2(ax,g,p)
    end
end
xticks(ax,plotx)
xlim(ax,[min(plotx) - vlscale *1.2, max(plotx) + vlscale *1.2])
set(ax, 'YGrid', 'on', 'XGrid', 'off')
% % xlim(ax,[0 length(d)+1])
% setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

