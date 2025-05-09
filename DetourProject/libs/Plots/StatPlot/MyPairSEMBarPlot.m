function  hd = MyPairSEMBarPlot(d,ax,xplot,plotallsample,varargin)
%MyPairSEMPlot  plot pairs comparison with SEM for each category

% input   d  is a matrix with rows being samples and columns being
% different categories to compare
%         ax axis handle to plot

% output  hd  is the handle with length of 2 matching d
% example hd = MyPairSEMPlot(rand(20,3),gca)
% % colors = {[0.83 0.42 0.4],[0.52 0.45 0.83],[0.68 0.88 0.43],[0.5 0.73 0.93],...
% %     [0.75 0.42 0.92],[0.94 0.64 0.94],[0.74 0.74 0.22],[0.22 0.74 0.74],[0.2 0.32 0.62]}; % colors used for plotting
% % colors = Color_Budapest;

%% check inputs
if nargin < 3
   xplot = 1:size(d,2); 
end
if nargin < 4
   plotallsample = 1;
end

args.markersize = 3; % marker size
args.linewidth = 2; % linewidth
args.samplelinewidth = 0.5; % linewidth
args.CapSize = 4;
args.sigplot = 0;
args.barwidth = 0.4;
args.emptyface = 0;
args.colors = Color_FrenchDisp_modify;
args = parseArgs(varargin, args);
nsample = size(d,1);

if ~iscell(args.colors)
    tmpcolors = args.colors;
    clear colors
    args.colors = cell(1,2);
    args.colors{1} = tmpcolors;
    args.colors{2} = tmpcolors;
end
%% make plot
hold(ax,'on')

meannow = nan(size(xplot));
sem = nan(size(xplot));
for icat = 1:size(d,2)
    meannow(icat) = nanmean(d(:,icat));
    stdnow = nanstd(d(:,icat));
    sem(icat) = stdnow./sqrt(nsample);
end

for icat = 1:size(d,2)
    colori = rem(icat,length(args.colors));
    colori(colori == 0) = length(args.colors);
    
    if args.emptyface
        bar(ax,xplot(icat),meannow(icat),args.barwidth,'edgecolor',0.8*args.colors{colori},'facecolor',[1,1,1])
    else
        bar(ax,xplot(icat),meannow(icat),args.barwidth,'facecolor',args.colors{colori},'edgecolor',0.8*args.colors{colori})
    end
end


if plotallsample
    for isample = 1:size(d,1)
        xoffset = (rand(size(xplot))-0.5)./5;
        xplotnow = xplot+xoffset;
        yplotnow = d(isample,:);
        goodidx = ~isnan(yplotnow);
        
        plot(ax,xplotnow(goodidx),yplotnow(goodidx),':','color',[0.8 0.8 0.8],'linewidth',args.samplelinewidth)

        for icat = 1:size(d,2)
            colori = rem(icat,length(args.colors));
            colori(colori == 0) = length(args.colors);
            
            plot(ax,xplotnow(icat),d(isample,icat),'o',...
                'markerfacecolor',1-0.3*(1-args.colors{colori}),...
                'markeredgecolor',0.8*(args.colors{colori}),...
                'markersize',args.markersize)
        end
    end
end

for icat = 1:size(d,2)
    colori = rem(icat,length(args.colors));
    colori(colori == 0) = length(args.colors);
    plot(ax,xplot(icat),meannow(icat),'-','linewidth',args.linewidth,'color',0.5*args.colors{colori})
    hd = errorbar(ax,xplot(icat),meannow(icat),sem(icat),'-',...
        'linewidth',args.linewidth,'CapSize',args.CapSize,'color',0.5*args.colors{colori});
end
% setplot(ax)

if args.sigplot == 1
   [g,p] = AllGroupTtest(d);
   myoverlapsigstar_v2(ax,g,p/2)
end

end

