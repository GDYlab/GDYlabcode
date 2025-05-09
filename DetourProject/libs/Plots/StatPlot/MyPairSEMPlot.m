function  hd = MyPairSEMPlot(d,ax,xplot,plotallsample,varargin)
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
args.CapSize = 4;
args.sigplot = 0;
args = parseArgs(varargin, args);
nsample = size(d,1);
%% make plot
hold(ax,'on')
if plotallsample
    for isample = 1:size(d,1)
        xoffset = (rand(size(xplot))-0.5)./8;
        xplotnow = xplot+xoffset;
        yplotnow = d(isample,:);
        goodidx = ~isnan(yplotnow);
        
        plot(ax,xplotnow(goodidx),yplotnow(goodidx),'--','color',[0.7 0.7 0.7])

        for icat = 1:size(d,2)
            plot(ax,xplotnow(icat),d(isample,icat),'o','color',[0.7 0.7 0.7],'markersize',args.markersize)
        end
    end
end

meannow = nan(size(xplot));
sem = nan(size(xplot));
for icat = 1:size(d,2)
    meannow(icat) = nanmean(d(:,icat));
    stdnow = nanstd(d(:,icat));
    sem(icat) = stdnow./sqrt(nsample);
end
plot(ax,xplot,meannow,'-k','linewidth',args.linewidth)
hd = errorbar(ax,xplot,meannow,sem,'k-','linewidth',args.linewidth,'CapSize',args.CapSize);
% setplot(ax)

if args.sigplot == 1
   [g,p] = AllGroupTtest(d);
   myoverlapsigstar_v2(ax,g,p)
end
end

