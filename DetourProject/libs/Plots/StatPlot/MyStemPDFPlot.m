function MyStemPDFPlot(d,shf,ax,varargin)
% 
%MyStemPDFPlot plot data stem along with it's shuffle control in ax
% input   d, data, it must be a vector 
%         shf, shuffle, must be cell with same length with data
%         each element contain the shuffle distribution
%
%         ax: axis handle to plot
%

% example:  




%% check inputs

args.facealpha = 0.4; 
args.linewidth = 2;
args.plotx = 1:length(d);
args.pdfbin = 40;
args.colors = Color_FrenchDisp_modify;
args.stargap = 0.0015;

args = parseArgs(varargin, args);

linewidth = args.linewidth;
plotx = args.plotx;
colors = args.colors;
facealpha = args.facealpha;
pdfbin = args.pdfbin;
stargap = args.stargap;


if length(d) ~= length(shf)
   error('Input data and shfpdf must have the same length') 
end

if ~iscell(colors)
    tmpcolors = colors;
    clear colors
    colors = cell(1,2);
    colors{1} = tmpcolors;
    colors{2} = tmpcolors;
end
% % colors = Color_Budapest;

hold(ax,'on')
nplot = length(d);
for iplot = 1:nplot
    % plot data
    %         stem(ax,ifact-0.15,betaori(ifact),'filled','color',colors{ifact},'linewidth',1.5)
    stem(ax,plotx(iplot)-0.15,d(iplot),'filled','color',colors{iplot},'linewidth',linewidth)
    % plot shuffle distribution
    shfd = shf{iplot};
    shfd = shfd(~isnan(shfd));
    [Nc,xbin] = histcounts(shfd,pdfbin);
    xbinc = (xbin(1:end-1) + xbin(2:end))/2;
    xx = berow(Nc)./max(Nc)*0.55;
    yy = berow(xbinc);
    xx = cat(2,xx,zeros(size(xx)));
    yy = cat(2,yy,fliplr(yy));
    fill(ax,xx+plotx(iplot)+0.05,yy,1-0.6*(1-colors{iplot}),'LineStyle','none','FaceAlpha',facealpha);
    pct = invprctile(shfd,d(iplot));
    pval = (100-pct)/100;
    myoverlapsigstar(ax,{[plotx(iplot)-0.15]},pval,9,d(iplot)*1.02,stargap)
    
end



end