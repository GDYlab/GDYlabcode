function  [hd,per] = My1PDF1VComparePlot_v2(distri,v1,ax,varargin)
%My1PDF1VComparePlot  plot one PDF against one single value

% input   distri, data of one distribution, should be vector or matrix
%         v1, another single value to compare
%         ax, axis handle to plot
%         range,  if specify, its the bin edges to get pdf
%         invflag, if 1, will add text indicating the percentage of v1
%         compare against distri
% output  hd  is the handle with length of 2 matching d
% example hd = My1PDF1VComparePlot(rand(300,1),0.3,gca,[0:0.1:1])

%% check inputs
args.invflag = 1; % show inv pct
args.linewidth = 1.5; % linewidth
args.normalize = 'pdf';  % linewidth
args.linecolor = [218,161,143]/255; % color of the arrow
args.pdfcolor = [0.2 0.2 0.2]; % color of the arrow
args.range = [];
args.ylim = [];
args.txtfont = 12;
args = parseArgs(varargin, args);

distri = distri(:);
% if range is not specified, we will use 20 bins and cover all the datasets
if isempty(args.range)
    mind = nanmin(distri); maxd = nanmax(distri);
    dx = (maxd-mind)/20;
    args.range = mind:dx:maxd;
end
rangec = (args.range(1:end-1)+args.range(2:end))/2;
%% make plot
hold(ax,'on')
hd(1) = histogram(ax,distri,args.range,'LineStyle','none','Facecolor',args.pdfcolor,'Normalization',args.normalize,'FaceAlpha',0.8);
if isempty(args.ylim)
    yl = ylim(ax);
    ydiff = diff(yl);
    yl(2) = yl(2) + ydiff/10;
else
    yl = args.ylim;
end

ylim(ax,yl)
hd(2) = plot(ax,[v1,v1],yl,'-','linewidth',args.linewidth,'color',args.linecolor);
xrange = diff(xlim(ax));

if args.invflag
   per = invprctile(distri,v1);
   text(ax,v1+xrange/35,yl(2)*0.9,[num2str(per,'%.3f'),'%'],'fontsize',args.txtfont)
end

% % histogram(X,edges)
% % 'DisplayStyle','stairs'
% setplot(ax)
% set(gcf,'CurrentAxes',ax)
% ylabel('PDF')
end

