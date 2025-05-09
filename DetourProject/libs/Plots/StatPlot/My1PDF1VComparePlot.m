function  [hd,per] = My1PDF1VComparePlot(distri,v1,ax,range,invflag,pdfcolor,normalize,txtfont)
%My1PDF1VComparePlot  plot one PDF against one single value

% input   distri, data of one distribution, should be vector or matrix
%         v1, another single value to compare
%         ax, axis handle to plot
%         range,  if specify, its the bin edges to get pdf
%         invflag, if 1, will add text indicating the percentage of v1
%         compare against distri
% output  hd  is the handle with length of 2 matching d
% example hd = My1PDF1VComparePlot(rand(300,1),0.3,gca,[0:0.1:1])

% % colors = {[0.83 0.42 0.4],[0.52 0.45 0.83],[0.68 0.88 0.43],[0.5 0.73 0.93],...
% %     [0.75 0.42 0.92],[0.94 0.64 0.94],[0.74 0.74 0.22],[0.22 0.74 0.74],[0.2 0.32 0.62]}; % colors used for plotting
colors = Color_FrenchDisp;
% % colors = Color_Budapest;

%% check inputs
if nargin>=6
    colors{1} = pdfcolor; 
end

if nargin<8
    txtfont = 14; 
end

if nargin<7
    normalize = 'pdf'; 
end


if nargin<5
    invflag = 1; 
end

if nargin<4
    range = []; 
end
distri = distri(:);
% if range is not specified, we will use 20 bins and cover all the datasets
if isempty(range)
    mind = nanmin(distri); maxd = nanmax(distri);
    dx = (maxd-mind)/20;
    range = mind:dx:maxd;
end
rangec = (range(1:end-1)+range(2:end))/2;
%% make plot
% % set(gcf,'CurrentAxes',ax)
% % hold on
hold(ax,'on')
% % [h1,~] = histcounts(d1,range);
% % h1 = h1./sum(h1)./dx; % normalize counts to pdf
% % [h2,~] = histcounts(d2,range);

% hd(1) = bar(rangec,h1,1,'LineStyle','none','Facecolor',colors{1},'FaceAlpha',0.8);
hd(1) = histogram(ax,distri,range,'LineStyle','none','Facecolor',colors{1},'Normalization',normalize,'FaceAlpha',0.8);
yl = ylim(ax);
ydiff = diff(yl);
yl(2) = yl(2) + ydiff/10;
ylim(ax,yl)
hd(2) = plot(ax,[v1,v1],yl,'-','linewidth',3,'color',colors{2});
xrange = diff(xlim(ax));

if invflag
   per = invprctile(distri,v1);
   text(ax,v1+xrange/35,yl(2)*0.9,[num2str(per,'%.3f'),'%'],'fontsize',txtfont)
end

% % histogram(X,edges)
% % 'DisplayStyle','stairs'
% setplot(ax)
% set(gcf,'CurrentAxes',ax)
% ylabel('PDF')
end

