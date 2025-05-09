function MySEMPlot_v2(d,ax,colors,plotx,barlen,varargin)
%MySEMPlot_v2 plot mean and sem specified in d and plot in axis ax

% input   d: data to compare, matrix with rows being categories, columns
%         being mean and sem
%
%         ax: axis handle to plot
%
% example:  
% %          figure
% %          a = rand(100,1);b = rand(20,20); c = rand(1,15);
% %          d = [mean(a),sem(a);mean(b),sem(b),mean(c),sem(c)];
% %          [~,pab] = ttest2(a(:),b(:)); [~,pac] = ttest2(a,c); [~,pbc] = ttest2(b(:),c);
% %          MySEMPlot_v2(d,gca)
% %          sigstar({[1,2],[1,3],[2,3]},[pab,pac,pbc])

%% check inputs
if nargin < 5
   barlen = 0.1;
end

if nargin < 4
   plotx = 1:size(d,1); 
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

args.markersize = 3; % arrow size
args.linewidth = 1.5; % linewidth
args = parseArgs(varargin, args);
%% do the plot
hold(ax,'on')
for i = 1:size(d,1)
    colori = rem(i,length(colors));
    colori(colori == 0) = length(colors);
    if args.markersize > 0
        plot(ax,[plotx(i) plotx(i)],[d(i,1) d(i,1)],'s','markerfacecolor',colors{colori},...
            'markeredgecolor',colors{colori},'linewidth',args.linewidth,'markersize',args.markersize)
    end
    plot(ax,[plotx(i) plotx(i)],d(i,1)+[-d(i,2) d(i,2)],'-','color',colors{colori},'linewidth',args.linewidth)
    plot(ax,plotx(i)+[-barlen barlen]/2,[d(i,1)-d(i,2) d(i,1)-d(i,2)],...
        '-','color',colors{colori},'linewidth',args.linewidth)
    plot(ax,plotx(i)+[-barlen barlen]/2,[d(i,1)+d(i,2) d(i,1)+d(i,2)],...
        '-','color',colors{colori},'linewidth',args.linewidth)
    
end
try
xticks(ax,plotx)
end
% % xlim(ax,[0 length(d)+1])
% setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

