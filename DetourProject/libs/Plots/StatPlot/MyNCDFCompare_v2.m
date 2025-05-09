function MyNCDFCompare_v2(d,ax,colors,ldgstr)
%MyNCDFCompare compare CDF of several data samples

% input   d: data to compare, it must be a cell array with length over 2,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot

% example:  
% %          a = rand(100,1);b = rand(20,20);c = 0.2+rand(10,20)
% %          MyNCDFCompare({a,b,c},gca)

%% check inputs

if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
end

if length(d) < 2
   error('input data must have no less than 2 groups') 
end

if nargin < 4
ldgstr = [];
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
hplot = gobjects(1,length(d));  
for id = 1:length(d)
    dtmp = d{id}(:);
    dtmp(isnan(dtmp)) = [];
    [N,EDGES] = histcounts(dtmp,1e3);
    binc = (EDGES(1:end-1) + EDGES(2:end))/2;
    cumN = cumsum(N)/sum(N);
    hplot(id) = plot(ax,binc,cumN,'-','color',colors{id},'linewidth',1.5);
end
if ~isempty(ldgstr)
    legend(hplot,ldgstr)
end
% % set(gcf,'CurrentAxes',ax)
end

