function My2CDFCompare(d,ax,colors,ldgstr,ranksump)
%My2CDFCompare compare CDF of two data samples

% input   d: data to compare, it must be a cell array with length of 2,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot

% example:  
% %          a = rand(100,1);b = rand(20,20);
% %          My2CDFCompare({a,b},gca)

%% check inputs

if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
end

if length(d) ~= 2
   error('input data must have 2 groups') 
end

if nargin < 4
ldgstr = [];
end

if nargin < 5
ranksump = 1;
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
    [cdfy,cdfx] = ecdf(dtmp);
    hplot(id) = plot(ax,cdfx,cdfy,'-','color',colors{id},'linewidth',1.5);
end
if ~isempty(ldgstr)
    legend(hplot,ldgstr)
end
if ranksump 
    p = ranksum(d{1},d{2});
    title(ax,['p',num2str(p)],'fontweight','normal')
end
% % set(gcf,'CurrentAxes',ax)
end

