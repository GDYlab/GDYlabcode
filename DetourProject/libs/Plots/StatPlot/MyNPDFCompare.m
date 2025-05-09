function MyNPDFCompare(d,ax,varargin)
%MyNCDFCompare compare PDF of several data samples

% input   d: data to compare, it must be a cell array with length over 2,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot

% example:  
% %          a = rand(100,1);b = rand(20,20);c = 0.2+rand(10,20)
% %          MyNPDFCompare({a,b,c},gca)

%% check inputs

if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
end

if length(d) < 2
   error('input data must have no less than 2 groups') 
end

args.colors = Color_FrenchDisp_modify;
args.ldgstr = [];
args.lwid = 1;
args = parseArgs(varargin, args);


colors = args.colors;
ldgstr = args.ldgstr;
lwid = args.lwid;
%% do the plot
hold(ax,'on')
hplot = gobjects(1,length(d));  
for id = 1:length(d)
    dtmp = d{id}(:);
    dtmp(isnan(dtmp)) = [];
    [N,EDGES] = histcounts(dtmp,100);
    binc = (EDGES(1:end-1) + EDGES(2:end))/2;
    %     cumN = cumsum(N)/sum(N);
    histogram(ax,dtmp,15,'DisplayStyle','stairs','linewidth',lwid,'Edgecolor',colors{id},'Normalization','pdf');
     hplot(id) = plot(ax,nan,nan,'-','color',colors{id},'linewidth',1.5);
end
if ~isempty(ldgstr)
    [lgd,icon] = legend(hplot,ldgstr);
    for id = 1:length(d)
       icon(id).Position(1) = icon(id).Position(1) - 0.35;
       icon(1+length(d) + 2*(id-1)).XData(2) = icon(1+length(d) + 2*(id-1)).XData(2) - 0.35;
    end
    lgd.Position(3) = lgd.Position(3) -0.05;
    lgd.Position(4) = lgd.Position(4) -0.03;
    lgd.Location = 'northwest';
    lgd.FontSize = 6;
    
    legend boxoff
end
% % set(gcf,'CurrentAxes',ax)
end

