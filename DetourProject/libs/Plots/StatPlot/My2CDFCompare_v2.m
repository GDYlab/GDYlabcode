function p = My2CDFCompare_v2(d,ax,colors,ldgstr,lgdfont,tail)
%My2CDFCompare compare CDF of two data samples, run the Two-sample Kolmogorov-Smirnov test
% and plot the results

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
    ldgstr = {'Data1','Data2'};
end

if nargin < 5
    lgdfont = 10;
end

if nargin < 6
    tail = 'unequal';
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
% ylabel(ax,'CDF')
ylim(ax,[0 1])
axis tight
xl = xlim(ax);

% do test between groups
[~,p] = kstest2(d{1}(:),d{2}(:),'tail',tail);
% p = ranksum(d{1}(:),d{2}(:),'tail',tail);
if p >= 0.05
   sigstr = 'n.s.'; 
end
if p < 0.05 && p >= 0.01
   sigstr = '*'; 
end
if p < 0.01 && p >= 0.001
   sigstr = '**'; 
end
if p < 0.001
   sigstr = '***'; 
end

xstart = 0.65;
xtend = 0.80;
y1 = 0.2;
y2 = 0.1;
xlgend = [xl(2)*xstart+xl(1)*(1-xstart),xl(2)*xtend+xl(1)*(1-xtend)];
% plot legend at corner
plot(ax,xlgend,[y1 y1],'-','color',colors{1},'linewidth',1.5);
text(ax,xl(2)*(xtend+0.02) + xl(1)*(0.98-xtend),y1,ldgstr{1},'FontSize',lgdfont)
plot(ax,xlgend,[y2 y2],'-','color',colors{2},'linewidth',1.5);
text(ax,xl(2)*(xtend+0.02) + xl(1)*(0.98-xtend),y2,ldgstr{2},'FontSize',lgdfont)

% plot sigbar
plot(ax,[xl(2)*(xstart-0.02)+xl(1)*(1.02-xstart) xl(2)*(xstart-0.02)+xl(1)*(1.02-xstart)],[y1 y2],'-k')
text(ax,xl(2)*(xstart-0.03)+xl(1)*(1.03-xstart),mean([y1,y2]),sigstr,'HorizontalAlignment','right','FontSize',lgdfont)
% % set(gcf,'CurrentAxes',ax)
end

