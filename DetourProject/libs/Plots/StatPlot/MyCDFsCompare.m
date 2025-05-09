function MyCDFsCompare(d,ax,colors,ldgstr,lgdfont,tail,ppairs)
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

% if length(d) < 2
%    error('input data must have more than groups') 
% end

if nargin < 4
    ldgstr = cell(1,length(d));
    for id = 1:length(d)
    ldgstr{id} = ['Data',num2str(id)];
    end
end

if nargin < 5
    lgdfont = 10;
end

if nargin < 6
    tail = 'unequal';
end

if nargin < 7
    ppairs = [];
end

if nargin < 3
colors = Color_FrenchDisp;
end

if ~iscell(colors)
    tmpcolors = colors;
    clear colors
    colors = cell(1,length(d));
    for id = 1:length(d)
        colors{id} = tmpcolors;
    end
end
% % colors = Color_Budapest;
%% do the plot
hold(ax,'on')
hplot = gobjects(1,length(d));  
for id = 1:length(d)
    dtmp = d{id}(:);
    dtmp(isnan(dtmp)) = [];
    if length(dtmp) < 5
        hplot(id) = plot(ax,nan,nan,'-','color',colors{id},'linewidth',1.5);
       continue 
    end
    [cdfy,cdfx] = ecdf(dtmp);
    hplot(id) = plot(ax,cdfx,cdfy,'-','color',colors{id},'linewidth',1.5);
end
% ylabel(ax,'CDF')
ylim(ax,[0 1])
axis tight
xl = xlim(ax);


xstart = 0.65;
xtend = 0.80;

ydis = 0.1;
ally = nan(1,length(d));
for id = 1:length(d)
    ally(id) = ydis*id;
end

xlgend = [xl(2)*xstart+xl(1)*(1-xstart),xl(2)*xtend+xl(1)*(1-xtend)];
% plot legend at corner
for id = 1:length(d)
    plot(ax,xlgend,[ally(id) ally(id)],'-','color',colors{id},'linewidth',1.5);
    text(ax,xl(2)*(xtend+0.02) + xl(1)*(0.98-xtend),ally(id),ldgstr{id},'FontSize',lgdfont)
end



% do test between groups
if ~isempty(ppairs)
    xdiff = 0.015;
    xoffset = 0.01;
    for ip = 1:size(ppairs,2)
        id = ppairs(1,ip);
        jd = ppairs(2,ip);
        [~,p(ip)] = kstest2(d{id}(:),d{jd}(:),'tail',tail);
        sigstr = sigstrfromp(p);
        % plot sigbar
        xoffset = xoffset + xdiff;
        plot(ax,[xl(2)*(xstart-xoffset)+xl(1)*(1+xoffset-xstart) xl(2)*(xstart-xoffset)+xl(1)*(1+xoffset-xstart)],[ally(id) ally(jd)],'-k')
        text(ax,xl(2)*(xstart-xoffset-0.01)+xl(1)*(1+xoffset+0.01-xstart),mean([ally(id),ally(jd)]),sigstr,'HorizontalAlignment','right','FontSize',lgdfont)
    end
end
end


function sigstr = sigstrfromp(p)

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

end
