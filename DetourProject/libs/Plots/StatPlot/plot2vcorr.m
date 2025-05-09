function [r,p] = plot2vcorr(ax,d,varargin)
% this function plot the correlation bewteen pair of data
% the pearson correlation will be computed and printed
% inputs:  ax, axis to plot
%           d, n*2 matrix, rows are realizations, correlation are plot and
%           computed between columns

r = nan;
p = nan;
%% default inputs
args.rpval = 1; % show R and p val
args.rpsigstaronly = 0; % show p val by sigstar
args.linefit = 1; % do linefit and plot
args.markersize = 2; % size of data point
args.linewidth = 1.5; % linewidth
args.color = [218,161,143]/255; % color of the plot
args.marker = 'o';
args.txtfont = 12;
args.zeroaxis = 0;
args.tail = 'both';

args = parseArgs(varargin, args);

%% make the plot
if isempty(d)
   return 
end
px = d(:,1);
py = d(:,2);
goodind = ~isnan(px) & ~isnan(py);
px = px(goodind); py = py(goodind);



hold(ax,'on')

c = polyfit(px,py,1);
y_est = polyval(c,px);
plot(ax,px,py,'o','markerfacecolor',args.color,'markeredgecolor',args.color,'markersize',args.markersize);




if args.linefit
    plot(ax,px,y_est,'-','linewidth',args.linewidth,'color',0.7*args.color)
end
if args.rpval
    xl = xlim(ax);
    yl = ylim(ax);
    [r,p] = corr(becolumn(px),becolumn(py),'tail',args.tail);
    if args.rpsigstaronly
        stars = getsigstarfromp(p);
        text(ax,0.95*xl(1)+0.05*xl(2),0.1*yl(1)+0.9*yl(2),stars,...
            'FontSize',args.txtfont)
    else
        text(ax,0.95*xl(1)+0.05*xl(2),0.1*yl(1)+0.9*yl(2),{['R =',num2str(r,'%.3f')],['p = ',num2str(p,'%.4f')]},...
            'FontSize',args.txtfont)
    end
end

if args.zeroaxis
    hy = plot(ax,[0 0],ylim(ax),'-','linewidth',0.5,'color',[0.4 0.4 0.4]);
    hx = plot(ax,xlim(ax),[0 0],'-','linewidth',0.5,'color',[0.4 0.4 0.4]);
    uistack(hx, 'bottom');
    uistack(hy, 'bottom');
end

end

function stars = getsigstarfromp(p)

    if p<=1E-3
        stars='***'; 
    elseif p<=1E-2
        stars='**';
    elseif p<=0.05
        stars='*';
    elseif isnan(p)
        stars='n.s.';
    else
        stars='n.s.';
    end
end
