function addSEMBartoPlot(ax,x,d,varargin)
% this function add error bar computed as SEM error to a given plot
% inputs ax  axis to plot
%         x  x-values to plot error bars
%         d  cell array of data, length should match x, SEM computed for
%            each element 


%% default inputs
args.capsize = 0.3; % length of bar cap
args.linewidth = 1.5; % linewidth
args.color = [0 0 0]; % color of the errorbar
args = parseArgs(varargin, args);


%% check inputs
if length(x) ~=  length(d)
   error('Length of group should match x-coordinate of plot') 
end

%% plot sem error bar
hold(ax,'on')
for iele = 1:length(d)
    dnow = d{iele};
    dnow(isnan(dnow)) = [];
    meannow = mean(dnow(:));
    semnow = sem(dnow(:));
    plot(ax,[x(iele),x(iele)],[meannow-semnow,meannow+semnow],'-','color',args.color,'LineWidth',args.linewidth)
    plot(ax,[x(iele)-args.capsize/2,x(iele)+args.capsize/2],[meannow-semnow,meannow-semnow],'-','color',args.color,'LineWidth',args.linewidth)
    plot(ax,[x(iele)-args.capsize/2,x(iele)+args.capsize/2],[meannow+semnow,meannow+semnow],'-','color',args.color,'LineWidth',args.linewidth)
    
end

end