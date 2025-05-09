function  MySampleComparePlot(d,ax,range,datalabel,varargin)
%MySampleComparePlot compare group of observations specified in d and plot in axis ax

% input   d data to compare, it must be a cell with length larger than 1,
%         in each cell, it's a vector of observations (matrix will convert to vector)
%         ax axis handle to plot
%         range  if specify, its the range of pdf (edge of histcounts)
%         datalabel is the label of data, the size or dimension should
%         match d, you can also set labels for some subgroups as empty
% % colors = {[0.83 0.42 0.4],[0.52 0.45 0.83],[0.68 0.88 0.43],[0.5 0.73 0.93],...
% %     [0.75 0.42 0.92],[0.94 0.64 0.94],[0.74 0.74 0.22],[0.22 0.74 0.74],[0.2 0.32 0.62]}; % colors used for plotting
colors = Color_FrenchDisp_modify;
% % colors = Color_Budapest;

%% check inputs
if nargin <4
   datalabel = []; 
end

if nargin<3
    range = []; 
end

if ~iscell(d)
   error('input data must be cell with elements being groups to compare') 
end

args.markersize = 2; % arrow size
args.linewidth = 1.5; % linewidth
args.median75 = 1; % linewidth
args = parseArgs(varargin, args);

%%
set(gcf,'CurrentAxes',ax)
hold on
for i = 1:length(d)
    ptemp = d{i};
    ptemp = ptemp(:);
    if ~isempty(datalabel)
        labtemp = datalabel{i};
        labtemp = labtemp(:);
        if ~isempty(labtemp) && length(labtemp)~=length(ptemp)
            error('label doesnt match data, please check')
        end
        if ~isempty(labtemp)
            labtemp(isnan(ptemp)) = [];
        end
    end
    ptemp(isnan(ptemp)) = [];
    if isempty(ptemp)
        continue
    end
    if length(ptemp) < 3
        plotx = i-0.2 + (rand(size(ptemp))-0.5)/5.*0.2;
        hdata = plot(plotx,ptemp,'o','color',colors{i},'markersize',2,'linewidth',2);
        if ~isempty(datalabel) && ~isempty(labtemp)
            customDataCursor(hdata,labtemp);
        end
        continue
    end
    
    ptemp = berow(ptemp);
    sp = sort(ptemp,'ascend');
    p25 = sp(round(length(sp)*0.25));
    p75 = sp(round(length(sp)*0.75));
    if  args.median75
        plot(i,nanmedian(ptemp),'s','color',colors{i},'markersize',4,'linewidth',4);
        plot([i i],[p25 p75],'-','color',colors{i},'linewidth',args.linewidth);
    end
    if isempty(range)
        rangetemp = [min(ptemp),max(ptemp)];
        pdfedge = rangetemp(1):diff(rangetemp)/20:rangetemp(2);
    else
        pdfedge = range;
    end
    [hv,hedge] = histcounts(ptemp,pdfedge);
    hbin = (hedge(1:end-1) + hedge(2:end))/2;
    hv = hv./max(hv)*0.6;
    xx = [berow(hv),zeros(size(berow(hv)))];
    xx = i+xx+0.05;
    ytemp = berow(hbin);
    yy = [ytemp,fliplr(ytemp)];
    fill(xx, yy, 1-1*(1-colors{i}),'LineStyle','none','FaceAlpha',0.6);
    
    hbintemp = hedge;
    hbintemp(1) = -inf;
    hbintemp(end) = inf;
    cind = discretize(ptemp,hbintemp); % prepare the x offset scale
    hlen = hv(cind);
    hlen = berow(hlen);
    hscale = (hlen./max(hv)).^(1)*1.5+0.05;
    plotx = i-0.2 + (rand(size(ptemp))-0.5)/5.*hscale;
    hdata = plot(plotx,ptemp,'o','color',colors{i},...
        'markersize',args.markersize ,'linewidth',args.linewidth);
    if ~isempty(datalabel) && ~isempty(labtemp) 
        customDataCursor(hdata,labtemp);
    end
    
end
xticks(1:length(d))
xlim([0 length(d)+1])
setplot
set(gcf,'CurrentAxes',ax)
end

