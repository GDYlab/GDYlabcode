function MyBootstrapMean_v2(d,ax,baseline,varargin)
%MyBootstrapMean compare group of observations specified in d and plot in axis ax

% input   d: data to compare, it must be a cell array with length equal or larger than 1,
%         each element is a vector of observations (matrix will be converted to vector)
%
%         ax: axis handle to plot
%
%         baseline: the baseline data of bootstrap compare, set to empty if
%         don't want the dash line
%
%         range:  if specify, it's the range of pdf
%
%         its: bootstrap datasamples, 1e3 as default, if set to 0, then it
%         equals to the size of datasample
%
%         datalabel: is the label of data, the size or dimension should
%         match d, you can also set labels for some groups as empty
%
%         morperc: plot mean+-sem (1) or median+-25%75% (0), mean+-sem as default
%
% example:  
% %          a = rand(100,1);b = rand(20,20); c = rand(1,15);
% %          [~,pab] = ttest2(a,b(:)); [~,pac] = ttest2(a,c); [~,pbc] = ttest2(b(:),c);
% %          MyBootstrapMean({a,b,c},gca,1)
% %          sigstar({[1,2],[1,3],[2,3]},[pab,pac,pbc])

%% check inputs



args.markersize = 4; % arrow size
args.morperc = 1; 
args.datalabel = [];
args.its = 1e3;
args.range = [];
args.colors = Color_FrenchDisp; % color of the arrow
args = parseArgs(varargin, args);

colors = args.colors;
markersize = args.markersize;
morperc = args.morperc;
range = args.range;
its = args.its;
datalabel = args.datalabel;

if ~iscell(colors)
    tmpcolors = colors;
    clear colors
    colors = cell(1,2);
    colors{1} = tmpcolors;
    colors{2} = tmpcolors;
end


if ~iscell(d)
   error('input data must be cell array with elements being groups to compare') 
end

% % if length(d) < 2
% %    error('input data must have more than 2 groups') 
% % end

if baseline>length(d)
   error('Baseline dataset exceed the size of datagroups') 
end

%% do the plot
hold(ax,'on')

for i = 1:length(d)
    colori = rem(i,length(colors));
    colori(colori == 0) = length(colors);
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
    % if too less datapoints
    if length(ptemp) <= 2
        hdata = plot(ax,i*ones(size(ptemp)),ptemp,'o','color',colors{colori},'markersize',2,'linewidth',2);
        if ~isempty(datalabel) && ~isempty(labtemp)
            customDataCursor(hdata,labtemp);
        end
        continue
    end
    
    if morperc
    % plot mean and sem (standard error measure)
    ptemp = berow(ptemp);
    nancount = sum(~isnan(ptemp));
    plot(ax,i,nanmean(ptemp),'s','color',colors{colori},'markersize',markersize,'linewidth',markersize);
    plot(ax,[i i],[nanmean(ptemp)-nanstd(ptemp)./sqrt(nancount),nanmean(ptemp)+nanstd(ptemp)./sqrt(nancount)],'-','color',colors{colori},'linewidth',1.5);
    else
    % plot median and 25% 75%
    ptemp = berow(ptemp);
    [ptemp,~] = sort(ptemp,'ascend');
    plot(ax,i,nanmedian(ptemp),'s','color',colors{colori},'markersize',markersize,'linewidth',markersize);
    plot(ax,[i i],[ptemp(round(length(ptemp)*0.25)),ptemp(round(length(ptemp)*0.75))],'-','color',colors{colori},'linewidth',1.5);
    end
    
    if i == baseline
        plot(ax,[0 length(d)+1],[nanmean(ptemp) nanmean(ptemp)],'--k','linewidth',1.5)
    end
    
    
    % do the violin plot of raw data
    if isempty(range)
        rangetemp = [min(ptemp),max(ptemp)];
    else
        rangetemp = range;
    end
    if diff(rangetemp) == 0
        rangetemp(2) = rangetemp(2)+1;
        rangetemp(1) = rangetemp(1)-1;
    end
    [hv,hbin] = hist(ptemp,rangetemp(1):diff(rangetemp)/20:rangetemp(2));
    hv = hv./max(hv)*0.6;
% %     xx = [berow(hv),zeros(size(berow(hv)))];
% %     xx = i+xx+0.05;
% %     ytemp = berow(hbin);
% %     yy = [ytemp,fliplr(ytemp)];
% %     fill(xx, yy, 1-1*(1-colors{i}),'LineStyle','none','FaceAlpha',0.6);
    
    hbintemp = hbin;
    hbintemp(1) = -inf;
    hbintemp(end) = inf;
    cind = discretize(ptemp,hbintemp); % prepare the x offset scale
    hlen = hv(cind);
    hlen = berow(hlen);
    hscale = (hlen./max(hv)).^(1/2)*0.7+0.3;
    plotx = i-0.25 + (rand(size(ptemp))-0.5)/5.*hscale;
    hdata = plot(ax,plotx,ptemp,'o','color',colors{colori},'markersize',2,'linewidth',1);
    if ~isempty(datalabel) && ~isempty(labtemp)
        customDataCursor(hdata,labtemp);
    end
    
    % do the bootstrap, plot the mean distribution
    if its == 0
       its = length(ptemp);
    end
    mtemp = bootstrp(its,@mean,ptemp);
    
    [hm,hmbin] = hist(mtemp,min(mtemp):(max(mtemp)-min(mtemp))/20:max(mtemp));
    if max(hm)~=length(mtemp)
        hm = hm./max(hm)*0.5;
        xx = [berow(hm),zeros(size(berow(hm)))];
        xx = i+xx+0.1;
        ytemp = berow(hmbin);
        yy = [ytemp,fliplr(ytemp)];
        fill(ax,xx, yy, 1-1*(1-colors{colori}),'LineStyle','none','FaceAlpha',0.6);
    end
    
end
xticks(ax,1:length(d))
xlim(ax,[0 length(d)+1])
setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

