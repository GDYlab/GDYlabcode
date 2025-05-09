function My2samplesCompare(d1,d2,ax,range,morperc,le,pointp,plb1,plb2)
%My2samplesCompare compare group of observations specified in d and plot in axis ax

% input   d1&d2: two data samples to compare, they must be a cell array with length equal or larger than 1,
%         each element is a vector of observations (matrix will be converted to vector)
%         their length must match
%
%         ax: axis handle to plot
%
%         range:  if specify, it's the range of pdf
%
%         morperc: plot mean+-std (1) or median+-25%75% (0), mean+-std as default
%        
%         le: legend of two groups
%         
%         pointp: plot data sample, default 0
%         plb: labels of data sample
% example:  
% %          a1 = rand(100,1);b1 = rand(20,20); c1 = rand(1,15);
% %          a2 = rand(100,1)+0.5;b2 = rand(20,20)+0.2; c2 = rand(1,15)+0.3;
% %          My2samplesCompare({a1,b1,c1},{a2,b2,c2},gca,[],0,{'g1','g2'})

%% check inputs
if nargin <9
   plb2 = [];
end

if nargin <8
   plb1 = []; 
   plb2 = [];
end

if nargin <7
   pointp = 0; 
end

if nargin <6
   le = []; 
end

if nargin <5
   morperc = 1; 
end

if nargin<4
    range = []; 
end

if ~iscell(d1) || ~iscell(d2)
   error('input data must be cell array with elements being groups to compare') 
end


if length(d1) ~= length(d2)
   error('data to compare must have the same length') 
end

if ~pointp
   plb1 = []; 
   plb2 = [];
end

colors = Color_FrenchDisp;
% % colors = Color_Budapest;
%% do the plot
hold(ax,'on')
for i = 1:length(d1)
    %% plot group 1
    ptemp = d1{i};
    ptemp = ptemp(:);
    if ~isempty(plb1)
        labtemp = plb1{i};
        labtemp = labtemp(:);
        if ~isempty(labtemp) && length(labtemp)~=length(ptemp)
            error('label doesnt match data, please check')
        end
        if ~isempty(labtemp)
            labtemp(isnan(ptemp)) = [];
        end
    end
    ptemp(isnan(ptemp)) = [];
    
      
    if length(ptemp) < 2
        if morperc
            % plot mean
            ptemp = berow(ptemp);
            hh(1) = plot(ax,i+0.03,nanmean(ptemp),'s','color',colors{1},'markersize',4,'linewidth',5);
        else
            % plot median
            ptemp = berow(ptemp);
            hh(1) = plot(ax,i+0.03,nanmedian(ptemp),'s','color',colors{1},'markersize',4,'linewidth',5);
        end
    end
    
    if length(ptemp) >= 2
        if morperc
            % plot mean and std
            ptemp = berow(ptemp);
            hh(1) = plot(ax,i+0.03,nanmean(ptemp),'s','color',colors{1},'markersize',4,'linewidth',4);
            plot(ax,[i+0.03 i+0.03],[nanmean(ptemp)-nanstd(ptemp),nanmean(ptemp)+nanstd(ptemp)],'-','color',colors{1},'linewidth',1.5);
        else
            % plot median and 25% 75%
            ptemp = berow(ptemp);
            [ptemp,~] = sort(ptemp,'ascend');
            hh(1) = plot(ax,i+0.03,nanmedian(ptemp),'s','color',colors{1},'markersize',4,'linewidth',4);
            plot(ax,[i+0.03 i+0.03],[ptemp(max(1,round(length(ptemp)*0.25))),ptemp(max(1,round(length(ptemp)*0.75)))],'-','color',colors{1},'linewidth',1.5);
        end
        if isempty(range)
            rangetemp = [min(ptemp),max(ptemp)];
        else
            rangetemp = range;
        end
        [ev,eb] = hist(ptemp,rangetemp(1):diff(rangetemp)/20:rangetemp(2));
        ev = ev./max(ev)*0.45; ev = berow(ev);
        xx = [ev,zeros(size(ev))];
        yy = [berow(eb),fliplr(berow(eb))];
        xx = xx + i + 0.07;
        try
           fill(ax,xx, yy, 1-0.6*(1-colors{1}),'LineStyle','none','FaceAlpha',0.7);
        end
    end
    
    if  pointp
        ptemp = berow(ptemp);
        px = rand(size(ptemp));
        px = (px-0.5)/20 + i+ 0.15;
        hdata = plot(ax,px,ptemp,'s','markerfacecolor',0.75*colors{1},...
            'markeredgecolor',0.75*colors{1},'markersize',4);
        if ~isempty(plb1) && ~isempty(labtemp)
            customDataCursor(hdata,labtemp);
        end
    end
    
    %% plot group 2
    ptemp = d2{i};
    ptemp = ptemp(:);
    if ~isempty(plb1)
        labtemp = plb2{i};
        labtemp = labtemp(:);
        if ~isempty(labtemp) && length(labtemp)~=length(ptemp)
            error('label doesnt match data, please check')
        end
        if ~isempty(labtemp)
            labtemp(isnan(ptemp)) = [];
        end
    end
    ptemp(isnan(ptemp)) = [];
    
    
    if length(ptemp) < 2
        if morperc
            % plot mean
            ptemp = berow(ptemp);
            hh(2) = plot(ax,i-0.03,nanmean(ptemp),'o','color',colors{2},'markersize',4,'linewidth',5);
        else
            % plot median
            ptemp = berow(ptemp);
            [ptemp,~] = sort(ptemp,'ascend');
            hh(2) = plot(ax,i-0.03,nanmedian(ptemp),'o','color',colors{2},'markersize',4,'linewidth',5);
        end
    end
    
    if length(ptemp) >= 2
        if morperc
            % plot mean and std
            ptemp = berow(ptemp);
            hh(2) = plot(ax,i-0.03,nanmean(ptemp),'o','color',colors{2},'markersize',4,'linewidth',4);
            plot(ax,[i-0.03 i-0.03],[nanmean(ptemp)-nanstd(ptemp),nanmean(ptemp)+nanstd(ptemp)],'-','color',colors{2},'linewidth',1.5);
        else
            % plot median and 25% 75%
            ptemp = berow(ptemp);
            [ptemp,~] = sort(ptemp,'ascend');
            hh(2) = plot(ax,i-0.03,nanmedian(ptemp),'o','color',colors{2},'markersize',4,'linewidth',4);
            plot(ax,[i-0.03 i-0.03],[ptemp(round(length(ptemp)*0.25)),ptemp(round(length(ptemp)*0.75))],'-','color',colors{2},'linewidth',1.5);
        end
        if isempty(range)
            rangetemp = [min(ptemp),max(ptemp)];
        else
            rangetemp = range;
        end
        [ev,eb] = hist(ptemp,rangetemp(1):diff(rangetemp)/20:rangetemp(2));
        ev = ev./max(ev)*0.45; ev = berow(ev);
        xx = [-ev,zeros(size(ev))];
        yy = [berow(eb),fliplr(berow(eb))];
        xx = xx + i - 0.07;
        try
            fill(ax,xx, yy, 1-0.6*(1-colors{2}),'LineStyle','none','FaceAlpha',0.7);
        end
    end
    
    if  pointp
        ptemp = berow(ptemp);
        px = rand(size(ptemp));
        px = (px-0.5)/20 + i - 0.15;
        hdata = plot(ax,px,ptemp,'o','markerfacecolor',0.75*colors{2},...
            'markeredgecolor',0.75*colors{2},'markersize',3);
        if ~isempty(plb2) && ~isempty(labtemp)
            customDataCursor(hdata,labtemp);
        end
    end
    
end

if ~isempty(le) && ~isempty(hh)
   legend(hh,le) 
end
xticks(ax,1:length(d1))
xlim(ax,[0 length(d1)+1])
setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

