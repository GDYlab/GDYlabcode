function  [hd,kldiv,symkldiv] = My2PDFComparePlot(d,ax,range,normalize,lwid,colors)
%My2PDFComparePlot  plot two PDFs to compare

% input   d data to compare, it must be a cell array with length of 2,
%         in each cell, it's a vector of observations (matrix will convert to vector)
%         ax axis handle to plot
%         range  if specify, its the bin edges to get pdf
% output  hd  is the handle with length of 2 matching sum pdf1*log(pdf1/pdf2)
%         kldiv  is the kl divergence of 2 pdfs, as 
% example hd = My2PDFComparePlot({rand(300,1),0.3+rand(100,1)/2},gca,[0:0.1:1])

% % colors = {[0.83 0.42 0.4],[0.52 0.45 0.83],[0.68 0.88 0.43],[0.5 0.73 0.93],...
% %     [0.75 0.42 0.92],[0.94 0.64 0.94],[0.74 0.74 0.22],[0.22 0.74 0.74],[0.2 0.32 0.62]}; % colors used for plotting
% % colors = Color_Budapest;

%% check inputs
if nargin < 6
   colors = Color_FrenchDisp; 
end

if nargin<5
    lwid = 5; 
end

if nargin<4
    normalize = 'pdf'; 
end

if nargin<3
    range = []; 
end

if ~iscell(d) 
   error('input data must be cell array with length of 2 and elements being groups to compare') 
end

if length(d)~=2 
   error('input data must be cell array with length of 2 and elements being groups to compare') 
end

% if range is not specified, we will use 20 bins and cover all the datasets
if isempty(range)
    mind = nan; maxd = nan;
    for i = 1:length(d)
        dtemp = d{i};
        dtemp(isinf(dtemp)) = nan;
        mind = nanmin([mind,nanmin(dtemp(:))]);
        maxd = nanmax([maxd,nanmax(dtemp(:))]);
        d{i} = dtemp;
    end
    dx = (maxd-mind)/20;
    range = mind:dx:maxd;
end
rangec = (range(1:end-1)+range(2:end))/2;
%% make plot
hold(ax,'on')
d1 = d{1}; d1 = d1(:);
d2 = d{2}; d2 = d2(:);

% % [h1,~] = histcounts(d1,range);
% % h1 = h1./sum(h1)./dx; % normalize counts to pdf
% % [h2,~] = histcounts(d2,range);

% hd(1) = bar(rangec,h1,1,'LineStyle','none','Facecolor',colors{1},'FaceAlpha',0.8);
hd(1) = histogram(ax,d1,range,'LineStyle','none','Facecolor',colors{1},'Normalization',normalize,'FaceAlpha',0.8);
hd(2) = histogram(ax,d2,range,'DisplayStyle','stairs','linewidth',lwid,'Edgecolor',colors{2},'Normalization',normalize);
% % histogram(X,edges)
% % 'DisplayStyle','stairs'
% % setplot
ylabel(ax,'PDF')

hc1 = histcounts(d1,range,'Normalization','pdf');
hc2 = histcounts(d2,range,'Normalization','pdf');
hc1 = hc1 + eps;
hc2 = hc2 + eps;
kldiv = sum(hc1.*log(hc1./hc2));
kldiv2 = sum(hc2.*log(hc2./hc1));
% l2 = norm(hc1-hc2);

symkldiv = (kldiv + kldiv2)/2;
end

