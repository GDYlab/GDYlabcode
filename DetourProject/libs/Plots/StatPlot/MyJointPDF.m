function MyJointPDF(px,py,ax,datalabel,morperc)
%MyJointPDF plot 2d scatter plot and the marginal pdf for each dimension
% input   px py: data to plot, must be vectors and must match in length
%
%         ax: axis handle to plot, default be gca
%
%         datalabel: is the label of data, the size or dimension should
%         match px and py, you can also set labels for some groups as empty
%
%         morperc: plot mean or median for both direction, if set empty,
%         then won't plot these, empty as default
%
% example:  
% %          figure();a = rand(30,1)-1;b = rand(30,1)+2.5; MyJointPDF(a,b,gca,[],'mean')
        

%% check inputs
if nargin < 5
   morperc = []; 
end
if nargin < 4
    datalabel = []; 
end
if nargin < 3
    ax = gca; 
end
if ~isvector(px) || ~isvector(py) 
   error('Input data must be vector') 
end
if length(px) ~= length(py)
    error('Input vectors must match in length')
end
if ~isempty(datalabel) && length(px) ~= length(datalabel)
    error('length of datalabel should match length of data')
end

colors = Color_FrenchDisp;
% % colors = Color_Budapest;
%% do the plot
hold(ax,'on')

% find the range in x and y direction
rangex = [nanmin(px),nanmax(px)];
xbin = rangex(1):diff(rangex)/20:rangex(2);
rangey = [nanmin(py),nanmax(py)];
ybin = rangey(1):diff(rangey)/20:rangey(2);

[hx,hbinx] = hist(px,xbin);
hx = hx./max(hx)*nanmean(diff(ybin))*5;

[hy,hbiny] = hist(py,ybin);
hy = hy./max(hy)*nanmean(diff(xbin))*5;

xx = [berow(hy),zeros(size(berow(hy)))];
xx = xx + rangex(1) - nanmean(diff(xbin));
ytemp = berow(hbiny);
yy = [ytemp,fliplr(ytemp)];
fill(ax,xx, yy, 1-1*(1-colors{1}),'LineStyle','none','FaceAlpha',0.6);

xx = [berow(hx),zeros(size(berow(hx)))];
xx = xx + rangey(1) - nanmean(diff(ybin));
ytemp = berow(hbinx);
yy = [ytemp,fliplr(ytemp)];
fill(ax,yy, xx, 1-1*(1-colors{2}),'LineStyle','none','FaceAlpha',0.6);

hdata = plot(ax,px,py,'o','color',[0.4 0.4 0.4],'markersize',2,'linewidth',2);
if ~isempty(datalabel)
    customDataCursor(hdata,datalabel);
end

if ~isempty(morperc)
    switch morperc
        case 'mean'
            plot(ax,nanmean(px),nanmean(py),'+','color',colors{6},'markersize',30,'linewidth',4);
% %             [~,p1] = ttest(px);
% %             [~,p2] = ttest(py);
            
        case 'median'
            plot(ax,nanmedian(px),nanmedian(py),'+','color',colors{6},'markersize',30,'linewidth',4);
            
    end
end
xlim(ax,[rangex(1) - nanmean(diff(xbin)) , rangex(2)])
ylim(ax,[rangey(1) - nanmean(diff(ybin)) , rangey(2)])
setplot(ax)
% % set(gcf,'CurrentAxes',ax)
end

