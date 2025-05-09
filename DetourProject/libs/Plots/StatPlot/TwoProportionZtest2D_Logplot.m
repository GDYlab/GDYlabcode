function h = TwoProportionZtest2D_Logplot(s1,s2,m1,m2,ax,s1add,s2add,m3)
% this function give the 2D proportion difference plot (See Usman Science 2019)


% inputs:   s1 is one population with rows being realizations, columns
%           being two measures of interest
%           s2 is another population with rows being realizations, columns
%           being two measures of interest
%           m1 is the vector specifying the step size of m1 measure
%           m2 is the vector specifying the step size of m2 measure
%           ax is the axis to plot
%           s3 and m3 is the additional criteria that the data is supposed
%           to meet. This is usually not necessary. One example of using
%           this is in the decoding, in addition to the linefit and jump
%           criteria, we require the spatial coverage to be larger than 50
%
%           note that in m1 and m2, the large value means 'good' samples,
%           and the small value means trivial samples, thus we will always
%           count the ratio larger than the criteria, you should adjust
%           measures in s1 and s2 accrodingly
%           note that we will test if rate in s1 is significant higher than s2
% output    h is the handle of the plot
cmap=getNCLColmap('MPL_rainbow.rgb',300);
cmap = flip(cmap,1);
crange = [-10,0];
lowlimit = 10^(crange(1));
plog = log10(0.05);
plogidx = round(size(cmap,1)/diff(crange)*(plog-crange(1)));
cmap(plogidx:end,:) = ones(size(cmap,1)-plogidx+1,3);
%% check inputs
if size(s1,2) ~= 2 || size(s2,2) ~= 2
   error('Input population should have two columns corresponding to two measures of interest') 
end
% if there is additional criteria
addcri = 0;

if nargin > 5
   addcri = 1;
end
%% create bin for two measures

xx = sort(m1,'ascend');
yy = sort(m2,'ascend');
ratdiffz = nan(length(yy),length(xx));  % rows corresponding to measure 2, columns corresponding to measure1
xxgrid = nanmean(diff(xx));
yygrid = nanmean(diff(yy));
%% count the ratio with step criteria
for ir = 1:length(yy)
    for ic = 1:length(xx)
        indr1 = s1(:,2) >= yy(ir);
        indc1 = s1(:,1) >= xx(ic);
        ind1 = indr1 & indc1;
        
        indr2 = s2(:,2) >= yy(ir);
        indc2 = s2(:,1) >= xx(ic);
        ind2 = indr2 & indc2;
        
        if addcri
            add1 = s1add >= m3;
            add2 = s2add >= m3;
            ind1 = ind1 & add1;
            ind2 = ind2 & add2;
        end
        
        ratdiffz(ir,ic) = myztest(sum(ind1),sum(ind2),size(s1,1),size(s2,1));
    end
end
% convert zscore to one tail p value see if S1 significant larger than S2
ratdiffp = 1-normcdf(ratdiffz);
%% plot the p value
set(gcf,'CurrentAxes',ax)
hold on
ratdiffp(ratdiffp<=lowlimit) = lowlimit;
ratdiffp = log10(ratdiffp);
h = imagesc([xx(1)+xxgrid/2 xx(end)+xxgrid/2],[yy(1)+yygrid/2 yy(end)+yygrid/2],ratdiffp);
ax = gca;
ax.YDir = 'normal';        
shading flat
colormap(cmap)         
caxis(crange)
setplot
axis tight
cb = addcolorbar;
cb.Label.String = 'log10(p)';
ticknow = cb.Ticks;
ticknow = cat(2,ticknow,plog);
[ticknow,sind] = sort(ticknow,'ascend');
cb.Ticks = ticknow;
cb.TickLabels{sind(end)} = 'p=0.05';

end