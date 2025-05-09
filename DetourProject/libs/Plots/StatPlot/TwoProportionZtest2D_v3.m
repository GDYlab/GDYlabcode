function h = TwoProportionZtest2D_v3(s1,s2,m1,m2,m1dir,m2dir,ax,s1add,s2add,m3,cmap)
% this function give the 2D proportion difference plot (See Usman Science 2019)


% inputs:   s1 is one population with rows being realizations, columns
%           being two measures of interest
%           s2 is another population with rows being realizations, columns
%           being two measures of interest
%           m1 is the vector specifying the step size of m1 measure
%           m2 is the vector specifying the step size of m2 measure
%           m1dir and m2sdir specify the "direction" of the measure, if 1
%           then we will accumulate values larger than the current pixel,
%           if 0 we will accumulate values smaller than the current pixel
%
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
if nargin < 9
%     cmap=getNCLColmap('WhiteBlueGreenYellowRed.rgb',200);
    cmap=getNCLColmap('MPL_rainbow.rgb',200);

    cmap = flip(cmap,1);
end


%% check inputs
if size(s1,2) ~= 2 || size(s2,2) ~= 2
   error('Input population should have two columns corresponding to two measures of interest') 
end
% if there is additional criteria
addcri = 0;

if nargin > 7
   addcri = 1;
end

%% create bin for two measures

xx = sort(m1,'ascend');
yy = sort(m2,'ascend');
ratdiffz = nan(length(yy),length(xx));  % rows corresponding to measure 2, columns corresponding to measure1
alphC = zeros(length(yy),length(xx)); % if there are data represented by pixel in first measure
xxgrid = nanmean(diff(xx));
yygrid = nanmean(diff(yy));
%% count the ratio with step criteria
for ir = 1:length(yy)
    for ic = 1:length(xx)
        if m2dir == 1
            indr1 = s1(:,2) >= yy(ir);
        else
            indr1 = s1(:,2) <= yy(ir);
        end
        
        if m1dir == 1
            indc1 = s1(:,1) >= xx(ic);
        else
            indc1 = s1(:,1) <= xx(ic);
        end
        ind1 = indr1 & indc1;
        
        
        if m2dir == 1
            indr2 = s2(:,2) >= yy(ir);
        else
            indr2 = s2(:,2) <= yy(ir);
        end
        
        if m1dir == 1
            indc2 = s2(:,1) >= xx(ic);
        else
            indc2 = s2(:,1) <= xx(ic);
        end
        ind2 = indr2 & indc2;
        
        if addcri
            add1 = s1add >= m3;
            add2 = s2add >= m3;
            ind1 = ind1 & add1;
            ind2 = ind2 & add2;
        end
        if sum(ind1) > 0
            alphC(ir,ic) = 1;
        end
        
        
        ratdiffz(ir,ic) = myztest(sum(ind1),sum(ind2),size(s1,1),size(s2,1));
    end
end
% convert zscore to one tail p value see if S1 significant larger than S2
ratdiffp = 1-normcdf(ratdiffz);
%% plot the p value
set(gcf,'CurrentAxes',ax)
hold on
h = imagesc([xx(1) xx(end)],[yy(1) yy(end)],ratdiffp);
set(h, 'AlphaData', alphC);
ax = gca;
ax.YDir = 'normal';        
shading flat
xlim(ax,[xx(1)-xxgrid/2,xx(end)+xxgrid/2])
ylim(ax,[yy(1)-yygrid/2,yy(end)+yygrid/2])
colormap(ax,cmap)         
caxis([0 0.05])
box on
set(ax,'layer','top')
% setplot
axis square
axis tight
% cb = addcolorbar(ax,0.0005);
% cb.Label.String = 'p';
end