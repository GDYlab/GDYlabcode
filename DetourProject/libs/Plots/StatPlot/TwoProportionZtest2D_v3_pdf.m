function h = TwoProportionZtest2D_v3_pdf(s1,s2,m1,m2,ax,s1add,s2add,m3,cmap)
% this function give the 2D proportion difference plot (See Usman Science 2019)


% inputs:   s1 is one population with rows being realizations, columns
%           being two measures of interest
%           s2 is another population with rows being realizations, columns
%           being two measures of interest
%           m1 is the vector specifying the step size of m1 measure, its
%           the edge of bins
%           m2 is the vector specifying the step size of m2 measure
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
xxcenter = (xx(1:end-1) + xx(2:end))/2;
yycenter = (yy(1:end-1) + yy(2:end))/2;

ratdiffz = nan(length(yycenter),length(xxcenter));  % rows corresponding to measure 2, columns corresponding to measure1
alphC = zeros(length(yycenter),length(xxcenter)); % if there are data represented by pixel in first measure
xxgrid = nanmean(diff(xx));
yygrid = nanmean(diff(yy));
%% count the ratio with step criteria
for ir = 1:length(yycenter)
    for ic = 1:length(xxcenter)

        indr1 = s1(:,2) >= yy(ir) & s1(:,2) < yy(ir+1);

        indc1 = s1(:,1) >= xx(ic) & s1(:,1) < xx(ic+1);

        ind1 = indr1 & indc1;
        
        
        indr2 = s2(:,2) >= yy(ir) & s2(:,2) < yy(ir+1);
        
        indc2 = s2(:,1) >= xx(ic) & s2(:,1) < xx(ic+1);

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
h = imagesc([xxcenter(1) xxcenter(end)],[yycenter(1) yycenter(end)],ratdiffp);
set(h, 'AlphaData', alphC);
ax = gca;
ax.YDir = 'normal';        
shading flat
xlim(ax,[xx(1),xx(end)])
ylim(ax,[yy(1),yy(end)])
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