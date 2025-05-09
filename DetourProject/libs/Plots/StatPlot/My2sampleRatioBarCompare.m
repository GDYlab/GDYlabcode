function My2sampleRatioBarCompare(d1,d2,ax,color1,color2,pztest,d1num,d2num)
%My2sampleRatioBarCompare compare two group of ratios specified in d and plot in axis ax

% input   d1&d2: two ratios to compare, they must have the same length 
%
%         ax: axis handle to plot
%
%         color1 & color 2:  color for each bar plot, if specified, they
%         should be cell array with length equal to data, set to 0 will
%         use default color
%
%         pztest: whether do two proportional z test, if 1, total number of
%         trails must be specified in d1num and d2num
%        
% example:  
%           a1 = rand(4,1);b1 = rand(4,1);
% 
%           My2sampleRatioBarCompare(a1,b1,gca,0,0,1,20,40)

%% check inputs
if nargin < 4
    color1 = 0;
end
if nargin < 5
    color2 = 0;
end
if nargin < 6
    pztest = 0;
end

if pztest == 1 && nargin ~= 8
   error('Must specify number of trails for two proportional z test') 
end

if length(d1) ~= length(d2)
   error('Input ratios must have same length ')
end

if pztest == 1 && (length(d1) ~= length(d1num) || length(d2) ~= length(d2num))
   error('Input ratios must match number of trails')
end
    
if color1 == 0
    color1 = cell(1,length(d1));
    for ic = 1:length(d1)
        color1{ic} = [0.83 0.42 0.4];
    end
end

if color2 == 0
    color2 = cell(1,length(d2));
    for ic = 1:length(d2)
        color2{ic} = [0.6 0.6 0.6];
    end
end

if length(color1) ~= length(d1) || length(color2) ~= length(d2)
   error('bar color must match input ratio')
end
%% do the plot
axes(ax)
hold on
px1 = 1:length(d1);
px2 = 1:length(d1);
px1 = px1 + 0.2;
px2 = px2 - 0.2;
for is = 1:length(d1)
    bar(px1(is),d1(is)*100,0.4,'FaceColor',color1{is},'EdgeColor',color1{is})
    bar(px2(is),d2(is)*100,0.4,'FaceColor',color2{is},'EdgeColor',color2{is})
    maxr = max([d1(is)*100,d2(is)*100]);
    
    if pztest
        ratdiffz = myztest(d1(is)*d1num(is),d2(is)*d2num(is),d1num(is),d2num(is));
        ratdiffpr = 1-normcdf(ratdiffz);
        ratdiffpl = normcdf(ratdiffz);
        psmall = min([ratdiffpr, ratdiffpl]);
        text(px2(is)-0.3,maxr*1.2,['p',num2str(psmall,'%.3f')],'Fontsize',12)
    end
end
setplot
end