function Plot2SetDisjointVenn(ax,set1,set2,overlap,allele,varargin)
% this function will plot the venn plot for 2 sets of disjoint sets to
% visualize their overlap
% inputs:ax,    axis to make the plot
%        set 1, vector with length of 2, number of elements set1a and set1b
%        set 2, vector with length of 2, number of elements set2a and set2b
%        note that set1a, set1b must be disjoint
%        note that set2a, set2b must be disjoint
%        overlap, a 2*2 matrix with elements being number of intersect of 
%        [set1a^set2a   set1b^set2a]
%        [set1a^set2b   set1b^set2b]
%        allele, total number of elements

args.set1color = {[0.8 0.42 0.2],[0.2 0.42 0.8]}; % color of the set1
args.set2color = {[0.76,0.76,0.3],[0.3,0.76,0.76]}; % color of the set2
args.alpha = 0.5; % color of the arrow
args.lgdposition = 'northoutside';
args.ldgstr = {'1a','1b','2a','2b'};
args = parseArgs(varargin, args);

rset1color{1} = [args.set1color{1},args.alpha];
rset1color{2} = [args.set1color{2},args.alpha];
rset2color{1} = [args.set2color{1},args.alpha];
rset2color{2} = [args.set2color{2},args.alpha];
%% normalize based on allele
set1 = set1./(allele);
set2 = set2./(allele);
overlap = overlap./(allele);

%% plot set1a and set1b as rectangles
hold(ax,'on')
rectangle(ax,'Position',[0, 0, set1(1), 1],'Curvature',0.02,'linestyle','none','facecolor',rset1color{1})
rectangle(ax,'Position',[1-set1(2), 0, set1(2), 1],'Curvature',0.02,'linestyle','none','facecolor',rset1color{2})
%% plot set2a and set2b based on overlap
rectangle(ax,'Position',[0, 0, set1(1), overlap(1,1)/set1(1)],'Curvature',0.02,'linestyle','none','facecolor',rset2color{1})
rectangle(ax,'Position',[1-set1(2), 0, set1(2), overlap(1,2)/set1(2)],'Curvature',0.02,'linestyle','none','facecolor',rset2color{1})
rectangle(ax,'Position',[set1(1), 0, 1-set1(2)-set1(1), (set2(1)-overlap(1,1)-overlap(1,2))/(1-set1(2)-set1(1))],...
    'Curvature',0.02,'linestyle','none','facecolor',rset2color{1})

rectangle(ax,'Position',[0, 1-overlap(2,1)/set1(1), set1(1), overlap(2,1)/set1(1)],'Curvature',0.02,'linestyle','none','facecolor',rset2color{2})
rectangle(ax,'Position',[1-set1(2), 1-overlap(2,2)/set1(2), set1(2), overlap(2,2)/set1(2)],'Curvature',0.02,'linestyle','none','facecolor',rset2color{2})
rectangle(ax,'Position',[set1(1), 1-(set2(2)-overlap(2,1)-overlap(2,2))/(1-set1(2)-set1(1)), 1-set1(2)-set1(1), (set2(2)-overlap(2,1)-overlap(2,2))/(1-set1(2)-set1(1))],...
    'Curvature',0.02,'linestyle','none','facecolor',rset2color{2})

hd(1) = scatter(ax,nan,nan,'s','markeredgecolor',args.set1color{1},'markerfacecolor',args.set1color{1},'MarkerFaceAlpha',args.alpha);
hd(2) = scatter(ax,nan,nan,'s','markeredgecolor',args.set1color{2},'markerfacecolor',args.set1color{2},'MarkerFaceAlpha',args.alpha);
hd(3) = scatter(ax,nan,nan,'s','markeredgecolor',args.set2color{1},'markerfacecolor',args.set2color{1},'MarkerFaceAlpha',args.alpha);
hd(4) = scatter(ax,nan,nan,'s','markeredgecolor',args.set2color{2},'markerfacecolor',args.set2color{2},'MarkerFaceAlpha',args.alpha);
legend(hd,args.ldgstr,'Location',args.lgdposition)
legend boxoff
axis square
axis off
end