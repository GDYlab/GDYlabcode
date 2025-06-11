%% Fig 4a
% Change to appropriate directory using files
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('TrackProbLearnEffectOnSleepFramesHPCRipp.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:4
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[1,2,4,5],[21,22,24,25],[6,7,9,10,11,12,14,15],[16,17,19,20]};
for Days=1:length(Sess)
for sess=1:length(Sess{Days})
    for jj=1:size(Probability1{Days,sess},2)
     x=(([(Probability1{Days,sess}{2,jj}(:,1)), (Probability1{Days,sess}{2,jj}(:,2))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,1)), (Probability1{Days,sess}{1,jj}(:,2))]'))';
     indices{1,Days}{jj}=[indices{1,Days}{jj}; (x-y)'];
     indices1{1,Days}{jj}=[indices1{1,Days}{jj}; (x)']; 
     indices1{2,Days}{jj}= [indices1{2,Days}{jj}; (y)'];         
     
end
end

end

col1={rgb('Black'),rgb('Gray'),rgb('Red'),...
    rgb('Blue'),rgb('Gray'),rgb('LightBlue')}
fig1=figure(1)
left=0.1; bottom=0.75; width=0.20; height=0.20;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'PreBa-1','PreNo-cue','PreBa-2'};
labels2={'PreBa-1','PreAsn','PreAsd'};
P=[];
for Days=1:4
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:3
    ddt1=[ddt1, movmean(indices{1,Days}{jj},[1,1],2)];
    ddt2=[ddt2, indices1{1,Days}{jj}];
    ddt3=[ddt3, indices1{2,Days}{jj}];
    end
    ddt_mean=nanmean(ddt1);
    ddt_sem=nansem(ddt1);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
    if Days==4
    plot([0 30],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[5 15 25],'xticklabel',{'PreBa-1','PreAsn/No-cue','PreAsd/Ba-2'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Ba-2/Asd - Ba1)','fontsize',5)
end

left=0.1; bottom=0.55; width=0.20; height=0.12;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};P=[];
for Days=1:4
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:3
    ddt1=[ddt1, indices{1,Days}{jj}];
    ddt2=[ddt2, indices1{1,Days}{jj}];
    ddt3=[ddt3, indices1{2,Days}{jj}];
    ddt_mean{Days}{jj}=nanmean(indices{1,Days}{jj}')';
    for ww=1:size(ddt_mean{Days}{jj}(:),1)
        x1=ddt_mean{Days}{jj}(:); hold on;
        plot(xpos(Days)+jj+(randn(1).*0.1),x1(ww),'marker','o','markersize',1,'color',col1{Days});
    end
    hh1(Days)=errorbar(xpos(Days)+jj,nanmean(nanmean(ddt_mean{Days}{jj})),nanmean(nansem(ddt_mean{Days}{jj})),'linestyle','-','marker','o','markersize',2,'linewidth',1,'color',col1{Days},'capsize',3); hold on;    
    end

    Groups=[1,2;1,3;2,3];h=[];pval=[];
    for ww1=1:size(Groups,1)
        [h pval(ww1) ci stats]=ttest(ddt_mean{Days}{Groups(ww1,1)},ddt_mean{Days}{Groups(ww1,2)});
        if pval(ww1) < 0.05
           hh2=sigstar(Groups(ww1,:)+xpos(Days),pval(ww1),0);
        end
        P{Days,2}{ww1,2}=pval(ww1);P{Days,2}{ww1,1}=Groups(ww1,:);
        P{Days,2}{ww1,3}=stats.df; P{Days,2}{ww1,4}=stats.tstat;
    end
    
    
    yl(Days,:)=get(gca,'ylim');
    Groups=[1:3];h=[];pval=[];
    for ww1=1:size(Groups,2)
        [h pval(ww1) ci stats]=ttest(ddt_mean{Days}{Groups(ww1)});
        P{Days,1}{ww1,2}=pval(ww1);P{Days,1}{ww1,1}=Groups(ww1);
        P{Days,1}{ww1,3}=stats.df; P{Days,1}{ww1,4}=stats.tstat;
        P{Days,1}{ww1,5}=nanmean(ddt_mean{Days}{Groups(ww1)});
        P{Days,1}{ww1,6}=nansem(ddt_mean{Days}{Groups(ww1)});
        if pval(ww1) < 0.05 & pval(ww1) > 0.01
          text(Groups(ww1)+xpos(Days),yl(Days,end)+0.001,siglabels{1},'fontsize',4,'color',col1{Days});
        elseif  pval(ww1) < 0.01 & pval(ww1) > 0.001
          text(Groups(ww1)+xpos(Days),yl(Days,end)+0.001,siglabels{2},'fontsize',4,'color',col1{Days});
        elseif  pval(ww1) < 0.001 
          text(Groups(ww1)+xpos(Days),yl(Days,end)+0.001,siglabels{3},'fontsize',4,'color',col1{Days});
        end
    end
    if Days==4
    plot([0 30],[0 0],'--k');
    xlim([0 xpos(Days)+4]);
    ylim([nanmin(yl(:,1))-0.001 nanmax(yl(:,2))+0.001])
    end
   
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[1 2 3 5 6 7 9 10 11 13 14 15],'xticklabel',[labels1,labels1,labels2,labels2],'xticklabelrotation',45);
    
    
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Ba-2/Asd - Ba1)','fontsize',5)
end

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300',strcat('TrackProbLearningEffectsOnSleepFramesRipp'));
print(fig1,'-painters','-depsc','-r300',strcat('TrackProbLearningEffectsOnSleepFramesRipp')); 
close all;

% Normalizing for baseline experience

% Change to appropriate directory using files
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('TrackProbLearnEffectOnSleepFramesHPCRipp.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:4
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[1,2,4,5],[21,22,24,25],[6,7,9,10,11,12,14,15],[16,17,19,20]};
for Days=1:length(Sess)
for sess=1:length(Sess{Days})
    for jj=2:size(Probability1{Days,sess},2)
     x=(([ (Probability1{Days,sess}{2,jj}(:,1)), (Probability1{Days,sess}{2,jj}(:,2))]'))'; 
     y=(([ (Probability1{Days,sess}{1,jj}(:,1)), (Probability1{Days,sess}{1,jj}(:,2))]'))';
     x1=(([(Probability1{Days,sess}{2,jj-1}(:,1)), (Probability1{Days,sess}{2,jj-1}(:,2))]'))'; 
     y1=(([(Probability1{Days,sess}{1,jj-1}(:,1)), (Probability1{Days,sess}{1,jj-1}(:,2))]'))';
     indices1{1,Days}{jj-1}=[indices1{1,Days}{jj-1}; (x-x1)']; 
     indices1{2,Days}{jj-1}= [indices1{2,Days}{jj-1}; (y-y1)'];         
          
end
end

end

col1={rgb('Black'),rgb('Gray'),rgb('Red'),...
    rgb('Blue'),rgb('Gray'),rgb('LightBlue')}
fig1=figure(1)
left=0.1; bottom=0.65; width=0.20; height=0.12;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];P=[];
xpos=[0,3,6,9];ddt_mean=[];siglabels={'*','**','***'};
for Days=1:4
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:2
    if jj==1
    ddt_mean{Days}{jj}=nanmean([indices1{2,Days}{2}-indices1{2,Days}{1}]') ;
    else
    ddt_mean{Days}{jj}=nanmean([indices1{1,Days}{2}-indices1{1,Days}{1}]') ;    
    end
    for ww=1:size(ddt_mean{Days}{jj}(:),1)
        x1=ddt_mean{Days}{jj}(:); hold on;
        plot(xpos(Days)+jj+(randn(1).*0.1),x1(ww),'marker','o','markersize',1,'color',col1{Days});
    end
    hh1(Days)=errorbar(xpos(Days)+jj,nanmean(nanmean(ddt_mean{Days}{jj})),nanmean(nansem(ddt_mean{Days}{jj})),'linestyle','-','marker','o','markersize',2,'linewidth',0.5,'color',col1{Days},'capsize',3); hold on;    
    end

    Groups=[1,2];h=[];pval=[];
    for ww1=1:size(Groups,1)
        [h pval(ww1) ci stats]=ttest(ddt_mean{Days}{Groups(ww1,1)}(:),ddt_mean{Days}{Groups(ww1,2)}(:));
        P{Days,2}{ww1,2}=pval(ww1);P{Days,2}{ww1,1}=Groups(ww1,:);
        P{Days,2}{ww1,3}=stats.df; P{Days,2}{ww1,4}=stats.tstat;
        if pval(ww1) < 0.05
           hh2=sigstar(Groups(ww1,:)+xpos(Days),pval(ww1),0);
        end
    end
    yl(Days,:)=get(gca,'ylim');
    Groups=[1:2];h=[];pval=[];
    for ww1=1:size(Groups,2)
        [h pval(ww1) ci stats]=ttest(ddt_mean{Days}{Groups(ww1)}(:));
        P{Days,1}{ww1,2}=pval(ww1);P{Days,1}{ww1,1}=Groups(ww1);
        P{Days,1}{ww1,3}=stats.df; P{Days,1}{ww1,4}=stats.tstat;
        P{Days,1}{ww1,5}=nanmean(ddt_mean{Days}{Groups(ww1)});
        P{Days,1}{ww1,6}=nansem(ddt_mean{Days}{Groups(ww1)});
        if pval(ww1) < 0.05 & pval(ww1) > 0.01
          text(Groups(ww1)+xpos(Days),yl(Days,end)+0.001,siglabels{1},'fontsize',4,'color',col1{Days});
        elseif  pval(ww1) < 0.01 & pval(ww1) > 0.001
          text(Groups(ww1)+xpos(Days),yl(Days,end)+0.001,siglabels{2},'fontsize',4,'color',col1{Days});
        elseif  pval(ww1) < 0.001 
          text(Groups(ww1)+xpos(Days),yl(Days,end)+0.001,siglabels{3},'fontsize',4,'color',col1{Days});
        end
    end
    if Days==4
    plot([0 30],[0 0],'--k');
    xlim([0 xpos(Days)+4]);
    ylim([nanmin(yl(:,1))-0.01 nanmax(yl(:,2))+0.01])
    end
   
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[1 2 4 5 7 8 10 11],'xticklabel',{'Ba-1','Asd'},'xticklabelrotation',45);
    
    
    xlabel('Awake sessions','fontsize',5);
    ylabel('\Delta prob Normalized','fontsize',5)
end

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('TrackProbLearningEffectsIncreaseNormalizedOnSleepFramesRippV2'));
print(fig1,'-painters','-depsc','-r300',strcat('TrackProbLearningEffectsIncreaseNormalizedOnSleepFramesRippV2')); 
close all;

%% Fig 4b
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('TrackProbLearnEffectOnSleepFramesHPCNewVsOldRipp.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,7,9,10,11,12,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=1:size(Probability1{Days,sess},2)
     x=(nanmean([(Probability1{Days,sess}{1,jj}(:,1)), (Probability1{Days,sess}{1,jj}(:,2))]'))'; 
     y=(nanmean([(Probability1{Days,sess}{1,jj}(:,3)), (Probability1{Days,sess}{1,jj}(:,4))]'))';
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
     indices1{1,1}{jj}=[indices1{1,1}{jj}; (x)']; 
     indices1{2,1}{jj}= [indices1{2,1}{jj}; (y)']; 
     x=(nanmean([(Probability1{Days,sess}{3,jj}(:,1)), (Probability1{Days,sess}{3,jj}(:,2))]'))'; 
     y=(nanmean([(Probability1{Days,sess}{3,jj}(:,3)), (Probability1{Days,sess}{3,jj}(:,4))]'))';
     indices{1,2}{jj}=[indices{1,2}{jj}; ((x-y))'];
     indices1{2,1}{jj}=[indices1{1,2}{jj}; (x)']; 
     indices1{2,2}{jj}= [indices1{2,2}{jj}; (y)'];   
     
end
end

end

col1={rgb('Red'),...
    rgb('Blue'),rgb('Green'),rgb('LightBlue')}
fig1=figure(1)
left=0.1; bottom=0.75; width=0.20; height=0.20;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'PreBa-1','PreNo-cue','PreBa-2'};
labels2={'PreBa-1','PreAsn','PreAsd'};
for Days=1:2
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:3
    ddt1=[ddt1, (indices{1,Days}{jj})];
    ddt2=[ddt2, indices1{1,Days}{jj}];
    ddt3=[ddt3, indices1{2,Days}{jj}];
    end
    ddt_mean=movmean(nanmean(ddt1),[1 1],2);
    ddt_sem=movmean(nansem(ddt1),[1 1],2);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
    if Days==2
    plot([0 30],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[5 15 25],'xticklabel',{'PreBa-1','PreAsn','PreAsd'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Changed-Unchanged)','fontsize',5)
end


pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};


pval=[]; P=[];   
  for jj=1:3
    ddt1=[];ddt2=[];ddt3=[];    
left=0.35; bottom=0.9-(jj-1).*0.08; width=0.08; height=0.05;
ax1(jj)=axes('Position',[left, bottom, width, height])      
    for Days=1:2    
    ddt_mean{Days}{jj}=(indices{1,Days}{jj}(:));
    [h1 stats]=ecdf(ddt_mean{Days}{jj});hold on
    plot(stats,h1,'color',col1{Days});
    plot([nanmedian(ddt_mean{Days}{jj}) nanmedian(ddt_mean{Days}{jj})],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{Days});
    
    
    end
    plot([0 0],[0 1],'--k');
     [pval(jj) h stats]=signrank(ddt_mean{1}{jj},ddt_mean{2}{jj});
     P{1,3}{jj,2}=pval(jj); P{1,3}{jj,3}=stats.zval; 
        if pval(jj) < 0.05 & pval(jj) > 0.01
          text(-0.02,0.8,siglabels{1},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.01 & pval(jj) > 0.001
          text(-0.02,0.8,siglabels{2},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.001 
          text(-0.02,0.8,siglabels{3},'fontsize',6,'color','k');
        end
       
    P{1,1}{jj,1}=1;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj});  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj});
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}); 
    
        [pval1(jj) h stats]=signrank(ddt_mean{1}{jj});  
        P{1,1}{jj,2}=pval1(jj);   P{1,1}{jj,3}=stats.zval;
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95,siglabels{1},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95,siglabels{2},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95,siglabels{3},'fontsize',6,'color',col1{1});
        end
    P{1,2}{jj,1}=2;
    P{1,2}{jj,4}=length(ddt_mean{2}{jj});  P{1,2}{jj,5}=nanmean(ddt_mean{2}{jj});
    P{1,2}{jj,6}=nansem(ddt_mean{2}{jj});   
        
        
        
         [pval2(jj) h stats]=signrank(ddt_mean{2}{jj});
         P{1,2}{jj,2}=pval2(jj);   P{1,2}{jj,3}=stats.zval;
         if pval2(jj) < 0.05 & pval2(jj) > 0.01
          text(-0.01,0.9,siglabels{1},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.01 & pval2(jj) > 0.001
          text(-0.01,0.9,siglabels{2},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.001 
          text(-0.01,0.9,siglabels{3},'fontsize',6,'color',col1{2});
        end    
    xl=get(gca,'xtick');
    xl1=[nanmin(xl) nanmax(xl)];
    set(gca,'xtick',xl1,'xticklabel',xl1,'fontsize',5);    
    set(gca,'box','off','tickdir','out','fontsize',5);
    if jj==3
    xlabel('Frame prob (Changed-Unchanged)','fontsize',5)
    end
    if jj==1
    ylabel('cdf','fontsize',5);    
    end
end


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('TrackProbLearningEffectsOnSleepFramesNewVsOldRipp'));
print(fig1,'-painters','-depsc','-r300',strcat('TrackProbLearningEffectsOnSleepFramesNewVsOldRipp')); 
close all;

%% Fig 4c

% Change to appropriate directory using files
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SleepFramesCorrMat10binsHPC');
load('IndicesPlaceArmsHPC')
load('IndicesFlavorHPC')
load('IndicesOdorHPC')

sd=2;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);gg=0;
Coact=[];MembersZcontOld=[];
rr=0;
for ii=[6:15]
   
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
 
 [c1 c2] = textread('celltype','%s %s');
elecreg=importdata('elecreg');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
ind=a12(:,1)==1 & a12(:,2)==1 & (elecreg(:,2)==1 | elecreg(:,2)==2);
ind_all=elecreg;
% ind_all=elecreg(ind_HPC | ind_PFC,2);
 
  
    tt=cell2mat(cellfun(@(x)(~isempty(x)),Indices_Place(ii,:),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    gg=0;rr=rr+1;
    if any(ii==[1:5,21:25,13,18,8]) ii==6 | ii==15; dd=0; else dd=1;end
    
   
    MembersIndices1=[];comm1=[];ll=0;
    MembersIndices2=[];comm2=[];
    MembersIndices3=[];comm3=[];
    MembersIndices4=[];comm4=[];mm1=0;
     for jj=[1,2,length(tt1)-dd,length(tt1)]   
      ll=ll+1; 
     MembersIndices1{ll}=(squeeze(Indices_Place{ii,jj}.Pvalues1(:,:,1,1)))<0.05;
     comm1{ll}=~(any(MembersIndices1{ll}(:,[2,3,5])')' & any(MembersIndices1{ll}(:,[1,4,6])')');
     
     MembersIndices2{ll}=(squeeze(Indices_Place{ii,jj}.Pvalues1(:,:,2,1)))<0.05;
     comm2{ll}=~(any(MembersIndices2{ll}(:,[2,3,5])')' & any(MembersIndices2{ll}(:,[1,4,6])')');
     
     MembersIndices3{ll}=(squeeze(Indices_Flavor{ii,jj}.Pvalues(:,:,1)))<0.05;
     comm3{ll}=~(any(MembersIndices3{ll}(:,[2,3,5])')' & any(MembersIndices3{ll}(:,[1,4,6])')');
     
     MembersIndices4{ll}=(squeeze(Indices_Odor{ii,jj}.Pvalues(:,:,1)))<0.05;
     comm4{ll}=~(any(MembersIndices4{ll}(:,[2,3,5])')' & any(MembersIndices4{ll}(:,[1,4,6])')');
     end
     
     comb=[1,3];
     if ii < 11; nn=[2,3,5];
     elseif ii > 10; nn=[1,4,6];
     end
     id11=[];id22=[];id1=[];id2=[];id33=[];id44=[];
     id3=[];id4=[];id5=[];id6=[];id7=[];id8=[];gg=0;
     for jj=1:size(comb,2);
     gg=gg+1; mm1=0;     
     for aa=[1,2,length(tt1)-dd,length(tt1)]  
     mm1=mm1+1;    
     zassemblycont=CorrMat{ii,aa};     
     for bb=1:length(nn)
         id11{jj}{bb}=find((MembersIndices1{comb(jj)}(:,nn(bb))  & comm1{comb(jj)}) |(MembersIndices2{comb(jj)}(:,nn(bb))  & comm2{comb(jj)}));
         id22{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)});
 
         id1{1,jj}{bb}=(id11{jj}{bb});
         id1{2,jj}{bb}=id22{jj}{bb};
                 
      end
     comb1=[1 2];
     for ww=1:size(comb1,1)
      
     for cc=1:length(zassemblycont)   
     xx=[];    
     for bb=1:length(nn)
     indices=[];    
     for s=1:length(id1{comb1(ww,1),jj}{bb})
     for s1=1:length(id1{comb1(ww,2),jj}{bb})
     if id1{comb1(ww,1),jj}{bb}(s)== id1{comb1(ww,2),jj}{bb}(s1) ; continue; end   
     indices=[indices; (zassemblycont{cc}(id1{comb1(ww,1),jj}{bb}(s),id1{comb1(ww,2),jj}{bb}(s1)))];
     end
     end
%      if size(zassemblycont{cc},2) > 10; zz=1:10; else zz=1:size(zassemblycont{cc},2); end
%      indices=nanmean(indices')';
     indices=indices(:);
     xx=[xx; nanmean(indices)];      
     end  
     Coact{1,rr}{gg,mm1}(cc,:)=nanmean(xx);
     end
     end
    
     end
     end
         
     if ii < 11; nn=[1,4,6];
     elseif ii > 10; nn=[2,3,5];
     end
     id11=[];id22=[];id1=[];id2=[];id33=[];id44=[];
     id3=[];id4=[];id5=[];id6=[];id7=[];id8=[];gg=0;
     for jj=1:size(comb,2);
     gg=gg+1; mm1=0;     
     for aa=[1,2,length(tt1)-dd,length(tt1)]  
     mm1=mm1+1;    
     zassemblycont=CorrMat{ii,aa};     
     for bb=1:length(nn)
         id11{jj}{bb}=find((MembersIndices1{comb(jj)}(:,nn(bb)) & comm1{comb(jj)}) | (MembersIndices2{comb(jj)}(:,nn(bb)) & comm2{comb(jj)}));
         id22{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)});
 
         id1{1,jj}{bb}=id11{jj}{bb};
         id1{2,jj}{bb}=id22{jj}{bb};
                 
        end
     comb1=[1 2];
     for ww=1:size(comb1,1)
      
     for cc=1:length(zassemblycont)   
     xx=[];    
     for bb=1:length(nn)
     indices=[];    
     for s=1:length(id1{comb1(ww,1),jj}{bb})
     for s1=1:length(id1{comb1(ww,2),jj}{bb})
     if id1{comb1(ww,1),jj}{bb}(s)== id1{comb1(ww,2),jj}{bb}(s1) ; continue; end   
     indices=[indices; (zassemblycont{cc}(id1{comb1(ww,1),jj}{bb}(s),id1{comb1(ww,2),jj}{bb}(s1)))];
     end
     end
%      if size(zassemblycont{cc},2) > 10; zz=1:10; else zz=1:size(zassemblycont{cc},2); end
%      indices=nanmean(indices')';
     indices=indices(:);
     xx=[xx; nanmean(indices)];      
     end  
     Coact{2,rr}{gg,mm1}(cc,:)=nanmean(xx);
     end
     end
    
     end
     end
        
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('CuePlaceCoactivityLearnEffectOnSleepFramesHPCNewVsOld.mat','Coact')

% plot 
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CuePlaceCoactivityLearnEffectOnSleepFramesHPCNewVsOld.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,7,9,10,11,12,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=1:size(Coact{Days,sess},2)
     x=(([(Coact{1,sess}{1,jj})]'))'; 
     y=(([(Coact{2,sess}{1,jj})]'))';
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
     indices1{1,1}{jj}=[indices1{1,1}{jj}; (x)']; 
     indices1{2,1}{jj}= [indices1{2,1}{jj}; (y)']; 
     x=(([(Coact{1,sess}{2,jj})]'))'; 
     y=(([(Coact{2,sess}{2,jj})]'))';
      x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,2}{jj}=[indices{1,2}{jj}; ((x-y))'];
     indices1{2,1}{jj}=[indices1{1,2}{jj}; (x)']; 
     indices1{2,2}{jj}= [indices1{2,2}{jj}; (y)'];   
     
end
end

end

col1={rgb('Red'),...
    rgb('Blue'),rgb('Green'),rgb('LightBlue')}
fig1=figure(1)
left=0.1; bottom=0.75; width=0.20; height=0.20;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'PreBa-1','PreNo-cue','PreBa-2'};
labels2={'PreBa-1','PreAsn','PreAsd'};
for Days=1:2
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:3
    ddt1=[ddt1, (indices{1,Days}{jj})];
    ddt2=[ddt2, indices1{1,Days}{jj}];
    ddt3=[ddt3, indices1{2,Days}{jj}];
    end
   ddt_mean=movmean(nanmean(ddt1),[1 1],2);
    ddt_sem=movmean(nansem(ddt1),[1 1],2);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
    if Days==2
    plot([0 30],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[5 15 25],'xticklabel',{'PreBa-1','PreAsn','PreAsd'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Changed-Unchanged)','fontsize',5)
end


pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};

binns=[-0.05:0.001:0.05];
pval=[]; P=[];   
  for jj=1:3
    ddt1=[];ddt2=[];ddt3=[];    
left=0.35; bottom=0.9-(jj-1).*0.08; width=0.08; height=0.05;
ax1(jj)=axes('Position',[left, bottom, width, height])      
    for Days=1:2    
    ddt_mean{Days}{jj}=((indices{1,Days}{jj}(:)')');
    [h1 stats]=ecdf(ddt_mean{Days}{jj});hold on
    plot(stats,h1,'color',col1{Days});
    plot([nanmedian(ddt_mean{Days}{jj}) nanmedian(ddt_mean{Days}{jj})],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{Days});
    
    
    end
    plot([0 0],[0 1],'--k');
     [pval(jj) h stats]=signrank(ddt_mean{1}{jj},ddt_mean{2}{jj});
     P{1,3}{jj,2}=pval(jj); P{1,3}{jj,3}=stats.zval; 
        if pval(jj) < 0.05 & pval(jj) > 0.01
          text(-0.02,0.8,siglabels{1},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.01 & pval(jj) > 0.001
          text(-0.02,0.8,siglabels{2},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.001 
          text(-0.02,0.8,siglabels{3},'fontsize',6,'color','k');
        end
       
    P{1,1}{jj,1}=1;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj});  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj});
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}); 
    
        [pval1(jj) h stats]=signrank(ddt_mean{1}{jj});  
        P{1,1}{jj,2}=pval1(jj);   P{1,1}{jj,3}=stats.zval;
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95,siglabels{1},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95,siglabels{2},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95,siglabels{3},'fontsize',6,'color',col1{1});
        end
    P{1,2}{jj,1}=2;
    P{1,2}{jj,4}=length(ddt_mean{2}{jj});  P{1,2}{jj,5}=nanmean(ddt_mean{2}{jj});
    P{1,2}{jj,6}=nansem(ddt_mean{2}{jj});   
        
        
        
         [pval2(jj) h stats]=signrank(ddt_mean{2}{jj});
         P{1,2}{jj,2}=pval2(jj);   P{1,2}{jj,3}=stats.zval;
         if pval2(jj) < 0.05 & pval2(jj) > 0.01
          text(-0.01,0.9,siglabels{1},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.01 & pval2(jj) > 0.001
          text(-0.01,0.9,siglabels{2},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.001 
          text(-0.01,0.9,siglabels{3},'fontsize',6,'color',col1{2});
        end    
    xl=get(gca,'xtick');
    xl1=[nanmin(xl) nanmax(xl)];
    set(gca,'xtick',xl1,'xticklabel',xl1,'fontsize',5);    
    set(gca,'box','off','tickdir','out','fontsize',5);
    if jj==3
    xlabel('Frame prob (Changed-Unchanged)','fontsize',5)
    end
    if jj==1
    ylabel('cdf','fontsize',5);    
    end
end



set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('PlaceCueCoactivityLearningEffectsOnSleepNewVsOld'));
print(fig1,'-painters','-depsc','-r300',strcat('PlaceCueCoactivityLearningEffectsOnSleepNewVsOld')); 
close all;

%% Fig 4d
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('GoalProbLearnEffectOnSleepFramesHPCNewVsOldRippPFC.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,7,9,10,11,12,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=1:size(Probability1{Days,sess},2)
     x=(([(Probability1{Days,sess}{1,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,2))]'))';
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
     indices1{1,1}{jj}=[indices1{1,1}{jj}; (x)']; 
     indices1{2,1}{jj}= [indices1{2,1}{jj}; (y)']; 
     x=(([(Probability1{Days,sess}{3,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{3,jj}(:,2))]'))';
     indices{1,2}{jj}=[indices{1,2}{jj}; ((x-y))'];
     indices1{2,1}{jj}=[indices1{1,2}{jj}; (x)']; 
     indices1{2,2}{jj}= [indices1{2,2}{jj}; (y)'];   
     
end
end

end


col1={rgb('Red'),...
    rgb('Blue'),rgb('Green'),rgb('LightBlue')}
fig1=figure(1)
left=0.1; bottom=0.75; width=0.20; height=0.20;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'PreBa-1','PreNo-cue','PreBa-2'};
labels2={'PreBa-1','PreAsn','PreAsd'};
for Days=1:2
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:3
    ddt1=[ddt1,(indices{1,Days}{jj})];
    ddt2=[ddt2, indices1{1,Days}{jj}];
    ddt3=[ddt3, indices1{2,Days}{jj}];
    end
    ddt_mean=movmean(nanmean(ddt1),[1 1],2);
    ddt_sem=movmean(nansem(ddt1),[1 1],2);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
    if Days==2
    plot([0 30],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[5 15 25],'xticklabel',{'PreBa-1','PreAsn','PreAsd'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Changed-Unchanged)','fontsize',5)
end


pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};

binns=[-0.1:0.01:0.1];
pval=[]; P=[];   
  for jj=1:3
    ddt1=[];ddt2=[];ddt3=[];    
left=0.35; bottom=0.9-(jj-1).*0.08; width=0.08; height=0.05;
ax1(jj)=axes('Position',[left, bottom, width, height])      
    for Days=1:2    
    ddt_mean{Days}{jj}=((indices{1,Days}{jj}(:)')');
    [h1 stats]=ecdf(ddt_mean{Days}{jj});hold on
    plot(stats,h1,'color',col1{Days});
    plot([nanmedian(ddt_mean{Days}{jj}) nanmedian(ddt_mean{Days}{jj})],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{Days});
    
    
    end
    plot([0 0],[0 1],'--k');
     [pval(jj) h stats]=signrank(ddt_mean{1}{jj},ddt_mean{2}{jj});
     P{1,3}{jj,2}=pval(jj); P{1,3}{jj,3}=stats.zval; 
        if pval(jj) < 0.05 & pval(jj) > 0.01
          text(-0.02,0.8,siglabels{1},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.01 & pval(jj) > 0.001
          text(-0.02,0.8,siglabels{2},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.001 
          text(-0.02,0.8,siglabels{3},'fontsize',6,'color','k');
        end
       
    P{1,1}{jj,1}=1;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj});  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj});
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}); 
    
        [pval1(jj) h stats]=signrank(ddt_mean{1}{jj});  
        P{1,1}{jj,2}=pval1(jj);   P{1,1}{jj,3}=stats.zval;
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95,siglabels{1},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95,siglabels{2},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95,siglabels{3},'fontsize',6,'color',col1{1});
        end
    P{1,2}{jj,1}=2;
    P{1,2}{jj,4}=length(ddt_mean{2}{jj});  P{1,2}{jj,5}=nanmean(ddt_mean{2}{jj});
    P{1,2}{jj,6}=nansem(ddt_mean{2}{jj});   
        
        
        
         [pval2(jj) h stats]=signrank(ddt_mean{2}{jj});
         P{1,2}{jj,2}=pval2(jj);   P{1,2}{jj,3}=stats.zval;
         if pval2(jj) < 0.05 & pval2(jj) > 0.01
          text(-0.01,0.9,siglabels{1},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.01 & pval2(jj) > 0.001
          text(-0.01,0.9,siglabels{2},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.001 
          text(-0.01,0.9,siglabels{3},'fontsize',6,'color',col1{2});
        end    
    xl=get(gca,'xtick');
    xl1=[nanmin(xl) nanmax(xl)];
    set(gca,'xtick',xl1,'xticklabel',xl1,'fontsize',5);    
    set(gca,'box','off','tickdir','out','fontsize',5);
    if jj==3
    xlabel('Frame prob (Changed-Unchanged)','fontsize',5)
    end
    if jj==1
    ylabel('cdf','fontsize',5);    
    end
end
  

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('GoalProbLearningEffectsOnSleepFramesNewVsOldRippPFC'));
print(fig1,'-painters','-depsc','-r300',strcat('GoalProbLearningEffectsOnSleepFramesNewVsOldRippPFC')); 
close all;

%% Fig 4e
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SleepFramesCorrMat10binsAllSleep20MultiHPCPFC');
load('IndicesPlaceArmsHPCPFC')
load('IndicesFlavorHPCPFC')
load('IndicesOdorHPCPFC')

sd=2;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);gg=0;
Coact=[];MembersZcontOld=[];
rr=0;
for ii=[6:15]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
 
 [c1 c2] = textread('celltype','%s %s');
elecreg=importdata('elecreg');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
ind=a12(:,1)==1 & a12(:,2)==1 & (elecreg(:,2)==1 | elecreg(:,2)==2);
ind_all=elecreg;
ind_all=elecreg(ind_HPC | ind_PFC,2);
if sum(ind_HPC) < 10 | sum(ind_PFC)<10; continue; end 
  
    tt=cell2mat(cellfun(@(x)(~isempty(x)),Indices_Place(ii,:),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    gg=0;rr=rr+1;
    if any(ii==[1:5,21:25,13,18,8]) ii==6 | ii==15; dd=0; else dd=1;end
    
   
    MembersIndices1=[];comm1=[];ll=0;
    MembersIndices2=[];comm2=[];
    MembersIndices3=[];comm3=[];
    MembersIndices4=[];comm4=[];mm1=0;
     for jj=[1,2,length(tt1)-dd,length(tt1)]   
      ll=ll+1; 
     MembersIndices1{ll}=(squeeze(Indices_Place{ii,jj}.Pvalues1(:,:,1,1)))<0.05;
     comm1{ll}=~(any(MembersIndices1{ll}(:,[2,3,5])')' & any(MembersIndices1{ll}(:,[1,4,6])')');
     
     MembersIndices2{ll}=(squeeze(Indices_Place{ii,jj}.Pvalues1(:,:,2,1)))<0.05;
     comm2{ll}=~(any(MembersIndices2{ll}(:,[2,3,5])')' & any(MembersIndices2{ll}(:,[1,4,6])')');
     
     MembersIndices3{ll}=(squeeze(Indices_Flavor{ii,jj}.Pvalues(:,:,1)))<0.05;
     comm3{ll}=~(any(MembersIndices3{ll}(:,[2,3,5])')' & any(MembersIndices3{ll}(:,[1,4,6])')');
     
     MembersIndices4{ll}=(squeeze(Indices_Odor{ii,jj}.Pvalues(:,:,1)))<0.05;
     comm4{ll}=~(any(MembersIndices4{ll}(:,[2,3,5])')' & any(MembersIndices4{ll}(:,[1,4,6])')');
     end
     
     comb=[1,3];
     if ii < 11; nn=[2,3,5];
     elseif ii > 10; nn=[1,4,6];
     end
     id11=[];id22=[];id1=[];id2=[];id33=[];id44=[];
     id3=[];id4=[];id5=[];id6=[];id7=[];id8=[];gg=0;
     for jj=1:size(comb,2);
     gg=gg+1; mm1=0;     
     for aa=[1,2,length(tt1)-dd,length(tt1)]  
     mm1=mm1+1;    
     zassemblycont=CorrMat{ii,aa};     
     for bb=1:length(nn)
         id11{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==1);
         id22{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==2);
 
         id1{1,jj}{bb}=(id11{jj}{bb});
         id1{2,jj}{bb}=id22{jj}{bb};
                 
      end
     comb1=[1 2];
     for ww=1:size(comb1,1)
      
     for cc=1:length(zassemblycont)   
     xx=[];    
     for bb=1:length(nn)
     indices=[];    
     for s=1:length(id1{comb1(ww,1),jj}{bb})
     for s1=1:length(id1{comb1(ww,2),jj}{bb})
     if id1{comb1(ww,1),jj}{bb}(s)== id1{comb1(ww,2),jj}{bb}(s1) ; continue; end   
     indices=[indices; (zassemblycont{cc}(id1{comb1(ww,1),jj}{bb}(s),id1{comb1(ww,2),jj}{bb}(s1)))];
     end
     end
%      if size(zassemblycont{cc},2) > 10; zz=1:10; else zz=1:size(zassemblycont{cc},2); end
%      indices=nanmean(indices')';
     indices=indices(:);
     xx=[xx; nanmean(indices)];      
     end  
     Coact{1,rr}{gg,mm1}(cc,:)=nanmean(xx);
     end
     end
    
     end
     end
         
     if ii < 11; nn=[1,4,6];
     elseif ii > 10; nn=[2,3,5];
     end
     id11=[];id22=[];id1=[];id2=[];id33=[];id44=[];
     id3=[];id4=[];id5=[];id6=[];id7=[];id8=[];gg=0;
     for jj=1:size(comb,2);
     gg=gg+1; mm1=0;     
     for aa=[1,2,length(tt1)-dd,length(tt1)]  
     mm1=mm1+1;    
     zassemblycont=CorrMat{ii,aa};     
     for bb=1:length(nn)
         id11{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==1);
         id22{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==2);
 
         id1{1,jj}{bb}=id11{jj}{bb};
         id1{2,jj}{bb}=id22{jj}{bb};
                 
        end
     comb1=[1 2];
     for ww=1:size(comb1,1)
      
     for cc=1:length(zassemblycont)   
     xx=[];    
     for bb=1:length(nn)
     indices=[];    
     for s=1:length(id1{comb1(ww,1),jj}{bb})
     for s1=1:length(id1{comb1(ww,2),jj}{bb})
     if id1{comb1(ww,1),jj}{bb}(s)== id1{comb1(ww,2),jj}{bb}(s1) ; continue; end   
     indices=[indices; (zassemblycont{cc}(id1{comb1(ww,1),jj}{bb}(s),id1{comb1(ww,2),jj}{bb}(s1)))];
     end
     end
%      if size(zassemblycont{cc},2) > 10; zz=1:10; else zz=1:size(zassemblycont{cc},2); end
%      indices=nanmean(indices')';
     indices=indices(:);
     xx=[xx; nanmean(indices)];      
     end  
     Coact{2,rr}{gg,mm1}(cc,:)=nanmean(xx);
     end
     end
    
     end
     end
        
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath);
save('CueCueCoactivityLearnEffectOnSleepFramesHPCPFCNewVsOld.mat','Coact')

% plot
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CueCueCoactivityLearnEffectOnSleepFramesHPCPFCNewVsOld.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,9,10,11,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=1:size(Coact{Days,sess},2)
     x=(([(Coact{1,sess}{1,jj})]'))'; 
     y=(([(Coact{2,sess}{1,jj})]'))';
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
     indices1{1,1}{jj}=[indices1{1,1}{jj}; (x)']; 
     indices1{2,1}{jj}= [indices1{2,1}{jj}; (y)']; 
     x=(([(Coact{1,sess}{2,jj})]'))'; 
     y=(([(Coact{2,sess}{2,jj})]'))';
      x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,2}{jj}=[indices{1,2}{jj}; ((x-y))'];
     indices1{2,1}{jj}=[indices1{1,2}{jj}; (x)']; 
     indices1{2,2}{jj}= [indices1{2,2}{jj}; (y)'];   
     
end
end

end

col1={rgb('Red'),...
    rgb('Blue'),rgb('Green'),rgb('LightBlue')}
fig1=figure(1)
left=0.1; bottom=0.75; width=0.20; height=0.20;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'PreBa-1','PreNo-cue','PreBa-2'};
labels2={'PreBa-1','PreAsn','PreAsd'};
for Days=1:2
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:3
    ddt1=[ddt1, (indices{1,Days}{jj})];
    ddt2=[ddt2, indices1{1,Days}{jj}];
    ddt3=[ddt3, indices1{2,Days}{jj}];
    end
    
    ddt_mean=movmean(nanmean(ddt1),[1 1],2);
    ddt_sem=movmean(nansem(ddt1),[1 1],2);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
    if Days==2
    plot([0 30],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',4,'tickdir','out',...
        'xtick',[5 15 25],'xticklabel',{'PreBa-1','PreAsn','PreAsd'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Changed-Unchanged)','fontsize',5)
    ylim([-0.05 0.035]) 
end


pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};

binns=[-0.05:0.001:0.05];
pval=[]; P=[];   
  for jj=1:3
    ddt1=[];ddt2=[];ddt3=[];    
left=0.35; bottom=0.9-(jj-1).*0.08; width=0.08; height=0.05;
ax1(jj)=axes('Position',[left, bottom, width, height])      
    for Days=1:2    
    ddt_mean{Days}{jj}=((indices{1,Days}{jj}(:)')');
    [h1 stats]=ecdf(ddt_mean{Days}{jj});hold on
    plot(stats,h1,'color',col1{Days});
    plot([nanmedian(ddt_mean{Days}{jj}) nanmedian(ddt_mean{Days}{jj})],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{Days});
    
    
    end
    plot([0 0],[0 1],'--k');
     [pval(jj) h stats]=signrank(ddt_mean{1}{jj},ddt_mean{2}{jj});
     P{1,3}{jj,2}=pval(jj); P{1,3}{jj,3}=stats.zval; 
        if pval(jj) < 0.05 & pval(jj) > 0.01
          text(-0.02,0.8,siglabels{1},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.01 & pval(jj) > 0.001
          text(-0.02,0.8,siglabels{2},'fontsize',6,'color','k');
        elseif  pval(jj) < 0.001 
          text(-0.02,0.8,siglabels{3},'fontsize',6,'color','k');
        end
       
    P{1,1}{jj,1}=1;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj});  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj});
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}); 
    
        [pval1(jj) h stats]=signrank(ddt_mean{1}{jj});  
        P{1,1}{jj,2}=pval1(jj);   P{1,1}{jj,3}=stats.zval;
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95,siglabels{1},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95,siglabels{2},'fontsize',6,'color',col1{1});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95,siglabels{3},'fontsize',6,'color',col1{1});
        end
    P{1,2}{jj,1}=2;
    P{1,2}{jj,4}=length(ddt_mean{2}{jj});  P{1,2}{jj,5}=nanmean(ddt_mean{2}{jj});
    P{1,2}{jj,6}=nansem(ddt_mean{2}{jj});   
        
        
        
         [pval2(jj) h stats]=signrank(ddt_mean{2}{jj});
         P{1,2}{jj,2}=pval2(jj);   P{1,2}{jj,3}=stats.zval;
         if pval2(jj) < 0.05 & pval2(jj) > 0.01
          text(-0.01,0.9,siglabels{1},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.01 & pval2(jj) > 0.001
          text(-0.01,0.9,siglabels{2},'fontsize',6,'color',col1{2});
        elseif  pval2(jj) < 0.001 
          text(-0.01,0.9,siglabels{3},'fontsize',6,'color',col1{2});
        end    
    xl=get(gca,'xtick');
    xl1=[nanmin(xl) nanmax(xl)];
    set(gca,'xtick',xl1,'xticklabel',xl1,'fontsize',5);    
    set(gca,'box','off','tickdir','out','fontsize',5);
    if jj==3
    xlabel('Frame prob (Changed-Unchanged)','fontsize',5)
    end
    if jj==1
    ylabel('cdf','fontsize',5);    
    end
end

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('CueCueCoactivityLearningEffectsOnSleepNewVsOld'));
print(fig1,'-painters','-depsc','-r300',strcat('CueCueCoactivityLearningEffectsOnSleepNewVsOld')); 
close all;


%% Fig 4f
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('TrackProbLearnEffectOnSleepFramesHPCNewVsOldRipp.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,7,9,10,11,12,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=2:size(Probability1{Days,sess},2)
     x=(nanmean([(Probability1{Days,sess}{1,jj}(:,1)), (Probability1{Days,sess}{1,jj}(:,2))]'))'; 
     y=(nanmean([(Probability1{Days,sess}{1,jj}(:,3)), (Probability1{Days,sess}{1,jj}(:,4))]'))';
     
     x1=(nanmean([(Probability1{Days,sess}{1,1}(:,1)), (Probability1{Days,sess}{1,1}(:,2))]'))'; 
     y1=(nanmean([(Probability1{Days,sess}{1,1}(:,3)), (Probability1{Days,sess}{1,1}(:,4))]'))';
     indices1{1,1}{jj-1}=[indices1{1,1}{jj-1}; (x-x1)']; 
     indices1{1,2}{jj-1}= [indices1{1,2}{jj-1}; (y-y1)']; 
     
     x=(nanmean([(Probability1{Days,sess}{3,jj}(:,1)), (Probability1{Days,sess}{3,jj}(:,2))]'))'; 
     y=(nanmean([(Probability1{Days,sess}{3,jj}(:,3)), (Probability1{Days,sess}{3,jj}(:,4))]'))';
     x1=(nanmean([(Probability1{Days,sess}{3,1}(:,1)), (Probability1{Days,sess}{3,1}(:,2))]'))'; 
     y1=(nanmean([(Probability1{Days,sess}{3,1}(:,3)), (Probability1{Days,sess}{3,1}(:,4))]'))';
     indices1{2,1}{jj-1}=[indices1{2,1}{jj-1}; (x-x1)']; 
     indices1{2,2}{jj-1}= [indices1{2,2}{jj-1}; (y-y1)'];   
     
end
end

end
col1={rgb('Red'),...
    rgb('Blue'),rgb('Red'),rgb('Blue')}
fig1=figure(1)

left=0.1; bottom=0.65; width=0.20; height=0.12;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval1=[]; cons=[];nn=0; pvalInd=[]; consInd=[];P=[];
xpos=[0,3,6,9];ddt_mean=[];siglabels={'*','**','***'};rr=0;ddt_fortest=[];
for Days=1
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:2
    rr=rr+1;    
    if jj==1
    ddt_mean{Days}{jj}=([indices1{2,1}{2}]'-[indices1{1,1}{2}]') ;
    else
    ddt_mean{Days}{jj}=([indices1{2,2}{2}]'-[indices1{1,2}{2}]') ;
    end
    ddt_fortest{rr}=ddt_mean{Days}{jj}(:);
    [h1 stats]=ecdf(ddt_mean{Days}{jj}(:));hold on
    plot(stats,h1,'color',col1{rr});
    plot([nanmedian(ddt_mean{Days}{jj}(:)) nanmedian(ddt_mean{Days}{jj}(:))],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{rr});
     plot([0 0],[0 1],'--k');
       [ pval1(jj) h stats]=signrank(ddt_mean{Days}{jj}(:));    
       P{1}{jj,1}=pval1(jj); P{1}{jj,4}=nanmean(ddt_mean{Days}{jj}(:));P{1}{jj,5}=nansem(ddt_mean{Days}{jj}(:))
       P{1}{jj,2}=stats.zval;    P{1}{jj,3}=stats.signedrank;
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col1{rr});
        end
       
    
    
    end

    
   
    set(gca,'box','off','fontsize',4,'tickdir','out');
    
    
    xlabel('\Delta prob normalized','fontsize',5);
    ylabel('cdf','fontsize',5)
    xlim([-0.04 0.06])
end

Groups=[1,2]; h=[]; pval1=[];
col2={rgb('Gray'),rgb('Black')};
for jj=1:size(Groups,1)
       [pval1(jj) h stats]=signrank(ddt_fortest{Groups(jj,1)},ddt_fortest{Groups(jj,2)});
      
      P{2}{jj,2}=pval1(jj);    P{2}{jj,3}=stats.zval;    P{2}{jj,4}=stats.signedrank; 
        P{2}{jj,1}=Groups(jj,:); 
      if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.85-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.85-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.85-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col2{jj});
        end
end      

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,3])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('TrackProbNormalizedLearningEffectsOnSleepFramesNewVsOldRipp'));
print(fig1,'-painters','-depsc','-r300',strcat('TrackProbNormalizedLearningEffectsOnSleepFramesNewVsOldRipp')); 
close all;

%% Fig 4g
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CuePlaceCoactivityLearnEffectOnSleepFramesHPCNewVsOld.mat');


%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,9,10,11,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=2:size(Coact{Days,sess},2)
     x=(([(Coact{1,sess}{1,jj})])); 
     y=(([(Coact{2,sess}{1,jj})]));
%        x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     x1=(([(Coact{1,sess}{1,1})])); 
     y1=(([(Coact{2,sess}{1,1})]));
%      x1(:,all(isnan(x1)))=0; y1(:,all(isnan(y1)))=0;
     indices1{1,1}{jj-1}=[indices1{1,1}{jj-1}; (x-x1)']; 
     indices1{1,2}{jj-1}= [indices1{1,2}{jj-1}; (y-y1)']; 
     
     x=(([(Coact{1,sess}{2,jj})])); 
     y=(([(Coact{2,sess}{2,jj})]));
%      x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     x1=(([(Coact{1,sess}{2,1})])); 
     y1=(([(Coact{2,sess}{2,1})]));
%      x1(:,all(isnan(x1)))=0; y1(:,all(isnan(y1)))=0;
     indices1{2,1}{jj-1}=[indices1{2,1}{jj-1}; (x-x1)']; 
     indices1{2,2}{jj-1}= [indices1{2,2}{jj-1}; (y-y1)'];   
     
end
end

end
col1={rgb('Red'),...
    rgb('Blue'),rgb('Red'),rgb('Blue')}
fig1=figure(1)

left=0.1; bottom=0.65; width=0.20; height=0.12;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval1=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,3,6,9];ddt_mean=[];siglabels={'*','**','***'};rr=0;ddt_fortest=[];
P=[];
for Days=1
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:2
    rr=rr+1;    
   if jj==1
    ddt_mean{Days}{jj}=([indices1{2,1}{2}]'-[indices1{1,1}{1}]') ;
    else
    ddt_mean{Days}{jj}=([indices1{2,2}{2}]'-[indices1{1,2}{1}]') ;
    end
    ddt_fortest{rr}=ddt_mean{Days}{jj}(:);
    [h1 stats]=ecdf(ddt_mean{Days}{jj}(:));hold on
    plot(stats,h1,'color',col1{rr});
    plot([nanmedian(ddt_mean{Days}{jj}(:)) nanmedian(ddt_mean{Days}{jj}(:))],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{rr});
    
       [ pval1(jj) h stats]=signrank(ddt_mean{Days}{jj}(:));    
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col1{rr});
        end
    P{1,1}{jj,1}=1; P{1,1}{jj,2}=pval1(jj);   P{1,1}{jj,3}=stats.zval;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj}(:));  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj}(:));
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}(:));    
    
    end   
   
    set(gca,'box','off','fontsize',4,'tickdir','out');
    xlabel('\Delta prob normalized','fontsize',5);
    ylabel('cdf','fontsize',5)
    xlim([-0.06 0.06])
end

Groups=[1,2]; h=[]; pval1=[];
col2={rgb('Gray'),rgb('Black')};
for jj=1:size(Groups,1)
       [pval1(jj) h stats]=ranksum(ddt_fortest{Groups(jj,1)},ddt_fortest{Groups(jj,2)});
        P{1,2}{jj,2}=pval1(jj); P{1,2}{jj,3}=stats.zval; 
         
       if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.85-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.85-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.85-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col2{jj});
        end
end      

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,3])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('PlaceCueCoactivityNormalizedLearningEffectsOnSleepFramesNewVsOld'));
print(fig1,'-painters','-depsc','-r300',strcat('PlaceCueCoactivityNormalizedLearningEffectsOnSleepFramesNewVsOld')); 
close all;

%% Fig 4h
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('GoalProbLearnEffectOnSleepFramesPFCNewVsOld.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,7,9,10,11,12,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=2:size(Probability1{Days,sess},2)
     x=(([(Probability1{Days,sess}{1,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,2))]'))';
     
     x1=(([(Probability1{Days,sess}{1,1}(:,1))]'))'; 
     y1=(([(Probability1{Days,sess}{1,1}(:,2))]'))';
     indices1{1,1}{jj-1}=[indices1{1,1}{jj-1}; (x-x1)']; 
     indices1{1,2}{jj-1}= [indices1{1,2}{jj-1}; (y-y1)']; 
     
     x=(([(Probability1{Days,sess}{3,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{3,jj}(:,2))]'))';
     x1=(([(Probability1{Days,sess}{3,1}(:,1))]'))'; 
     y1=(([(Probability1{Days,sess}{3,1}(:,2))]'))';
     indices1{2,1}{jj-1}=[indices1{2,1}{jj-1}; (x-x1)']; 
     indices1{2,2}{jj-1}= [indices1{2,2}{jj-1}; (y-y1)'];   
     
end
end

end
col1={rgb('Red'),...
    rgb('Blue'),rgb('Red'),rgb('Blue')}
fig1=figure(1)

left=0.05; bottom=0.65; width=0.20; height=0.12;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,3,6,9];ddt_mean=[];siglabels={'*','**','***'};rr=0;ddt_fortest=[];
for Days=1
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:2
    rr=rr+1;    
    if jj==1
    ddt_mean{Days}{jj}=([indices1{2,1}{2}]'-[indices1{1,1}{2}]') ;
    else
    ddt_mean{Days}{jj}=([indices1{2,2}{2}]'-[indices1{1,2}{2}]') ;
    end
    ddt_fortest{rr}=ddt_mean{Days}{jj}(:);
    [h1 stats]=ecdf(ddt_mean{Days}{jj}(:));hold on
    plot(stats,h1,'color',col1{rr});
    plot([nanmedian(ddt_mean{Days}{jj}(:)) nanmedian(ddt_mean{Days}{jj}(:))],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{rr});
    
       [ pval1(jj) h]=signrank(ddt_mean{Days}{jj}(:));    
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col1{rr});
        end
       
    
    
    end

    
   
    set(gca,'box','off','fontsize',4,'tickdir','out');
    
    
    xlabel('\Delta prob normalized','fontsize',5);
    ylabel('cdf','fontsize',5)
    xlim([-0.04 0.06])
end

Groups=[1,2]; h=[]; pval1=[];
col2={rgb('Gray'),rgb('Black')};
for jj=1:size(Groups,1)
       [pval1(jj) h ]=signrank(ddt_fortest{Groups(jj,1)},ddt_fortest{Groups(jj,2)});
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.85-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.85-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.85-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col2{jj});
        end
end      



set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('GoalProbNormalizedLearningEffectsOnSleepFramesNewVsOldPFCRipp'));
print(fig1,'-painters','-depsc','-r300',strcat('GoalProbNormalizedLearningEffectsOnSleepFramesNewVsOldPFCRipp')); 
close all;

%% Fig 4i

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CueCueCoactivityLearnEffectOnSleepFramesHPCPFCNewVsOld.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:4
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,9,10,11,14,15]};
for Days=1
for sess=1:length(Sess{Days})
    for jj=2:size(Coact{Days,sess},2)
     x=(([(Coact{1,sess}{1,jj})])); 
     y=(([(Coact{2,sess}{1,jj})]));
       x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     x1=(([(Coact{1,sess}{1,1})])); 
     y1=(([(Coact{2,sess}{1,1})]));
     x1(:,all(isnan(x1)))=0; y1(:,all(isnan(y1)))=0;
     indices1{1,1}{jj-1}=[indices1{1,1}{jj-1}; (x-x1)']; 
     indices1{1,2}{jj-1}= [indices1{1,2}{jj-1}; (y-y1)']; 
     
     x=(([(Coact{1,sess}{2,jj})])); 
     y=(([(Coact{2,sess}{2,jj})]));
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     x1=(([(Coact{1,sess}{2,1})])); 
     y1=(([(Coact{2,sess}{2,1})]));
     x1(:,all(isnan(x1)))=0; y1(:,all(isnan(y1)))=0;
     indices1{2,1}{jj-1}=[indices1{2,1}{jj-1}; (x-x1)']; 
     indices1{2,2}{jj-1}= [indices1{2,2}{jj-1}; (y-y1)'];   
     
end
end

end
col1={rgb('Red'),...
    rgb('Blue'),rgb('Red'),rgb('Blue')}
fig1=figure(1)

left=0.1; bottom=0.65; width=0.20; height=0.12;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval1=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,3,6,9];ddt_mean=[];siglabels={'*','**','***'};rr=0;ddt_fortest=[];
P=[];
for Days=1
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:2
    rr=rr+1;    
    if jj==1
    ddt_mean{Days}{jj}=([indices1{2,1}{2}]'-[indices1{1,1}{2}]') ;
    else
    ddt_mean{Days}{jj}=([indices1{2,2}{2}]'-[indices1{1,2}{2}]') ;
    end
    ddt_fortest{rr}=ddt_mean{Days}{jj}(:);
    [h1 stats]=ecdf(ddt_mean{Days}{jj}(:));hold on
    plot(stats,h1,'color',col1{rr});
    plot([nanmedian(ddt_mean{Days}{jj}(:)) nanmedian(ddt_mean{Days}{jj}(:))],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{rr});
    plot([0 0],[0 1],'--k');
       [ pval1(jj) h stats]=signrank(ddt_mean{Days}{jj}(:));    
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.95-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.95-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col1{rr});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.95-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col1{rr});
        end
    P{1,1}{jj,1}=1; P{1,1}{jj,2}=pval1(jj);   P{1,1}{jj,3}=stats.zval;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj}(:));  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj}(:));
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}(:));     
    
    end
   set(gca,'box','off','fontsize',4,'tickdir','out');    
    xlabel('\Delta prob normalized','fontsize',5);
    ylabel('cdf','fontsize',5)
    xlim([-0.06 0.06])
end

Groups=[1,2]; h=[]; pval1=[];
col2={rgb('Gray'),rgb('Black')};
for jj=1:size(Groups,1)
       [pval1(jj) h stats]=ranksum(ddt_fortest{Groups(jj,1)},ddt_fortest{Groups(jj,2)});
        P{1,2}{jj,2}=pval1(jj); P{1,2}{jj,3}=stats.zval; 
         
       if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.01,0.85-(rr-1).*0.05,siglabels{1},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.01,0.85-(rr-1).*0.05,siglabels{2},'fontsize',6,'color',col2{jj});
        elseif  pval1(jj) < 0.001 
          text(-0.01,0.85-(rr-1).*0.05,siglabels{3},'fontsize',6,'color',col2{jj});
        end
end      

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,4,3])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig4/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('CueCueCoactivityNormalizedLearningEffectsOnSleepFramesHPCPFCNewVsOld'));
print(fig1,'-painters','-depsc','-r300',strcat('CueCueCoactivityNormalizedLearningEffectsOnSleepFramesHPCPFCNewVsOld')); 
close all;




