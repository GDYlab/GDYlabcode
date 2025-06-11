

%% Learning Curves panel 1D

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load Target_arms
load Chosen_arms

performance{1}=Chosen_arms1;
performance{2}=Chosen_arms2;
performance{3}=Chosen_arms3;
performance{4}=Chosen_arms4;
performance{5}=Chosen_arms5;
performance{6}=Chosen_arms6;
performance{7}=Chosen_arms7;
% performance{8}=(cellfun(@length,Chosen_arms8)-1)
performance{5}=[performance{5}, num2cell(ones(18,1).*NaN)];
performance{6}=[performance{6}, num2cell(ones(18,1).*NaN)];
performance{7}=[performance{7}, num2cell(ones(18,1).*NaN)];


f1=figure
c={'r','g','b','m','r','g','b','m'};
c1='k';
performance13=[];
for i=1:7
    startvisits=cellfun(@(x) sum(x==3 | x==6),performance{i},'UniformOutput',false); 
    performance13=[performance13; nanmean(cell2mat(startvisits))]
end
performance13(5:7,36)=NaN;
performance12=[];
for i=1:7
    startvisits=cellfun(@(x) length(x)-1,performance{i},'UniformOutput',false); 
    performance12=[performance12; nanmean(cell2mat(startvisits))];
end
performance12(5:7,36)=NaN;
yyaxis left
errorbar(1:size(performance12,2), nanmean(performance12),nansem(performance12),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color','k','capsize',2); hold on;
Lyax =gca; 
Lyact= get(Lyax, 'YTick');     
p=[];stats=[];
for a=1:size(performance12,2)
    [h(a) p(a) CI stats{a}] = ttest(performance12(:,a),performance13(:,a).*0.875+2.5,'Tail','left');
    yl=get(gca,'ylim'); 
    if p(a) < 0.05 & p(a)>0.01 & ~any(a==[1,2,3,4,5,6])
    t= text(a-0.1,min(yl),'*','fontsize',5,'color','r');
      t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 & ~any(a==[1,2,3,4,5,6])
     t=text(a-0.1,min(yl),'**','fontsize',5,'color','r');  
     t.Rotation=270;
    elseif p(a) < 0.001 & ~any(a==[1,2,3,4,5,6])
     t=text(a-0.1,min(yl),'***','fontsize',5,'color','r');  
     t.Rotation=270;
    end
end
p=[];
for a=1:size(performance12,2)
    [h(a) p(a) CI stats] = ttest(performance12(:,a),ones(size(performance12(:,a),1),1).*2.5,'Alpha',0.01,'Tail','left');
    yl=get(gca,'ylim'); 
    if p(a) < 0.05 & p(a)>0.01 & ~any(a==[1,2,3,4,5,6])
    t= text(a-0.1,min(yl)+0.4,'*','fontsize',5,'color','k');
    t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 & ~any(a==[1,2,3,4,5,6])
    t=text(a-0.1,min(yl)+0.4,'**','fontsize',5,'fontweight','bold','color','k'); 
    t.Rotation=270;
    elseif p(a) < 0.001 & ~any(a==[1,2,3,4,5,6])
    t=text(a-0.1,min(yl)+0.4,'***','fontsize',5,'color','k');  
    t.Rotation=270;
    end
end

plot(1:size(performance13,2),nanmean(performance13).*0.875+2.5,'--','color','r');
plot([0 51],[2.5 2.5],'--','color',c1);
% errorbar(1:size(performance12,2), nanmean(performance12),nansem(performance12),'linewidth', 1, 'marker', 'o', 'markersize', 3,'color','m');hold on;
set(gca,'box','off','fontsize',6);
xlabel('Training Day','fontsize',6);
ylabel('Average Number of Errors','fontsize',6);
xlim([0 37]);
set(gca,'xtick',[0 10 20 30 36],'xticklabel',[0 10 20 30 36],'fontsize',6,'tickdir','out')
set(gca,'ytick',[0:5],'yticklabel',[0,1,2,3,4,5],'fontsize',6,'tickdir','out');
ylim([0 5]);
set(gca,'ydir','reverse');

performance112=nanmean(performance12);
sem112=nansem(performance12);
performance112=(6-(performance112))/6.*100;
sem=zeros(size(performance112,2));
c={'r','g','b','m','r','g','b','m'};
c1='k';
hold on
yyaxis right
errorbar(1:size(performance12,2), nanmean(performance12),nansem(performance12),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color','k','capsize',2); hold on;
ylim([0 5]);

Ryax = gca; 
                           % Right axis handle
Ryaxt = get(Ryax, 'YTick');                 % Get Left axis ticks
RyaxDegC = ((6-(Ryaxt-1))/6).*100;     % Convert to °C (or whatever you want)

RyaxDegC = round(((6-([5,4,3,2,1]))/6).*100); 
set(Ryax, 'YTick',[1:5], 'YTickLabel', flip(RyaxDegC) ) ;   % New Y-tick values
set(gca,'ydir','reverse')

set(gca,'box','off','fontsize',6);
xlabel('Training Day','fontsize',6);
ylabel('performance (%)','fontsize',6);
xlim([0 37]);
% Set the color of each axis to black
Ryax .YAxis(1).Color = [0 0 0];
Ryax .YAxis(2).Color = [0 0 0];

GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
set(gcf,'paperunits','centimeters')
set(gcf,'papertype','a4')
set(gcf,'paperposition',[1,1,3.5,3]) 
print(f1,'-djpeg','-r300','Behavior_details_DynamicErrorTogethorRev');
print(f1,'-painters','-depsc','-r300','-r300','Behavior_details_DynamicErrorTogethorRev');
close(f1);


%% Strategies panel 1e
%% compare cumulative improvements with reward history vs against chance levels
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load ChoicesAllAnimals.mat
load Chosen_arms

performance_learn{1}=Chosen_arms1;
performance_learn{2}=Chosen_arms2;
performance_learn{3}=Chosen_arms3;
performance_learn{4}=Chosen_arms4;
performance_learn{5}=Chosen_arms5;
performance_learn{6}=Chosen_arms6;
performance_learn{7}=Chosen_arms7;
Days=[1,10; 11,20;21,30];
NCDays={[31,33,35],[32,34,36]};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];
fig1=figure

col1={rgb('Red');rgb('Salmon');rgb('Pink');rgb('LightBlue');rgb('DeepSkyBlue');rgb('Blue')};


% col1={rgb('Red'),rgb('Black'),rgb('Red');...
%     rgb('Red'),rgb('Blue'),rgb('Green');...
%     rgb('Red'),rgb('Gray'),rgb('Red');...
%     rgb('Red'),rgb('Blue'),rgb('Green')}
IntraGroup=[];BetweenGroup=[];Means=[];P1=[];Tstat1=[];Sems=[];df1=[];P2=[];df2=[];Tstat2=[];
xplot=[1,2;4,5;7,8;10,11];
for jj=1:3
 left=0.1; bot=0.75; width=0.15; height=0.15;    
%  ax1(jj)=axes('Position',[left bot width height]);
 perf1=[];T11=[];  
for ii=1:7   
p1=[];
Trials=[1:6];T1=[];
for kk=1:3
x11=(performance_learn{ii}(blocks(kk,1):blocks(kk,2),Days(jj,1):Days(jj,2)));
x12=cellfun(@(x) length(x)-1,x11,'UniformOutput',false); 
x1=100-((cell2mat(x12)/5)'.*100);
if size(x1,1) > 1
p1=[p1, (x1)];
else
p1=[p1, (x1)];
end
T1=[T1, Trials];
end
perf1=[perf1; nanmean(p1)];
T11=[T11; (T1)];
hold on;
% plot([1:18]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',col1{ss,jj},'linestyle','none'); 
end
incl1=[1:2;7,8;13:14];incl2=[5:6;11,12;17:18];Strategy1=[];
for kk=1:3
Strategy1{2}(kk,:)=nanmean((perf1(:,incl2(kk,:))-perf1(:,incl1(kk,:)))');
Strategy1{1}(kk,:)=nanmean((perf1(:,[incl1(kk,1):incl2(kk,end)])')-50);
end
Strategy{1}=Strategy1{1}(:);
Strategy{2}=Strategy1{2}(:);
[hh pval2 ci stats]=ttest(Strategy{1},Strategy{2});
P2=[P2; pval2];Tstat2=[Tstat2; stats.tstat];df2=[df2; stats.df];
if pval2<0.05
hh1=sigstar(xplot(jj,[1,2]),pval2,0)
end
col2={rgb('red'),rgb('Black')};pval1=[];
for ss1=1:length(Strategy)
    hold on;
h1(xplot(jj,ss1))=bar(xplot(jj,ss1),nanmean(Strategy{ss1}));
h1(xplot(jj,ss1)).FaceColor='none'; h1(xplot(jj,ss1)).EdgeColor=col2{ss1};

for ss=1:size(Strategy{ss1},1)
    plot(xplot(jj,ss1)+randn(1).*0.1,Strategy{ss1}(ss),'marker','o','markersize',0.5,'color',col2{ss1});
end
errorbar(xplot(jj,ss1),nanmean(Strategy{ss1}),nansem(Strategy{ss1}),'capsize',3,'color',col2{ss1});
Means=[Means; nanmean(Strategy{ss1})];Sems=[Sems; nansem(Strategy{ss1})];
[hh pval1(ss1) ci stats]=ttest(Strategy{ss1});
P1=[P1; pval1(ss1)];Tstat1=[Tstat1; stats.tstat]; df1=[df1; stats.df];
if nanmean(Strategy{ss1}) > 0
  if pval1(ss1) < 0.05 & pval1(ss1) > 0.01
        t=text(xplot(jj,ss1), 40, '*','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.01 & pval1(ss1) > 0.001
        t=text(xplot(jj,ss1), 40, '**','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.001 
        t=text(xplot(jj,ss1), 40, '***','fontsize',6,'color',col2{ss1},'rotation',90);
  end
end

end
yl(jj,:)=get(gca,'ylim');



end

xlim([0, 8.5]);
ylim([min(min(yl))-40 55]);

for jj=4
 left=0.1; bot=0.75; width=0.15; height=0.15;    
%  ax1(jj)=axes('Position',[left bot width height]);
 perf1=[];T11=[];  
for ii=1:7   
p1=[];
Trials=[1:6];T1=[];
for kk=1:3
x11=performance_learn{ii}(blocks(kk,1):blocks(kk,2),NCDays{1,1});
x12=cellfun(@(x) length(x)-1,x11,'UniformOutput',false); 
x1=100-((cell2mat(x12)/5)'.*100);
if size(x1,1) > 1
p1=[p1, (x1)];
else
p1=[p1, (x1)];
end
T1=[T1, Trials];
end
perf1=[perf1; nanmean(p1)];
T11=[T11; (T1)];
hold on;
% plot([1:18]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',col1{ss,jj},'linestyle','none'); 
end
incl1=[1:2;7,8;13:14];incl2=[5:6;11,12;17:18];Strategy1=[];
for kk=1:3
Strategy1{2}(kk,:)=nanmean((perf1(:,incl2(kk,:))-perf1(:,incl1(kk,:)))');
Strategy1{1}(kk,:)=nanmean((perf1(:,[incl1(kk,1):incl2(kk,end)])')-50);
end
Strategy{1}=Strategy1{1}(:);
Strategy{2}=Strategy1{2}(:);
[hh pval2 ci stats]=ttest(Strategy{1},Strategy{2});
P2=[P2; pval2];Tstat2=[Tstat2; stats.tstat];df2=[df2; stats.df];
if pval2<0.05
hh1=sigstar(xplot(jj,[1,2]),pval2,0)
end
col2={rgb('red'),rgb('Black')};pval1=[];
for ss1=1:length(Strategy)
    hold on;
h1(xplot(jj,ss1))=bar(xplot(jj,ss1),nanmean(Strategy{ss1}));
h1(xplot(jj,ss1)).FaceColor='none'; h1(xplot(jj,ss1)).EdgeColor=col2{ss1};

for ss=1:size(Strategy{ss1},1)
    plot(xplot(jj,ss1)+randn(1).*0.1,Strategy{ss1}(ss),'marker','o','markersize',0.5,'color',col2{ss1});
end
errorbar(xplot(jj,ss1),nanmean(Strategy{ss1}),nansem(Strategy{ss1}),'capsize',3,'color',col2{ss1});
Means=[Means; nanmean(Strategy{ss1})];Sems=[Sems; nansem(Strategy{ss1})];
[hh pval1(ss1) ci stats]=ttest(Strategy{ss1});
P1=[P1; pval1(ss1)];Tstat1=[Tstat1; stats.tstat]; df1=[df1; stats.df];
if nanmean(Strategy{ss1}) > 0
  if pval1(ss1) < 0.05 & pval1(ss1) > 0.01
        t=text(xplot(jj,ss1), 50, '*','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.01 & pval1(ss1) > 0.001
        t=text(xplot(jj,ss1), 50, '**','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.001 
        t=text(xplot(jj,ss1), 50, '***','fontsize',6,'color',col2{ss1},'rotation',90);
  end
end

end
yl(jj,:)=get(gca,'ylim');


end

IntraGroup=[Means, Sems, df1, Tstat1, P1];
BetweenGroup=[df2, Tstat2, P2];

xlim([0 12]);
ylim([-40 70]);
set(gca,'xtick',[1.5,4.5,7.5,10.5],'xticklabel',...
    {'Early', 'Middle', 'Late','No cue'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'fontsize',5,'tickdir','out');
set(gca,'box','off','fontsize',5);
ylabel('Strategy Index','fontsize',5);

set(gcf,'paperunits','centimeters')
set(gcf,'papertype','a4')
set(gcf,'paperposition',[1,1,3.5,4])   

GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r600','StrategyIndex_Training_Performance'); 
print(gcf,'-painters','-depsc','-r600','StrategyIndex_Training_Performance'); 
close(fig1);


%% Panel 1h and J
clear all; close all;

cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BehavError.mat','Err','per_new1','per_old1')
c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
perf1=[perf1; nanmean(Err{ii}{1})];
hold on;
plot([1:3]+rand(1).*0.2, nanmean(Err{ii}{1}),'marker', 'o', 'markersize', 2,'color',rgb('Gray'),'linestyle','--'); 
end
errorbar([1:3],nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',c{6},'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)    
    [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; 
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{1}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),2.5);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 3.25, '*','fontsize',8);
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 3.25, '**','fontsize',8);
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 3.25, '***','fontsize',8);
    end
end
IntraGroup{1}=[G,df,tstat,p,means,sems]; 

pval=[];df=[];G=[];p=[];tstat=[];

newLim = ylim();
plot([0 51],[2.5 2.5],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
ylim([0 3.5]);
set(gca,'ytick',[0:4],'yticklabel',[0:4],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Baseline', 'No cues', 'Baseline'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Number of errors','fontsize',5);

% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.75,1])  

GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','No_Cue_Test_Behavior'); 
print(gcf,'-painters','-depsc','-r300','No_Cue_Test_Behavior'); 
close(fig1)

c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
perf1=[perf1; nanmean(Err{ii}{3})];
hold on;
plot([1:3]+rand(1).*0.2, nanmean(Err{ii}{3}),'marker', 'o', 'markersize', 2,'color',rgb('Gray'),'linestyle','--'); 
end
errorbar([1:3],nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',c{6},'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)
   [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat;
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{4}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),2.5);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 3.25, '*','fontsize',8);
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 3.25, '**','fontsize',8);
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 3.25, '***','fontsize',8);
    end
end
IntraGroup{4}=[G,df,tstat,p,means,sems]; 

newLim = ylim();
plot([0 51],[2.5 2.5],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
ylim([0 3.5]);
set(gca,'ytick',[0:4],'yticklabel',[0:4],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Baseline', 'No cues', 'Baseline'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Number of errors','fontsize',5);

% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.75,1])  

GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Blocking_Test_Behavior'); 
print(gcf,'-painters','-depsc','-r300','Blocking_Test_Behavior'); 
close(fig1)

sess1={1,[2:4],5; 1,2,[3:4];1,2,3;1,2,[3:4];1,2,[3:4]};

c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
p1=[];
for jj=1:3
x1=nanmean(nanmean(Err{ii}{2}(:,sess1{ii,jj})'));
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 2,'color',rgb('Gray'),'linestyle','--'); 
end
errorbar([1:3],nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',c{6},'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)
   [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat;
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{2}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),2.5);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 2.75, '*','fontsize',8);
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 2.75, '**','fontsize',8);
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 2.75, '***','fontsize',8);
    end
end
IntraGroup{2}=[G,df,tstat,p,means,sems]; 
newLim = ylim();
plot([0 51],[2.5 2.5],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
ylim([0 3]);
set(gca,'ytick',[0:4],'yticklabel',[0:4],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Baseline', 'Assn', 'Assd'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Number of errors','fontsize',5);

% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.75,1]) 

GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim1_Behavior'); 
print(gcf,'-painters','-depsc','-r300','Assim1_Behavior'); 
close(fig1)


sess1={6,[7:8],9; 5,6,[7:8];4,5,6;5,6,[7:8];5,6,[7]};

c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
p1=[];
for jj=1:3
x1=nanmean(nanmean(Err{ii}{2}(:,sess1{ii,jj})'));
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 2,'color',rgb('Gray'),'linestyle','--'); 
end
errorbar([1:3],nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',c{6},'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)
   [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat;
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{3}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),2.5);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 2.75, '*','fontsize',8);
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 2.75, '**','fontsize',8);
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 2.75, '***','fontsize',8);
    end
end
IntraGroup{3}=[G,df,tstat,p,means,sems];  

newLim = ylim();
plot([0 51],[2.5 2.5],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
ylim([0 3]);
set(gca,'ytick',[0:4],'yticklabel',[0:4],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Baseline', 'Assn', 'Assd'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Number of errors','fontsize',5);

% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.75,1])  

GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim2_Behavior'); 
print(gcf,'-painters','-depsc','-r300','Assim2_Behavior'); 
close(fig1)


sess1={1,[2],[3,4]; 1,[2],[3,4];1,2,3;1,[2],[3,4];1,[2],[3,4];};

c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
p1=[];
for jj=1:3
x1=nanmean(nanmean(Err{ii}{4}(:,sess1{ii,jj})));
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 2,'color',rgb('Gray'),'linestyle','--'); 
end
errorbar([1:3],nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',c{6},'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)
   [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat;
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{5}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),2.5);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 3.25, '*','fontsize',8);
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 3.25, '**','fontsize',8);
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 3.25, '***','fontsize',8);
    end
end
IntraGroup{5}=[G,df,tstat,p,means,sems]; 


newLim = ylim();
plot([0 51],[2.5 2.5],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
ylim([0 3.5]);
set(gca,'ytick',[0:4],'yticklabel',[0:4],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Baseline', 'Assn', 'Assd'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Number of errors','fontsize',5);

% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.75,1])  
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim3_Behavior'); 
print(gcf,'-painters','-depsc','-r300','Assim3_Behavior'); 
close(fig1)


%% Panel 1j V2 changed Unchanged Separated

%% New and old divided
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BehavError.mat','Err','per_new1','per_old1')
sess1={1,[2:4],5; 1,2,[3:4];1,2,3;1,2,[3:4];1,2,[3:4]};
c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
p1=[];per_new1{ii}(per_new1{ii}==-1)=NaN;
for jj=1:3
x1=nanmean(nanmean(1-(per_new1{ii}(:,sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',rgb('Red'),'linestyle','--'); 
end
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Red'),'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)    
    [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; 
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{2}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),50);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 45, '*','fontsize',6,'color','r');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 45, '**','fontsize',6,'color','r');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 45, '***','fontsize',6,'color','r');
    end
end
IntraGroup{1}=[G,df,tstat,p,means,sems]; 
c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];

for ii=1:5
p1=[];per_old1{ii}(per_old1{ii}==-1)=NaN;
for jj=1:3
x1=nanmean(nanmean(1-(per_old1{ii}(:,sess1{ii,jj})/5)').*100);
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',rgb('Blue'),'linestyle','--'); 
end
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Blue'),'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)    
    [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; 
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{2}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),50);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 50, '*','fontsize',6,'color','b');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 50, '**','fontsize',6,'color','b');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 50, '***','fontsize',6,'color','b');
    end
end
IntraGroup{2}=[G,df,tstat,p,means,sems]; 
pval=[];df=[];G=[];p=[];tstat=[];
ylim([40 100]);
newLim = ylim();
plot([0 4],[50 50],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
% % ylim([40 100]);
set(gca,'fontsize',5,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Ba-1', 'Asn', 'Asd'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('performance','fontsize',5);
Ryax = gca;
yyaxis right  
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Blue'),'capsize',3,'linewidth',1); hold on;  % Right axis handle
ylim(newLim);
Ryaxt = [40:10:100];                 % Get Left axis ticks
RyaxDegC = 5-(([40:10:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
% set(gca,'ydir','reverse')
% 
% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.85,1])  
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim1_Behavior_performanceV2'); 
print(gcf,'-painters','-depsc','-r300','Assim1_Behavior_performanceV2'); 
close(fig1)


sess1={6,[7:8],9; 5,6,[7:8];4,5,6;5,6,[7:8];5,6,[7]};


c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
p1=[];per_new1{ii}(per_new1{ii}==-1)=NaN;
for jj=1:3
x1=nanmean(nanmean(1-(per_new1{ii}(:,sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',rgb('Red'),'linestyle','--'); 
end
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Red'),'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)    
    [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; 
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{3}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),50);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 45, '*','fontsize',6,'color','r');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 45, '**','fontsize',6,'color','r');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 45, '***','fontsize',6,'color','r');
    end
end
IntraGroup{3}=[G,df,tstat,p,means,sems]; 
c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];

for ii=1:5
p1=[];per_old1{ii}(per_old1{ii}==-1)=NaN;
for jj=1:3
x1=nanmean(nanmean(1-(per_old1{ii}(:,sess1{ii,jj})/5)').*100);
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',rgb('Blue'),'linestyle','--'); 
end
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Blue'),'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)    
    [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; 
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{4}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),50);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 50, '*','fontsize',6,'color','b');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 50, '**','fontsize',6,'color','b');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 50, '***','fontsize',6,'color','b');
    end
end
IntraGroup{4}=[G,df,tstat,p,means,sems]; 
pval=[];df=[];G=[];p=[];tstat=[];
ylim([40 100]);
newLim = ylim();
plot([0 4],[50 50],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
% % ylim([40 100]);
set(gca,'fontsize',5,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Ba-1', 'Asn', 'Asd'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('performance','fontsize',5);
Ryax = gca;
yyaxis right  
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Blue'),'capsize',3,'linewidth',1); hold on;  % Right axis handle
ylim(newLim);
Ryaxt = [40:10:100];                 % Get Left axis ticks
RyaxDegC = 5-(([40:10:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
% set(gca,'ydir','reverse')
% 
% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.85,1])  
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim2_Behavior_performanceV2'); 
print(gcf,'-painters','-depsc','-r300','Assim2_Behavior_performanceV2'); 
close(fig1)


sess1={1,[2],[3,4]; 1,[2],[3,4];1,2,3;1,[2],[3,4];1,[2],[3,4];};

c={'r','m','b','g','c','k'};
c1=[rand rand rand];perf1=[];
fig1=figure
for ii=1:5
p1=[];
for jj=1:3
x1=nanmean(nanmean(1-(Err{ii}{4}(:,sess1{ii,jj})/5)').*100);
p1=[p1, x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',rgb('Red'),'linestyle','--'); 
end
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Red'),'capsize',3,'linewidth',1); hold on;  
Groups=[1,2;1,3;2,3];
pval=[];df=[];G=[];p=[];tstat=[];
for ii=1:size(Groups,1)    
    [h pval(ii) ci stat]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
    G(ii,:)=Groups(ii,:); p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; 
    if pval(ii) < 0.05
      h1=sigstar(Groups(ii,:),pval(ii),0);
    end
end
BetweenGroups{5}=[G,df,tstat,p]; 
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ii=1:3
    [h pval(ii) ci stat]=ttest(perf1(:,ii),50);
    G(ii,:)=ii; p(ii,1)=pval(ii); df(ii,1)=stat.df; tstat(ii,1)=stat.tstat; means(ii,1)=nanmean(perf1(:,ii)); sems(ii,1)=nansem(perf1(:,ii)); 
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii, 45, '*','fontsize',6);
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii, 45, '**','fontsize',6);
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii, 45, '***','fontsize',6);
    end
end
IntraGroup{5}=[G,df,tstat,p,means,sems]; 

pval=[];df=[];G=[];p=[];tstat=[];
ylim([35 90]);
newLim = ylim();
plot([0 4],[50 50],'--','color','k');   
xlim([0.5 3.5]);
% if k ~=3
% % ylim([40 100]);
set(gca,'fontsize',5,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:3],'xticklabel',...
    {'Ba-1', 'Asn', 'Asd'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('performance','fontsize',5);
Ryax = gca;
yyaxis right  
errorbar([1:3],(nanmean(perf1)),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color',rgb('Red'),'capsize',3,'linewidth',1); hold on;  % Right axis handle
ylim(newLim);
Ryaxt = [40:10:90];                 % Get Left axis ticks
RyaxDegC = 5-(([40:10:90]*5)./100);     

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,0.85,1])  
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim3_Behavior_PerformanceV2'); 
print(gcf,'-painters','-depsc','-r300','Assim3_Behavior_PerformanceV2'); 
close(fig1)

%% Panel 1k

clear all; close all;

cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BehavError.mat','Err','per_new1','per_old1')



%% Block Wise New & Old
sess1={1,[2:4],5; 1,2,[3:4];1,2,3;1,2,[3:4];1,2,[3:4]};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];
fig1=figure

c={'r','m','b','g','c','k'};
c1=[rand rand rand];

ss=[1,4,7,2,5,8,3,6,9];behav=[];
for ss1=1:3
for jj=1:3
    for kk=1:3
    behav{ss1,jj,kk}=[];
    end
end
end
for jj=1:3 
perf1=[];  
for ii=1:5    
p1=[];per_new1{ii}(per_new1{ii}==-1)=NaN;
for kk=1:3
x1=nanmean(nanmean(1-(per_new1{ii}(blocks(kk,1):blocks(kk,2),sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
behav{1,jj,kk}=[behav{1,jj,kk}; x1];
end
perf1=[perf1; p1];

hold on;
plot([1:3]+plott(jj)+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color','r','linestyle','none'); 
end

% errorbar([1:3]+plott(jj),nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','r','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii) CI Stats]=ttest(perf1(:,ii),50);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 15, '*','fontsize',6,'color','r');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 15, '**','fontsize',6,'color','r');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 15, '***','fontsize',6,'color','r');
    end
end


perf2=[];  
for ii=1:5    
p1=[];per_old1{ii}(per_old1{ii}==-1)=NaN;
for kk=1:3
x1=nanmean(nanmean(1-(per_old1{ii}(blocks(kk,1):blocks(kk,2),sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
behav{2,jj,kk}=[behav{2,jj,kk}; x1];
end
perf2=[perf2; p1];

hold on;
plot([1:3]+plott(jj)+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color','b','linestyle','none'); 
end

% errorbar([1:3]+plott(jj),nanmean(perf2),nansem(perf2),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','b','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii)]=ttest(perf2(:,ii),50);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 10, '*','fontsize',6,'color','b');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 10, '**','fontsize',6,'color','b');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 10, '***','fontsize',6,'color','b');
    end
end


perf=[perf1; perf2];
% errorbar([1:3]+plott(jj),nanmean(perf),nansem(perf),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','k','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii)]=ttest(perf(:,ii),50);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 5, '*','fontsize',6,'color','k');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 5, '**','fontsize',6,'color','k');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 5, '***','fontsize',6,'color','k');
    end
end


end


xplot=[1:9];
col1={rgb('Red'),rgb('Blue'),rgb('Black')};
for ss1=1
    behav1=[];
for jj=1:3
    for kk=1:3
        behav1=[behav1, [behav{1,jj,kk}; behav{2,jj,kk}]];
    end
end
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{3},'capsize',3,'linewidth',1); hold on;  
   
end
xplot=[1:9];
col1={rgb('Red'),rgb('Blue'),rgb('Black')};
for ss1=1:2
    behav1=[];
for jj=1:3
    for kk=1:3
        behav1=[behav1, behav{ss1,jj,kk}];
    end
end
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{ss1},'capsize',3,'linewidth',1); hold on;  
   
end
hold on
ylim([10 100]);
newLim = ylim();
plot([0 51],[50 50],'--','color','k');   
plot([3.5 3.5],[10 100],'--k');
plot([6.5 6.5],[10 100],'--k');

xlim([0.5 9.5]);
% if k ~=3

set(gca,'ytick',[10:10:100],'yticklabel',[10:10:100],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:9],'xticklabel',...
    {'Block 1', 'Block 2', 'Block 3'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Performance','fontsize',5);
Ryax = gca;
yyaxis right  
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{ss1},'capsize',3,'linewidth',1); hold on;  
ylim(newLim);
Ryaxt = [10:10:100];                 % Get Left axis ticks
RyaxDegC = 5-(([10:10:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,1.5,1])    
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim1_NewOld_Behavior_Performance'); 
print(gcf,'-painters','-depsc','-r300','Assim1_NewOld_Behavior_Performance'); 
close(fig1)

%% Day 3 Block Wise New & Old
sess1={6,[7:8],9; 5,6,[7:8];4,5,6;5,6,[7:8];5,6,[7]};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];
fig1=figure


c={'r','m','b','g','c','k'};
c1=[rand rand rand];

ss=[1,4,7,2,5,8,3,6,9];behav=[];
for ss1=1:3
for jj=1:3
    for kk=1:3
    behav{ss1,jj,kk}=[];
    end
end
end
for jj=1:3 
perf1=[];  
for ii=1:5    
p1=[];per_new1{ii}(per_new1{ii}==-1)=NaN;
for kk=1:3
x1=nanmean(nanmean(1-(per_new1{ii}(blocks(kk,1):blocks(kk,2),sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
behav{1,jj,kk}=[behav{1,jj,kk}; x1];
end
perf1=[perf1; p1];

hold on;
plot([1:3]+plott(jj)+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color','r','linestyle','none'); 
end

% errorbar([1:3]+plott(jj),nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','r','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii) CI Stats]=ttest(perf1(:,ii),50);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 15, '*','fontsize',6,'color','r');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 15, '**','fontsize',6,'color','r');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 15, '***','fontsize',6,'color','r');
    end
end


perf2=[];  
for ii=1:5    
p1=[];per_old1{ii}(per_old1{ii}==-1)=NaN;
for kk=1:3
x1=nanmean(nanmean(1-(per_old1{ii}(blocks(kk,1):blocks(kk,2),sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
behav{2,jj,kk}=[behav{2,jj,kk}; x1];
end
perf2=[perf2; p1];

hold on;
plot([1:3]+plott(jj)+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color','b','linestyle','none'); 
end

% errorbar([1:3]+plott(jj),nanmean(perf2),nansem(perf2),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','b','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii)]=ttest(perf2(:,ii),50);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 10, '*','fontsize',6,'color','b');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 10, '**','fontsize',6,'color','b');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 10, '***','fontsize',6,'color','b');
    end
end


perf=[perf1; perf2];
% errorbar([1:3]+plott(jj),nanmean(perf),nansem(perf),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','k','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii)]=ttest(perf(:,ii),50);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 5, '*','fontsize',6,'color','k');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 5, '**','fontsize',6,'color','k');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 5, '***','fontsize',6,'color','k');
    end
end


end


xplot=[1:9];
col1={rgb('Red'),rgb('Blue'),rgb('Black')};
for ss1=1
    behav1=[];
for jj=1:3
    for kk=1:3
        behav1=[behav1, [behav{1,jj,kk}; behav{2,jj,kk}]];
    end
end
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{3},'capsize',3,'linewidth',1); hold on;  
   
end
xplot=[1:9];
col1={rgb('Red'),rgb('Blue'),rgb('Black')};
for ss1=1:2
    behav1=[];
for jj=1:3
    for kk=1:3
        behav1=[behav1, behav{ss1,jj,kk}];
    end
end
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{ss1},'capsize',3,'linewidth',1); hold on;  
   
end
hold on
ylim([10 100]);
newLim = ylim();
plot([0 51],[50 50],'--','color','k');   
plot([3.5 3.5],[10 100],'--k');
plot([6.5 6.5],[10 100],'--k');

xlim([0.5 9.5]);
% if k ~=3

set(gca,'ytick',[10:10:100],'yticklabel',[10:10:100],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:9],'xticklabel',...
    {'Block 1', 'Block 2', 'Block 3'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Performance','fontsize',5);
Ryax = gca;
yyaxis right  
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{ss1},'capsize',3,'linewidth',1); hold on;  
ylim(newLim);
Ryaxt = [10:10:100];                 % Get Left axis ticks
RyaxDegC = 5-(([10:10:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
% 
% labels={'Rat 1','Rat 2','Rat 3','Rat 4','Rat 5','Average'};
% l1=legend(gca,labels)
% set(l1,'box','off','fontsize',4,'location','northeast')

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,1.5,1])   
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r300','Assim2_NewOld_Behavior_Performance'); 
print(gcf,'-painters','-depsc','-r300','Assim2_NewOld_Behavior_Performance'); 
close(fig1)


%% Block Wise Day 5
sess1={1,[2],[3,4]; 1,[2],[3,4];1,2,3;1,[2],[3,4];1,[2],[3,4];};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];
fig1=figure

c={'r','m','b','g','c','k'};
c1=[rand rand rand];
behav=[];
for ss1=1
for jj=1:3
    for kk=1:3
    behav{ss1,jj,kk}=[];
    end
end
end
ss=[1,4,7,2,5,8,3,6,9];
for jj=1:3 
 perf1=[];  
for ii=1:5    
p1=[];Err{ii}{4}(Err{ii}{4}==-1)=NaN;
for kk=1:3
x1=nanmean(nanmean(1-(Err{ii}{4}(blocks(kk,1):blocks(kk,2),sess1{ii,jj}))/5)'.*100);
p1=[p1, x1];
behav{1,jj,kk}=[behav{1,jj,kk}; x1];
end
perf1=[perf1; p1];
hold on;
plot([1:3]+plott(jj)+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color','r','linestyle','none'); 
end

% errorbar([1:3]+plott(jj),nanmean(perf1),nansem(perf1),'linestyle', '-', 'marker', 'o', 'markersize', 2,'color','r','capsize',3,'linewidth',1); hold on;  

% Groups=[1,2;1,3;2,3];pval=[];
% for ii=1:size(Groups,1)
%     [h pval(ii)]=ttest(perf1(:,Groups(ii,1)), perf1(:,Groups(ii,2)));
%     if pval(ii) < 0.05
%       h1=sigstar(Groups(ii,:),pval(ii),0);
%     end
% end


pval=[];
for ii=1:3
    [h pval(ii)]=ttest(perf1(:,ii),50,'Alpha',0.05);
    if pval(ii) < 0.05 & pval(ii) > 0.01
        text(ii+plott(jj), 15, '*','fontsize',6,'color','r');
    elseif pval(ii) < 0.01 & pval(ii) > 0.001
        text(ii+plott(jj), 15, '**','fontsize',6,'color','r');
     elseif pval(ii) < 0.001 & pval(ii) > 0
        text(ii+plott(jj), 15, '***','fontsize',6,'color','r');
    end
end

end
xplot=[1:9];
col1={rgb('Red'),rgb('Blue'),rgb('Black')};
for ss1=1
    behav1=[];
for jj=1:3
    for kk=1:3
        behav1=[behav1, [behav{1,jj,kk}]];
    end
end
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{1},'capsize',3,'linewidth',1); hold on;  
   
end
hold on
ylim([10 100]);
newLim = ylim();
plot([0 51],[50 50],'--','color','k');   
plot([3.5 3.5],[10 100],'--k');
plot([6.5 6.5],[10 100],'--k');

xlim([0.5 9.5]);
% if k ~=3

set(gca,'ytick',[10:10:100],'yticklabel',[10:10:100],'fontsize',4,'tickdir','out');
% title('No cue test','fontsize',6)
set(gca,'xtick',[1:9],'xticklabel',...
    {'Block 1', 'Block 2', 'Block 3'},...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'box','off','fontsize',5);
ylabel('Performance','fontsize',5);
Ryax = gca;
yyaxis right  
errorbar(xplot,nanmean(behav1),nansem(behav1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color',col1{ss1},'capsize',3,'linewidth',1); hold on;  
ylim(newLim);
Ryaxt = [10:10:100];                 % Get Left axis ticks
RyaxDegC = 5-(([10:10:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter') 
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,1.5,1])   
print(gcf,'-djpeg','-r300','Assim3_NewOld_Behavior_Performance'); 
print(gcf,'-painters','-depsc','-r300','Assim3_NewOld_Behavior_Performance'); 
close(fig1)




%% Panel l,m and n: Strategy Indices for learning days
%% for each block during experiment
clear all; close all;

cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BehavError.mat','Err','per_new1','per_old1')


%% Each Block Reset in Errors
sess1{1}={1,[2],[3]; 1,[2],[3];1,2,3;1,[2],[3];1,[2],[3]};
sess1{2}={[1,6],[2:4,7],[5,8,9]; [1,5],[2,6],[3,4,7,8];[1,4],[2,5],[3,6]; [1,5],[2,6],[3,4,7,8]; [1,5],[2,6],[3,4,7]};
sess1{3}={1,[2],[3]; 1,[2],[3];1,2,3;1,[2],[3];1,[2],[3]};
sess1{4}={1,[2],[3,4]; 1,[2],[3,4];1,2,3;1,[2],[3,4];1,[2],[3,4]};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];
fig1=figure
label1={'Ba 1','Nocue','Ba 2';'Ba 1','Assn','Assd';'Ba 1','Blkg','Ba 2';'Ba 1','Assn','Assd'};

col1={rgb('Red'),rgb('Black'),rgb('Red');...
    rgb('Red'),rgb('Blue'),rgb('Green');...
    rgb('Red'),rgb('Gray'),rgb('Red');...
    rgb('Red'),rgb('Blue'),rgb('Green')}
xplot=[1,2;4,5;7,8;10,11];nn1=[1,2,3,4];
IntraGroup=[];BetweenGroup=[];Means=[];P1=[];Tstat1=[];Sems=[];df1=[];P2=[];df2=[];Tstat2=[];M1=[]; S1=[];
for ss=1:4
     left=0.1; bot=0.80-((ss-1).*0.2); width=0.15; height=0.15;    
 ax1(ss)=axes('Position',[left bot width height]);
for jj=1:3 

 perf1=[];T11=[];  
for ii=1:5    
p1=[];Err{ii}{nn1(ss)}(Err{ii}{nn1(ss)}==-1)=NaN;
Trials=[1:6];T1=[];
for kk=1:3
x1=((Err{ii}{nn1(ss)}(blocks(kk,1):blocks(kk,2),sess1{ss}{ii,jj})'));
x1=100-((x1/5).*100);
if size(x1,1) > 1
p1=[p1, nanmean(x1)];
else
p1=[p1, (x1)];
end
T1=[T1, Trials];
end
perf1=[perf1; p1];
T11=[T11; T1];
hold on;
% plot([1:18]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',col1{ss,jj},'linestyle','none'); 
end
incl1=[1:2;7,8;13:14];incl2=[5:6;11,12;17:18];Strategy1=[];
for kk=1:3
Strategy1{2}(kk,:)=nanmean((perf1(:,incl2(kk,:))-perf1(:,incl1(kk,:)))');
Strategy1{1}(kk,:)=nanmean((perf1(:,[incl1(kk,1):incl2(kk,end)])')-50);
end
Strategy{1}=Strategy1{1}(:);
Strategy{2}=Strategy1{2}(:);
[hh pval2 ci stats]=ttest(Strategy{1},Strategy{2});
P2=[P2; pval2];Tstat2=[Tstat2; stats.tstat];df2=[df2; stats.df];
if pval2<0.05
hh1=sigstar(xplot(jj,[1,2]),pval2,0)
end
col2={rgb('red'),rgb('Black')};pval1=[];
for ss1=1:length(Strategy)
    hold on;
h1(xplot(jj,ss1))=bar(xplot(jj,ss1),nanmean(Strategy{ss1}));
h1(xplot(jj,ss1)).FaceColor='none'; h1(xplot(jj,ss1)).EdgeColor=col2{ss1};

for ss2=1:size(Strategy{ss1},1)
    plot(xplot(jj,ss1)+randn(1).*0.1,Strategy{ss1}(ss2),'marker','o','markersize',0.5,'color',col2{ss1});
end
errorbar(xplot(jj,ss1),nanmean(Strategy{ss1}),nansem(Strategy{ss1}),'capsize',3,'color',col2{ss1});
[hh pval1(ss1) ci stats]=ttest(Strategy{ss1});
P1=[P1; pval1(ss1)];Tstat1=[Tstat1; stats.tstat]; df1=[df1; stats.df];
M1=[M1; nanmean(Strategy{ss1})]; S1=[S1; nansem(Strategy{ss1})]; 
if nanmean(Strategy{ss1}) > 0
  if pval1(ss1) < 0.05 & pval1(ss1) > 0.01
        t=text(xplot(jj,ss1), 40, '*','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.01 & pval1(ss1) > 0.001
        t=text(xplot(jj,ss1), 40, '**','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.001 
        t=text(xplot(jj,ss1), 40, '***','fontsize',6,'color',col2{ss1},'rotation',90);
  end
end

end
yl(jj,:)=get(gca,'ylim');
end
ylim([-70 100]);
xlim([0 9.5]);
set(gca,'xtick',[1.5,4.5,7.5],'xticklabel',...
    label1(ss,:),...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'fontsize',5,'tickdir','out');
end


xlim([0, 9]);


set(gca,'fontsize',5,'tickdir','out');
set(gca,'box','off','fontsize',5);
ylabel('Strategy Index','fontsize',5);

IntraGroup=[M1, S1, df1, Tstat1, P1];
BetweenGroup=[df2, Tstat2, P2];

set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter');
set(gcf,'paperposition',[0.1,0.1,3,4]) ;   
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r600','StrategyIndex_Behavior_Performance'); 
print(gcf,'-painters','-depsc','-r600','StrategyIndex_Behavior_Performance'); 
close(fig1)



%% Each Block Reset in Errors
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BehavError.mat','Err','per_new1','per_old1')
sess1{1}={1,[2],[3]; 1,[2],[3];1,2,3;1,[2],[3];1,[2],[3]};
sess1{2}={[1,6],[2:4,7],[5,8,9]; [1,5],[2,6],[3,4,7,8];[1,4],[2,5],[3,6]; [1,5],[2,6],[3,4,7,8]; [1,5],[2,6],[3,4,7]};
sess1{3}={1,[2],[3]; 1,[2],[3];1,2,3;1,[2],[3];1,[2],[3]};
sess1{4}={1,[2],[3,4]; 1,[2],[3,4];1,2,3;1,[2],[3,4];1,[2],[3,4]};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];
fig1=figure
label1={'Ba 1','Nocue','Ba 2';'Ba 1','Assn','Assd';'Ba 1','Blkg','Ba 2';'Ba 1','Assn','Assd'};

col1={rgb('Red'),rgb('Black'),rgb('Red');...
    rgb('Red'),rgb('Blue'),rgb('Green');...
    rgb('Red'),rgb('Gray'),rgb('Red');...
    rgb('Red'),rgb('Blue'),rgb('Green')}
xplot=[1,2;4,5;7,8;10,11];nn1=[1,2,3,4];
IntraGroup=[];BetweenGroup=[];Means=[];P1=[];Tstat1=[];Sems=[];df1=[];P2=[];df2=[];Tstat2=[];
for ss=1:2
 if ss==1; 
    Err=per_new1;
else
   Err=per_old1; 
end 
     left=0.1; bot=0.80-((ss-1).*0.2); width=0.15; height=0.15;    
 ax1(ss)=axes('Position',[left bot width height]);
for jj=2:3 

 perf1=[];T11=[];  
for ii=1:5    
p1=[];Err{ii}(Err{ii}==-1)=NaN;
Trials=[1:6];T1=[];
for kk=1:3
x1=((Err{ii}(blocks(kk,1):blocks(kk,2),sess1{2}{ii,jj})'));
x1=100-((x1/5).*100);
if size(x1,1) > 1
p1=[p1, nanmean(x1)];
else
p1=[p1, (x1)];
end
T1=[T1, Trials];
end
perf1=[perf1; p1];
T11=[T11; T1];
hold on;
% plot([1:18]+rand(1).*0.2, (p1),'marker', 'o', 'markersize', 1,'color',col1{ss,jj},'linestyle','none'); 
end
incl1=[1:2;7,8;13:14];incl2=[5:6;11,12;17:18];Strategy1=[];
for kk=1:3
Strategy1{2}(kk,:)=nanmean((perf1(:,incl2(kk,:))-perf1(:,incl1(kk,:)))');
Strategy1{1}(kk,:)=nanmean((perf1(:,[incl1(kk,1):incl2(kk,end)])')-50);
end
Strategy{1}=Strategy1{1}(:);
Strategy{2}=Strategy1{2}(:);
[hh pval2 ci stats]=ttest(Strategy{1},Strategy{2});
P2=[P2; pval2];Tstat2=[Tstat2; stats.tstat];df2=[df2; stats.df];
if pval2<0.05
hh1=sigstar(xplot(jj,[1,2]),pval2,0)
end
col2={rgb('red'),rgb('Black')};pval1=[];
for ss1=1:length(Strategy)
    hold on;
h1(xplot(jj,ss1))=bar(xplot(jj,ss1),nanmean(Strategy{ss1}));
h1(xplot(jj,ss1)).FaceColor='none'; h1(xplot(jj,ss1)).EdgeColor=col2{ss1};

for ss2=1:size(Strategy{ss1},1)
    plot(xplot(jj,ss1)+randn(1).*0.1,Strategy{ss1}(ss2),'marker','o','markersize',0.5,'color',col2{ss1});
end
errorbar(xplot(jj,ss1),nanmean(Strategy{ss1}),nansem(Strategy{ss1}),'capsize',3,'color',col2{ss1});
[hh pval1(ss1) ci stats]=ttest(Strategy{ss1});
P1=[P1; pval1(ss1)];Tstat1=[Tstat1; stats.tstat]; df1=[df1; stats.df];
Means=[Means; nanmean(Strategy{ss1})]; Sems=[Sems; nansem(Strategy{ss1})]; 
if nanmean(Strategy{ss1}) > 0
  if pval1(ss1) < 0.05 & pval1(ss1) > 0.01
        t=text(xplot(jj,ss1), 50, '*','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.01 & pval1(ss1) > 0.001
        t=text(xplot(jj,ss1), 50, '**','fontsize',6,'color',col2{ss1},'rotation',90);
  elseif  pval1(ss1) < 0.001 
        t=text(xplot(jj,ss1), 50, '***','fontsize',6,'color',col2{ss1},'rotation',90);
  end
end

end
yl(jj,:)=get(gca,'ylim');
ylim([-70 100]);
xlim([3 9]);
end

set(gca,'xtick',[1.5,4.5,7.5],'xticklabel',...
    label1(2,:),...
    'xticklabelrotation',45,'fontsize',5)
set(gca,'fontsize',5,'tickdir','out');
set(gca,'box','off','fontsize',5);
ylabel('Strategy Index','fontsize',5);
end


IntraGroup=[Means, Sems, df1, Tstat1, P1];
BetweenGroup=[df2, Tstat2, P2];


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,2.5,4])    
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r600','StrategyIndex_NewOld_Performance'); 
print(gcf,'-painters','-depsc','-r600','StrategyIndex_NewOld_Performance'); 
close(fig1)




%% Panel 1o; simulation of n trials and effects on strategies in various types of sessions
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BehavError.mat','Err','per_new1','per_old1')

for ss1=1:5
    for ss2=1:4
        perf1{ss1,ss2}=[];
        T11{ss1,ss2}=[];
    end
end


%% Each Block Reset in Errors
sess1{1}={1,[2],[3]; 1,[2],[3];1,2,3;1,[2],[3];1,[2],[3]};
sess1{2}={[1,6],[2:4,7],[5,8,9]; [1,5],[2,6],[3,4,7,8];[1,4],[2,5],[3,6]; [1,5],[2,6],[3,4,7,8]; [1,5],[2,6],[3,4,7]};
sess1{3}={1,[2],[3]; 1,[2],[3];1,2,3;1,[2],[3];1,[2],[3]};
sess1{4}={1,[2],[3,4]; 1,[2],[3,4];1,2,3;1,[2],[3,4];1,[2],[3,4]};
blocks=[1,6;7,12;13,18];
plott=[0,3,6];

for ss=1:4
for jj=1:3 

 
for ii=1:5    
p1=[];Err{ii}{ss}(Err{ii}{ss}==-1)=NaN;
Trials=[1:6];T1=[];
for kk=1:3
x1=((Err{ii}{ss}(blocks(kk,1):blocks(kk,2),sess1{ss}{ii,jj})'));
x1=100-((x1/5).*100);
if size(x1,1) > 1
p1=[p1, nanmean(x1)];
else
p1=[p1, (x1)];
end
T1=[T1, Trials];
end
perf1{ss,jj}=[perf1{ss,jj}; p1];
T11{ss,jj}=[T11{ss,jj}; T1];

end


end
end


Perf{1}=[perf1{1,1}; perf1{2,1}; perf1{3,1}; perf1{4,1}; perf1{5,1}; perf1{1,3}; perf1{3,3}];
Perf{2}=[perf1{1,2}];
Perf{3}=[perf1{3,2}];
Perf{4}=[perf1{4,2}];
Perf{5}=[perf1{2,3}; perf1{4,3}];

Bins=[1:6; 7:12; 13:18];
Bins=flip(Bins,2);
for ss1=1:4
    for ii=1:6
TrialPerf{ss1,ii}=[];
    end
end
for ss1=1:5
for ii=1:size(Bins,2)
    xx1=Perf{ss1}(:,Bins(:,1:ii));
    TrialPerf{ss1,ii}=xx1(:);
end
end
fig1=figure
set(fig1,'defaultAxesColorOrder',[rgb('Black'); rgb('Black')]);
col1={rgb('Red'),rgb('Black'),rgb('Gray'),...
    rgb('Blue'),rgb('Magenta'),rgb('Orange')};
left=0.1; bot=0.70; width=0.2; height=0.2;
ax1(1)=axes('position',[left bot width height]);
% 
for ss1=1:5
xx1=[];    
for ii=2:6    
    xx1{ii-1}=TrialPerf{ss1,ii};   
end
    plot([1:size(xx1,2)],cellfun(@nanmean,xx1),'marker','o','markersize',2,'linewidth',1,'linestyle','-','color',col1{ss1});
    pp=patch([[1:size(xx1,2)] fliplr([1:size(xx1,2)])], ...
    [cellfun(@nanmean,xx1)+cellfun(@nansem,xx1) fliplr(cellfun(@nanmean,xx1)-cellfun(@nansem,xx1))], col1{ss1})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on
    hold on;    
   
end

plot([1:5],ones(1,5).*50,'marker','o','markersize',2,'linewidth',1,'linestyle','-','color',col1{6});
RewardHistErr=[0,0.25,0.5,0.75,1,1.25];RewardHist=[];
for ss1=2:length(RewardHistErr)
    RewardHist(ss1-1)=100-(RewardHistErr(ss1)/(ss1-1)).*100;
end
% plot([0 6],[75 75],'linestyle','--','marker','none','color',col1{6});

RewardHistErr=[90,80,70,60,50];
plot(1:5,RewardHistErr,'marker','o','markersize',2,'linewidth',1,'linestyle','--','color','c');
xlim([0.5 5.5])
set(gca,'box','off','tickdir','out','fontsize',5);
xlabel('Number of trials in the task','fontsize',5)
ylabel('Performance','fontsize',5)
newLim=get(gca,'ylim');
set(gca,'xtick',[1:5],'xticklabel',[2:6],'fontsize',5);
ylim([40 100]);
set(gca,'fontsize',5,'tickdir','out');

% comparison of strategies
left=0.4; bot=0.70; width=0.2; height=0.2;
ax1(2)=axes('position',[left bot width height]);
% 
for ss1=1:5
    for ss2=1:3
        xx1{ss1,ss2}=[];
    end
end
xbins=[0,1,0,2,0,3]
for ss1=1:5
    
for ii=[2,4,6]   
    xx1{ss1,xbins(ii)}=TrialPerf{ss1,ii};   
end
       
end
RewardHistErr=[90,70,50];
xplot=[1,2,3,4,5;7,8,9,10,11;13,14,15,16,17];
for ss1=1:5
    
for ii=1:3   
    h1(ss1,ii)=bar(xplot(ii,ss1),nanmean(xx1{ss1,ii})); hold on
    h1(ss1,ii).FaceColor=col1{ss1};   h1(ss1,ii).EdgeColor='none';
    errorbar(xplot(ii,ss1),nanmean(xx1{ss1,ii}),nansem(xx1{ss1,ii}),'capsize',1,'color',col1{ss1});
    
    [h11, pval]=ttest(xx1{ss1,ii},50);
  if pval < 0.05 & pval > 0.01
        t=text(xplot(ii,ss1)-0.2, 95, '*','fontsize',6,'color',col1{6});
  elseif  pval < 0.01 & pval > 0.001
        t=text(xplot(ii,ss1)-0.2, 95, '**','fontsize',6,'color',col1{6});
  elseif  pval < 0.001 
        t=text(xplot(ii,ss1)-0.2, 95, '***','fontsize',6,'color',col1{6});
  end
  t.Rotation=270;
%      [h11, pval]=ttest(xx1{ss1,ii},75);
%   if pval < 0.05 & pval > 0.01
%         t=text(xplot(ii,ss1)-0.2, 105, '*','fontsize',6,'color',col1{6});
%   elseif  pval < 0.01 & pval > 0.001
%         t=text(xplot(ii,ss1)-0.2, 105, '**','fontsize',6,'color',col1{6});
%   elseif  pval < 0.001 
%         t=text(xplot(ii,ss1)-0.2, 105, '***','fontsize',6,'color',col1{6});
%   end
%   t.Rotation=270;

  [h11, pval]=ttest(xx1{ss1,ii},RewardHistErr(ii));
  if pval < 0.05 & pval > 0.01
        t=text(xplot(ii,ss1)-0.2, 105, '*','fontsize',6,'color','c');
  elseif  pval < 0.01 & pval > 0.001
        t=text(xplot(ii,ss1)-0.2, 105, '**','fontsize',6,'color','c');
  elseif  pval < 0.001 
        t=text(xplot(ii,ss1)-0.2, 105, '***','fontsize',6,'color','c');
  end
  t.Rotation=270;
  ylim([40 110]);
end
       
end

plot([3.5,9.5,15.5],ones(1,3).*50,'marker','o','markersize',2,'linewidth',1,'linestyle','--','color',col1{6});

RewardHistErr=[90,70,50];
plot([3.5,9.5,15.5],RewardHistErr,'marker','o','markersize',2,'linewidth',1,'linestyle','--','color','c');
set(gca,'box','off','tickdir','out','fontsize',5);
xlabel('Number of trials in the task','fontsize',5)
ylabel('Performance','fontsize',5)
newLim=get(gca,'ylim');
set(gca,'xtick',[3,9,15],'xticklabel',[2,4,6],'fontsize',5);
set(gca,'fontsize',5,'tickdir','out');

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,4,4])    
GenPath=strcat(paths{2},'Results/PlotFigures/Fig1/Figures');
cd(GenPath)
print(gcf,'-djpeg','-r600','RewardHistoryWorkingMemoryFunctionTrialNumber'); 
print(gcf,'-painters','-depsc','-r600','RewardHistoryWorkingMemoryFunctionTrialNumber'); 
close(fig1)








