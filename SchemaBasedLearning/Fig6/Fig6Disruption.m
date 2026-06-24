%% Panel a
clear all
xxlim=[-1 37];
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('Paths')
GenPath=strcat(paths{2},'Results/PlotFigures/Database');
cd(GenPath)

load TrainingData

performance{1}=Chosen_arms10;
performance{2}=Chosen_arms11;
performance{3}=Chosen_arms12;
performance{4}=Chosen_arms13;



performance13=[];
for i=1:4
    startvisits=cellfun(@(x) sum(x==3 | x==6),performance{i},'UniformOutput',false); 
    performance13=[performance13; nanmean(cell2mat(startvisits))];
    
end
% performance13(5:7,36)=NaN;
c={'r','g','b','m','r','g','b','m'};
c1='k';

f1=figure
errorbar(1:size(performance13,2), nanmean(performance13),nansem(performance13),'linestyle', '-', 'marker', 'none', 'markersize', 2,'color','k','capsize',2); hold on;

% plot([0 51],[1.75 1.75],'--','color',c1);


% for k=1:7
% errorbar(1:size(performance{k},2), nanmean(performance {k}),nansem(performance{k}),'linestyle', 'none', 'marker', 'o', 'markersize', 1,'color',c{k},'capsize',2); hold on;
% end
% plot([0 51],[2.5 2.5],'--','color',c1);
% errorbar(1:size(performance12,2), nanmean(performance12),nansem(performance12),'linewidth', 1, 'marker', 'o', 'markersize', 3,'color','m');hold on;
set(gca,'box','off','fontsize',8);
xlabel('Training Day','fontsize',8);
ylabel('visits per trial','fontsize',8);
xlim(xxlim+1);
set(gca,'xtick',[0 10 20 30 ],'xticklabel',[0 10 20 30 36],'fontsize',6,'tickdir','out')
set(gca,'ytick',[0:0.2:1],'yticklabel',[0:0.2:1],'fontsize',6,'tickdir','out');
ylim([0 1.1]);
% if k ~=3
% end
% title(strcat(['Average']),'fontsize',8)
GenPath=strcat(paths{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath)
set(gcf,'paperunits','centimeters')
set(gcf,'papertype','a4')
set(gcf,'paperposition',[1,1,4,4]) 
print('-djpeg','-r300','ErrorsVisitingStartArmsRD');
print(f1,'-painters','-depsc','-r300','-r300','ErrorsVisitingStartArmsRD');
close(f1);



f1=figure

performance13=[];
for i=1:4
    startvisits=cellfun(@(x) sum(x==3 | x==6),performance{i},'UniformOutput',false); 
    performance13=[performance13; nanmean(cell2mat(startvisits))]
end
% performance13(5:7,36)=NaN;
performance12=[];
for i=1:4
    startvisits=cellfun(@(x) length(x)-1,performance{i},'UniformOutput',false); 
    performance12=[performance12; nanmean(cell2mat(startvisits))];
end

ErrPerf=performance12(:,[31,33,35]);
ErrPerf1=ErrPerf(:);
ErrPerf=(100-((ErrPerf1)/5.*100));
% performance12(5:7,36)=NaN;
yyaxis left
x1=100-((performance12/5).*100);
errorbar(1:size(performance12,2), nanmean(x1),nansem(x1),'linestyle', '-', 'marker', 'o', 'markersize', 1,'color','k','capsize',2); hold on;
ddt_mean=ones(1,size(performance12,2)).*nanmean(ErrPerf);
ddt_sem=ones(1,size(performance12,2)).*nansem(ErrPerf);
plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color','k'); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], rgb('Gray'))
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
Lyax =gca; 
Lyact= get(Lyax, 'YTick');     
p=[];stats=[];
for a=1:size(performance12,2)
    [h(a) p(a) CI stats{a}] = ttest2(performance12(:,a),ErrPerf1);
    yl=get(gca,'ylim'); 
    if p(a) < 0.05 & p(a)>0.01 & nanmean(performance12(:,a)) < nanmean(ErrPerf1)
    t= text(a-0.3,max(yl),'*','fontsize',5,'color','r');
      t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 & nanmean(performance12(:,a)) < nanmean(ErrPerf1)
     t=text(a-0.3,max(yl),'**','fontsize',5,'color','r');  
     t.Rotation=270;
    elseif p(a) < 0.001 & nanmean(performance12(:,a)) < nanmean(ErrPerf1)
     t=text(a-0.3,max(yl),'***','fontsize',5,'color','r');  
     t.Rotation=270;
    end
end

Lyax.YColor = 'k';

set(gca,'box','off','fontsize',5);
xlabel('Training Day','fontsize',5);
ylabel('Performance','fontsize',5);
xlim([0 37]);
set(gca,'xtick',[0 10 20 30 36],'xticklabel',[0 10 20 30 36],'fontsize',5,'tickdir','out')
set(gca,'ytick',[0:20:100],'fontsize',5,'tickdir','out');
ylim([0 100]);

Ryax = gca;
yyaxis right  
ylim([0 100]);
Ryaxt = [0:20:100];                 % Get Left axis ticks
RyaxDegC = 5-(([0:20:100]*5)./100);     % Convert to °C (or whatever you want)
Ryax.YColor = 'k';
% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
GenPath=strcat(paths{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath)
set(gcf,'paperunits','centimeters')
set(gcf,'papertype','a4')
set(gcf,'paperposition',[1,1,5,4]) 
print(f1,'-djpeg','-r300','Behavior_details_DynamicErrorTogethorRD');
print(f1,'-painters','-depsc','-r300','-r300','Behavior_details_DynamicErrorTogethorRD');
close(f1);



%% Panel e and f

clear all; close all;
%% Test2

clear all

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
files= importdata('filesDisruption');
FolderNumber=[1,4,6,11,13,15,18,21,23];

sd=4;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
col1={rgb('red'),rgb('blue'),rgb('Black'),rgb('Gray')};
MultiSum1=[];MultiAct3=[];MultiAct4=[];

for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
MultiRates=[];
dn1=pwd;
load RipplePowerSpectrum3ZV2.mat
if ii==21 
    sess1=[3]; 
else
    sess1=[3,4];
end
for is=sess1
  
   if any(is==[3,4])

   sample1=MultiRates{1,is}.Hist;
sample2=MultiRates{1,is}.SumHist;
    [spectrum1, f] = pwelch(sample1, hamming(40), 10, 64, 1000);
    [spectrum2, f] = pwelch(sample2, hamming(40), 10, 64, 1000);
    spectrum1=nanmean(spectrum1')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum1(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum1 = log_pxx - log_pxx_fit;
  
    spectrum1(spectrum1==inf | spectrum1==-inf)=NaN;
%  
     spectrum2=nanmean(spectrum2')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum2(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum2 = log_pxx - log_pxx_fit;
    
    
    
    
   MultiAct3=[MultiAct3; ((spectrum1'))];   
   HistMulti=(MultiRates{1, is}.SumHist(10:end,:)); 
   MultiAct4=[MultiAct4; (spectrum2')];
% plot(f,zscore(nanmean(spectrum1)),'-r'); hold on; plot(f,zscore(nanmean(spectrum2)),'-k');
   end
   
end
end

  
fig1=figure;   
left=0.1; bottom=0.5; width=0.35; height=0.4;
ax1(1)=axes('Position',[left bottom width height])  
plot(f(2:end),(nanmean(MultiAct3)),'-r'); 
hold on; plot(f(2:end),(nanmean(MultiAct4)),'-k');
ddt_sem1=nansem(MultiAct3);
ddt_mean1=nanmean(MultiAct3);
ddt_sem2=nansem(MultiAct4);
ddt_mean2=nanmean(MultiAct4);
pp1=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem1+ddt_mean1 fliplr(ddt_mean1-ddt_sem1)], col1{1})
pp2=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem2+ddt_mean2 fliplr(ddt_mean2-ddt_sem2)], col1{2})
    pp1.EdgeColor = 'red';
    alpha(0.4); hold on ;
    pp2.EdgeColor = 'blue';
    alpha(0.4); hold on  
    yl=get(gca,'ylim');
xlim([50 400])
plot([140 140],[min(yl) 0.5],'--k');
plot([250 250],[min(yl) 0.5],'--k');
set(gca, 'XScale', 'log');
set(gca,'tickdir','out','box','off','fontsize',5);
set(gca, 'xtick',[50,100,200, 400])
xlabel('frequency (Hz)','fontsize',5);
ylabel('log(Power) whitened','fontsize',5);

mmm1{1}=ddt_mean1([5:26])';
mmm1{2}=ddt_mean2([5:26])';
mmm2{1}=ddt_sem1([5:26])';
mmm2{2}=ddt_sem2([5:26])';
mmm3=f([5:26]);


f1=f(2:end);p=[];
ind=find(f1>30 & f1 < 35);

for a=ind:size(MultiAct3,2)-10
    [p(a) h]=signrank(MultiAct3(:,a), MultiAct4(:,a));
       if p(a) < 0.05 & p(a)>0.01 
    t= text(f1(a)-0.1,max(yl)-0.2,'*','fontsize',5,'color','r');
      t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 
     t=text(f1(a)-0.1,max(yl)-0.2,'**','fontsize',5,'color','r');  
     t.Rotation=270;
    elseif p(a) < 0.001 
     t=text(f1(a)-0.1,max(yl)-0.2,'***','fontsize',5,'color','r');  
     t.Rotation=270;
    end
end
ylim([min(yl) max(yl)-0.1])


FolderNumber=[3,9,10,17,20,22];

sd=4;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
col1={rgb('red'),rgb('blue'),rgb('Black'),rgb('Gray')};
MultiSum1=[];MultiAct3=[];MultiAct4=[];

for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
MultiRates=[];
dn1=pwd;
load RipplePowerSpectrum3ZV2.mat
if ii==22 
    sess1=[3]; 
else
    sess1=[3,4];
end
for is=sess1
  
   if any(is==[3,4])

   sample1=MultiRates{1,is}.Hist;
   sample2=MultiRates{1,is}.SumHist;
   sample1=sample1(:,~any(isnan(sample1)));
   sample2=sample2(:,~any(isnan(sample2)));
    [spectrum1, f] = pwelch(sample1, hamming(40), 15, 64, 1000);
    [spectrum2, f] = pwelch(sample2, hamming(40), 15, 64, 1000);
    spectrum1=nanmean(spectrum1')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum1(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum1 = log_pxx - log_pxx_fit;
  
    spectrum1(spectrum1==inf | spectrum1==-inf)=NaN;
%  
     spectrum2=nanmean(spectrum2')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum2(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum2 = log_pxx - log_pxx_fit;
    
    
    
    
   MultiAct3=[MultiAct3; ((spectrum1'))];   
   HistMulti=(MultiRates{1, is}.SumHist(10:end,:)); 
   MultiAct4=[MultiAct4; (spectrum2')];
% plot(f,zscore(nanmean(spectrum1)),'-r'); hold on; plot(f,zscore(nanmean(spectrum2)),'-k');
   end
   
end
end

  
 
left=0.6; bottom=0.5; width=0.35; height=0.4;

ax1(2)=axes('Position',[left bottom width height])  
plot(f(2:end),(nanmean(MultiAct3)),'-r'); 
hold on; plot(f(2:end),(nanmean(MultiAct4)),'-k');
ddt_sem1=nansem(MultiAct3);
ddt_mean1=nanmean(MultiAct3);
ddt_sem2=nansem(MultiAct4);
ddt_mean2=nanmean(MultiAct4);
pp1=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem1+ddt_mean1 fliplr(ddt_mean1-ddt_sem1)], col1{1})
pp2=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem2+ddt_mean2 fliplr(ddt_mean2-ddt_sem2)], col1{2})
    pp1.EdgeColor = 'red';
    alpha(0.4); hold on ;
    pp2.EdgeColor = 'blue';
    alpha(0.4); hold on 
      yl=get(gca,'ylim');
xlim([50 400])
plot([140 140],[min(yl) 0.5],'--k');
plot([250 250],[min(yl) 0.5],'--k');
set(gca, 'XScale', 'log');
set(gca,'tickdir','out','box','off','fontsize',5);
set(gca, 'xtick',[50,100,200, 400])
xlabel('frequency (Hz)','fontsize',5);
ylabel('log(Power) whitened','fontsize',5);

mmm1{1}=ddt_mean1([5:26])';
mmm1{2}=ddt_mean2([5:26])';
mmm2{1}=ddt_sem1([5:26])';
mmm2{2}=ddt_sem2([5:26])';
mmm3=f([5:26]);
f1=f(2:end);p=[];
ind=find(f1>30 & f1 < 35);

for a=ind:size(MultiAct3,2)-10
    [p(a) h]=signrank(MultiAct3(:,a), MultiAct4(:,a));
       if p(a) < 0.05 & p(a)>0.01 
    t= text(f1(a)-0.1,max(yl)-0.2,'*','fontsize',5,'color','r');
      t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 
     t=text(f1(a)-0.1,max(yl)-0.2,'**','fontsize',5,'color','r');  
     t.Rotation=270;
    elseif p(a) < 0.001 
     t=text(f1(a)-0.1,max(yl)-0.2,'***','fontsize',5,'color','r');  
     t.Rotation=270;
    end
end
ylim([min(yl) max(yl)-0.1])
  set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter') 
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,3,2])   
print(gcf,'-djpeg','-r300','RipplePowerSpectrumNorm3ZV2InclR11'); 
print(gcf,'-painters','-depsc','-r300','RipplePowerSpectrumNorm3ZV2InclR11'); 
close(fig1)


%% zscored power non normalized
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
files= importdata('filesDisruption');
FolderNumber=[1,4,6,11,13,15,18,21,23];

sd=4;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
col1={rgb('red'),rgb('blue'),rgb('Black'),rgb('Gray')};
MultiSum1=[];MultiAct3=[];MultiAct4=[];

for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
MultiRates=[];
dn1=pwd;
load RipplePowerSpectrum3ZV2.mat
if ii==21 
    sess1=[3]; 
else
    sess1=[3,4];
end
for is=sess1
  
   if any(is==[3,4])

   sample1=MultiRates{1,is}.Hist;
sample2=MultiRates{1,is}.SumHist;
nfft = 2^nextpow2(40);
    [spectrum1, f] = pwelch(sample1, hamming(40), 10, 64, 1000);
    [spectrum2, f] = pwelch(sample2, hamming(40), 10, 64, 1000);
    spectrum1=nanmean(spectrum1')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum1(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum1 = (10*log_pxx);
  
    spectrum1(spectrum1==inf | spectrum1==-inf)=NaN;
%  
     spectrum2=nanmean(spectrum2')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum2(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum2 = (10*log_pxx);
    
    
    
    
     
   MultiAct3=[MultiAct3; zscore((spectrum1'))];   

   MultiAct4=[MultiAct4; zscore(spectrum2')];
% plot(f,zscore(nanmean(spectrum1)),'-r'); hold on; plot(f,zscore(nanmean(spectrum2)),'-k');
   end
   
end
end

  
fig1=figure;   
left=0.1; bottom=0.5; width=0.35; height=0.4;
ax1(1)=axes('Position',[left bottom width height])  
plot(f(2:end),(nanmean(MultiAct3)),'-r'); 
hold on; plot(f(2:end),(nanmean(MultiAct4)),'-k');
ddt_sem1=nansem(MultiAct3);
ddt_mean1=nanmean(MultiAct3);
ddt_sem2=nansem(MultiAct4);
ddt_mean2=nanmean(MultiAct4);
pp1=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem1+ddt_mean1 fliplr(ddt_mean1-ddt_sem1)], col1{1})
pp2=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem2+ddt_mean2 fliplr(ddt_mean2-ddt_sem2)], col1{2})
    pp1.EdgeColor = 'red';
    alpha(0.4); hold on ;
    pp2.EdgeColor = 'blue';
    alpha(0.4); hold on  
    yl=get(gca,'ylim');
xlim([50 400])
plot([140 140],[min(yl) max(yl)-3],'--k');
plot([250 250],[min(yl) max(yl)-3],'--k');
set(gca, 'XScale', 'log');
set(gca,'tickdir','out','box','off','fontsize',5);
set(gca, 'xtick',[50,100,200, 400])
xlabel('frequency (Hz)','fontsize',5);
ylabel('log(Power) zscored','fontsize',5);
mmm1{1}=ddt_mean1([5:26])';
mmm1{2}=ddt_mean2([5:26])';
mmm2{1}=ddt_sem1([5:26])';
mmm2{2}=ddt_sem2([5:26])';
mmm3=f([5:26]);
f1=f(2:end);p=[];
ind=find(f1>30 & f1 < 35);

for a=ind:size(MultiAct3,2)-10
    [p(a) h]=signrank(MultiAct3(:,a), MultiAct4(:,a));
       if p(a) < 0.05 & p(a)>0.01 
    t= text(f1(a)-0.1,max(yl)-1.5,'*','fontsize',5,'color','r');
      t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 
     t=text(f1(a)-0.1,max(yl)-1.5,'**','fontsize',5,'color','r');  
     t.Rotation=270;
    elseif p(a) < 0.001 
     t=text(f1(a)-0.1,max(yl)-1.5,'***','fontsize',5,'color','r');  
     t.Rotation=270;
    end
end
ylim([min(yl) max(yl)-1.5])


FolderNumber=[3,9,10,17,20,22];

sd=4;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
col1={rgb('red'),rgb('blue'),rgb('Black'),rgb('Gray')};
MultiSum1=[];MultiAct3=[];MultiAct4=[];

for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
MultiRates=[];
dn1=pwd;
load RipplePowerSpectrum3ZV2.mat
if ii==22 
    sess1=[3]; 
else
    sess1=[3,4];
end
for is=sess1
  
   if any(is==[3,4])

   sample1=MultiRates{1,is}.Hist;
   sample2=MultiRates{1,is}.SumHist;
   sample1=sample1(:,~any(isnan(sample1)));
   sample2=sample2(:,~any(isnan(sample2)));
   nfft = 2^nextpow2(40);
   [spectrum1, f] = pwelch(sample1, hamming(40), 15, 64, 1000);
    [spectrum2, f] = pwelch(sample2, hamming(40), 15, 64, 1000);
    spectrum1=nanmean(spectrum1')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum1(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum1 = 10*log_pxx;
  
    spectrum1(spectrum1==inf | spectrum1==-inf)=NaN;
%  
     spectrum2=nanmean(spectrum2')';
    log_f = log10(f(f > 0));  % Exclude 0 Hz
log_pxx = log10(spectrum2(f > 0));
    % Linear regression
p = polyfit(log_f, log_pxx, 1);  % p(1): slope, p(2): intercept

% Predicted 1/f trend
log_pxx_fit = polyval(p, log_f);

% Subtract the fit (this is the whitening)
spectrum2 = 10*log_pxx;
    
    
    
    
   MultiAct3=[MultiAct3; zscore((spectrum1'))];   

   MultiAct4=[MultiAct4; zscore(spectrum2')];
% plot(f,zscore(nanmean(spectrum1)),'-r'); hold on; plot(f,zscore(nanmean(spectrum2)),'-k');
   end
   
end
end

  
 
left=0.6; bottom=0.5; width=0.35; height=0.4;
ax1(1)=axes('Position',[left bottom width height])  
plot(f(2:end),(nanmean(MultiAct3)),'-r'); 
hold on; plot(f(2:end),(nanmean(MultiAct4)),'-k');
ddt_sem1=nansem(MultiAct3);
ddt_mean1=nanmean(MultiAct3);
ddt_sem2=nansem(MultiAct4);
ddt_mean2=nanmean(MultiAct4);
pp1=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem1+ddt_mean1 fliplr(ddt_mean1-ddt_sem1)], col1{1})
pp2=patch([[f(2:end)'] fliplr([f(2:end)'])], ...
       [ddt_sem2+ddt_mean2 fliplr(ddt_mean2-ddt_sem2)], col1{2})
    pp1.EdgeColor = 'red';
    alpha(0.4); hold on ;
    pp2.EdgeColor = 'blue';
    alpha(0.4); hold on  
    yl=get(gca,'ylim');
xlim([50 400])
plot([140 140],[min(yl) max(yl)],'--k');
plot([250 250],[min(yl) max(yl)],'--k');
set(gca, 'XScale', 'log');
set(gca,'tickdir','out','box','off','fontsize',5);
set(gca, 'xtick',[50,100,200, 400])
xlabel('frequency (Hz)','fontsize',5);
ylabel('log(Power) z-scored','fontsize',5);
mmm1{1}=ddt_mean1([5:26])';
mmm1{2}=ddt_mean2([5:26])';
mmm2{1}=ddt_sem1([5:26])';
mmm2{2}=ddt_sem2([5:26])';
mmm3=f([5:26]);
f1=f(2:end);p=[];
ind=find(f1>30 & f1 < 35);

for a=ind:size(MultiAct3,2)-10
    [p(a) h]=signrank(MultiAct3(:,a), MultiAct4(:,a));
       if p(a) < 0.05 & p(a)>0.01 
    t= text(f1(a)-0.1,max(yl)-1.5,'*','fontsize',5,'color','r');
      t.Rotation=270;
    elseif p(a) < 0.01 & p(a)>0.001 
     t=text(f1(a)-0.1,max(yl)-1.5,'**','fontsize',5,'color','r');  
     t.Rotation=270;
    elseif p(a) < 0.001 
     t=text(f1(a)-0.1,max(yl)-1.5,'***','fontsize',5,'color','r');  
     t.Rotation=270;
    end
end
ylim([min(yl) max(yl)-1.5]);
  set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter'); 
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath);
set(gcf,'paperposition',[0.1,0.1,3,2]);   
print(gcf,'-djpeg','-r300','RipplePowerSpectrumZscored3ZV2InclR11'); 
print(gcf,'-painters','-depsc','-r300','RipplePowerSpectrumZscored3ZV2InclR11'); 
close(fig1)

%% Spectrograms

%% zscored power non normalized
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
files= importdata('filesDisruption');
FolderNumber=[1,4,6,11,13,15];ss2=0;
fig1=figure
for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
MultiRates=[];
dn1=pwd;
load ('ripplespectrogramV2.mat')
allave=[];rr1=0; ff = (wlff);
for nn1=3:4
  wltt1=-0.15:0.001:0.15;
 wltt = -0.1:0.001:0.1;
%  aveawl=allwl{nn1}-allwl{nn1}(:,:,1);
 aveawl=allwl{nn1};
 avewl=nanmean(aveawl.^2,1);
 avewl=(squeeze(avewl(1,:,:)));
 avewl=(avewl(ff>100 & ff < 350,(wltt1>= -0.11 & wltt1 <0.0905)));
 rr1=rr1+1;
 allave(rr1,:,:)=avewl;
end
 avewl=nanmean(allave,1);
 avewl=(squeeze(avewl(1,:,:)));
 clim1=[nanmin(nanmin(avewl(6:12,1:100))) 5*10^5];
%         for ss1=1:size(avewl,1)
%             avewl(ss1,:)=zscore(avewl(ss1,:));
%         end
%         for ss1=1:size(avewl,2)
%             avewl(:,ss1)=zscore(avewl(:,ss1));
%         end

avewl(:,[102:126])=NaN;
avewl=avewl;
        ff=ff(ff>100 & ff < 350);
        ss2=ss2+1;
        subplot(3,2,ss2)
        pcolor(wltt,ff,(avewl))
        colorbar
%                 caxis([nanmin(nanmin(avewl)) nanmean(nanmean(avewl))])
%         caxis([nanmin(nanmin(avewl)) nanmean(nanmean(avewl))-(nanmean(nanmean(avewl)))/2])
          caxis(clim1)
%         cmap=getNCLColmap('MPL_jet.rgb',128)
%         colormap(gca,cmap);

        ytk = [100,150,200,250,300];
        ytklabels = cell(1,length(ytk));
        for ll = 1:length(ytk) 
            ytklabels{ll} = num2str(ytk(ll));
        end
       shading interp;
%         colormap(ax,cmap)
        % %             relat = [min(laptime)-min(ts) max(laptime)-min(ts)];
%         ylim(gca,([20 300]))
        yticks(gca,(ytk))
        yticklabels(gca,ytklabels)
        hold on;
        plot(zeros(size(ff)),ff,'-k');
        plot(wltt,ones(size(wltt)).*140,'-k');
        plot(wltt,ones(size(wltt)).*250,'-k');
        set(gca,'box','off','tickdir','out','fontsize',5);
        xlabel('Time (s)','fontsize',5);
         ylabel('Frequency (Hz)','fontsize',5);
%         set(gca, 'YScale', 'log')
end
set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter') ;
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath);
set(gcf,'paperposition',[0.1,0.1,5,5]);   
print(gcf,'-djpeg','-r300','DisruptionSpectrogramAllSess'); 
print(gcf,'-painters','-depsc','-r300','DisruptionSpectrogramAllSess'); 
close(fig1)



FolderNumber=[3,9,10,17];ss2=0;
fig1=figure
for ii=FolderNumber
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
MultiRates=[];
dn1=pwd;
load ('ripplespectrogramV2.mat')
allave=[];rr1=0; ff = (wlff);
for nn1=1:2
wltt1=-0.15:0.001:0.15;
 wltt = -0.1:0.001:0.1;
 aveawl=allwl{nn1};
 avewl=nanmean(aveawl.^2,1);
 avewl=(squeeze(avewl(1,:,:)));
 avewl=(avewl(ff>100 & ff < 350,(wltt1>= -0.11 & wltt1 <0.0905)));
 rr1=rr1+1;
 allave(rr1,:,:)=avewl;
end
 avewl=nanmean(allave,1);
 avewl=(squeeze(avewl(1,:,:)));
%         for ss1=1:size(avewl,1)
%             avewl(ss1,:)=zscore(avewl(ss1,:));
%         end
%         for ss1=1:size(avewl,2)
%             avewl(:,ss1)=zscore(avewl(:,ss1));
%         end


        ff=ff(ff>100 & ff < 350);
        ss2=ss2+1;
        subplot(3,2,ss2)
        pcolor(wltt,ff,(avewl))
        colorbar
%                 caxis([nanmin(nanmin(avewl)) nanmean(nanmean(avewl))])
       clim1=[nanmin(nanmin(avewl(6:12,1:100))) nanmean(nanmax(avewl))+(nanmean(nanmax(avewl)))/2];
       caxis(clim1);
%         cmap=getNCLColmap('MPL_jet.rgb',128)
%         colormap(gca,cmap);
%         caxis([-2 2])
        ytk = [100,150,200,250,300];
        ytklabels = cell(1,length(ytk));
        for ll = 1:length(ytk) 
            ytklabels{ll} = num2str(ytk(ll));
        end
       shading interp;
%         colormap(ax,cmap)
        % %             relat = [min(laptime)-min(ts) max(laptime)-min(ts)];
%         ylim(gca,([20 300]))
        yticks(gca,(ytk))
        yticklabels(gca,ytklabels)
        hold on;
        plot(zeros(size(ff)),ff,'-k');
        plot(wltt,ones(size(wltt)).*140,'-k');
        plot(wltt,ones(size(wltt)).*250,'-k');
        set(gca,'box','off','tickdir','out','fontsize',5);
        xlabel('Time (s)','fontsize',5);
         ylabel('Frequency (Hz)','fontsize',5);
%         set(gca, 'YScale', 'log')
end
set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter') ;
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,5,5]) ; 
print(gcf,'-djpeg','-r300','CtrlDisruptionSpectrogramAllSess'); 
print(gcf,'-painters','-depsc','-r300','CtrlDisruptionSpectrogramAllSess'); 
close(fig1)


%% panels i and j

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
files= importdata('filesDisruption');
load RatsDisruptionBehavior 'DisruptionData1' 'DisruptionData2'

load TrainingData

performance{1}=Chosen_arms10;
performance{2}=Chosen_arms11;
performance{3}=Chosen_arms12;
performance{4}=Chosen_arms13;

performance12=[];
for i=1:4
    startvisits=cellfun(@(x) length(x)-1,performance{i},'UniformOutput',false); 
    performance12=[performance12; nanmean(cell2mat(startvisits))];
end

ErrPerf=performance12(:,[31,33,35]);
ErrPerf1=ErrPerf(:);
ErrPerf=(100-((ErrPerf1)/5.*100));

col1={rgb('red'),rgb('blue'),rgb('black')};
Sess=[1,2];

Perf=[];Target=[];


for ii=1:6
    perf1=cellfun(@(x) length(x)-1,DisruptionData1{ii},'uniformoutput',false);
    Perf{ii}=cell2mat(perf1);
    Target{ii}=cell2mat(cellfun(@(x) x(end),DisruptionData1{ii},'uniformoutput',false));
end
for ii=1:5
    perf1=cellfun(@(x) length(x)-1,DisruptionData2{ii},'uniformoutput',false);
    Perf{6+ii}=cell2mat(perf1);
    Target{6+ii}=cell2mat(cellfun(@(x) x(end),DisruptionData2{ii},'uniformoutput',false));
end  

NewOldInfo={[1,5,8],[2,4,7];[2,4,7],[1,5,8];[1,2,4,5,7,8],[];[2,5,7],[1,4,8];[1,4,8],[2,5,7];[2,5,7],[1,4,8];...
    [1,5,8],[2,4,7];[2,4,7],[1,5,8];[1,2,4,5,7,8],[];[1,4,8],[2,5,7];[2,5,7],[1,4,8]};

SessInfo={[1,3,4,7,9,10],[2,6,8,11]};


labels1={'Pre-Ex 1','Pre-Ex 2','Pre-Ex 3',...
    'Ex 1','Ex 2','Ex 3',...
    'Ex 4','Ex 5','Ex 6',...
    'Ex 7','Ex 8','Ex 9',...
    'Ex 10','Ex 11','Ex 12'};

labels2={'Pre-Ex 1','Pre-Ex 2','Pre-Ex 3',...
    'Ex 1','Ex 2','Ex 3',...
    'Ex 4','Ex 5','Ex 6',...
    'Ex 7','Ex 8','Ex 9',...
    'Ex 10','Ex 11','Ex 12'};
perfnew=[];perfold=[];
for ii=1:2
    
 perfnew{ii}=[];perfold{ii}=[];
end
for sess=1:size(SessInfo,2)
rr1=0;    
for ii=SessInfo{sess}
    rr1=rr1+1;
    if ~any(ii==3 | ii==9)
    NewTrialIds=(Target{ii}==NewOldInfo{ii,1}(1) | Target{ii}==NewOldInfo{ii,1}(2) | Target{ii}==NewOldInfo{ii,1}(3));
    else
    NewTrialIds=logical(ones(size(Perf{ii})));
    end
  
   pp1=[];pp2=[];
    for jj=1:size(Perf{ii},2)
        perfnew1=Perf{1,ii}(NewTrialIds(:,jj),jj);
        perfold1=Perf{1,ii}(~NewTrialIds(:,jj),jj); 
        pp1(:,jj)=((1-(perfnew1/5)).*100);
        pp2(:,jj)=((1-( perfold1/5)).*100);
    end
    
    if size(pp1,1)==9
    pp11=reshape(pp1,[3],[]);
    pp22=reshape(pp2,[3],[]);
    else
        pp11=reshape(pp1,[6],[]);
        pp22=[];
    end
    
        perfnew{sess}=[perfnew{sess}; nanmean(pp11)];
        if ~isempty(perfold1)
        perfold{sess}=[perfold{sess}; nanmean(pp22)];
        end
    end
end

fig1=figure  
xlim1=[0.5,15.5;0.5,12.5];
for ii=1:2
     if any(ii==[1,3,4])
    left=0.1; bottom=0.85-(ii-1).*0.15; width=0.5; height=0.1;
     else
    left=0.1; bottom=0.85-(ii-1).*0.15; width=0.42; height=0.1;
     end
    ax1(ii)=axes('Position',[left bottom width height])
     for ss1=1:size(perfnew{ii},2)
        plot(ss1+randn(size(perfnew{ii},1),1).*0.1,perfnew{ii}(:,ss1),'marker','o','markersize',1,'linestyle','none','color',col1{1}); hold on;
    end
       for ss1=1:size(perfold{ii},2)
        plot([ss1+randn(size(perfold{ii},1),1).*0.1]+0.125,perfold{ii}(:,ss1),'marker','o','markersize',1,'linestyle','none','color',col1{2}); hold on;
    end
        hh1(ii)=errorbar(1:size(perfnew{ii},2),nanmean(perfnew{ii}),nansem(perfnew{ii}),'color',col1{1},'capsize',2,'marker','o','markersize',2); 
         mmm1{ii,1}=(nanmean(perfnew{ii}))';
 mmm2{ii,1}=(nansem(perfnew{ii}))';
        hold on;        
        hh2(ii)=errorbar([1:size(perfold{ii},2)]+0.125,nanmean(perfold{ii}),nansem(perfold{ii}),'color',col1{2},'capsize',2,'marker','o','markersize',2);
       mmm1{ii,2}=(nanmean(perfold{ii}))';
 mmm2{ii,2}=(nansem(perfold{ii}))';
        ddt_mean=ones(1,size(perfnew{ii},2)).*nanmean(ErrPerf);
ddt_sem=ones(1,size(perfnew{ii},2)).*nansem(ErrPerf);
plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color','k'); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], rgb('Gray'))
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
        if any(ii==[1,3,4])
        labels=labels1;    
        text(12.6, 20,'1 day later','fontsize',4);
        text(9.25, 20,{ 'SWR', 'disruption'},'fontsize',4);
        text(6.25, 20,{'SWR', 'disruption'},'fontsize',4);
        plot([12.5 12.5],[0 100],'--k');
        plot([3.5 3.5],[0 100],'--k');
        plot([6.5 6.5],[0 100],'--k');
        plot([9.5 9.5],[0 100],'--k');
        elseif any(ii==[2,6])
        labels=labels2;    
        text(9.25, 20,{'control', 'SWR', 'disruption'},'fontsize',4);
        text(6.25, 20,{'control', 'SWR', 'disruption'},'fontsize',4); 
        plot([3.5 3.5],[0 100],'--k');
        plot([6.5 6.5],[0 100],'--k');
        plot([9.5 9.5],[0 100],'--k');
        else
        labels=labels2;    
        text(3.25, 20,{'SWR', 'disruption'},'fontsize',4);      
        end
         set(gca,'box','off','tickdir','out','xlim',[0 size(perfnew{ii},2)+1],'ylim',[0 100],...
            'xtick',[1:size(perfnew{ii},2)],'xticklabel',labels,'xticklabelrotation',45,'fontsize',5);
%         plot([0 size(perfnew{ii},2)+1],[50 50],'--k');
        newLim=get(gca,'ylim');
        Ryax = gca;
yyaxis right  
ylim(newLim);
Ryaxt = [0:20:100];                 % Get Left axis ticks
RyaxDegC = 5-(([0:20:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
ax1(ii).YColor='k';
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ss1=1:size(perfnew{ii},2)
    [h pval(ss1) ci stat]=ttest(perfnew{ii}(:,ss1),nanmean(ErrPerf));
    G(ss1,:)=ss1; p(ss1,1)=pval(ss1); df(ss1,1)=stat.df; tstat(ss1,1)=stat.tstat; 
    if pval(ss1) < 0.05 & pval(ss1) > 0.01 & nanmean(perfnew{ii}(:,ss1)) > 50
        text(ss1-0.2, 5, '*','fontsize',6,'color',col1{1}) ;
    elseif pval(ss1) < 0.01 & pval(ss1) > 0.001 & nanmean(perfnew{ii}(:,ss1)) > 50
        text(ss1-0.2, 5, '**','fontsize',6 ,'color',col1{1}) ;
     elseif pval(ss1) < 0.001 & pval(ss1) > 0 & nanmean(perfnew{ii}(:,ss1)) > 50
        text(ss1-0.2, 5, '***','fontsize',6 ,'color',col1{1}) ;
         else
         text(ss1-0.2, 5, 'n.s.','fontsize',4 ,'color',col1{1}) ;
    end
end
IntraGroup{1}=[G,df,tstat,p,means,sems]; 


pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ss1=1:size(perfold{ii},2)
    [h pval(ss1) ci stat]=ttest(perfold{ii}(:,ss1),nanmean(ErrPerf));
    G(ss1,:)=ss1; p(ss1,1)=pval(ss1); df(ss1,1)=stat.df; tstat(ss1,1)=stat.tstat; 
    if pval(ss1) < 0.05 & pval(ss1) > 0.01 & nanmean(perfold{ii}(:,ss1)) > 50
        text(ss1-0.2, 15, '*','fontsize',6,'color',col1{2}) ;
    elseif pval(ss1) < 0.01 & pval(ss1) > 0.001 & nanmean(perfold{ii}(:,ss1)) > 50
        text(ss1-0.2, 15, '**','fontsize',6 ,'color',col1{2}) ;
     elseif pval(ss1) < 0.001 & pval(ss1) > 0 & nanmean(perfold{ii}(:,ss1)) > 50
        text(ss1-0.2, 15, '***','fontsize',6 ,'color',col1{2}) ;
    else
         text(ss1-0.2, 15, 'n.s.','fontsize',4 ,'color',col1{2}) ;
    end
end
IntraGroup{2}=[G,df,tstat,p,means,sems]; 
xlim(xlim1(ii,:));
end
set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter') ;
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,4,6]) ;  
print(gcf,'-djpeg','-r300','DisruptionPerformanceBlockWiseChance'); 
print(gcf,'-painters','-depsc','-r300','DisruptionPerformanceBlockWiseChance'); 
close(fig1)

%% Averaged Data

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load RatsDisruptionBehavior 'DisruptionData1' 'DisruptionData2' 'DisruptionData3'
load TrainingData

performance{1}=Chosen_arms10;
performance{2}=Chosen_arms11;
performance{3}=Chosen_arms12;
performance{4}=Chosen_arms13;

performance12=[];
for i=1:4
    startvisits=cellfun(@(x) length(x)-1,performance{i},'UniformOutput',false); 
    performance12=[performance12; nanmean(cell2mat(startvisits))];
end

ErrPerf=performance12(:,[31,33,35]);
ErrPerf1=ErrPerf(:);
ErrPerf=(100-((ErrPerf1)/5.*100));
col1={rgb('red'),rgb('blue'),rgb('black')};
Sess=[1,2];

Perf=[];Target=[];


for ii=1:6
    perf1=cellfun(@(x) length(x)-1,DisruptionData1{ii},'uniformoutput',false);
    Perf{ii}=cell2mat(perf1);
    Target{ii}=cell2mat(cellfun(@(x) x(end),DisruptionData1{ii},'uniformoutput',false));
end
for ii=1:5
    perf1=cellfun(@(x) length(x)-1,DisruptionData2{ii},'uniformoutput',false);
    Perf{6+ii}=cell2mat(perf1);
    Target{6+ii}=cell2mat(cellfun(@(x) x(end),DisruptionData2{ii},'uniformoutput',false));
end  
% MissingTargets=[4,1,5,8,2,7,5,1,4,2,7,8,2,8,7,5,4,1]';
% for ii=1:5
%     perf1=cellfun(@(x) length(x)-1,DisruptionData3{ii},'uniformoutput',false);
%     Perf{ii+11}=cell2mat(perf1); 
%     for ss1=1:size(DisruptionData3{ii},1)
%         for ss2=1:size(DisruptionData3{ii},2)
%             if ~isempty(DisruptionData3{ii}{ss1,ss2})
%           Target{ii+11}(ss1,ss2)=  DisruptionData3{ii}{ss1,ss2}(end);
%             else
%            Target{ii+11}(ss1,ss2)=MissingTargets(ss1,1);
%             end
%         end
%     end
% %     Target{ii}=cell2mat(cellfun(@(x) x(end),DisruptionData3{ii},'uniformoutput',false));
% end


NewOldInfo={[1,5,8],[2,4,7];[2,4,7],[1,5,8];[1,2,4,5,7,8],[];[2,5,7],[1,4,8];[1,4,8],[2,5,7];[2,5,7],[1,4,8];...
    [1,5,8],[2,4,7];[2,4,7],[1,5,8];[1,2,4,5,7,8],[];[1,4,8],[2,5,7];[2,5,7],[1,4,8]};

SessInfo={[1,3,4,7,9,10],[2,6,8,11]};



labels1={'Pre-Ex-sess','Ex-sess 1','Ex-sess 2',...
    'Ex-sess 3','Ex-sess 4'};
labels2={'Pre-Ex-sess','Ex-sess 1','Ex-sess 2',...
    'Ex-sess 3','Ex-sess 4'};
perfnew=[];perfold=[];
for ii=1:2
    
 perfnew{ii}=[];perfold{ii}=[];
end
for sess=1:size(SessInfo,2)
rr1=0;    
for ii=SessInfo{sess}
    rr1=rr1+1;
    if ~any(ii==3 | ii==9)
    NewTrialIds=(Target{ii}==NewOldInfo{ii,1}(1) | Target{ii}==NewOldInfo{ii,1}(2) | Target{ii}==NewOldInfo{ii,1}(3));
    else
    NewTrialIds=logical(ones(size(Perf{ii})));
    end
  
   pp1=[];pp2=[];
    for jj=1:size(Perf{ii},2)
        perfnew1=Perf{1,ii}(NewTrialIds(:,jj),jj);
        perfold1=Perf{1,ii}(~NewTrialIds(:,jj),jj); 
        pp1(:,jj)=((1-(perfnew1/5)).*100);
        pp2(:,jj)=((1-( perfold1/5)).*100);
    end
    
    if size(pp1,1)==9
    pp11=reshape(pp1,[3],[]);
    pp22=reshape(pp2,[3],[]);
    else
        pp11=reshape(pp1,[6],[]);
        pp22=[];
    end
    pp111=[];rr1=0;
   for ss1=1:3:size(pp11,2)
       rr1=rr1+1;
       pp111(:,rr1)=nanmean(pp11(:,[ss1:ss1+2]));
   end
     pp222=[];rr1=0;
     if ~isempty(pp22)
   for ss1=1:3:size(pp11,2)
       rr1=rr1+1;
       pp222(:,rr1)=nanmean(pp22(:,[ss1:ss1+2]));
   end
     end    
        perfnew{sess}=[perfnew{sess}; (pp111)];
        if ~isempty(perfold1)
        perfold{sess}=[perfold{sess}; (pp222)];
        end
    end
end

fig1=figure  
xlim1=[0.5,5.5;0.5,4.5];
for ii=1:2
     if any(ii==[1,3,4])
    left=0.1; bottom=0.85-(ii-1).*0.15; width=0.5; height=0.1;
     else
    left=0.1; bottom=0.85-(ii-1).*0.15; width=0.42; height=0.1;
     end
    ax1(ii)=axes('Position',[left bottom width height])
    for ss1=1:size(perfnew{ii},2)
        plot(ss1+randn(size(perfnew{ii},1),1).*0.1,perfnew{ii}(:,ss1),'marker','o','markersize',1,'linestyle','none','color',col1{1}); hold on;
    end
       for ss1=1:size(perfold{ii},2)
        plot([ss1+randn(size(perfold{ii},1),1).*0.1]+0.125,perfold{ii}(:,ss1),'marker','o','markersize',1,'linestyle','none','color',col1{2}); hold on;
    end
        hh1(ii)=errorbar(1:size(perfnew{ii},2),nanmean(perfnew{ii}),nansem(perfnew{ii}),'color',col1{1},'capsize',2,'marker','o','markersize',2); 
 mmm1{ii,1}=(nanmean(perfnew{ii}))';
 mmm2{ii,1}=(nansem(perfnew{ii}))';
        
        hold on;        
        hh2(ii)=errorbar([1:size(perfold{ii},2)]+0.125,nanmean(perfold{ii}),nansem(perfold{ii}),'color',col1{2},'capsize',2,'marker','o','markersize',2);
 mmm1{ii,2}=(nanmean(perfold{ii}))';
 mmm2{ii,2}=(nansem(perfold{ii}))';      
ddt_mean=ones(1,size(perfnew{ii},2)).*nanmean(ErrPerf);
ddt_sem=ones(1,size(perfnew{ii},2)).*nansem(ErrPerf);
plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color','k'); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], rgb('Gray'))
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
 
        
        if any(ii==[1,3,4])
        labels=labels1;    
        text(4.6, 20,'1 day later','fontsize',4);
        text(3.25, 20,{ 'SWR', 'disruption'},'fontsize',4);
        text(2.25, 20,{'SWR', 'disruption'},'fontsize',4);
        plot([4.5 4.5],[0 100],'--k');
        plot([1.5 1.5],[0 100],'--k');
        elseif any(ii==[2,6])
        labels=labels2;    
        text(3.25, 20,{'control', 'SWR', 'disruption'},'fontsize',4);
        text(2.25, 20,{'control', 'SWR', 'disruption'},'fontsize',4);   
        else
        labels=labels2;    
        text(1.25, 20,{'SWR', 'disruption'},'fontsize',4);      
        end
         set(gca,'box','off','tickdir','out','xlim',[0 size(perfnew{ii},2)+1],'ylim',[0 100],...
            'xtick',[1:size(perfnew{ii},2)],'xticklabel',labels,'xticklabelrotation',45,'fontsize',5);
%         plot([0 size(perfnew{ii},2)+1],[50 50],'--k');
        newLim=get(gca,'ylim');
        Ryax = gca;
yyaxis right  
ylim(newLim);
Ryaxt = [0:20:100];                 % Get Left axis ticks
RyaxDegC = 5-(([0:20:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
ax1(ii).YColor='k'
pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ss1=1:size(perfnew{ii},2)
    [h pval(ss1) ci stat]=ttest(perfnew{ii}(:,ss1),nanmean(ErrPerf));
    G(ss1,:)=ss1; p(ss1,1)=pval(ss1); df(ss1,1)=stat.df; tstat(ss1,1)=stat.tstat; 
    if pval(ss1) < 0.025 & pval(ss1) > 0.01 & nanmean(perfnew{ii}(:,ss1)) > 50
        text(ss1-0.2, 5, '*','fontsize',6,'color',col1{1}) ;
    elseif pval(ss1) < 0.01 & pval(ss1) > 0.001 & nanmean(perfnew{ii}(:,ss1)) > 50
        text(ss1-0.2, 5, '**','fontsize',6 ,'color',col1{1}) ;
     elseif pval(ss1) < 0.001 & pval(ss1) > 0 & nanmean(perfnew{ii}(:,ss1)) > 50
        text(ss1-0.2, 5, '***','fontsize',6 ,'color',col1{1}) ;
         else
         text(ss1-0.2, 5, 'n.s.','fontsize',4 ,'color',col1{1}) ;
    end
end
IntraGroup{1}=[G,df,tstat,p,means,sems]; 


pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
for ss1=1:size(perfold{ii},2)
    [h pval(ss1) ci stat]=ttest(perfold{ii}(:,ss1),nanmean(ErrPerf));
    G(ss1,:)=ss1; p(ss1,1)=pval(ss1); df(ss1,1)=stat.df; tstat(ss1,1)=stat.tstat; 
    if pval(ss1) < 0.05 & pval(ss1) > 0.01 & nanmean(perfold{ii}(:,ss1)) > 50
        text(ss1-0.2, 15, '*','fontsize',6,'color',col1{2}) ;
    elseif pval(ss1) < 0.01 & pval(ss1) > 0.001 & nanmean(perfold{ii}(:,ss1)) > 50
        text(ss1-0.2, 15, '**','fontsize',6 ,'color',col1{2}) ;
     elseif pval(ss1) < 0.001  & nanmean(perfold{ii}(:,ss1)) > 50
        text(ss1-0.2, 15, '***','fontsize',6 ,'color',col1{2}) ;
    else
         text(ss1-0.2, 15, 'n.s.','fontsize',4 ,'color',col1{2}) ;
    end
end
IntraGroup{2}=[G,df,tstat,p,means,sems]; 
xlim(xlim1(ii,:));
end
set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter') ;
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,3,5])   
print(gcf,'-djpeg','-r300','DisruptionPerformanceOverallChance'); 
print(gcf,'-painters','-depsc','-r300','DisruptionPerformanceOverallChance'); 
close(fig1)


%% Block transition panel k



%% Averaged Data Blockwise Performance

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load RatsDisruptionBehavior 'DisruptionData1' 'DisruptionData2'

col1={rgb('red'),rgb('blue'),rgb('black')};
Sess=[1,2];

Perf=[];Target=[];


for ii=1:6
    perf1=cellfun(@(x) length(x)-1,DisruptionData1{ii},'uniformoutput',false);
    Perf{ii}=cell2mat(perf1);
    Target{ii}=cell2mat(cellfun(@(x) x(end),DisruptionData1{ii},'uniformoutput',false));
end
for ii=1:5
    perf1=cellfun(@(x) length(x)-1,DisruptionData2{ii},'uniformoutput',false);
    Perf{6+ii}=cell2mat(perf1);
    Target{6+ii}=cell2mat(cellfun(@(x) x(end),DisruptionData2{ii},'uniformoutput',false));
end  

NewOldInfo={[1,5,8],[2,4,7];[2,4,7],[1,5,8];[1,2,4,5,7,8],[];[2,5,7],[1,4,8];[1,4,8],[2,5,7];[2,5,7],[1,4,8];...
    [1,5,8],[2,4,7];[2,4,7],[1,5,8];[1,2,4,5,7,8],[];[1,4,8],[2,5,7];[2,5,7],[1,4,8]};

SessInfo={[1,3,4,7,9,10],[2,6,8,11]};


labels1={'Pre-Ex 1','Pre-Ex 2','Pre-Ex 3',...
    'Ex 1','Ex 2','Ex 3',...
    'Ex 4','Ex 5','Ex 6',...
    'Ex 7','Ex 8','Ex 9',...
    'Ex 10','Ex 11','Ex 12'};

labels2={'Pre-Ex 1','Pre-Ex 2','Pre-Ex 3',...
    'Ex 1','Ex 2','Ex 3',...
    'Ex 4','Ex 5','Ex 6',...
    'Ex 7','Ex 8','Ex 9',...
    'Ex 10','Ex 11','Ex 12'};
perfnew=[];perfold=[];
for ii=1:2
    
 perfnew{ii}=[];perfold{ii}=[];
end
for sess=1:size(SessInfo,2)
rr1=0;    
for ii=SessInfo{sess}
    rr1=rr1+1;
    if ~any(ii==3 | ii==9)
    NewTrialIds=(Target{ii}==NewOldInfo{ii,1}(1) | Target{ii}==NewOldInfo{ii,1}(2) | Target{ii}==NewOldInfo{ii,1}(3));
    else
    NewTrialIds=logical(ones(size(Perf{ii})));
    end
  
   pp1=[];pp2=[];
    for jj=1:size(Perf{ii},2)
        perfnew1=Perf{1,ii}(NewTrialIds(:,jj),jj);
        perfold1=Perf{1,ii}(~NewTrialIds(:,jj),jj); 
        pp1(:,jj)=((1-(perfnew1/5)).*100);
        pp2(:,jj)=((1-( perfold1/5)).*100);
    end
    
    if size(pp1,1)==9
    pp11=reshape(pp1,[3],[]);
    pp22=reshape(pp2,[3],[]);
    else
        pp11=reshape(pp1,[6],[]);
        pp22=[];
    end
    
        perfnew{sess}=[perfnew{sess}; nanmean(pp11)];
        if ~isempty(perfold1)
        perfold{sess}=[perfold{sess}; nanmean(pp22)];
        end
    end
end


BlockCrossNew{1}=perfnew{1}(:,[6,7,9,10]);
BlockCrossNew{2}=perfnew{2}(:,[6,7,9,10]);
BlockCrossOld{1}=perfold{1}(:,[6,7,9,10]);
BlockCrossOld{2}=perfold{2}(:,[6,7,9,10]);

col2={rgb('red'),rgb('blue'),rgb('red'),rgb('blue')};
fig1=figure  
xlim1=[0.5,15.5;0.5,12.5];
for ii=1:2
    
    left=0.1; bottom=0.85-(ii-1).*0.15; width=0.4; height=0.1;
    
    ax1(ii)=axes('Position',[left bottom width height])
       for ss1=1:size(BlockCrossNew{ii},1)
        plot([1,2],[BlockCrossNew{ii}(ss1,1) BlockCrossNew{ii}(ss1,2)],'marker','o','markersize',1,'linestyle','--','color',rgb('red')); hold on;
       end
        for ss1=1:size(BlockCrossNew{ii},1)
        plot([3,4],[BlockCrossNew{ii}(ss1,3) BlockCrossNew{ii}(ss1,4)],'marker','o','markersize',1,'linestyle','--','color',rgb('red')); hold on;
       end
        hh1(ii)=errorbar(1:2,nanmean(BlockCrossNew{ii}(:,1:2)),nansem(BlockCrossNew{ii}(:,1:2)),'color',col1{1},'capsize',2,'marker','o','markersize',2); 
        hold on;        
        hh1(ii)=errorbar(3:4,nanmean(BlockCrossNew{ii}(:,3:4)),nansem(BlockCrossNew{ii}(:,3:4)),'color',col1{1},'capsize',2,'marker','o','markersize',2); 
        
       for ss1=1:size(BlockCrossOld{ii},1)
        plot([1,2]+0.25,[BlockCrossOld{ii}(ss1,1) BlockCrossOld{ii}(ss1,2)],'marker','o','markersize',1,'linestyle','--','color',rgb('blue')); hold on;
       end
        for ss1=1:size(BlockCrossOld{ii},1)
        plot([3,4]+0.25,[BlockCrossOld{ii}(ss1,3) BlockCrossOld{ii}(ss1,4)],'marker','o','markersize',1,'linestyle','--','color',rgb('blue')); hold on;
       end
        hh1(ii)=errorbar([1:2]+0.25,nanmean(BlockCrossOld{ii}(:,1:2)),nansem(BlockCrossOld{ii}(:,1:2)),'color',col1{2},'capsize',2,'marker','o','markersize',2); 
        hold on;        
        hh1(ii)=errorbar([3:4]+0.25,nanmean(BlockCrossOld{ii}(:,3:4)),nansem(BlockCrossOld{ii}(:,3:4)),'color',col1{2},'capsize',2,'marker','o','markersize',2); 
 newLim=get(gca,'ylim');
Ryax = gca; 
yyaxis right  
ylim(newLim);
Ryaxt = [0:20:100];                 % Get Left axis ticks
RyaxDegC = 5-(([0:20:100]*5)./100);     % Convert to °C (or whatever you want)

% RyaxDegC = round(5-(([3,2,1,0]*5)./100)); 
set(Ryax, 'YTick',Ryaxt, 'YTickLabel', (RyaxDegC) )    % New Y-tick values
set(gca,'fontsize',5,'tickdir','out');
ax1(ii).YColor='k';      
xlim([0 5])

end
close all

BlockCrossNew{1}=perfnew{1}(:,[6,7,9,10]);
BlockCrossNew{2}=perfnew{2}(:,[6,7,9,10]);
BlockCrossOld{1}=perfold{1}(:,[6,7,9,10]);
BlockCrossOld{2}=perfold{2}(:,[6,7,9,10]);


fig1=figure  
xlim1=[0.5,15.5;0.5,12.5];
for ii=1:2
    
    left=0.1; bottom=0.85-(ii-1).*0.15; width=0.5; height=0.1;
    
    ax1(ii)=axes('Position',[left bottom width height])
    if ii==1; xx2=BlockCrossNew; else xx2=BlockCrossOld; end
        xx1{1}= [xx2{1}(:,2)-xx2{1}(:,1)]; 
        xx1{2}= [xx2{2}(:,2)-xx2{2}(:,1)]; 
        xx1{3}= [xx2{1}(:,4)-xx2{1}(:,3)]; 
        xx1{4}= [xx2{2}(:,4)-xx2{2}(:,3)];   
       for ss1=1:size(xx1,2)        
        plot([ss1]+randn(size(xx1{ss1},1),1).*0.1,[xx1{ss1}],'marker','o','markersize',2,'linestyle','none','color',col2{ss1}); hold on;
       end
       for ss1=1:size(xx1,2) 
        hh1(ii)=errorbar(ss1,nanmean([xx1{ss1}]),nansem(xx1{ss1}),'color',col2{ss1},'capsize',2,'marker','o','markersize',2); 
       mmm1{ii}(ss1,1)=nanmean([xx1{ss1}]);
       mmm2{ii}(ss1,1)=nansem([xx1{ss1}]);
       end 
 plot([0 5],[0 0],'--k');     
 newLim=get(gca,'ylim');
 ylim(newLim);
set(gca, 'box','off','tickdir','out','xtick',[1.5 3.5],'xticklabel',{'First','Second'},'fontsize',5);
xlabel('Disruption','fontsize',5);
ylabel('Performance difference','fontsize',5);
ax1(ii).YColor='k';      
xlim([0 5])


pval=[];df=[];G=[];p=[];tstat=[];means=[]; sems=[];
test1={'left','right','left','right'}
for ss1=1:size(xx1,2)
    [h, pval(ss1) ci stat]=ttest(xx1{ss1},0, 'tail',test1{ss1});
    G(ss1,:)=ss1; p(ss1,1)=pval(ss1); df(ss1,1)=stat.df;
    if pval(ss1) < 0.05 & pval(ss1) > 0.01 
        text(ss1-0.2, 65, '*','fontsize',6,'color',col1{2}) ;
    elseif pval(ss1) < 0.01 & pval(ss1) > 0.001 
        text(ss1-0.2, 65, '**','fontsize',6 ,'color',col1{2}) ;
     elseif pval(ss1) < 0.001  
        text(ss1-0.2, 65, '***','fontsize',6 ,'color',col1{2}) ;
    else
         text(ss1-0.2, 65, 'n.s.','fontsize',4 ,'color',col1{2}) ;
    end
end

Groups=[1,2;3,4];
for ss1=1:size(Groups,1)
    [p h]=ranksum(xx1{Groups(ss1,1)},xx1{Groups(ss1,2)})
    if p<0.05
      hh4=sigstar(Groups(ss1,:),p,0);
    end
end


 ylim([-40 75]);

end

set(fig1,'paperunits','inches') ;
set(fig1,'papertype','usletter') ;
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig6/Figures');
cd(GenPath);
set(gcf,'paperposition',[0.1,0.1,2,8])   
print(gcf,'-djpeg','-r300','BlockTransitionDisruptionCompared'); 
print(gcf,'-painters','-depsc','-r300','BlockTransitionDisruptionCompared'); 
close(fig1)

