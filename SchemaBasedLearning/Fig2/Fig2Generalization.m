

%% Fig 2b

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load PlFieldTempNAll
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
for i=1:size(PlFieldTemp,2)
     plf{i,1}=PlFieldTemp(i).PlaceFields;
     animal(i,1)=PlFieldTemp(i).Animal;
     day(i,1)=PlFieldTemp(i).Day;
     session(i,1)=PlFieldTemp(i).session;
     track(i,1)=PlFieldTemp(i).Track;
     direction(i,1)=PlFieldTemp(i).Direction;
end
Num_startframes=[];  rNum_startframes=[];
for ii=[14]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load mazel
load LFPDataV2
dn1=pwd;
clusters = importdata('clusters');
elecreg = importdata('elecreg');
% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1; 
LapInfo1=LapInfo;
Replayfile1=strcat(paths,Replayfile{ii});
cd(Replayfile1{1});
dn=pwd;
fstruct = dir('*FrameHPCShuff*.mat');
fstruct1=array2table(fstruct);
if ~isempty(fstruct1);
for i=[2]
     dn=pwd;
     filename=fstruct(i).name;
     load(filename);
     dotLocations = find(filename == '.')
     filename1 = filename(1:dotLocations(1)-1)
     assignin ( 'caller', 'framefile',eval(filename1));
     clear((filename1))
     sess=i/2;
     
     indx1=LapInfo1(:,5)==sess;
     LapInfo=LapInfo1(indx1,:);
     indx=~isnan(LapInfo(:,10)) & LapInfo(:,9)==0;     
     time_start=LapInfo(indx,1:4);
     start_arm=LapInfo(indx,6);
     target_arm=LapInfo(indx,10);
     [pos b]=find(indx);next_visited_arm=[];
     for a=1:sum(indx)
         next_visited_arm(a)=LapInfo(pos(a)+1,6);
     end
     visited_arms=[];
     for a=1:sum(indx)-1
         visited_arms{a}=LapInfo(pos(a)+1:pos(a+1)-1,6);
     end
     visited_arms{a+1}=LapInfo(pos(a+1)+1:length(indx)-1,6);
     frame_times=[];
     for j=1:length(framefile)
         t1=squeeze(framefile(j).t(sess,1,1,:,:));
         frame_times(j,:)=[t1(1), t1(end)];
     end
     frameids=[];
     for j=1:length(time_start);
         if ~isnan(time_start(j,1))
          indx=(frame_times(:,1) >= time_start(j,2) & frame_times(:,1) <= time_start(j,3));     
          frameids(:,j)=indx;
         else
          indx=(frame_times(:,1) >= time_start(j,2) & frame_times(:,1) <= time_start(j,3));     
          frameids(:,j)=indx;
         end
     end
     start_Pvalues=[];start_SeqScore=[];start_rz=[];start_jd=[];sig_arms=[];trial_number=[];
     current=[];next=[];target=[];visit=[];tt=1;start_t=[];start_pdf=[];
   
     
 for j=12;
        fig1=figure(1)
        tvec=[time_start(j,2) time_start(j,3)];
        indx=LFPData(i).t >=tvec(1) & LFPData(i).t <=tvec(end); 
        LFPamp1=LFPData(i).LFP(indx);
        LFPtvec=LFPData(i).t(indx);
        LFPamp2=LFPData(i).firLFP{1, 3}(indx);
        LFPtheta=LFPData(i).firLFP{1, 2}(indx);
        ax=axes('position',[0.1 0.85 0.8 0.08])  
        plot(LFPtvec,LFPamp1,'-k','linewidth',0.1);hold on;
        set(gca,'box','off','XTickLabel', [],'TickLength', [0 0],'xtick',[],'XColor',get(gcf,'color'),...
            'fontsize',5);
        xlim([LFPtvec(1) LFPtvec(end)]);
        set(gca,'ytick',[nanmin(LFPamp1)+100 0 nanmax(LFPamp1)-100],'fontsize',5)
        ax=axes('position',[0.1 0.70 0.8 0.08])  
        plot(LFPtvec,LFPamp2,'-k','linewidth',0.1);hold on;
        set(gca,'box','off','XTickLabel', [],'TickLength', [0 0],'xtick',[],'XColor',get(gcf,'color'),...
            'fontsize',5);
        xlim([LFPtvec(1) LFPtvec(end)]);
        set(gca,'ytick',[nanmin(LFPamp2)+100 0 nanmax(LFPamp2)-100],'fontsize',5)
%         ax=axes('position',[0.1 0.45 0.30 0.08]) 
cd(dn1)              
cltime1=[];cltime2=[];
for nn=1:2
if nn ==1    
kkk=1:length(clusters);
kk=kkk(ind_HPC);
cl_interest=clusters(ind_HPC);
if sum(ind_HPC) < 10; continue; end
for a=1:length(cl_interest)
clG=cl2mat(cl_interest{a});
% get time, velocity and ID
cltime1= [cltime1; clG.featuredata((clG.featuredata(:,8)>=tvec(1) & clG.featuredata(:,8)<=tvec(end)),[8,11]) ...
           repmat(a,size(clG.featuredata((clG.featuredata(:,8)>=tvec(1) & clG.featuredata(:,8)<=tvec(end)),8),1),1)];
clear clG
end
else
kkk=1:length(clusters);
kk=kkk(ind_PFC);
cl_interest=clusters(ind_PFC);
if sum(ind_PFC) < 10; continue; end
for a=1:length(cl_interest)
clG=cl2mat(cl_interest{a});
% get time, velocity and ID
cltime2= [cltime2; clG.featuredata((clG.featuredata(:,8)>=tvec(1) & clG.featuredata(:,8)<=tvec(end)),[8,11]) ...
           repmat(a,size(clG.featuredata((clG.featuredata(:,8)>=tvec(1) & clG.featuredata(:,8)<=tvec(end)),8),1),1)];
clear clG
end
end
end
cd(dn)
fig1=figure(1)
ax=axes('position',[0.1 0.1 0.80 0.55])
for a=1:length(cltime1)
line([cltime1(a,1) cltime1(a,1)],[cltime1(a,3)-0.5 cltime1(a,3)+0.5],'linestyle','-','color','r');hold on;
end
for b=1:length(cltime2)
line([cltime2(b,1) cltime2(b,1)],[cltime2(b,3)-0.5 cltime2(b,3)+0.5]+cltime1(a,3),'linestyle','-','color','k');hold on;
end
xlim([LFPtvec(1) LFPtvec(end)]);
ylim([0.5 cltime2(b,3)+0.5+cltime1(a,3)])

xx=get(gca,'xlim');

set(gca,'box','off','tickdir','out','fontsize',6); xlabel('time (sec)','fontsize',6);
ylabel('Cell ID','fontsize',6);
set(gca,'ytick',[1 sum(ind_HPC) sum(ind_PFC)+sum(ind_HPC)],'fontsize',5);
set(gca,'xtick',[tvec(1)+0.1 tvec(2)-0.1],'fontsize',5);
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.3,0.3,1.5,1.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-depsc','-painters','-r300',strcat('AnimalFiringExample',num2str(ii),'Session',num2str(sess),'TrialNumber',num2str(j)));
print(fig1,'-djpeg','-r300',strcat('AnimalFiringExample',num2str(ii),'Session',num2str(sess),'TrialNumber',num2str(j)));        
close(fig1)
close all
cd(dn)   
       
end
end
end       
end

%% Fig 2C and 2D
% Group Asn TSNE plots
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SessionWiseRDMMatHPCV2','spk3');
spk{1,1}=spk3;
load('SessionWiseRDMMatPFCV2','spk3');
spk{1,2}=spk3;
corrMat=[];
%% Part Days

corrMat=[];

    for j=1:2
     
        corrMat{j}=[];
        Block{j}=[];
        Arms{j}=[];
        TrialTypes{j}=[];
        TimeBins{j}=[];
    end
sess={[6,7,9,10,11,12,14,15],[6,8,9,10,11,13,14,15]}
for Reg=1:2
for ii=sess{Reg}
x1=[]; x2=[]; x3=[]; x4=[];x5=[];   
for ll=1:3     
for kk=1:2 
for jj=1:6
for rr=1:25 
x1=[x1;  spk{1,Reg}{ii,ll}{jj,kk}(:,rr)'];   
x2=[x2; ll];
x3=[x3; jj];
x4=[x4; kk];
x5=[x5; rr];
end
end
end
end
Block{Reg}=[Block{Reg},x2];
Arms{Reg}=[Arms{Reg},x3];
TrialTypes{Reg}=[TrialTypes{Reg}, x4];
TimeBins{Reg}=[TimeBins{Reg}, x5];
corrMat{Reg}=[corrMat{Reg},x1];
end
end

%% 3D
col1{1}={rgb('Red'),rgb('IndianRed'),rgb('LightSalmon')};
col1{2}={rgb('Blue'),rgb('DodgerBlue'),rgb('LightBlue')};
col1{3}={rgb('Purple'),rgb('MediumPurple'),rgb('Plum')};
col1{4}={rgb('Green'),rgb('LimeGreen'),rgb('LightGreen')};
col1{5}={rgb('Brown'),rgb('Chocolate'),rgb('SandyBrown')};
col1{6}={rgb('Black'),rgb('Gray'),rgb('DarkGray')};
Mar={'o','.'};

corrMat{1,1}(isnan(corrMat{1,1}))=0;
goodrows = not(any(isnan(corrMat{1,1}),2));
rng(31)
Y = tsne(corrMat{1,1},'Algorithm','barneshut','Distance','correlation','Standardize',true,'NumDimensions',2);
fig1=figure   
ss=0;sizes=[2,2]
for ii=1:6
    
        for kk=1:3
            for jj=1:2
goodrows1 = Arms{1,1}(goodrows,1) ==ii &  TrialTypes{1,1}(goodrows,1) ==jj & Block{1,1}(goodrows,1) ==kk;
if isempty(Y(goodrows1 ,1)); continue; end
hold on
plot(Y(goodrows1 ,1),Y(goodrows1 ,2),'Marker',Mar{jj},'color',col1{ii}{kk},'MarkerSize',sizes(jj),'linestyle','none','LineWidth',0.2)
        end
  labels2={'Ba','Assn','Assd'};
ss=ss+1;
text(+45,50-(ss-1).*3,strcat('Goal',num2str(ii),'in',labels2(kk)),'Fontsize',4,'Color',col1{ii}{kk});
        
    end
end

text(-50,50,strcat('o: positive outcome'),'Fontsize',6,'Color','k');
text(-50,47,strcat('.: negative outcome'),'Fontsize',6,'Color','k');
xlim([-60 60]);
ylim([-60 60]);
set(gca,'tickdir','out','box','off','fontsize',5)
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,1.5,1.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','PartDaysHPCTSNEEmbeddingsGoal');
print(fig1,'-painters','-depsc','-r600','PartDaysHPCTSNEEmbeddingsGoal'); 
close all;

% PFC Group Asn
%% 24
corrMat{1,2}(isnan(corrMat{1,2}))=0;
goodrows = not(any(isnan(corrMat{1,2}),2));
rng(24)
Y = tsne(corrMat{1,2},'Algorithm','barneshut','Distance','correlation','Standardize',true,'NumDimensions',2);
fig1=figure 
ss=0;sizes=[2,2]
for ii=1:6
    
        for kk=1:3
            for jj=1:2
goodrows1 = Arms{1,2}(goodrows) ==ii &  TrialTypes{1,2}(goodrows) ==jj & Block{1,2}(goodrows) ==kk;
if isempty(Y(goodrows1 ,1)); continue; end
hold on
plot(Y(goodrows1 ,1),Y(goodrows1 ,2),'Marker',Mar{jj},'color',col1{ii}{kk},'MarkerSize',sizes(jj),'linestyle','none','LineWidth',0.2)
        end
    labels2={'Ba','Assn','Assd'};
ss=ss+1;
text(+45,50-(ss-1).*3,strcat('Goal',num2str(ii),'in',labels2(kk)),'Fontsize',4,'Color',col1{ii}{kk});
        
    end
end

text(-50,50,strcat('o: positive outcome'),'Fontsize',6,'Color','k');
text(-50,47,strcat('x: negative outcome'),'Fontsize',6,'Color','k');
xlim([-50 50]);
ylim([-50 50])
set(gca,'tickdir','out','box','off','fontsize',5)
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,1.5,1.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','PartDaysPFCTSNEEmbeddingsGoal');
print(fig1,'-painters','-depsc','-r300','PartDaysPFCTSNEEmbeddingsGoal'); 
close all;

% Ctrl 1 TSNE
%% Ctrl1

corrMat=[];

    for j=1:2
     
        corrMat{j}=[];
        Block{j}=[];
        Arms{j}=[];
        TrialTypes{j}=[];
        TimeBins{j}=[];
    end
sess={[1,2,4,5],[1,3,4,5]}
for Reg=1:2
for ii=sess{Reg}
x1=[];    
for ll=1:3     
for kk=1:2 
for jj=1:6
for rr=1:25 
x1=[x1;  spk{1,Reg}{ii,ll}{jj,kk}(:,rr)'];   

Block{Reg}=[Block{Reg}; ll];
Arms{Reg}=[Arms{Reg}; jj];
TrialTypes{Reg}=[TrialTypes{Reg}; kk];
TimeBins{Reg}=[TimeBins{Reg}; rr];
end
end
end
end
corrMat{Reg}=[corrMat{Reg},x1];
end
end


Mar={'o','.'};
corrMat{1,1}(isnan(corrMat{1,1}))=0;
goodrows = not(any(isnan(corrMat{1,1}),2));
% rng(35)
rng(29)
Y = tsne(corrMat{1,1},'Algorithm','barneshut','Distance','correlation','Standardize',true,'NumDimensions',2);
fig1=figure   
ss=0;sizes=[2,2]
for ii=1:6
    
        for kk=1:3
            for jj=1:2
goodrows1 = Arms{1,1}(goodrows) ==ii &  TrialTypes{1,1}(goodrows) ==jj & Block{1,1}(goodrows) ==kk;
if isempty(Y(goodrows1 ,1)); continue; end
hold on
plot(Y(goodrows1 ,1),Y(goodrows1 ,2),'Marker',Mar{jj},'color',col1{ii}{kk},'MarkerSize',sizes(jj),'linestyle','none','LineWidth',0.2)
        end
  labels2={'Ba','Assn','Assd'};
ss=ss+1;
text(+45,50-(ss-1).*2,strcat('Goal',num2str(ii),'in',labels2(kk)),'Fontsize',4,'Color',col1{ii}{kk});
        
    end
end

text(-50,50,strcat('o: positive outcome'),'Fontsize',6,'Color','k');
text(-50,47,strcat('x: negative outcome'),'Fontsize',6,'Color','k');
xlim([-40 50]);
ylim([-50 50])
set(gca,'tickdir','out','box','off','fontsize',5)
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,1.5,1.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','Ctrl1HPCTSNEEmbeddingsGoal');
print(fig1,'-painters','-depsc','-r300','Ctrl1HPCTSNEEmbeddingsGoal'); 
close all;

corrMat{1,2}(isnan(corrMat{1,2}))=0;
goodrows = not(any(isnan(corrMat{1,2}),2));
rng(29)
Y = tsne(corrMat{1,2},'Algorithm','barneshut','Distance','correlation','Standardize',true,'NumDimensions',2);
fig1=figure 
ss=0;
for ii=1:6
    
        for kk=1:3
            for jj=1:2
goodrows1 = Arms{1,2}(goodrows) ==ii &  TrialTypes{1,2}(goodrows) ==jj & Block{1,2}(goodrows) ==kk;
if isempty(Y(goodrows1 ,1)); continue; end
hold on
plot(Y(goodrows1 ,1),Y(goodrows1 ,2),'Marker',Mar{jj},'color',col1{ii}{kk},'MarkerSize',sizes(jj),'linestyle','none','LineWidth',0.2)
        end
    labels2={'Ba','Assn','Assd'};
ss=ss+1;
text(+45,50-(ss-1).*2,strcat('Goal',num2str(ii),'in',labels2(kk)),'Fontsize',4,'Color',col1{ii}{kk});
        
    end
end

text(-50,50,strcat('o: positive outcome'),'Fontsize',6,'Color','k');
text(-50,47,strcat('x: negative outcome'),'Fontsize',6,'Color','k');
xlim([-40 40]);
ylim([-40 40])
set(gca,'tickdir','out','box','off','fontsize',5)
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,1.5,1.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','Ctrl1PFCTSNEEmbeddingsGoal');
print(fig1,'-painters','-depsc','-r300','Ctrl1PFCTSNEEmbeddingsGoal'); 
close all;

%% Fig 2e
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('PlotExampleRSMDataAllDays','corrMat1');
Days=3; 
fig1=figure
for ee=1:3
for Reg=1:2    
left=0.1+(ee-1).*0.305; bot=0.7-(Reg-1).*0.30; width=0.25; height=0.20;
ax(ee,Reg)=axes('position',[left bot width height])
mm=squeeze(corrMat1(ee,Reg,Days,:,:));
mm(logical(eye(size(mm))))=NaN;
corrMat2=squeeze(corrMat1(ee,Reg,Days,:,:));
corrMat2(logical(eye(size(mm))))=NaN;
x1=corrMat2;
if Reg==2
imagesc(x1,[nanmin(nanmin(x1))+nanmax(nanmax(x1))/6 ...
    nanmax(nanmax(x1))])
else
imagesc(x1,[nanmin(nanmin(x1))+nanmax(nanmax(x1))/6 ...
    nanmax(nanmax(x1))])    
end
% colormap jet
cmap=getNCLColmap('MPL_Greens.rgb',128)
colormap(gca,cmap)
colorbar
cb=colorbar;
set(cb,'Position',[left+width+0.005 bot 0.01 0.05])
cb.FontSize=4;
hold on;
plot([0 36],[12.5 12.5],'linestyle','--','Marker','none','linewidth',1,'color','b');
plot([0 36],[24.5 24.5],'linestyle','--','Marker','none','linewidth',1,'color','b');
plot([12.5 12.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','b');
plot([24.5 24.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','b');

plot([0 36],[6.5 6.5],'linestyle','--','Marker','none','linewidth',1,'color','r');
plot([0 36],[18.5 18.5],'linestyle','--','Marker','none','linewidth',1,'color','r');
plot([0 36],[30.5 30.5],'linestyle','--','Marker','none','linewidth',1,'color','r');

plot([6.5 6.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','r');
plot([18.5 18.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','r');
plot([30.5 30.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','r');

plot([0 36],[3.5 3.5],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([0 36],[9.5 9.5],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([0 36],[15.5 15.5],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([0 36],[21.5 21.5],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([0 36],[27.5 27.5],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([0 36],[33.5 33.5],'linestyle','--','Marker','none','linewidth',1,'color','c');

plot([3.5 3.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([9.5 9.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([15.5 15.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([21.5 21.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([27.5 27.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','c');
plot([33.5 33.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','c');

% set(gca,'box','off','tickdir','out','fontsize',3)
% set(gca,'xtick',[1:36],'xticklabel',labels,'fontsize',3,'xticklabelrotation',90);
% set(gca,'ytick',[1:36],'yticklabel',labels,'fontsize',3);



set(gca,'box','off','tickdir','out','fontsize',3)
set(gca,'xtick',[],'xticklabel',[],'fontsize',3,'xticklabelrotation',90);
set(gca,'ytick',[],'yticklabel',[],'fontsize',3,'xticklabelrotation',90);
% set(gca,'ytick',[1:36],'yticklabel',labels,'fontsize',3);
end
end
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.05,0.1,4,5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','ExampleRDMPartialAssn');
print(fig1,'-painters','-depsc','-r300','ExampleRDMPartialAssn'); 
close all;

%% Fig 2F

SpecificInfo=[];SessInfo=[];Outcome=[];NovelInfo=[];SessID=[];
outcomei=[ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1];
outcomej=[ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1];

Noveli=[ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,...
    ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1];
Novelj=[ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,...
    ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1];
Across=zeros(36,36);
SessIDi=[ones(1,12),ones(1,12).*3,ones(1,12).*5];
SessID=[];
for ii=1:36
    for jj=1:36
       
        if ii==jj | jj==ii+6 | jj==ii+12 | jj==ii+18 | jj==ii+24 | jj==ii+30 | jj==ii+36 ...
         | jj==ii | jj==ii-6 | jj==ii-12 | jj==ii-18 | jj==ii-24 | jj==ii-30 | jj==ii-36          
        SpecificInfo(ii,jj)=1;
       else
         SpecificInfo(ii,jj)=-1;  
       end
       
       if (ii <=12 & jj <= 12) | ((ii > 12 & ii <= 24) &  ((jj > 12 & jj <= 24)) ) | (ii > 24 & jj > 24) 
       SessInfo(ii,jj)=1;
       else
       SessInfo(ii,jj)=-1; 
       end
       Outcome(ii,jj)=outcomei(ii)*outcomej(jj);     
       NovelInfo(ii,jj)=Noveli(ii)*Novelj(jj); 
       if jj>ii
           Across(ii,jj)=logical(1);
       end
       
       SessID(ii,jj)=SessIDi(ii)*SessIDi(jj);
     end
end

SpecificInfo1=SpecificInfo;
SpecificInfo11=SpecificInfo1;

SpecificInfo11(SessInfo==1 | Across==0)=0.4;

Outcome11=Outcome;

Outcome11(SessInfo==1 | Across==0)=0.4;

SpecificInfo1(SpecificInfo1==-1)=0; 
Interact=SpecificInfo1.*Outcome;
Interact11=Interact;
% Interact11(Interact11==0)=0.4;
Interact11(SessInfo==1 | Across==0)=0.4;
fig1=figure
subplot(3,3,1)
imagesc(Outcome11)
set(gca,'box','off','xtick','','ytick','');
hold on;
plot([0 36],[12.5 12.5],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([0 36],[24.5 24.5],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([12.5 12.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([24.5 24.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','b');

plot([0 36],[6.5 6.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 36],[18.5 18.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 36],[30.5 30.5],'linestyle','-','Marker','none','linewidth',1,'color','r');

plot([6.5 6.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([18.5 18.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([30.5 30.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');

subplot(3,3,4)
imagesc(SpecificInfo11)
set(gca,'box','off','xtick','','ytick','');
hold on;
plot([0 36],[12.5 12.5],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([0 36],[24.5 24.5],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([12.5 12.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([24.5 24.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','b');

plot([0 36],[6.5 6.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 36],[18.5 18.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 36],[30.5 30.5],'linestyle','-','Marker','none','linewidth',1,'color','r');

plot([6.5 6.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([18.5 18.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([30.5 30.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');

subplot(3,3,7)
imagesc(Interact11)
set(gca,'box','off','xtick','','ytick','');hold on;
plot([0 36],[12.5 12.5],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([0 36],[24.5 24.5],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([12.5 12.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','b');
plot([24.5 24.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','b');

plot([0 36],[6.5 6.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 36],[18.5 18.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 36],[30.5 30.5],'linestyle','-','Marker','none','linewidth',1,'color','r');

plot([6.5 6.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([18.5 18.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([30.5 30.5],[0 36],'linestyle','-','Marker','none','linewidth',1,'color','r');

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,4,4.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','RSMSimilarityCartoons');
print(fig1,'-painters','-depsc','-r300','RSMSimilarityCartoons'); 
close all;


%% Fig 2g
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CPDRSMSpatialV2.mat','CPD1','CPDS')
CPDSpatial=CPD1;
CPDSpatialShuff=CPDS;

%spatial mat dir days reg timebins shuffledata factors
load('CPDRSMTemporalV2.mat','CPD1','CPDS')
CPDTemporal=CPD1;
CPDTemporalShuff=CPDS;

%spatial mat days reg timebins shuffledata factors

col1={rgb('Red'),rgb('Blue'),rgb('Black'),rgb('Black')};
labels={'1','4','2&3','5'};
titlesReg={'PFC','HPC'};
siglabels={'*','**','***'};
ax1=[];x1=[];x1S=[];x2=[];x2S=[];
fig1=figure
for nn=1:3
    
for ee=1:3
   
left=0.1+(ee-1).*0.2; bot=0.80-(nn-1).*0.2; width=0.15; height=0.15;
ax1(nn,ee)=axes('position',[left bot width height]);
if ee==1
x1{nn,ee}=squeeze(CPDSpatial(1,3,1,:,nn)).*100;
x1S{nn,ee}=squeeze(CPDSpatialShuff(1,3,1,:,:,nn)).*100;

x2{nn,ee}=squeeze(CPDSpatial(1,3,2,:,nn)).*100;
x2S{nn,ee}=squeeze(CPDSpatialShuff(1,3,2,:,:,nn)).*100;
elseif ee==3
    
x1{nn,ee}=flip(squeeze(CPDSpatial(2,3,1,:,nn))).*100;
x1S{nn,ee}=flip(squeeze(CPDSpatialShuff(2,3,1,:,:,nn))).*100; 

x2{nn,ee}=flip(squeeze(CPDSpatial(2,3,2,:,nn))).*100;
x2S{nn,ee}=flip(squeeze(CPDSpatialShuff(2,3,2,:,:,nn))).*100;
elseif ee==2
x1{nn,ee}=squeeze(CPDTemporal(3,1,:,nn)).*100;
x1S{nn,ee}=squeeze(CPDTemporalShuff(3,1,:,:,nn)).*100; 

x2{nn,ee}=squeeze(CPDTemporal(3,2,:,nn)).*100;
x2S{nn,ee}=squeeze(CPDTemporalShuff(3,2,:,:,nn)).*100; 
end
y1=movmean(x1{nn,ee},[2,2],1);
plot(1:size(y1,1),y1,'color',col1{1},'linestyle','-','linewidth',2);
hold on;
y1=movmean(x2{nn,ee},[2,2],1);
plot(1:size(y1,1),y1,'color',col1{2},'linestyle','-','linewidth',2);
set(gca,'box','off','tickdir','out','fontsize',5);
if ee==1
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Center to goal)','fontsize',5);
ylabel('Percent variance explained','fontsize',5);
elseif ee==2
set(gca,'xtick',[0 5 25],'xticklabel',[-1 0 5],'fontsize',5);
xlabel('Time at goal (s)','fontsize',5); 
else
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Goal to center)','fontsize',5);    
end
    
end
end


linkaxes([ax1(1,1) ax1(1,2) ax1(1,3) ],'y');

linkaxes([ax1(2,1) ax1(2,2) ax1(2,3) ],'y');

linkaxes([ax1(3,1) ax1(3,2) ax1(3,3) ],'y');
for nn=1:3
    for ee=1:3
    yl(nn,:) = (get(ax1(nn,3), 'Ylim')) ;
    
     Pval=[];sig=[];
     for ss=1:size(x1{nn,ee},1)    
         Pval(ss)=1-sum(x1{nn,ee}(ss,:)> x1S{nn,ee}(ss,:))./size(x1S{nn,ee}(ss,:),2);
         if Pval(ss) <0.01 & Pval(ss) >= 0.001 ;
         sig{ss}=siglabels{1};
         elseif Pval(ss) <0.001 & Pval(ss) >= 0.0001; 
         sig{ss}=siglabels{2};  
         elseif Pval(ss) <0.0001 & Pval(ss) >= 0;
         sig{ss}=siglabels{3}; 
         else
         sig{ss}=NaN;     
         end
     if ~isnan(sig{ss})    
     t=text('Parent', ax1(nn,ee), 'Position', [ss yl(nn,end)+yl(nn,end)*25/100], 'String', sig{ss});
     t.FontSize=3; t.Color=col1{1};
     t.Rotation=270;
     end
     end
     
       Pval=[];sig=[];
     for ss=1:size(x2{nn,ee},1)   
         
             
         Pval(ss)=1-sum(x2{nn,ee}(ss,:)> x2S{nn,ee}(ss,:))./size(x2S{nn,ee}(ss,:),2);
         if Pval(ss) <0.01 & Pval(ss) >= 0.001 ;
         sig{ss}=siglabels{1};
         elseif Pval(ss) <0.001 & Pval(ss) >= 0.0001; 
         sig{ss}=siglabels{2};  
         elseif Pval(ss) <0.0001 & Pval(ss) >= 0;
         sig{ss}=siglabels{3}; 
         else
         sig{ss}=NaN;     
         end
     if ~isnan(sig{ss})    
     t=text('Parent', ax1(nn,ee), 'Position', [ss yl(nn,end)+yl(nn,end)*15/100], 'String', sig{ss});
     t.FontSize=3; t.Color=col1{2};
     t.Rotation=270;
     end
     end
   
    
       Pval=[];sig=[];
     for ss=1:size(x2{nn,ee},1) 
         if nn==1
             diff=1-sum((x2{nn,ee}(ss,:)-x1{nn,ee}(ss,:))> (x2S{nn,ee}(ss,:)-x1S{nn,ee}(ss,:)))./size(x2S{nn,ee}(ss,:),2);
         else
             diff=1-sum((x1{nn,ee}(ss,:)-x2{nn,ee}(ss,:))> (x1S{nn,ee}(ss,:)-x2S{nn,ee}(ss,:)))./size(x2S{nn,ee}(ss,:),2);
         end
         Pval(ss)=diff
         if Pval(ss) <0.01 & Pval(ss) >= 0.001 ;
         sig{ss}=siglabels{1};
         elseif Pval(ss) <0.001 & Pval(ss) >= 0.0001; 
         sig{ss}=siglabels{2};  
         elseif Pval(ss) <0.0001 & Pval(ss) >= 0;
         sig{ss}=siglabels{3}; 
         else
         sig{ss}=NaN;     
         end
     if ~isnan(sig{ss})    
     t=text('Parent', ax1(nn,ee), 'Position', [ss yl(nn,end)+yl(nn,end)*5/100], 'String', sig{ss});
     t.FontSize=3; t.Color=col1{4};
     t.Rotation=270;
     end
     end
    
    end    
    
end
for nn=1:3
set(ax1(nn,1),'ylim',[0 yl(nn,end)+yl(nn,end)*28/100]);     
end


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,4,4.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','ComparisionTimelinesDay2&3AcrossSessCode');
print(fig1,'-painters','-depsc','-r300','ComparisionTimelinesDay2&3AcrossSessCode'); 
close all;


%% Fig 2h
% PFC Outcome
P=[];
col1={rgb('Black'),rgb('Blue'),rgb('Green'),rgb('Red')};
labels={'1','4','2&3','5'};
titlesReg={'PFC','HPC'};
siglabels={'*','**','***'};
ax1=[];x1=[];x1S=[];x2=[];x2S=[];Pval=[];
fig1=figure
ss=0;Mmeans=[];
for nn=1
ss=ss+1; 
for ee=1:3
   
left=0.1+(ee-1).*0.12; bot=0.80-(nn-1).*0.2; width=0.10; height=0.15;
ax1(ss,ee)=axes('position',[left bot width height]);
    for zz=1:4

if ee==1
x1{nn,ee,zz}=movmean(squeeze(CPDSpatial(1,zz,2,:,nn)).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDSpatialShuff(1,zz,2,:,:,nn)).*100,[5 5],1);

elseif ee==3
    
x1{nn,ee,zz}=movmean(flip(squeeze(CPDSpatial(2,zz,2,:,nn))).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(flip(squeeze(CPDSpatialShuff(2,zz,2,:,:,nn))).*100,[5 5],1); 


elseif ee==2
x1{nn,ee,zz}=movmean(squeeze(CPDTemporal(zz,2,:,nn)).*100,[1 1],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDTemporalShuff(zz,2,:,:,nn)).*100,[1 1],1); 


end
y1=movmean(x1{nn,ee,zz},[1,1],1);
plot(1:size(y1,1),y1,'color',col1{zz},'linestyle','-','linewidth',2);

Mmeans{nn,ee}(zz,1)=nanmean(y1);
Mmeans{nn,ee}(zz,2)=nansem(y1);
hold on;

set(gca,'box','off','tickdir','out','fontsize',5);
if ee==1
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Center to goal)','fontsize',4);
ylabel('Percent variance explained','fontsize',4);
elseif ee==2
set(gca,'xtick',[0 5 25],'xticklabel',[-1 0 5],'fontsize',5);
xlabel('Time at goal (s)','fontsize',4); 
else
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Goal to center)','fontsize',4);    
end
    
    end
   Groups=[1,2;1,3;1,4;2,3;2,4;3,4];
for zz= 1:size(Groups,1)
    [ h Pval1]=ttest2(x1{nn,ee,Groups(zz,2)},x1{nn,ee,Groups(zz,1)});
    Pval(nn,ee,zz)=round(Pval1,4);
    
end
   
end
end
P{2,1}= squeeze(Pval(ss,:,:));
Means{2,1}=Mmeans;
linkaxes([ax1(ss,1) ax1(ss,2) ax1(ss,3) ],'y');

% linkaxes([ax1(2,1) ax1(2,2) ax1(2,3) ],'y');
% 
% linkaxes([ax1(3,1) ax1(3,2) ax1(3,3) ],'y');
for nn=ss
    for ee=1:3
    yl(nn,:) = (get(ax1(nn,3), 'Ylim')) ;
     str=num2str(squeeze(Pval(nn,ee,:)));
     t=text('Parent', ax1(nn,ee), 'Position', [1 yl(nn,end)+yl(nn,end)*25/100], 'String', str);
     t.FontSize=4; t.Color=col1{1};
         
    end    
    
end
for nn=ss
set(ax1(nn,1),'ylim',[0 yl(nn,end)+yl(nn,end)*28/100]);     
end


% PFC Outcome at specific arms
col1={rgb('Black'),rgb('Blue'),rgb('Green'),rgb('Red')};
labels={'1','4','2&3','5'};
titlesReg={'PFC','HPC'};
siglabels={'*','**','***'};
x1=[];x1S=[];x2=[];x2S=[];Mmeans=[];
for nn=1
 ss=ss+1   
for ee=1:3
   
left=0.5+(ee-1).*0.12; bot=0.80-(nn-1).*0.2; width=0.10; height=0.15;
ax1(ss,ee)=axes('position',[left bot width height]);
    for zz=1:4

if ee==1
x1{nn,ee,zz}=movmean(squeeze(CPDSpatial(1,zz,2,:,3)).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDSpatialShuff(1,zz,2,:,:,3)).*100,[5 5],1);

elseif ee==3
    
x1{nn,ee,zz}=movmean(flip(squeeze(CPDSpatial(2,zz,2,:,3))).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(flip(squeeze(CPDSpatialShuff(2,zz,2,:,:,3))).*100,[5 5],1); 


elseif ee==2
x1{nn,ee,zz}=movmean(squeeze(CPDTemporal(zz,2,:,3)).*100,[1 1],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDTemporalShuff(zz,2,:,:,3)).*100,[1 1],1); 


end
y1=movmean(x1{nn,ee,zz},[1,1],1);
plot(1:size(y1,1),y1,'color',col1{zz},'linestyle','-','linewidth',2);
Mmeans{nn,ee}(zz,1)=nanmean(y1);
Mmeans{nn,ee}(zz,2)=nansem(y1);
hold on;

set(gca,'box','off','tickdir','out','fontsize',5);
if ee==1
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Center to goal)','fontsize',4);
ylabel('Percent variance explained','fontsize',4); 
elseif ee==2
set(gca,'xtick',[0 5 25],'xticklabel',[-1 0 5],'fontsize',5);
xlabel('Time at goal (s)','fontsize',4); 
else
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Goal to center)','fontsize',4);    
end
   
    end
   Groups=[1,2;1,3;1,4;2,3;2,4;3,4];
for zz= 1:size(Groups,1)
    [ h Pval1]=ttest2(x1{nn,ee,Groups(zz,2)},x1{nn,ee,Groups(zz,1)});
    Pval(ss,ee,zz)=round(Pval1,4);
end
   
end
end
P{2,2}= squeeze(Pval(ss,:,:));
Means{2,2}=Mmeans;
linkaxes([ax1(ss,1) ax1(ss,2) ax1(ss,3) ],'y');

% linkaxes([ax1(2,1) ax1(2,2) ax1(2,3) ],'y');
% 
% linkaxes([ax1(3,1) ax1(3,2) ax1(3,3) ],'y');
for nn=ss
    for ee=1:3
    yl(nn,:) = (get(ax1(nn,3), 'Ylim')) ;
     str=num2str(squeeze(Pval(nn,ee,:)));
     t=text('Parent', ax1(nn,ee), 'Position', [1 yl(nn,end)+yl(nn,end)*25/100], 'String', str);
     t.FontSize=4; t.Color=col1{1};
         
    end    
    
end
for nn=ss
set(ax1(nn,1),'ylim',[0 yl(nn,end)+yl(nn,end)*28/100]);     
end



% HPC Outcome 

col1={rgb('Black'),rgb('Blue'),rgb('Green'),rgb('Red')};
labels={'1','4','2&3','5'};
titlesReg={'PFC','HPC'};
siglabels={'*','**','***'};
ax1=[];x1=[];x1S=[];x2=[];x2S=[];Pval=[];
Mmeans=[];
for nn=1
ss=ss+1; 
for ee=1:3
   
left=0.1+(ee-1).*0.12; bot=0.60-(nn-1).*0.2; width=0.10; height=0.15;
ax1(ss,ee)=axes('position',[left bot width height]);
    for zz=1:4

if ee==1
x1{nn,ee,zz}=movmean(squeeze(CPDSpatial(1,zz,1,:,nn)).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDSpatialShuff(1,zz,1,:,:,nn)).*100,[5 5],1);

elseif ee==3
    
x1{nn,ee,zz}=movmean(flip(squeeze(CPDSpatial(2,zz,1,:,nn))).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(flip(squeeze(CPDSpatialShuff(2,zz,1,:,:,nn))).*100,[5 5],1); 


elseif ee==2
x1{nn,ee,zz}=movmean(squeeze(CPDTemporal(zz,1,:,nn)).*100,[1 1],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDTemporalShuff(zz,1,:,:,nn)).*100,[1 1],1); 


end
y1=movmean(x1{nn,ee,zz},[1,1],1);
plot(1:size(y1,1),y1,'color',col1{zz},'linestyle','-','linewidth',2);
hold on;
Mmeans{nn,ee}(zz,1)=nanmean(y1);
Mmeans{nn,ee}(zz,2)=nansem(y1);
set(gca,'box','off','tickdir','out','fontsize',5);
if ee==1
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
ylabel('Percent variance explained','fontsize',4);
xlabel('Distance (Goal to center)','fontsize',4);  
elseif ee==2
set(gca,'xtick',[0 5 25],'xticklabel',[-1 0 5],'fontsize',5);
xlabel('Time at goal (s)','fontsize',4); 
else
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Goal to center)','fontsize',4);    
end
  
    end
   Groups=[1,2;1,3;1,4;2,3;2,4;3,4];
for zz= 1:size(Groups,1)
    [ h Pval1]=ttest2(x1{nn,ee,Groups(zz,2)},x1{nn,ee,Groups(zz,1)});
    Pval(ss,ee,zz)=round(Pval1,4);
end

end
end
P{1,1}= squeeze(Pval(ss,:,:));
Means{1,1}=Mmeans;
linkaxes([ax1(ss,1) ax1(ss,2) ax1(ss,3) ax1(1,1) ax1(1,2) ax1(1,3)],'y');

% linkaxes([ax1(2,1) ax1(2,2) ax1(2,3) ],'y');
% 
% linkaxes([ax1(3,1) ax1(3,2) ax1(3,3) ],'y');
for nn=ss
    for ee=1:3
    yl(nn,:) = (get(ax1(nn,3), 'Ylim')) ;
     str=num2str(squeeze(Pval(ss,ee,:)));
     t=text('Parent', ax1(nn,ee), 'Position', [1 yl(nn,end)+yl(nn,end)*25/100], 'String', str);
     t.FontSize=4; t.Color=col1{1};
         
    end    
    
end
for nn=ss
set(ax1(nn,1),'ylim',[0 yl(nn,end)+yl(nn,end)*28/100]);     
end

% HPC Outcome at specific arms

col1={rgb('Black'),rgb('Blue'),rgb('Green'),rgb('Red')};
labels={'1','4','2&3','5'};
titlesReg={'PFC','HPC'};
siglabels={'*','**','***'};
x1=[];x1S=[];x2=[];x2S=[];Mmeans=[];
for nn=1
 ss=ss+1   
for ee=1:3
   
left=0.5+(ee-1).*0.12; bot=0.60-(nn-1).*0.2; width=0.10; height=0.15;
ax1(ss,ee)=axes('position',[left bot width height]);
    for zz=1:4

if ee==1
x1{nn,ee,zz}=movmean(squeeze(CPDSpatial(1,zz,1,:,3)).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDSpatialShuff(1,zz,1,:,:,3)).*100,[5 5],1);

elseif ee==3
    
x1{nn,ee,zz}=movmean(flip(squeeze(CPDSpatial(2,zz,1,:,3))).*100,[5 5],1);
x1S{nn,ee,zz}=movmean(flip(squeeze(CPDSpatialShuff(2,zz,1,:,:,3))).*100,[5 5],1); 


elseif ee==2
x1{nn,ee,zz}=movmean(squeeze(CPDTemporal(zz,1,:,3)).*100,[1 1],1);
x1S{nn,ee,zz}=movmean(squeeze(CPDTemporalShuff(zz,1,:,:,3)).*100,[1 1],1); 


end
y1=movmean(x1{nn,ee,zz},[1,1],1);
plot(1:size(y1,1),y1,'color',col1{zz},'linestyle','-','linewidth',2);
hold on;
Mmeans{nn,ee}(zz,1)=nanmean(y1);
Mmeans{nn,ee}(zz,2)=nansem(y1);
set(gca,'box','off','tickdir','out','fontsize',5);
if ee==1
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
ylabel('Percent variance explained','fontsize',4);
elseif ee==2
set(gca,'xtick',[0 5 25],'xticklabel',[-1 0 5],'fontsize',5);
xlabel('Time at goal (s)','fontsize',4); 
else
set(gca,'xtick',[0 71],'xticklabel',[0 70],'fontsize',5);
xlabel('Distance (Goal to center)','fontsize',4);    
end
    
    end
   Groups=[1,2;1,3;1,4;2,3;2,4;3,4];
for zz= 1:size(Groups,1)
    [ h Pval1]=ttest2(x1{nn,ee,Groups(zz,2)},x1{nn,ee,Groups(zz,1)});
    Pval(ss,ee,zz)=round(Pval1,4);
end

end
end
P{1,2}= squeeze(Pval(ss,:,:));
Means{1,2}=Mmeans;
linkaxes([ax1(ss,1) ax1(ss,2) ax1(ss,3) ax1(2,1) ax1(2,2) ax1(2,3)],'y');
linkaxes([ax1(2,1) ax1(4,1) ],'y');

% linkaxes([ax1(2,1) ax1(2,2) ax1(2,3) ],'y');
% 
% linkaxes([ax1(3,1) ax1(3,2) ax1(3,3) ],'y');
for nn=ss
    for ee=1:3
    yl(nn,:) = (get(ax1(nn,3), 'Ylim')) ;
     str=num2str(squeeze(Pval(nn,ee,:)));
     t=text('Parent', ax1(nn,ee), 'Position', [1 yl(nn,end)+yl(nn,end)*25/100], 'String', str);
     t.FontSize=4; t.Color=col1{1};
         
    end    
    
end
for nn=ss
set(ax1(nn,1),'ylim',[0 yl(nn,end)+yl(nn,end)*28/100]);     
end

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,5,4.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig2/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','DaysComparisionTimelinesDay2&3AcrossSessCode');
print(fig1,'-painters','-depsc','-r300','DaysComparisionTimelinesDay2&3AcrossSessCode'); 
close all;





