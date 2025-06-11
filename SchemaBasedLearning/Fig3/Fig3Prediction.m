%% Fig 3a and b

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)

for ii=[19]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load maze
dn1=pwd;

cluster=importdata('clusters');
elecreg=importdata('elecreg');

[c1 c2] = textread('celltype','%s %s');
elecreg=importdata('elecreg');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;

tic;

FiringRates1=FiringRates1(:,1);
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
dn=pwd;
fstruct = dir('*FrameHPCShuff*.mat');
if ii==6;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
fstruct1=array2table(fstruct);
LapInfo1=LapInfo;
if ~isempty(fstruct)
   
  
for ll=4
LapInfo = LapInfo1;    
plf1=[];plfConcat=[];  
arms=[1,2,4,5,7,8];
plf_forplot=[];
for nn=1:6
    plf_forplot{nn}=[];
end

for nn=1:length(FiringRates1) 
plf11=FiringRates1{nn};
plf12=[];
ind=[16:40];
for jj=1:6
    plf1=plf11(LapInfo(:,5)==ll & LapInfo(:,6) ==arms(jj) & LapInfo(:,9) ==1,ind);
    plf12=[plf12, nanmean(plf1)];
    plf_forplot{jj}=[plf_forplot{jj}; nanmean(plf1)];
    
end
plfConcat=[plfConcat; (plf12)];
end
   
                  
    for nn=3
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==1));
    ind1=find(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==1));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(LapInfo(:,6)==arms(nn) & idx,:);     
    
       
 for ss=1:6
     
    toplot= plf_forplot{ss}(ind_PFC,:);
    if ss==nn
     
      [maxidx,imax]=max(toplot,[],2);
      [dummy,idx]=sort(imax,'ascend') ;
    end
 end
 fig1=figure(1)
 for ss=1:6
    ax=axes('Position',[0.05+((ss-1).*0.1) 0.05 0.07 0.10]); 
    toplot= plf_forplot{ss}(ind_PFC,:);
    toplot1=toplot(idx,:);
    h=plot_plfield(ax,toplot1);
    if ss==1
    ylabel('cell ID');
    set(gca,'ytick',[1 size(toplot1,1)],'yticklabel',[1 size(toplot1,1)],'fontsize',4)
    else
    set(gca,'ytick',[])    
    end
    set(gca,'xtick',[1, 6, 25],'xticklabel',[-1 0 5],'fontsize',4)
    set(gca,'fontsize',5);
    xlabel('Time (s)');
 end
    
    
    for rr=1:size(LapInfo,1)     
   t_stimulus=LapInfo(rr,7);
   t_bin=0.05;
   tbin=0.02;
   tvec2=[t_stimulus-1   t_stimulus+5];       
   binns=nanmin(tvec2):tbin:nanmax(tvec2);
   binedges=(binns(1:end-1)+binns(2:end))/2; 
cd(dn1);   
kkk=1:length(cluster);
kk=kkk(ind_PFC);   
trackcltime=[];
b11=cluster;
b11=b11(ind_PFC);
for c=1:size(b11,1)
       clG=cl2mat(b11{c});              
       trackcltime{c}= [clG.featuredata((clG.featuredata(:,8)>=tvec2(1) & clG.featuredata(:,8)<=tvec2(2)),[8,11]) ...
           repmat(kk(c),size(clG.featuredata((clG.featuredata(:,8)>=tvec2(1) & clG.featuredata(:,8)<=tvec2(2)),8),1),1)];
       clear clG
end

  cltimevec=[];
   
   for c=1:size(trackcltime,2)
   cltimevec=[cltimevec; trackcltime{c}(:,[1,3])];
   end
   cltimemat=[];
   [~,ind]=sort(cltimevec(:,1));
   cltimemat=cltimevec(ind,:);

 [dprobp dprob binedges]=BayesDecode_FramesWithPlot(tvec2,plfConcat,cltimevec,tbin);
% [dprobp dprob binedges]=BayesDecode_FramesWithOverlap(tvec2,plfConcat,cltimevec,tbin,0.5);
 bins=[1:25:151];   
 
 toplot= dprobp(bins(rr):bins(rr+1)-1,:);
  
zz = toplot';
   
   
 
%         if isempty(zz)
%             continue; end
zz(isnan(zz))=0;        

zz=movmean(zz,[5 5],1,'omitnan');  
zz=movmean(zz,[0.5 1],2,'omitnan'); 
cax=[nanmin(nanmin(zz)) 0.1]; 
 
 
 for ss=1:6
     left=0.05+((ss-1).*0.1); bottom=0.05+0.15+((rr-1).*0.15); width= 0.07; height=0.10;
    ax=axes('Position',[left bottom width height]); 
    toplot= dprobp(bins(ss):bins(ss+1)-1,:);
  
    zz = toplot';
 
zz(isnan(zz))=0;        

zz=movmean(zz,[5 5],1,'omitnan');  
zz=movmean(zz,[0.5 1],2,'omitnan'); 

imagesc(zz,[cax])
cmap=getNCLColmap('WhiteGreen.rgb',100);
% cmap2=getNCLColmap('hlu_default.rgb',30);
shading interp   
colormap(gca,cmap)
% end
if ss==6

colorbar
cb=colorbar;
set(cb,'Position',[left+width+0.01 bottom 0.01 0.05])
cb.FontSize=4;
end
    set(gca,'ydir','normal');
    xlabel('Predicted time (s)');
    set(gca,'xtick',[1, 6, 25],'xticklabel',[-1 0 5],'fontsize',4)
    if ss==1
    ylabel('True time (s)');    
    set(gca,'ytick',[1, 50, 300],'yticklabel',[-1 0 5],'fontsize',4)
    else
    set(gca,'ytick',[]);    
    end
 end
    
         
           
    end 
a =axes('Position',[0.35 0.05+0.15+((rr-1).*0.15)+0.1 0.07 0.10]);
t1 = title(strcat('HPCTimeSeq','Animal',num2str(ii),'Session',num2str(ll),'Track',num2str(nn)),'fontsize',6);
a.Visible = 'off'; % set(a,'Visible','off');
t1.Visible = 'on'; % set(t1,'Visible','on');
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,4,4])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300',strcat('HPCTimeSeq','Animal',num2str(ii),'Session',num2str(ll),'Track',num2str(nn)));
print(fig1,'-painters','-depsc','-r300',strcat('HPCTimeSeq','Animal',num2str(ii),'Session',num2str(ll),'Track',num2str(nn))); 
close all;        
cd(dn1)          
 end           
 end
 end    
end

%% Fig 3 c 
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('InstantaneousTimeErrHPC')

Err=[]; gg=0;Prediction=[];
for aa=1:25
    if isempty(MedErr{aa,1}) ; continue; end
    gg=gg+1;tt=cell2mat(cellfun(@(x)(~isempty(x)),MedErr,'UniformOutput',false));
    tt1=tt(aa,tt(aa,:)~=0);nn=0;
    if any(aa==[1:5,21:25,6,8,15]); dd=0; else dd=1;end
    for bb=[1,2,length(tt1)-dd, length(tt1)]
    nn=nn+1;Err{aa,nn}=[];    
    for a=1:size(MedErr{aa,bb},2)  
    ss= 1:length(MedErr{aa,bb}{a});
    ppred=[];Err11=[];tt=0;
    for b=ss
       
    TotPred=Predictions{aa,bb}{a,b}; 
    if isempty(TotPred); continue; end
     tt=tt+1;
    ppred(tt,:,:) =TotPred;
    Err11=[Err11; MedErr{aa,bb}{a}(b)];
    end
    ppred1=nanmean(ppred,1);
    Prediction(gg,nn,a,:,:)=squeeze(ppred1(1,:,:));
    Err{aa,nn}=[Err{aa,nn}; (Err11)];
    end
    end
end


ErrS=[];  gg=0;Prediction=[];
for aa=1:25
    if isempty(MedErrShuff{aa,1}) ; continue; end
    gg=gg+1;tt=cell2mat(cellfun(@(x)(~isempty(x)),MedErr,'UniformOutput',false));
    tt1=tt(aa,tt(aa,:)~=0);nn=0;
    if any(aa==[1:5,21:25,6,8,15]); dd=0; else dd=1;end
    for bb=[1,2,length(tt1)-dd, length(tt1)]
    nn=nn+1;ErrS{aa,nn}=[];    
    for a=1:size(MedErrShuff{aa,bb},2)  
    ss= 1:size(MedErrShuff{aa,bb}{a},1);
    ppred=[];Err11=[];tt=0;
    for b=ss
       
    TotPred=Predictions{aa,bb}{a,b}; 
    if isempty(TotPred); continue; end
     tt=tt+1;
    ppred(tt,:,:) =TotPred;
    Err11=[Err11; MedErrShuff{aa,bb}{a}(b,:)'];
    end
    ppred1=nanmean(ppred,1);
    Prediction(gg,nn,a,:,:)=squeeze(ppred1(1,:,:));
    ErrS{aa,nn}=[ErrS{aa,nn}; (Err11)];
    end
    end
end

Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
for ii=1:5
    MPred{ii}=[];
    MPredS{ii}=[];
end

for ii=1:25
      
      for nn=1:length(Err(ii,:));
      MPred{Animals(ii)}=[MPred{Animals(ii)}; (Err{ii,nn})];
      MPredS{Animals(ii)}=[MPredS{Animals(ii)}; (ErrS{ii,nn})];
      end
            
end
col1=[rgb('Red'); rgb('Gray')];
fig1=figure
rr=0;subplot(1,2,1)      

violin(MPred([1,2,4,5]),'x',[1 2 3 4],'facecolor',rgb('White'),'edgecolor',col1(1,:),...
'bw',0.1,'mc','k','medc','r-.')
axis([0.5 5 -0.5 5]); hold on

[h L Mex Med bw]=violin(MPredS([1,2,4,5]),'x',[1 2 3 4]+0.25,'facecolor',rgb('White'),'edgecolor',col1(2,:),...
'bw',0.1,'mc','k','medc','r-.')
set(gca, 'xtick',[1,1.25,2,2.25,3,3.25,4,4.25],'xticklabel',{'Rat1','Rat1 shuff','Rat2', 'Rat2 shuff','Rat4','Rat4 shuff','Rat5','Rat5 shuff'},'fontsize',5,'xticklabelrotation',45)
set(gca,'tickdir','out','box','off','fontsize',5);
xlim([0 5])
ylim([0 4.5])
ylabel('Decoding error (s)','fontsize',5);
rr=0;
for ii=1:5
      if isempty(MPred{ii}); continue; end
    rr=rr+1;
    [pval h ]=ranksum(MPred{ii},MPredS{ii})
    if pval < 0.05
        h1=sigstar([rr rr+0.25],pval,0)
    end
end
L.FontSize=5;

cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('InstantaneousTimeErrPFC')

Err=[];  gg=0;Prediction=[];
for aa=1:25
    if isempty(MedErr{aa,1}) ; continue; end
    gg=gg+1;tt=cell2mat(cellfun(@(x)(~isempty(x)),MedErr,'UniformOutput',false));
    tt1=tt(aa,tt(aa,:)~=0);nn=0;
    if any(aa==[1:5,21:25,6,8,15]); dd=0; else dd=1;end
    for bb=[1,2,length(tt1)-dd, length(tt1)]
    nn=nn+1;Err{aa,nn}=[];    
    for a=1:size(MedErr{aa,bb},2)  
    ss= 1:length(MedErr{aa,bb}{a});
    ppred=[];Err11=[];tt=0;
    for b=ss
       
    TotPred=Predictions{aa,bb}{a,b}; 
    if isempty(TotPred); continue; end
     tt=tt+1;
    ppred(tt,:,:) =TotPred;
    Err11=[Err11; MedErr{aa,bb}{a}(b)];
    end
    ppred1=nanmean(ppred,1);
    Prediction(gg,nn,a,:,:)=squeeze(ppred1(1,:,:));
    Err{aa,nn}=[Err{aa,nn}; (Err11)];
    end
    end
end


ErrS=[];  gg=0;Prediction=[];
for aa=1:25
    if isempty(MedErrShuff{aa,1}) ; continue; end
    gg=gg+1;tt=cell2mat(cellfun(@(x)(~isempty(x)),MedErr,'UniformOutput',false));
    tt1=tt(aa,tt(aa,:)~=0);nn=0;
    if any(aa==[1:5,21:25,6,8,15]); dd=0; else dd=1;end
    for bb=[1,2,length(tt1)-dd, length(tt1)]
    nn=nn+1;ErrS{aa,nn}=[];    
    for a=1:size(MedErrShuff{aa,bb},2)  
    ss= 1:size(MedErrShuff{aa,bb}{a},1);
    ppred=[];Err11=[];tt=0;
    for b=ss
       
    TotPred=Predictions{aa,bb}{a,b}; 
    if isempty(TotPred); continue; end
     tt=tt+1;
    ppred(tt,:,:) =TotPred;
    Err11=[Err11; MedErrShuff{aa,bb}{a}(b,:)'];
    end
    ppred1=nanmean(ppred,1);
    Prediction(gg,nn,a,:,:)=squeeze(ppred1(1,:,:));
    ErrS{aa,nn}=[ErrS{aa,nn}; (Err11)];
    end
    end
end

Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
for ii=1:5
    MPred{ii}=[];
    MPredS{ii}=[];
end

for ii=1:25
      
      for nn=1:length(Err(ii,:));
      MPred{Animals(ii)}=[MPred{Animals(ii)}; (Err{ii,nn})];
      MPredS{Animals(ii)}=[MPredS{Animals(ii)}; (ErrS{ii,nn})];
      end
            
end
rr=0;subplot(1,2,2)      

violin(MPred([1,3,4,5]),'x',[1 2 3 4],'facecolor',rgb('White'),'edgecolor',col1(1,:),...
'bw',0.1,'mc','k','medc','r-.')
axis([0.5 5 -0.5 5]); hold on

[h L Mex Med bw]=violin(MPredS([1,3,4,5]),'x',[1 2 3 4]+0.25,'facecolor',rgb('White'),'edgecolor',col1(2,:),...
'bw',0.1,'mc','k','medc','r-.')
set(gca, 'xtick',[1,1.25,2,2.25,3,3.25,4,4.25],'xticklabel',{'Rat1','Rat1 shuff','Rat3', 'Rat3 shuff','Rat4','Rat4 shuff','Rat5','Rat5 shuff'},'fontsize',5,'xticklabelrotation',45)
set(gca,'tickdir','out','box','off','fontsize',5);
xlim([0 5])
ylim([0 4.5])
ylabel('Decoding error (s)','fontsize',5);
rr=0;
for ii=1:5
      if isempty(MPred{ii}); continue; end
    rr=rr+1;
    [pval h ]=ranksum(MPred{ii},MPredS{ii})
    if pval < 0.05
        h1=sigstar([rr rr+0.25],pval,0)
    end
end
L.FontSize=5;
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,2,1])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','DecodingErrorsComparison');
print(fig1,'-painters','-depsc','-r300','DecodingErrorsComparison'); 
close all;

%% Fig 3d
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('InstantaneousTimeErrHPC')

Err=[]; gg=0;Prediction=[];
 bins=[1:25:151];  
for aa=1:25
    if isempty(MedErr{aa,1}) ; continue; end
    gg=gg+1;tt=cell2mat(cellfun(@(x)(~isempty(x)),MedErr,'UniformOutput',false));
    tt1=tt(aa,tt(aa,:)~=0);nn=0;
    if any(aa==[1:5,21:25,6,8,15]); dd=0; else dd=1;end
    for bb=[1,2,length(tt1)-dd, length(tt1)]
    nn=nn+1;Err{aa,nn}=[];  ErrS{aa,nn}=[];  
    for a=1:size(MedErr{aa,bb},2)  
    ss= 1:length(MedErr{aa,bb}{a});
    ppred=[];Err11=[];tt=0;
    for b=ss
    binsofInterest=bins(a):bins(a+1)-1;    
    TotPred=nanmean(sum(TotPredicitions{aa,bb}{a,b}(binsofInterest,:))); 
    TotPredS=[];
    for ss1=1:500
 ntimebin=1:size(TotPredicitions{aa,bb}{a,b},1);   
 ntimebin = ntimebin(randperm(length(ntimebin)));
  TotPredS(ss1,1)=nanmean(sum(TotPredicitions{aa,bb}{a,b}(ntimebin(1:25),:))); 
    end  
    
    Err{aa,nn}=[Err{aa,nn}; (TotPred)];
    ErrS{aa,nn}=[ErrS{aa,nn}; (TotPredS)];
    end
    end
end
end


Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
for ii=1:5
    MPred{ii}=[];
    MPredS{ii}=[];
end

for ii=1:25
      
      for nn=1:length(Err(ii,:));
      MPred{Animals(ii)}=[MPred{Animals(ii)}; (Err{ii,nn})];
      MPredS{Animals(ii)}=[MPredS{Animals(ii)}; (ErrS{ii,nn})];
      end
            
end
col1=[rgb('Red'); rgb('Gray')];
fig1=figure
rr=0;subplot(1,2,1)      
violin(MPred([1,2,4,5]),'x',[1 2 3 4],'facecolor',rgb('White'),'edgecolor',col1(1,:),...
'mc','k','medc','r-.')
axis([0.5 5 -0.5 5]); hold on

[h L Mex Med bw]=violin(MPredS([1,2,4,5]),'x',[1 2 3 4]+0.25,'facecolor',rgb('White'),'edgecolor',col1(2,:),...
'mc','k','medc','r-.')
set(gca, 'xtick',[1,1.25,2,2.25,3,3.25,4,4.25],'xticklabel',{'Rat1','Rat1 shuff','Rat2', 'Rat2 shuff','Rat4','Rat4 shuff','Rat5','Rat5 shuff'},'fontsize',5,'xticklabelrotation',90)
set(gca,'tickdir','out','box','off','fontsize',5);
xlim([0 5])
ylim([0 1])
ylabel('Predictions','fontsize',5);
rr=0;
for ii=1:5
      if isempty(MPred{ii}); continue; end
    rr=rr+1;
    [pval h ]=ranksum(MPred{ii},MPredS{ii})
    if pval < 0.05
        h1=sigstar([rr rr+0.25],pval,0)
    end
end
L.FontSize=5;


cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('InstantaneousTimeErrPFC')

Err=[];  gg=0;Prediction=[];
 bins=[1:25:151];  
for aa=1:25
    if isempty(MedErr{aa,1}) ; continue; end
    gg=gg+1;tt=cell2mat(cellfun(@(x)(~isempty(x)),MedErr,'UniformOutput',false));
    tt1=tt(aa,tt(aa,:)~=0);nn=0;
    if any(aa==[1:5,21:25,6,8,15]); dd=0; else dd=1;end
    for bb=[1,2,length(tt1)-dd, length(tt1)]
    nn=nn+1;Err{aa,nn}=[];  ErrS{aa,nn}=[];  
    for a=1:size(MedErr{aa,bb},2)  
    ss= 1:length(MedErr{aa,bb}{a});
    ppred=[];Err11=[];tt=0;
    for b=ss
    binsofInterest=bins(a):bins(a+1)-1;    
     TotPred=nanmean(sum(TotPredicitions{aa,bb}{a,b}(binsofInterest,:))); 
    TotPredS=[];
    for ss1=1:500
 ntimebin=1:size(TotPredicitions{aa,bb}{a,b},1);   
 ntimebin = ntimebin(randperm(length(ntimebin)));
  TotPredS(ss1,1)=nanmean(sum(TotPredicitions{aa,bb}{a,b}(ntimebin(1:25),:))); 
    end  
    
    Err{aa,nn}=[Err{aa,nn}; (TotPred)];
    ErrS{aa,nn}=[ErrS{aa,nn}; (TotPredS)];
    end
    end
end
end
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
for ii=1:5
    MPred{ii}=[];
    MPredS{ii}=[];
end

for ii=1:25
      
      for nn=1:length(Err(ii,:));
      MPred{Animals(ii)}=[MPred{Animals(ii)}; (Err{ii,nn})];
      MPredS{Animals(ii)}=[MPredS{Animals(ii)}; (ErrS{ii,nn})];
      end
            
end
rr=0;subplot(1,2,2)      

violin(MPred([1,3,4,5]),'x',[1 2 3 4],'facecolor',rgb('White'),'edgecolor',col1(1,:),...
'mc','k','medc','r-.')
axis([0.5 5 -0.5 5]); hold on

[h L Mex Med bw]=violin(MPredS([1,3,4,5]),'x',[1 2 3 4]+0.25,'facecolor',rgb('White'),'edgecolor',col1(2,:),...
'mc','k','medc','r-.')
set(gca, 'xtick',[1,1.25,2,2.25,3,3.25,4,4.25],'xticklabel',{'Rat1','Rat1 shuff','Rat3', 'Rat3 shuff','Rat4','Rat4 shuff','Rat5','Rat5 shuff'},'fontsize',5,'xticklabelrotation',90)
set(gca,'tickdir','out','box','off','fontsize',5);
xlim([0 5])
ylim([0 1])
ylabel('Predictions','fontsize',5);
rr=0;
for ii=1:5
      if isempty(MPred{ii}); continue; end
    rr=rr+1;
    [pval h ]=ranksum(MPred{ii},MPredS{ii})
    if pval < 0.05
        h1=sigstar([rr rr+0.25],pval,0)
    end
end
L.FontSize=5;
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,2,1])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','DecodingErrorsCorrectVsShuffComparison');
print(fig1,'-painters','-depsc','-r300','DecodingErrorsCorrectVsShuffComparison'); 
close all;


%% Fig 3e-i
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,5,5,5,5,5,4,4,4,4,4]);
ll=1;glmstats=[];GoalProbsStacked=[];Goalbins=[];Goaldprobp=[];

Ex={[10,4,7],[2,3,5],[19,4,2]};
col1={rgb('Red'),rgb('Green'),rgb('Blue'),...
    rgb('Black'),rgb('Cyan'),rgb('Magenta'),rgb('Gray')};
for ii=Ex{3}(1)
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
dn1=pwd;
load LapInfo
load PlFieldTempNLaps
load LFPDataV2
dn=pwd;
LapInfo1=LapInfo;
load mazel
load CenterTimePos.mat

maze=[];
for i=1:length(mazel)
    x=Center_Pos{i}(:,1);y=Center_Pos{i}(:,2);
    T=Center_Time{i};
    vx = dxdt(T,x./2.5);
    vy = dxdt(T,y./2.5);
    speed1=sqrt(vx.^2+vy.^2)';
    Centerl=ones(size(speed1,1),size(mazel{i},2)).*NaN;
    Centerl(:,[1:2,6,5])=[x,y,speed1,T];
    maze{i}=[mazel{i}; Centerl];
end

cluster=importdata('clusters');
elecreg=importdata('elecreg');
tic;
ss=ii-5;
% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % good is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates10(ind_HPC,1);
sel=1:100;
if sum(ind_HPC) < 10; continue; end
sigma = 2; % pick sigma value for the gaussian
gaussFilter = gausswin(3*sigma + 1)';
gaussFilter = gaussFilter / sum(gaussFilter); % normalize

   
for ll=Ex{3}(2)
tic;  
LapInfo = LapInfo1;
plf1=[];plfConcat=[];
arms=[1,2,4,5,7,8];
for nn=1:length(rates) 
plf11=rates{nn};
plf12=[];
ind=[16:40];
for jj=1:6
    plf1=plf11(LapInfo(:,5)==ll & LapInfo(:,6) ==arms(jj) & LapInfo(:,9) ==1,ind);
%     plf1=conv2( plf1, gaussFilter, 'same');
    plf12=[plf12, conv2(nanmean(plf1), gaussFilter, 'same')];       
end
plfConcat=[plfConcat; (plf12)];
end

    LapInfo = LapInfo1;
    if ii==12 & ll==4
    LapInfo=[LapInfo; ones(1,size(LapInfo,2)).*NaN];
    LapInfo(end,9)=0; LapInfo(end,5)=ll;Center_Epochs{end+1}=[];
    LapInfo1=LapInfo;
    end
    
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==0));
    ind1=find(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    
    EpochsTrialID=find(LapInfo1(:,5)==ll & (LapInfo1(:,9)==0));    
   
    target_info= LapInfo(:,10);
    iid=target_info>=3 & target_info<=6;
    target_info(iid,:)=target_info(iid,:)-1;
     iid=target_info>6 & target_info <=9;
    target_info(iid,:)=target_info(iid,:)-2;
    
    current_info=LapInfo(:,6);
    next_info= LapInfo1(ind1+1,6);
     iid=next_info>=3 & next_info<=6;
    next_info(iid,:)=next_info(iid,:)-1;
     iid=next_info>6 & next_info <=9;
    next_info(iid,:)=next_info(iid,:)-2;
    
    
      last_info= [NaN; LapInfo1(ind1(2:end)-1,6)];
       iid=last_info>=3 & last_info<=6;
    last_info(iid,:)=last_info(iid,:)-1;
     iid=last_info>6 & last_info <=9;
    last_info(iid,:)=last_info(iid,:)-2;
      
       MeanProb=[];
    
       
    
 for rr=Ex{3}(3)   
   t_stimulus=LapInfo(rr,7);
   t_bin=0.05;
   tbin=0.02;
   tvec2=[t_stimulus-1   t_stimulus+5];       
   binns=nanmin(tvec2):tbin:nanmax(tvec2);
   binedges=(binns(1:end-1)+binns(2:end))/2; 
cd(dn1);   
kkk=1:length(rates);
% kk=kkk(ind_HPC);   
trackcltime=[];
b11=cluster;
b11=b11(ind_HPC);
for c=1:size(b11,1)
       clG=cl2mat(b11{c});              
       trackcltime{c}= [clG.featuredata((clG.featuredata(:,8)>=tvec2(1) & clG.featuredata(:,8)<=tvec2(2)),[8,11]) ...
           repmat(kkk(c),size(clG.featuredata((clG.featuredata(:,8)>=tvec2(1) & clG.featuredata(:,8)<=tvec2(2)),8),1),1)];
       clear clG
end

  cltimevec=[];
   
   for c=1:size(trackcltime,2)
   cltimevec=[cltimevec; trackcltime{c}(:,[1,3])];
   end
   cltimemat=[];
   [~,ind]=sort(cltimevec(:,1));
   cltimemat=cltimevec(ind,:);
%    
[dprobp dprob binedges]=BayesDecode_FramesWithOverlap(tvec2,plfConcat,cltimevec,tbin,0.5);

fig1=figure

left=0.3; bottom=0.24; width=0.3; height=0.3;
ax=axes('position',[left bottom width height])

zz1=dprobp;
zz1(isnan(zz1))=0;        
% create filter
sigma = 2; % pick sigma value for the gaussian
gaussFilter = gausswin(5*sigma + 1)';
gaussFilter = gaussFilter / sum(gaussFilter); % normalize
bins=[1:25:151];zz=[];
for yy=1:size(bins,2)-1
zz(bins(yy):bins(yy+1)-1,:)= conv2(zz1(bins(yy):bins(yy+1)-1,:)', gaussFilter, 'same')'; 
end


zz=movmean(zz,[2 4],1,'omitnan');  
zz=movmean(zz,[2 4],2,'omitnan'); 

imagesc(zz,[nanmax(nanmax(zz))/10 nanmax(nanmax(zz))/3])
set(gca,'YDir','normal')
cmap=getNCLColmap('MPL_Greys.rgb',100);
% cmap2=getNCLColmap('hlu_default.rgb',30);
shading interp
colormap(cmap); 
cb=colorbar;
set(cb,'Position',[0.6 0.22 0.01 0.05])
cb.FontSize=3;
hold on;

hold on;   
binn=[1:25:151];prob_means=[];
for yy=1:6    
x = [0 size(zz,2) size(zz,2) 0];
y = [binn(yy) binn(yy) binn(yy+1) binn(yy+1)];
if yy==(next_info(rr)) & yy==(target_info(rr))
pp=patch(x,y,col1{1});  
pp.EdgeColor = 'none';
alpha(pp,0.1); hold on
elseif yy==(target_info(rr)) 
pp1=patch(x,y,col1{2});  
pp1.EdgeColor = 'none';
alpha(pp1,0.1); hold on
elseif yy==(next_info(rr)) 
pp2=patch(x,y,col1{3});
pp2.EdgeColor = 'none';
alpha(pp2,0.1); hold on
elseif yy==(last_info(rr)) & rr>1
pp2=patch(x,y,col1{6});
pp2.EdgeColor = 'none';
alpha(pp2,0.1); hold on


end
if yy> 1 & yy<7
   plot([0 size(zz,2)],[binn(yy) binn(yy)],'--k'); hold on;
end
sigma = 4; % pick sigma value for the gaussian
gaussFilter = gausswin(6*sigma + 1)';
gaussFilter = gaussFilter / sum(gaussFilter); % normalize
prob_means(yy,:)=conv2(nansum(zz(binn(yy):binn(yy+1)-1,:)), gaussFilter, 'same');
end    
ylabel('Predicted time at goals','fontsize',6);
set(gca,'ytick',(binn(1:end-1)+binn(2:end))/2,'yticklabel',{'A','B','C','D','E','F'},'fontsize',4);
xlabel('Time at start (s)','fontsize',6);  
set(gca,'xtick',[1 100 200 300 400 500 599],'xticklabel',[-1 0 1 2 3 4 5],'fontsize',5);
ax=axes('position',[left bottom-0.12 width 0.10]);
if rr > 1
plot([1:size(prob_means,2)], prob_means(last_info(rr),:),'color',col1{6},'linewidth',2); hold on;
end
if next_info(rr) == target_info(rr)
plot([1:size(prob_means,2)], prob_means(target_info(rr),:),'color',col1{1},'linewidth',2); hold on;
else
plot([1:size(prob_means,2)], prob_means(next_info(rr),:),'color',col1{3},'linewidth',2); hold on;
plot([1:size(prob_means,2)], prob_means(target_info(rr),:),'color',col1{2},'linewidth',2); hold on;
end

sst=setdiff([1,2,3,4,5,6],[target_info(rr),last_info(rr)]);
for ee=1:size(sst,2)
plot([1:size(prob_means,2)], (prob_means(sst(ee),:)),'color',col1{7},'linewidth',1); hold on;
end
yl=get(gca,'ylim');
xt=get(gca,'xtick');
set(gca,'box','off','xtick',[1 100 200 300 400 500 599],'xticklabel',[-1 0 1 2 3 4 5],...
        'ytick',[0 yl(end)],'YTickLabel', [0 round(yl(end),2)],'TickLength', [0 0],...
            'fontsize',4);
        xlim([0 size(zz,2)]);
xlabel('Time (s)','fontsize',4)
ylabel('Prediction','fontsize',4);  

ax=axes('position',[left+width+0.05 bottom-0.12 0.05 0.10]);
hh(1)=bar([1],[nanmean(prob_means(target_info(rr),:))]); hold on
hh(2)=bar([2],[nanmean(prob_means(last_info(rr),:))]); 
hh(3)=bar([3],[nanmean(nanmean(prob_means(sst,:)))]); 
hh(1).FaceColor='r';hh(2).FaceColor=col1{6};hh(3).FaceColor=col1{7};
hh1=errorbar([1:3],[nanmean(prob_means(target_info(rr),:)),nanmean(prob_means(last_info(rr),:)),nanmean(nanmean(prob_means(sst,:)))],...
    [nansem(prob_means(target_info(rr),:)),nansem(prob_means(last_info(rr),:)),nansem(nanmean(prob_means(sst,:)))], 'capsize',3,'linestyle','none','color','k')
Groups=[1,2;1,3;2,3];pval=[];
for ee=1:size(Groups,1)
    if ee==1
    [pval(ee) h4 ]=ranksum(prob_means(target_info(rr),:),prob_means(last_info(rr),:));
    elseif ee==2
    [pval(ee) h4 ]=ranksum(prob_means(target_info(rr),:),nanmean(prob_means(sst,:)));
    else
    [pval(ee) h4 ]=ranksum(prob_means(last_info(rr),:),nanmean(prob_means(sst,:)));
    end

if pval<0.05
    hh5=sigstar(Groups(ee,:),pval(ee),0);
end
end
set(ax,'tickdir','out','box','off','fontsize',4);
ylabel('Probability','fontsize',4);
set(gca,'xtick',[]);


        i=ll*2;
        indx=LFPData(i).t >tvec2(1) & LFPData(i).t <tvec2(2); 
        LFPamp=LFPData(i).LFP(indx);
        LFPtvec=LFPData(i).t(indx);
%         LFPripple=LFPData(i).firLFP{1, 3}(indx);
        LFPtheta=LFPData(i).firLFP{1, 4}(indx);
        ax=axes('position',[left bottom-0.19 width 0.05]);
        plot(LFPtvec,LFPtheta,'-k');
        set(gca,'box','off','XTick',[tvec2(1)  tvec2(1)+1 tvec2(1)+2 tvec2(1)+3 tvec2(1)+4 tvec2(1)+5 tvec2(end) ],'TickLength', [0 0],'xticklabel',[-1:5]);
        set(gca,'YTick',[nanmin(LFPtheta) nanmax(LFPtheta)],'YTickLabel', [nanmin(LFPtheta) nanmax(LFPtheta)],...
            'fontsize',4);
        xlim([tvec2(1) tvec2(2)])
        ylabel('Amplitude','fontsize',4);
        
        indx=maze{ll}(:,5) >=tvec2(1) &  maze{ll}(:,5) <=tvec2(end);
        totalxy=maze{ll}(:,[1,2]);
        [totalxy,outside1]=set_outside_V2(totalxy);
        animalxy1=maze{ll}(indx,[1,2]);
        animalxy1=nanmedian(animalxy1);
        vel=maze{ll}(indx,6);        
        mazetvec=maze{ll}(indx,5);
        ax=axes('position',[left-0.25 0.3 0.2 0.2])       
        plot(totalxy(:,1),totalxy(:,2),'linestyle','none','marker','o','markersize',1,'color','k');hold on;
     

        set(gca,'xlim',[0 630],'ylim',[0 500],'box','off'...
            ,'XTickLabel', [],'TickLength', [0 0],'xtick',[],'XColor','none','YColor','none');
       
       
       tvec_next=[LapInfo1(EpochsTrialID(rr)+1,1); LapInfo1(EpochsTrialID(rr)+1,4)];
       tvec_target=[LapInfo1(EpochsTrialID(rr+1)-1,1); LapInfo1(EpochsTrialID(rr+1)-1,2)];
           
       if rr==1
       tvec_all=[LapInfo1(EpochsTrialID(rr),3); LapInfo1(EpochsTrialID(rr+1)-1,2)];
       else
       tvec_last=[LapInfo1(EpochsTrialID(rr)-1,3); LapInfo1(EpochsTrialID(rr)-1,4)];    
       tvec_all=[tvec_last(1); LapInfo1(EpochsTrialID(rr+1)-1,2)];
       end
       hold on;
       indx=maze{ll}(:,5) >=tvec_all(1) &  maze{ll}(:,5) <=tvec_all(end);
       animalxy=maze{ll}(indx,[1,2]);
       [animalxy,outside1]=set_outside_V2(animalxy);
       plot(animalxy(:,1),animalxy(:,2),'linestyle','none','marker','o','markersize',3,'color',col1{5})
       if (next_info(rr)) ==(target_info(rr))
       indx=maze{ll}(:,5) >=tvec_target(1) &  maze{ll}(:,5) <=tvec_target(end);
       animalxy=maze{ll}(indx,[1,2]);
       [animalxy,outside1]=set_outside_V2(animalxy);
       plot(animalxy(:,1),animalxy(:,2),'linestyle','none','marker','o','markersize',3,'color',col1{1})
       else
       indx=maze{ll}(:,5) >=tvec_next(1) &  maze{ll}(:,5) <=tvec_next(end);
       animalxy=maze{ll}(indx,[1,2]);
       [animalxy,outside1]=set_outside_V2(animalxy);
       plot(animalxy(:,1),animalxy(:,2),'linestyle','none','marker','o','markersize',3,'color',col1{3})
       indx=maze{ll}(:,5) >=tvec_target(1) &  maze{ll}(:,5) <=tvec_target(2);
        animalxy=maze{ll}(indx,[1,2]);
        [animalxy,outside1]=set_outside_V2(animalxy);
       plot(animalxy(:,1),animalxy(:,2),'linestyle','none','marker','o','markersize',3,'color',col1{2})    
       end 
       if rr>1
       indx=maze{ll}(:,5) >=tvec_last(1) &  maze{ll}(:,5) <=tvec_last(end);
       animalxy=maze{ll}(indx,[1,2]);
       [animalxy,outside1]=set_outside_V2(animalxy);
       plot(animalxy(:,1),animalxy(:,2),'linestyle','none','marker','o','markersize',3,'color',col1{6})
       end
       plot(animalxy1(:,1),animalxy1(:,2),'linestyle','none','marker','pentagram','markersize',10,'color','y','markerfacecolor','y');
        text(190,450,'S2','fontsize',10,'color','r','FontWeight','bold');
        text(370,450,'D','fontsize',10,'color','r','FontWeight','bold');        
        text(571,373,'C','fontsize',10,'color','r','FontWeight','bold');
        text(587,179,'S1','fontsize',10,'color','r','FontWeight','bold');
        text(448,30,'B','fontsize',10,'color','r','FontWeight','bold');
        text(237,32,'A','fontsize',10,'color','r','FontWeight','bold');
        text(74,144,'F','fontsize',10,'color','r','FontWeight','bold');
        text(51,332,'E','fontsize',10,'color','r','FontWeight','bold');
ax=axes('position',[left bottom+0.32 width 0.1])         
spktimes=cltimevec(:,1);cellid=cltimevec(:,2);
for r=1:length(spktimes)
line([spktimes(r) spktimes(r)],[cellid(r)-0.5 cellid(r)+0.5],'linestyle','-','color','r','linewidth',0.5)
end
ylim([0 max(cellid)+0.5]);
xlim([tvec2(1) tvec2(2)]);
set(gca,'box','off','XTick',[tvec2(1)  tvec2(1)+1 tvec2(1)+2 tvec2(1)+3 tvec2(1)+4 tvec2(1)+5 tvec2(end) ],'TickLength', [0 0],'xticklabel',[-1:5],...
    'YTickLabel', [],'TickLength', [0 0],'ytick',[]); hold on;
plot([tvec2(1)+1 tvec2(1)+1],[0 max(cellid)+0.5],'--k');
ylabel('cell IDs','fontsize',6);
set(gca,'fontsize',4);
ArmsAlpha={'A','B','C','D','E','F'};        
        
 if rr>1      
 annotation('textbox', [left-0.02 0.69, 0, 0], 'string', ...
    strcat(['StartArm',num2str(current_info(rr)),'GoingToArm',ArmsAlpha{next_info(rr,1)}, ...
    'TargetArm',ArmsAlpha{target_info(rr,1)},'LastArm',ArmsAlpha{last_info(rr,1)},'Animal',num2str(Animals(ii)),...
    'Day',num2str(Days(ii)),'Session',num2str(ll),...
    'Lap',num2str(rr)])...
    ,'fontsize',6);
 else
  annotation('textbox', [left-0.02 0.69, 0, 0], 'string', ...
    strcat(['StartArm',num2str(current_info(rr)),'GoingToArm',ArmsAlpha{next_info(rr,1)}, ...
    'TargetArm',ArmsAlpha{target_info(rr,1)},'LastArm','None','Animal',num2str(Animals(ii)),...
    'Day',num2str(Days(ii)),'Session',num2str(ll),...
    'Lap',num2str(rr)])...
    ,'fontsize',6);
 end

set(fig1,'paperunits','inches');
set(fig1,'papertype','usletter');
set(gcf,'paperposition',[0.5,0.5,5,5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300',strcat('Animal',num2str(Animals(ii)),'Day',num2str(Days(ii)),'Session',num2str(ll),'LapNumber',num2str(rr),'StartArm',num2str(current_info(rr)),'Prediction'));        
print(fig1,'-depsc','-painters','-r300',strcat('Animal',num2str(Animals(ii)),'Day',num2str(Days(ii)),'Session',num2str(ll),'LapNumber',num2str(rr),'StartArm',num2str(current_info(rr)),'Prediction'));        
close(fig1);
close all          
end
toc
end
end


%% Fig 3j
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('Spec_power100.mat','Spec_power')

Spectrum1=[];Spectrum2=[];

    for Days=1:2
        for jj=1:3
            for kk=1:3
                Spectrum1{Days,jj}=[];
                Spectrum2{Days,jj}=[];
                Spectrum3{Days,jj}=[];
                Spectrum4{Days,jj}=[];
                Spectrum5{Days,jj}=[];
                Spectrum6{Days,jj}=[];
            end
        end
    end


Sess={[1:5]; [6:10, 11:15, 16:20 ]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};

for Days=1:2
for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),(Spec_power(ii,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
         Err=Spec_power{ii,jj}{1,6};
            x2=[];x4=[];
        for kk=1:size(Err,1)
            
             if Err(kk) == 0
             x1=(Spec_power{ii,jj}{1,1}(kk,:));
             x2=[x2; x1]; 
             else
             x1=(Spec_power{ii,jj}{1,1}(kk,:));
             x4=[x4; x1];
             end
                
        end
       
        Spectrum1{Days,ll}=[Spectrum1{Days,ll}; (x2)]; 
        Spectrum2{Days,ll}=[Spectrum2{Days,ll}; (x4)]; 
        

        
        x2=[];x4=[];
        for kk=1:size(Spec_power{ii,jj}{1,2},1)
            
             
             x1=(Spec_power{ii,jj}{1,2}(kk,:));
             x2=[x2; x1]; 
            
                
        end
         for kk=1:size(Spec_power{ii,jj}{1,3},1)            
             
             x1=(Spec_power{ii,jj}{1,3}(kk,:));
             x4=[x4; x1]; 
            
                
        end
        
        
        Spectrum3{Days,ll}=[Spectrum3{Days,ll}; (x2)]; 
        Spectrum4{Days,ll}=[Spectrum4{Days,ll}; (x4)];
        
        
             x2=[];x4=[];
        for kk=1:size(Spec_power{ii,jj}{1,4},1)
            
             
             x1=(Spec_power{ii,jj}{1,4}(kk,:));
             x2=[x2; x1]; 
            
                
        end
         for kk=1:size(Spec_power{ii,jj}{1,5},1)            
             
             x1=(Spec_power{ii,jj}{1,5}(kk,:));
             x4=[x4; x1]; 
            
                
        end
        
        
        Spectrum5{Days,ll}=[Spectrum5{Days,ll}; (x2)]; 
        Spectrum6{Days,ll}=[Spectrum6{Days,ll}; (x4)];
       
        
     end
end

end
Spectrum_Pw1=[];
Spectrum_Pw1{1}=[Spectrum1{2,3}; Spectrum2{2,3}];
% Spectrum_Pw1{2}=[movmean(Spectrum2{2,3},[100 100],2)];
Spectrum_Pw1{2}=Spectrum3{2,3};
% Spectrum_Pw1{3}=[movmean(Spectrum4{2,3},[100 100],2)];
% Spectrum_Pw1{5}=[movmean(Spectrum5{2,3},[100 100],2)];
% Spectrum_Pw1{3}=Spectrum6{2,3};


col1={rgb('Red'),rgb('Blue'),rgb('Black'),rgb('Orange'),rgb('Green'),rgb('Orange')};
fig1=figure(1);
ax1(1)=axes('Position',[0.1 0.35 0.10 0.15]); 
for a=1:size(Spectrum_Pw1,2)
    mm=nanmean((Spectrum_Pw1{a})./[100]);
    ms=nansem((Spectrum_Pw1{a})./[100]);
       H1(a)=plot(1:size(mm,2),mm ,'linewidth',1,'marker','none','color',col1{a}); hold on;
    pp=patch([[1:length(mm)] fliplr([1:length(mm)])], ...
       [ms+mm fliplr(mm-ms)], col1{a})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl=get(gca,'ylim');
    set(gca,'Xscale','log','Yscale','log','box','off','tickdir','out');
     set(gca,'fontsize',6);
    xlabel('Frequency (Hz)','fontsize',6);
	ylabel('Power','fontsize',6);
    set(gca,'xtick',[1 6 10 100 250],'xticklabel',[1 6 10 50 250],'fontsize',5)
end
ylim([0 max(yl)])
plot([6 6],[0 max(yl)],'--k');
plot([10 10],[0 max(yl)],'--k');
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter') 
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,3,3])   
print(gcf,'-djpeg','-r300','PowerSpectrumCompared'); 
print(gcf,'-painters','-depsc','-r300','PowerSpectrumCompared'); 
close(fig1)

%% Fig 3k-l
%% Ripple power comparison
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('Profile_ThetaRipplepowerZscore.mat','Spec_power')

Spectrum1=[];Spectrum2=[];

    for Days=1:2
        for jj=1:3
            for kk=1:3
                Spectrum1{Days,jj}=[];
                Spectrum2{Days,jj}=[];
                Spectrum3{Days,jj}=[];
                Spectrum4{Days,jj}=[];
                Spectrum5{Days,jj}=[];
                Spectrum6{Days,jj}=[];
            end
        end
    end


Sess={[1:5]; [6:10, 11:15, 16:20 ]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};

for Days=1:2
for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),(Spec_power(ii,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
         Err=Spec_power{ii,jj}{1,6};
            x2=[];x4=[];
        for kk=1:size(Err,1)
            
             if Err(kk) == 0
             x1=(Spec_power{ii,jj}{1,1}{2}(kk,:));
             x2=[x2; x1]; 
             else
             x1=(Spec_power{ii,jj}{1,1}{2}(kk,:));
             x4=[x4; x1];
             end
                
        end
       
        Spectrum1{Days,ll}=[Spectrum1{Days,ll}; (x2)]; 
        Spectrum2{Days,ll}=[Spectrum2{Days,ll}; (x4)]; 
        

        
        x2=[];x4=[];
        for kk=1:size(Spec_power{ii,jj}{1,2}{2},1)
            
             
             x1=(Spec_power{ii,jj}{1,2}{2}(kk,:));
             x2=[x2; x1]; 
            
                
        end
         for kk=1:size(Spec_power{ii,jj}{1,3}{2},1)            
             
             x1=(Spec_power{ii,jj}{1,3}{2}(kk,:));
             x4=[x4; x1]; 
            
                
        end
        
        
        Spectrum3{Days,ll}=[Spectrum3{Days,ll}; (x2)]; 
        Spectrum4{Days,ll}=[Spectrum4{Days,ll}; (x4)];
        
        
             x2=[];x4=[];
        for kk=1:size(Spec_power{ii,jj}{1,4}{2},1)
            
             
             x1=(Spec_power{ii,jj}{1,4}{2}(kk,:));
             x2=[x2; x1]; 
            
                
        end
         for kk=1:size(Spec_power{ii,jj}{1,5}{2},1)            
             
             x1=(Spec_power{ii,jj}{1,5}{2}(kk,:));
             x4=[x4; x1]; 
            
                
        end
        
        
        Spectrum5{Days,ll}=[Spectrum5{Days,ll}; (x2)]; 
        Spectrum6{Days,ll}=[Spectrum6{Days,ll}; (x4)];
       
        
     end
end

end
Spectrum_Pw1=[];
Spectrum_Pw1{1}=[movmean([Spectrum1{2,3}; Spectrum2{2,3}],[100 100],2)];
% Spectrum_Pw1{2}=[movmean(Spectrum2{2,3},[100 100],2)];
Spectrum_Pw1{2}=[movmean(Spectrum3{2,3},[100 100],2)];
% Spectrum_Pw1{3}=[movmean(Spectrum4{2,3},[100 100],2)];
% Spectrum_Pw1{5}=[movmean(Spectrum5{2,3},[100 100],2)];
Spectrum_Pw1{3}=[movmean(Spectrum6{2,3},[100 100],2)];

col1={rgb('Red'),rgb('Blue'),rgb('Black'),rgb('Orange'),rgb('Green'),rgb('Orange')};
fig1=figure(1);
ax1(1)=axes('Position',[0.1 0.55 0.15 0.2]); 
data1=[];labels1=[];
for a=1:size(Spectrum_Pw1,2)
    mm=movmean(nanmean((Spectrum_Pw1{a})),[50 50],2);
    ms=movmean(nansem((Spectrum_Pw1{a})),[50 50],2);
    mm1=(Spectrum_Pw1{a});
    mm1=mm1;
    mm2{a}=mm1(:);
    data1=[data1; nanmean(mm1)'];
    labels1=[labels1; ones(size(mm1,2),1).*a];
    H1(a)=plot(1:size(mm,2),mm ,'linewidth',1,'marker','none','color',col1{a}); hold on;
    pp=patch([[1:length(mm)] fliplr([1:length(mm)])], ...
       [ms+mm fliplr(mm-ms)], col1{a})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl=get(gca,'ylim');
    xlabel('Peri-event time (s)','fontsize',6);
	ylabel('Ripple Power (z-score)','fontsize',6);
    set(gca,'fontsize',6);
    set(gca,'xtick',[0:200:1200],'box','off','xticklabel',[-1:1:5],'fontsize',5);
    
    
end
plot([201 201],[min(yl) max(yl)],'--k');
ax1(2)=axes('Position',[0.1 0.15 0.15 0.25]); 
H=boxplot(data1,labels1,'boxstyle','outline','whisker',2, 'colors',[col1{1};col1{2};col1{3}] );

Groups=[1,2;1,3;2,3];
for ss1=1:size(Groups,1)
    [p h]=ranksum(mm2{Groups(ss1,1)},mm2{Groups(ss1,2)});
    if p<0.05
        h1=sigstar(Groups(ss1,:),p,0)
    end
end
set(gca,'fontsize',5,'box','off','xtick',[1,2,3],'xticklabel',{'Cue','Goal','Movement'},'xticklabelrotation',90);
xlabel('Events','fontsize',5);
ylabel('Ripple power (z-score)','fontsize',5);


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter') 
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,3,3])   
print(gcf,'-djpeg','-r300','RipplePowerComparedZscore'); 
print(gcf,'-painters','-depsc','-r300','RipplePowerComparedZscore'); 
close(fig1)

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('Profile_ThetaRipplepowerZscore.mat','Spec_power')

Spectrum1=[];Spectrum2=[];

    for Days=1:2
        for jj=1:3
            for kk=1:3
                Spectrum1{Days,jj}=[];
                Spectrum2{Days,jj}=[];
                Spectrum3{Days,jj}=[];
                Spectrum4{Days,jj}=[];
                Spectrum5{Days,jj}=[];
                Spectrum6{Days,jj}=[];
            end
        end
    end


Sess={[1:5]; [6:10, 11:15, 16:20 ]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};

for Days=1:2
for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),(Spec_power(ii,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
         Err=Spec_power{ii,jj}{1,6};
            x2=[];x4=[];
        for kk=1:size(Err,1)
            
             if Err(kk) == 0
             x1=(Spec_power{ii,jj}{1,1}{1}(kk,:));
             x2=[x2; x1]; 
             else
             x1=(Spec_power{ii,jj}{1,1}{1}(kk,:));
             x4=[x4; x1];
             end
                
        end
        Spectrum1{Days,ll}=[Spectrum1{Days,ll}; (x2)]; 
        Spectrum2{Days,ll}=[Spectrum2{Days,ll}; (x4)]; 
        

        
        x2=[];x4=[];
        for kk=1:size(Spec_power{ii,jj}{1,2}{1},1)
            
             
             x1=(Spec_power{ii,jj}{1,2}{1}(kk,:));
             x2=[x2; x1]; 
            
                
        end
         for kk=1:size(Spec_power{ii,jj}{1,3}{1},1)            
             
             x1=(Spec_power{ii,jj}{1,3}{1}(kk,:));
             x4=[x4; x1]; 
            
                
        end
        
        
        Spectrum3{Days,ll}=[Spectrum3{Days,ll}; (x2)]; 
        Spectrum4{Days,ll}=[Spectrum4{Days,ll}; (x4)];
        
        
             x2=[];x4=[];
        for kk=1:size(Spec_power{ii,jj}{1,4}{1},1)
            
             
             x1=(Spec_power{ii,jj}{1,4}{1}(kk,:));
             x2=[x2; x1]; 
            
                
        end
         for kk=1:size(Spec_power{ii,jj}{1,5}{1},1)            
             
             x1=(Spec_power{ii,jj}{1,5}{1}(kk,:));
             x4=[x4; x1]; 
            
                
        end
        
        
        Spectrum5{Days,ll}=[Spectrum5{Days,ll}; (x2)]; 
        Spectrum6{Days,ll}=[Spectrum6{Days,ll}; (x4)];
       
        
     end
end

end

Spectrum_Pw1=[];
Spectrum_Pw1{1}=[movmean([Spectrum1{2,3}; Spectrum2{2,3}],[100 100],2)];
% Spectrum_Pw1{2}=[movmean(Spectrum2{2,3},[100 100],2)];
Spectrum_Pw1{2}=[movmean(Spectrum3{2,3},[100 100],2)];
% Spectrum_Pw1{3}=[movmean(Spectrum4{2,3},[100 100],2)];
% Spectrum_Pw1{5}=[movmean(Spectrum5{2,3},[100 100],2)];
Spectrum_Pw1{3}=[movmean(Spectrum6{2,3},[100 100],2)];


col1={rgb('Red'),rgb('Blue'),rgb('Black'),rgb('Orange'),rgb('Green'),rgb('Orange')};
fig1=figure(1);
ax1(1)=axes('Position',[0.1 0.55 0.15 0.2]); 
data1=[];labels1=[];
for a=1:size(Spectrum_Pw1,2)
    mm=movmean(nanmean((Spectrum_Pw1{a})),[50 50],2);
    ms=movmean(nansem((Spectrum_Pw1{a})),[50 50],2);
    mm1=(Spectrum_Pw1{a});
    mm1=mm1;
    mm2{a}=mm1(:);
    data1=[data1; nanmean(mm1)'];
    labels1=[labels1; ones(size(mm1,2),1).*a];
    H1(a)=plot(1:size(mm,2),mm ,'linewidth',1,'marker','none','color',col1{a}); hold on;
    pp=patch([[1:length(mm)] fliplr([1:length(mm)])], ...
       [ms+mm fliplr(mm-ms)], col1{a})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl=get(gca,'ylim');
    xlabel('Peri-event time (s)','fontsize',6);
	ylabel('Theta Power (z-score)','fontsize',6);
    set(gca,'fontsize',6);
    set(gca,'xtick',[0:200:1200],'box','off','xticklabel',[-1:1:5],'fontsize',5);
    
    
end
plot([201 201],[min(yl) max(yl)],'--k');
ax1(2)=axes('Position',[0.1 0.15 0.15 0.25]); 
H=boxplot(data1,labels1,'boxstyle','outline','whisker',2, 'colors',[col1{1};col1{2};col1{3}] );

Groups=[1,2;1,3;2,3];
for ss1=1:size(Groups,1)
    [p h]=ranksum(mm2{Groups(ss1,1)},mm2{Groups(ss1,2)});
    if p<0.05
        h1=sigstar(Groups(ss1,:),p,0)
    end
end
set(gca,'fontsize',5,'box','off','xtick',[1,2,3],'xticklabel',{'Cue','Goal','Movement'},'xticklabelrotation',90);
xlabel('Events','fontsize',5);
ylabel('Theta power (z-score)','fontsize',5);


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter') 
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
set(gcf,'paperposition',[0.1,0.1,3,3])   
print(gcf,'-djpeg','-r300','ThetaPowerComparedZscore'); 
print(gcf,'-painters','-depsc','-r300','ThetaPowerComparedZscore'); 
close(fig1)



%% Fig 3m and n

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('StartArmDecodingTemporal50Shuff')
glm1{1}=MeanProbs1;
glm1{2}=MeanProbs2;
glm2{1}=MeanProbsShuff1;
glm2{2}=MeanProbsShuff2;
Err1{1}=Errors1; Err1{2}=Errors2;
CPDMat1=[];CPDMat2=[];
for Reg=1:2
    for Days=1:2
        for jj=1:3
            for kk=1:3
                CPDMat1{Reg,Days,jj,kk}=[];
                CPDMat2{Reg,Days,jj,kk}=[];
            end
        end
    end
end

Sess={[1:5,21:25]; [6:10, 11:15 , 16:20]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};
for Reg=1:2
for Days=1:2

for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),squeeze(glm1{Reg}(ii,1,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
        for ss=1:3
            x2=[];x4=[];
        for kk=1:size(squeeze(glm1{Reg}{ii,jj,ss}),2)
            
             if Err1{Reg}{ii,jj}(kk) > 0; continue; end 
             if isempty(squeeze(glm1{Reg}{ii,jj,1}{kk}));
             x1=ones(1,121).*NaN;
             else
             x1=(squeeze(glm1{Reg}{ii,jj,ss}{kk}(1:121,:))');
             x3=(squeeze(glm2{Reg}{ii,jj,ss}{kk}(1:121)));
%              x1(~all(x1'),:)=NaN;
             x1=nansum(x1);
             end
         x2=[x2; x1]; 
         x4=[x4; x3];
        end
        if size(x2,1) > 1; 
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; nanmean(x2)]; 
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; nanmean(x4)];
        else
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; (x2)];
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; (x4)];
        end
        
        end
        
     end
end

end
end


col1={rgb('Black'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('Magenta'),rgb('Grey')};
labels={'1','4','2','3','5'};
titlesReg={'HPC','PFC'};


Means1=[];Means2=[];
for ss=1:3
for nn=1:3

for Reg=1
 
for Days=1:2    
    x11=nanmean(movmean(squeeze(CPDMat1{Reg,Days,ss,nn}),[1 6],2)'); 
    x22=nanmean(movmean(squeeze(CPDMat2{Reg,Days,ss,nn}),[1 6],2)');
    Means1{ss,Days,nn}=((x11));  
    Means2{ss,Days,nn}=((x22)); 
end
end
end
end
fig1=figure(1)
labels1={'Ba-1','Ctrl','Ba-2'};
labels2={'Ba1','Asn','Asd'};
labels3={'Ba1','Assn','Assd'};
ss=-1;x1=[];x2=[];P=[];
for i=1:2
    if i==1 | i==2; sess=[1:3]; else; sess=[1:3]; end
    left=0.1+(i-1).*0.25; bottom=0.75; width=0.2; height=0.2;
    ax1(i)=axes('position',[left bottom width height])
 ss=ss+1;rr=0;   
 for j=sess
    ss=ss+1; 
    x1{j}=squeeze(Means1{j,i,1});  
    x2{j}=squeeze(Means2{j,i,1});
    [p1 h1 stats]=signrank(x1{j},x2{j});
    P{i,1}(j,2)=p1;P{i,1}(j,1)=j;P{i,1}(j,3)=length(x1{j});
 
   hold on;loc1=[];
    for nn=1:size(x1{j},2)
     loc1(nn)=+j;   
      plot(loc1(nn), x1{j}(nn),'marker','o','markersize',1,'color',col1{j})
    end
   loc2=[];
    for nn=1:size(x1{j},2)
      loc2(nn)= j+0.5; 
      plot(loc2(nn), x2{j}(nn),'marker','o','markersize',1,'color',col1{6})
    end
    
     for nn=1:size(x1{j},2)
      plot([loc1(nn) loc2(nn)], [x1{j}(nn) x2{j}(nn)],'linestyle','--','marker','none','color',rgb('Gray'),'linewidth',0.5);
    end

    HH(i,j)=errorbar(j-0.1,nanmean(x1{j}),nansem(x1{j}),'color',col1{j},'marker','o','capsize',3,'markersize',2); hold on;
    HH1(i,j)=errorbar(j+0.5+0.1,nanmean(x2{j}),nansem(x2{j}),'color',col1{6},'marker','o','capsize',3,'markersize',2); hold on;
    
    
    
    if p1 < 0.05  
      sigstar([j j+0.5],p1,0);
   end
    P{i,1}(j,4)=nanmean(x1{j}); P{i,1}(j,5)=nansem(x1{j});
    P{i,1}(j,6)=nanmean(x2{j}); P{i,1}(j,7)=nansem(x2{j});
end
Groups=[1,2;1,3;2,3];
p=[];h=[];
for ii=1:size(Groups,1)
    [p(ii) h(ii), stats]=ranksum(x1{Groups(ii,1)},x1{Groups(ii,2)});
    if p(ii) < 0.05
    gg=sigstar(Groups(ii,:),p(ii),0);
    end
    P{i,2}{ii,2}=p(ii);P{i,2}{ii,1}=Groups(ii,:);P{i,2}{ii,3}=stats.ranksum;
end
if i==1
 set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels1,...
    'fontsize',5,'xticklabelrotation',45)  
elseif i==2
   set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels2,...
    'fontsize',5,'xticklabelrotation',45)    
else
     set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels3,...
    'fontsize',5,'xticklabelrotation',45)   
end
set(gca,'fontsize',5);
if i==1
 ylabel('Goal prediction','fontsize',5)   
end
ylim([0.09 0.27])
xlim([0.5 3.75])
end
linkaxes(ax1([1:2]),'xy');



cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('StartArmDecodingTemporal50Shuff')
glm1{1}=MeanProbs1;
glm1{2}=MeanProbs2;
glm2{1}=MeanProbsShuff1;
glm2{2}=MeanProbsShuff2;
Err1{1}=Errors1; Err1{2}=Errors2;
CPDMat1=[];CPDMat2=[];
for Reg=1:2
    for Days=1:2
        for jj=1:3
            for kk=1:3
                CPDMat1{Reg,Days,jj,kk}=[];
                CPDMat2{Reg,Days,jj,kk}=[];
            end
        end
    end
end

Sess={[1:5,21:25]; [6:10, 11:15 , 16:20]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};
for Reg=1:2
for Days=1:2

for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),squeeze(glm1{Reg}(ii,1,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
        for ss=1:3
            x2=[];x4=[];
        for kk=1:size(squeeze(glm1{Reg}{ii,jj,ss}),2)
            
             if Err1{Reg}{ii,jj}(kk) < 1; continue; end 
             if isempty(squeeze(glm1{Reg}{ii,jj,1}{kk}));
             x1=ones(1,121).*NaN;
             else
             x1=(squeeze(glm1{Reg}{ii,jj,ss}{kk}(1:121,:))');
             x3=(squeeze(glm2{Reg}{ii,jj,ss}{kk}(1:121)));
%              x1(~all(x1'),:)=NaN;
             x1=nansum(x1);
             end
         x2=[x2; x1]; 
         x4=[x4; x3];
        end
        if size(x2,1) > 1; 
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; nanmean(x2)]; 
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; nanmean(x4)];
        else
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; (x2)];
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; (x4)];
        end
        
        end
        
     end
end

end
end


col1={rgb('Black'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('Magenta'),rgb('Grey')};
labels={'1','4','2','3','5'};
titlesReg={'HPC','PFC'};


Means1=[];Means2=[];
for ss=1:3
for nn=1:3

for Reg=1
 
for Days=1:2    
    x11=nanmean(movmean(squeeze(CPDMat1{Reg,Days,ss,nn}),[1 6],2)'); 
    x22=nanmean(movmean(squeeze(CPDMat2{Reg,Days,ss,nn}),[1 6],2)');
    Means1{ss,Days,nn}=((x11));  
    Means2{ss,Days,nn}=((x22)); 
end
end
end
end
labels1={'Ba-1','Ctrl','Ba-2'};
labels2={'Ba1','Asn','Asd'};
ss=-1;x1=[];x2=[];P=[];
for i=1:2
    if i==1 | i==2; sess=[1:3]; else; sess=[1:3]; end
    left=0.1+(i-1).*0.25; bottom=0.5; width=0.2; height=0.2;
    ax1(i)=axes('position',[left bottom width height])
 ss=ss+1;rr=0;   
 for j=sess
    ss=ss+1; 
    x1{j}=squeeze(Means1{j,i,1});  
    x2{j}=squeeze(Means2{j,i,1});
[p1 h1 stats]=signrank(x1{j},x2{j});
    P{i,1}(j,2)=p1;P{i,1}(j,1)=j;P{i,1}(j,3)=length(x1{j});
   hold on;loc1=[];
    for nn=1:size(x1{j},2)
     loc1(nn)=+j;   
      plot(loc1(nn), x1{j}(nn),'marker','o','markersize',1,'color',col1{j})
    end
   loc2=[];
    for nn=1:size(x1{j},2)
      loc2(nn)= j+0.5; 
      plot(loc2(nn), x2{j}(nn),'marker','o','markersize',1,'color',col1{6})
    end
    
     for nn=1:size(x1{j},2)
      plot([loc1(nn) loc2(nn)], [x1{j}(nn) x2{j}(nn)],'linestyle','--','marker','none','color',rgb('Gray'),'linewidth',0.5);
    end
 
  HH(i,j)=errorbar(j-0.1,nanmean(x1{j}),nansem(x1{j}),'color',col1{j},'marker','o','capsize',3,'markersize',2); hold on;
    HH1(i,j)=errorbar(j+0.5+0.1,nanmean(x2{j}),nansem(x2{j}),'color',col1{6},'marker','o','capsize',3,'markersize',2); hold on;
    
    
    
    if p1 < 0.05  
      sigstar([j j+0.5],p1,0);
   end
    P{i,1}(j,4)=nanmean(x1{j}); P{i,1}(j,5)=nansem(x1{j});
    P{i,1}(j,6)=nanmean(x2{j}); P{i,1}(j,7)=nansem(x2{j});
end
Groups=[1,2;1,3;2,3];
p=[];h=[];
for ii=1:size(Groups,1)
    [p(ii) h(ii) stats]=ranksum(x1{Groups(ii,1)},x1{Groups(ii,2)},'method','approximate');
    if p(ii) < 0.05
    gg=sigstar(Groups(ii,:),p(ii),0);
    end
    P{i,2}{ii,2}=p(ii);P{i,2}{ii,1}=Groups(ii,:);P{i,2}{ii,3}=stats.zval;
end
if i==1
 set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels1,...
    'fontsize',5,'xticklabelrotation',45)  
elseif i==2
   set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels2,...
    'fontsize',5,'xticklabelrotation',45)    
else
     set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels3,...
    'fontsize',5,'xticklabelrotation',45)   
end
set(gca,'fontsize',5);
if i==1
 ylabel('Goal prediction','fontsize',5)   
end
ylim([0.09 0.27])
xlim([0.5 3.75])
end
linkaxes(ax1([1:2]),'xy');


labels1={'Ba-1','Ctrl','Ba-2'};
labels2={'Ba1','Asn','Asd'};
labels3={'Ba1','Assn','Assd'};
ss=-1;x1=[];x2=[];P=[];
for i=1:2
    if i==1 | i==2; sess=[1:3]; else; sess=[1:3]; end
    left=0.1+(i-1).*0.25; bottom=0.25; width=0.2; height=0.2;
    ax1(i)=axes('position',[left bottom width height])
 ss=ss+1;rr=0;   
 for j=sess
    ss=ss+1; 
    x1{j}=squeeze(Means1{j,i,2});  
    x2{j}=squeeze(Means2{j,i,2});
   [p1 h1 stats]=signrank(x1{j},x2{j});
    P{i,1}(j,2)=p1;P{i,1}(j,1)=j;P{i,1}(j,3)=length(x1{j});
 
   hold on;loc1=[];
    for nn=1:size(x1{j},2)
     loc1(nn)=+j;   
      plot(loc1(nn), x1{j}(nn),'marker','o','markersize',1,'color',col1{j})
    end
   loc2=[];
    for nn=1:size(x1{j},2)
      loc2(nn)= j+0.5; 
      plot(loc2(nn), x2{j}(nn),'marker','o','markersize',1,'color',col1{6})
    end
    
     for nn=1:size(x1{j},2)
      plot([loc1(nn) loc2(nn)], [x1{j}(nn) x2{j}(nn)],'linestyle','--','marker','none','color',rgb('Gray'),'linewidth',0.5);
    end
    
 HH(i,j)=errorbar(j-0.1,nanmean(x1{j}),nansem(x1{j}),'color',col1{j},'marker','o','capsize',3,'markersize',2); hold on;
    HH1(i,j)=errorbar(j+0.5+0.1,nanmean(x2{j}),nansem(x2{j}),'color',col1{6},'marker','o','capsize',3,'markersize',2); hold on;
    
    
    
    if p1 < 0.05  
      sigstar([j j+0.5],p1,0);
   end
    P{i,1}(j,4)=nanmean(x1{j}); P{i,1}(j,5)=nansem(x1{j});
    P{i,1}(j,6)=nanmean(x2{j}); P{i,1}(j,7)=nansem(x2{j});
end
Groups=[1,2;1,3;2,3];
p=[];h=[];
for ii=1:size(Groups,1)
    [p(ii) h(ii) stats]=ranksum(x1{Groups(ii,1)},x1{Groups(ii,2)},'method','approximate');
    if p(ii) < 0.05
    gg=sigstar(Groups(ii,:),p(ii),0);
    end
    P{i,2}{ii,2}=p(ii);P{i,2}{ii,1}=Groups(ii,:);P{i,2}{ii,3}=stats.zval;
end
if i==1
 set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels1,...
    'fontsize',5,'xticklabelrotation',45)  
elseif i==2
   set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels2,...
    'fontsize',5,'xticklabelrotation',45)    
else
     set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels3,...
    'fontsize',5,'xticklabelrotation',45)   
end
set(gca,'fontsize',5);
if i==1
 ylabel('Goal prediction','fontsize',5)   
end
ylim([0.09 0.27])
xlim([0.5 3.75])
end
linkaxes(ax1([1:2]),'xy');


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(fig1,'paperposition',[0.2,0.2,3,3])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','StartArmDecodingTemporalTargetVsShuffSessCompressed');
print(fig1,'-painters','-depsc','-r300','StartArmDecodingTemporalTargetVsShuffSessCompressed'); 
close all;

%% Fig 3o
% HPC
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('StartArmDecodingSpatial50Shuff')
glm1{1}=MeanProbs1;

glm2{1}=MeanProbsShuff1;

Err1{1}=Errors1; 
CPDMat1=[];CPDMat2=[];
for Reg=1
    for Days=1:2
        for jj=1:3
            for kk=1:3
                CPDMat1{Reg,Days,jj,kk}=[];
                CPDMat2{Reg,Days,jj,kk}=[];
            end
        end
    end
end

Sess={[1:5,21:25]; [6:10, 11:15 , 16:20]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};
for Reg=1
for Days=1:2

for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),squeeze(glm1{Reg}(ii,1,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
        for ss=1:3
            x2=[];x4=[];
        for kk=1:size(squeeze(glm1{Reg}{ii,jj,ss}),2)
            
             if Err1{Reg}{ii,jj}(kk) > 0; continue; end 
             if isempty(squeeze(glm1{Reg}{ii,jj,1}{kk}));
             x1=ones(1,121).*NaN;
             else
             x1=(squeeze(glm1{Reg}{ii,jj,ss}{kk}(1:121,:))');
             x3=(squeeze(glm2{Reg}{ii,jj,ss}{kk}(1:121)));
%              x1(~all(x1'),:)=NaN;
             x1=nansum(x1);
             end
         x2=[x2; x1]; 
         x4=[x4; x3];
        end
        if size(x2,1) > 1; 
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; nanmean(x2)]; 
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; nanmean(x4)];
        else
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; (x2)];
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; (x4)];
        end
        
        end
        
     end
end

end
end


col1={rgb('Black'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('Magenta'),rgb('Grey')};
labels={'1','4','2','3','5'};
titlesReg={'HPC','PFC'};


Means1=[];Means2=[];
for ss=1:3
for nn=1:3

for Reg=1
 
for Days=1:2    
    x11=nanmean(movmean(squeeze(CPDMat1{Reg,Days,ss,nn}),[1 6],2)'); 
    x22=nanmean(movmean(squeeze(CPDMat2{Reg,Days,ss,nn}),[1 6],2)');
    Means1{ss,Days,nn}=((x11));  
    Means2{ss,Days,nn}=((x22)); 
end
end
end
end
fig1=figure(1)
labels1={'Ba-1','Ctrl','Ba-2'};
labels2={'Ba1','Asn','Asd'};
labels3={'Ba1','Assn','Assd'};
ss=-1;x1=[];x2=[];P=[];
for i=1:2
    if i==1 | i==2; sess=[1:3]; else; sess=[1:3]; end
    left=0.1+(i-1).*0.25; bottom=0.75; width=0.2; height=0.2;
    ax1(i)=axes('position',[left bottom width height])
 ss=ss+1;rr=0;   
 for j=sess
    ss=ss+1; 
    x1{j}=squeeze(Means1{j,i,1});  
    x2{j}=squeeze(Means2{j,i,1});
    [p1 h1 stats]=signrank(x1{j},x2{j});
    P{i,1}(j,2)=p1;P{i,1}(j,1)=j;P{i,1}(j,3)=length(x1{j});
 
   hold on;loc1=[];
    for nn=1:size(x1{j},2)
     loc1(nn)=+j;   
      plot(loc1(nn), x1{j}(nn),'marker','o','markersize',1,'color',col1{j})
    end
   loc2=[];
    for nn=1:size(x1{j},2)
      loc2(nn)= j+0.5; 
      plot(loc2(nn), x2{j}(nn),'marker','o','markersize',1,'color',col1{6})
    end
    
     for nn=1:size(x1{j},2)
      plot([loc1(nn) loc2(nn)], [x1{j}(nn) x2{j}(nn)],'linestyle','--','marker','none','color',rgb('Gray'),'linewidth',0.5);
    end
    
 
    HH(i,j)=errorbar(j-0.1,nanmean(x1{j}),nansem(x1{j}),'color',col1{j},'marker','o','capsize',3,'markersize',2); hold on;
    HH1(i,j)=errorbar(j+0.5+0.1,nanmean(x2{j}),nansem(x2{j}),'color',col1{6},'marker','o','capsize',3,'markersize',2); hold on;
    
    
    
    if p1 < 0.05  
      sigstar([j j+0.5],p1,0);
   end
    P{i,1}(j,4)=nanmean(x1{j}); P{i,1}(j,5)=nansem(x1{j});
    P{i,1}(j,6)=nanmean(x2{j}); P{i,1}(j,7)=nansem(x2{j});
end
Groups=[1,2;1,3;2,3];
p=[];h=[];
for ii=1:size(Groups,1)
    [p(ii) h(ii), stats]=ranksum(x1{Groups(ii,1)},x1{Groups(ii,2)});
    if p(ii) < 0.05
    gg=sigstar(Groups(ii,:),p(ii),0);
    end
    P{i,2}{ii,2}=p(ii);P{i,2}{ii,1}=Groups(ii,:);P{i,2}{ii,3}=stats.ranksum;
end
if i==1
 set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels1,...
    'fontsize',5,'xticklabelrotation',45)  
elseif i==2
   set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels2,...
    'fontsize',5,'xticklabelrotation',45)    
else
     set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels3,...
    'fontsize',5,'xticklabelrotation',45)   
end
set(gca,'fontsize',5);
if i==1
 ylabel('Goal prediction','fontsize',5)   
end
ylim([0 0.14])
xlim([0.5 3.75])
end
linkaxes(ax1([1:2]),'xy');

cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('StartArmDecodingSpatial50Shuff')
glm1{1}=MeanProbs1;

glm2{1}=MeanProbsShuff1;

Err1{1}=Errors1; 
CPDMat1=[];CPDMat2=[];
for Reg=1
    for Days=1:2
        for jj=1:3
            for kk=1:3
                CPDMat1{Reg,Days,jj,kk}=[];
                CPDMat2{Reg,Days,jj,kk}=[];
            end
        end
    end
end

Sess={[1:5,21:25]; [6:10, 11:15 , 16:20]};
sess1={[1:3],[1:4],[1:4],[1:4],[1:4]};
for Reg=1
for Days=1:2

for ii=Sess{Days}
    
    tt=cell2mat(cellfun(@(x)(~isempty(x)),squeeze(glm1{Reg}(ii,1,:)),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    ll=0;
     for jj=sess1{Days}
         if jj<4
         ll=ll+1;
         end
        for ss=1:3
            x2=[];x4=[];
        for kk=1:size(squeeze(glm1{Reg}{ii,jj,ss}),2)
            
             if Err1{Reg}{ii,jj}(kk) < 1; continue; end 
             if isempty(squeeze(glm1{Reg}{ii,jj,1}{kk}));
             x1=ones(1,121).*NaN;
             else
             x1=(squeeze(glm1{Reg}{ii,jj,ss}{kk}(1:121,:))');
             x3=(squeeze(glm2{Reg}{ii,jj,ss}{kk}(1:121)));
%              x1(~all(x1'),:)=NaN;
             x1=nansum(x1);
             end
         x2=[x2; x1]; 
         x4=[x4; x3];
        end
        if size(x2,1) > 1; 
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; nanmean(x2)]; 
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; nanmean(x4)];
        else
        CPDMat1{Reg,Days,ll,ss}=[CPDMat1{Reg,Days,ll,ss}; (x2)];
        CPDMat2{Reg,Days,ll,ss}=[CPDMat2{Reg,Days,ll,ss}; (x4)];
        end
        
        end
        
     end
end

end
end


col1={rgb('Black'),rgb('Blue'),rgb('Red'),rgb('Red'),rgb('Magenta'),rgb('Grey')};
labels={'1','4','2','3','5'};
titlesReg={'HPC','PFC'};


Means1=[];Means2=[];
for ss=1:3
for nn=1:3

for Reg=1
 
for Days=1:2    
    x11=nanmean(movmean(squeeze(CPDMat1{Reg,Days,ss,nn}),[1 6],2)'); 
    x22=nanmean(movmean(squeeze(CPDMat2{Reg,Days,ss,nn}),[1 6],2)');
    Means1{ss,Days,nn}=((x11));  
    Means2{ss,Days,nn}=((x22)); 
end
end
end
end
labels1={'Ba-1','Ctrl','Ba-2'};
labels2={'Ba1','Asn','Asd'};
ss=-1;x1=[];x2=[];P=[];
for i=1:2
    if i==1 | i==2; sess=[1:3]; else; sess=[1:3]; end
    left=0.1+(i-1).*0.25; bottom=0.5; width=0.2; height=0.2;
    ax1(i)=axes('position',[left bottom width height])
 ss=ss+1;rr=0;   
 for j=sess
    ss=ss+1; 
    x1{j}=squeeze(Means1{j,i,1});  
    x2{j}=squeeze(Means2{j,i,1});
[p1 h1 stats]=signrank(x1{j},x2{j});
    P{i,1}(j,2)=p1;P{i,1}(j,1)=j;P{i,1}(j,3)=length(x1{j});
   hold on;loc1=[];
    for nn=1:size(x1{j},2)
     loc1(nn)=+j;   
      plot(loc1(nn), x1{j}(nn),'marker','o','markersize',1,'color',col1{j})
    end
   loc2=[];
    for nn=1:size(x1{j},2)
      loc2(nn)= j+0.5; 
      plot(loc2(nn), x2{j}(nn),'marker','o','markersize',1,'color',col1{6})
    end
    
     for nn=1:size(x1{j},2)
      plot([loc1(nn) loc2(nn)], [x1{j}(nn) x2{j}(nn)],'linestyle','--','marker','none','color',rgb('Gray'),'linewidth',0.5);
    end

  HH(i,j)=errorbar(j-0.1,nanmean(x1{j}),nansem(x1{j}),'color',col1{j},'marker','o','capsize',3,'markersize',2); hold on;
    HH1(i,j)=errorbar(j+0.5+0.1,nanmean(x2{j}),nansem(x2{j}),'color',col1{6},'marker','o','capsize',3,'markersize',2); hold on;
    
    
    
    if p1 < 0.05  
      sigstar([j j+0.5],p1,0);
   end
    P{i,1}(j,4)=nanmean(x1{j}); P{i,1}(j,5)=nansem(x1{j});
    P{i,1}(j,6)=nanmean(x2{j}); P{i,1}(j,7)=nansem(x2{j});
end
Groups=[1,2;1,3;2,3];
p=[];h=[];
for ii=1:size(Groups,1)
    [p(ii) h(ii) stats]=ranksum(x1{Groups(ii,1)},x1{Groups(ii,2)},'method','approximate');
    if p(ii) < 0.05
    gg=sigstar(Groups(ii,:),p(ii),0);
    end
    P{i,2}{ii,2}=p(ii);P{i,2}{ii,1}=Groups(ii,:);P{i,2}{ii,3}=stats.zval;
end
if i==1
 set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels1,...
    'fontsize',5,'xticklabelrotation',45)  
elseif i==2
   set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels2,...
    'fontsize',5,'xticklabelrotation',45)    
else
     set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels3,...
    'fontsize',5,'xticklabelrotation',45)   
end
set(gca,'fontsize',5);
if i==1
 ylabel('Goal prediction','fontsize',5)   
end
ylim([0.04 0.14])
xlim([0.5 3.75])
end
linkaxes(ax1([1:2]),'xy');

labels1={'Ba-1','Ctrl','Ba-2'};
labels2={'Ba1','Asn','Asd'};
labels3={'Ba1','Assn','Assd'};
ss=-1;x1=[];x2=[];P=[];
for i=1:2
    if i==1 | i==2; sess=[1:3]; else; sess=[1:3]; end
    left=0.1+(i-1).*0.25; bottom=0.25; width=0.2; height=0.2;
    ax1(i)=axes('position',[left bottom width height])
 ss=ss+1;rr=0;   
 for j=sess
    ss=ss+1; 
    x1{j}=squeeze(Means1{j,i,2});  
    x2{j}=squeeze(Means2{j,i,2});
    [p1 h1 stats]=signrank(x1{j},x2{j});
    P{i,1}(j,2)=p1;P{i,1}(j,1)=j;P{i,1}(j,3)=length(x1{j});
 
   hold on;loc1=[];
    for nn=1:size(x1{j},2)
     loc1(nn)=+j;   
      plot(loc1(nn), x1{j}(nn),'marker','o','markersize',1,'color',col1{j})
    end
   loc2=[];
    for nn=1:size(x1{j},2)
      loc2(nn)= j+0.5; 
      plot(loc2(nn), x2{j}(nn),'marker','o','markersize',1,'color',col1{6})
    end
    
     for nn=1:size(x1{j},2)
      plot([loc1(nn) loc2(nn)], [x1{j}(nn) x2{j}(nn)],'linestyle','--','marker','none','color',rgb('Gray'),'linewidth',0.5);
    end
    
 
  HH(i,j)=errorbar(j-0.1,nanmean(x1{j}),nansem(x1{j}),'color',col1{j},'marker','o','capsize',3,'markersize',2); hold on;
    HH1(i,j)=errorbar(j+0.5+0.1,nanmean(x2{j}),nansem(x2{j}),'color',col1{6},'marker','o','capsize',3,'markersize',2); hold on;
    
    
    
    if p1 < 0.05  
      sigstar([j j+0.5],p1,0);
   end
    P{i,1}(j,4)=nanmean(x1{j}); P{i,1}(j,5)=nansem(x1{j});
    P{i,1}(j,6)=nanmean(x2{j}); P{i,1}(j,7)=nansem(x2{j});
end
Groups=[1,2;1,3;2,3];
p=[];h=[];
for ii=1:size(Groups,1)
    [p(ii) h(ii) stats]=ranksum(x1{Groups(ii,1)},x1{Groups(ii,2)},'method','approximate');
    if p(ii) < 0.05
    gg=sigstar(Groups(ii,:),p(ii),0);
    end
    P{i,2}{ii,2}=p(ii);P{i,2}{ii,1}=Groups(ii,:);P{i,2}{ii,3}=stats.zval;
end
if i==1
 set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels1,...
    'fontsize',5,'xticklabelrotation',45)  
elseif i==2
   set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels2,...
    'fontsize',5,'xticklabelrotation',45)    
else
     set(gca,'box','off','tickdir','out',...
    'xtick',[1:3],'xticklabel',labels3,...
    'fontsize',5,'xticklabelrotation',45)   
end
set(gca,'fontsize',5);
if i==1
 ylabel('Goal prediction','fontsize',5)   
end
ylim([0.04 0.14])
xlim([0.5 3.75])
end
linkaxes(ax1([1:2]),'xy');

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(fig1,'paperposition',[0.2,0.2,3,3])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
print(fig1,'-djpeg','-r300','StartArmDecodingSpatialTargetVsShuffSessHPCCompressed');
print(fig1,'-painters','-depsc','-r300','StartArmDecodingSpatialTargetVsShuffSessHPCCompressed'); 
close all;


