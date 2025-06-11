
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
ll=1;glmstats=[];GoalProbsStacked=[];Goalbins=[];Goaldprobp=[];

for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
dn1=pwd;
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
load mazel
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

   
for ll=1:length(mazel)
tic;  
LapInfo = LapInfo1;
plf1=[];plfConcat=[];
arms=[1,2,4,5,7,8];
for nn=1:length(FiringRates1) 
plf11=FiringRates1{nn};
plf12=[];
ind=[16:40];
for jj=1:6
    plf1=plf11(LapInfo(:,5)==ll & LapInfo(:,6) ==arms(jj)  & (LapInfo(:,9)==1),ind);
    plf12=[plf12, nanmean(plf1)];
       
end
plfConcat=[plfConcat; (plf12)];
end

    LapInfo = LapInfo1;
    if ii==12 & ll==4
    LapInfo=[LapInfo; ones(1,size(LapInfo,2)).*NaN];
    LapInfo(end,9)=0; LapInfo(end,5)=ll;
    LapInfo1=LapInfo;
    end

    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==0));
    ind1=find(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    
        
   
    target_info= LapInfo(:,10);
    current_info=LapInfo(:,6);
    Next_info= LapInfo1(ind1+1,6);
    Last_Info= [current_info(1); LapInfo1(ind1(2:end)-1,6)];
    MeanProb=[];     
 for rr=1:size(LapInfo,1)     
   t_stimulus=LapInfo(rr,7);
   t_bin=0.05;
   tbin=0.05;
   tvec2=[t_stimulus-5   t_stimulus+10];       
   binns=nanmin(tvec2):tbin:nanmax(tvec2);
   binedges=(binns(1:end-1)+binns(2:end))/2; 
cd(dn1);   
kkk=1:length(cluster);
kk=kkk(ind_HPC);   
trackcltime=[];
b11=cluster;
b11=b11(ind_HPC);
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
%    
[dprobp dprob binedges]=BayesDecode_FramesWithPlot(tvec2,plfConcat,cltimevec,tbin);
GoalProbsStacked{ii,ll}{rr}=dprob;
Goalbins{ii,ll}{rr}=binedges;
Goaldprobp{ii,ll}{rr}=dprobp;
 end
toc
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('GoalProbsSessionsHPCRT50','GoalProbsStacked','Goalbins','Goaldprobp','-v7.3') 


clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
ll=1;glmstats=[];GoalProbsStacked=[];Goalbins=[];Goaldprobp=[];

for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
dn1=pwd;
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
load mazel
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
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates10(ind_HPC,1);
sel=1:100;
if sum(ind_HPC) < 10; continue; end

   
for ll=1:length(mazel)
tic;  
LapInfo = LapInfo1;
plf1=[];plfConcat=[];
arms=[1,2,4,5,7,8];
for nn=1:length(FiringRates1) 
plf11=FiringRates1{nn};
plf12=[];
ind=[16:40];
for jj=1:6
    plf1=plf11(LapInfo(:,5)==ll & LapInfo(:,6) ==arms(jj)  & (LapInfo(:,9)==1),ind);
    plf12=[plf12, nanmean(plf1)];
       
end
plfConcat=[plfConcat; (plf12)];
end

    LapInfo = LapInfo1;
    if ii==12 & ll==4
    LapInfo=[LapInfo; ones(1,size(LapInfo,2)).*NaN];
    LapInfo(end,9)=0; LapInfo(end,5)=ll;
    LapInfo1=LapInfo;
    end
    
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==0));
    ind1=find(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)==0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    
        
   
    target_info= LapInfo(:,10);
    current_info=LapInfo(:,6);
    Next_info= LapInfo1(ind1+1,6);
    Last_Info= [current_info(1); LapInfo1(ind1(2:end)-1,6)];
    MeanProb=[];     
 for rr=1:size(LapInfo,1)     
   t_stimulus=LapInfo(rr,7);
   t_bin=0.05;
   tbin=0.05;
   tvec2=[t_stimulus-5   t_stimulus+10];       
   binns=nanmin(tvec2):tbin:nanmax(tvec2);
   binedges=(binns(1:end-1)+binns(2:end))/2; 
cd(dn1);   
kkk=1:length(cluster);
kk=kkk(ind_HPC);   
trackcltime=[];
b11=cluster;
b11=b11(ind_HPC);
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
%    
 [dprobp dprob binedges]=BayesDecode_FramesWithPlot(tvec2,plfConcat,cltimevec,tbin);
GoalProbsStacked{ii,ll}{rr}=dprob;
Goalbins{ii,ll}{rr}=binedges;
Goaldprobp{ii,ll}{rr}=dprobp;
 end
toc
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('GoalProbsSessionsPFCRT50','GoalProbsStacked','Goalbins','Goaldprobp','-v7.3') 

