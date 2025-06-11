clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load PlFieldTempNAll
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,5,5,5,5,5,4,4,4,4,4]);
for i=1:size(PlFieldTemp,2)
     plf{i,1}=PlFieldTemp(i).PlaceFields;
     animal(i,1)=PlFieldTemp(i).Animal;
     day(i,1)=PlFieldTemp(i).Day;
     session(i,1)=PlFieldTemp(i).session;
     track(i,1)=PlFieldTemp(i).Track;
     direction(i,1)=PlFieldTemp(i).Direction;
end
ll=1;glmstats=[];
TrackProbsStacked=[];Trackbins=[];Trackdprobp=[];Trackbinsp=[];
gaussFilter = gausswin(3,2);
kernel = gaussFilter / sum(gaussFilter); % Normalize.
% cd('/media/baburam/DiskStorage1/Database/TimeSequence/Awake');
% load('TrackProbsSessionsAllTrials50','TrackProbsStacked','Trackbins','Trackbinsp','Trackdprobp') 
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
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
% rates1=FiringRates10(ind_HPC,1);
sel=1:100;
if sum(ind_HPC) < 10; continue; end

   
for ll=1:length(mazel)
tic;  
LapInfo = LapInfo1;
arms=[1,2,4,5,7,8];
plf1=[];plfConcat=[];  
for kk=1:2
for jj=1:6
plf2=[];    
plf1=plf{animal==Animals(ii) & day==Days(ii) & session==ll & track==arms(jj) & direction==kk};
plf2(:,1:70)=(plf1(:,1:end-1)+plf1(:,2:end))/2;
plfConcat=[plfConcat, plf2];
end
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
   tvec2=[t_stimulus-1   t_stimulus+6];       
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
 [dprobp dprob binedges binedgesp]=BayesDecode_FramesWithOverlap(tvec2,plfConcat,cltimevec,tbin,0.7);
TrackProbsStacked{ii,ll}{rr}=dprob;
Trackbins{ii,ll}{rr}=binedges;
Trackbinsp{ii,ll}{rr}=binedgesp;
Trackdprobp{ii,ll}{rr}=dprobp;
end
toc
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('TrackProbsSessionsAllTrials50','TrackProbsStacked','Trackbins','Trackbinsp','Trackdprobp', '-v7.3') 

