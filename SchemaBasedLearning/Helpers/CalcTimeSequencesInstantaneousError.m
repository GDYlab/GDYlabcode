%% HPC instantaneous decoding

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)

MedErr=[];Predictions=[];TotPredicitions=[];TotPredicitionsp=[];
TrueTime=[];MedErrShuff=[];
for i=1:25
 TrueTime=[TrueTime; ones(12,1).*i];
end
for ii=[1:25]
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
   
  
for ll=1:length(maze)
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
    plf1=plf11(LapInfo(:,5)==ll & LapInfo(:,6) ==arms(jj) & LapInfo(:,9)==1,ind);
    plf12=[plf12, nanmean(plf1)];
    plf_forplot{jj}=[plf_forplot{jj}; nanmean(plf1)];
    
end
plfConcat=[plfConcat; (plf12)];
end

ErrorMaxProb=[];BinsDecoded=[];ErrorShuff=[];
Decoded=[];                 
 for nn=1:length(arms)
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & LapInfo(:,9)==1);
    ind1=find(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & LapInfo(:,9)==1);
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(LapInfo(:,6)==arms(nn) & idx,:);  
    
   bins=[1:25:151];          
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

 [dprobp dprob, binedges]=BayesDecode_FramesWithPlot(tvec2,plfConcat,cltimevec,tbin);
 binsofInterest=bins(nn):bins(nn+1)-1;

 binsProb = dprob(binsofInterest,:); 
 [MaxProbs MaxIDs]=max(binsProb);
 binsProb = dprobp(binsofInterest,:);   
 [MaxProbsp MaxIDsp]=max(binsProb);
 MaxIDsp(:,isnan(MaxProbs))=NaN;
 idx=MaxProbsp>0 & ~isnan(MaxProbs);
 TrueTime1=TrueTime(idx);
 MaxIDsp=MaxIDsp(idx);
 
 TotBinsDecoded{nn,rr}=dprob;
 TotBinsDecodedp{nn,rr}=dprobp;
 BinsDecoded{nn,rr}=binsProb;
  if ~isempty(BinsDecoded{nn,rr})
 ErrorMaxProb{nn}(rr)=nanmedian(abs(MaxIDsp-TrueTime1')./4);
 else
 ErrorMaxProb{nn}(rr)=NaN;
  end
 for ss1=1:500
 dprobs=[];dprobp=[];ntimebin=size(dprob,1);   
 for j=1:size(dprob,2)
  dprobs(:,j)=circshift( dprob(:,j), round(rand(1)*ntimebin),1);
  dprobsp(:,j)=circshift( dprob(:,j), round(rand(1)*ntimebin),1);
 end    
 binsofInterest=bins(nn):bins(nn+1)-1;

 binsProb = dprobs(binsofInterest,:); 
 [MaxProbs MaxIDs]=max(binsProb);
 binsProb = dprobsp(binsofInterest,:);   
 [MaxProbsp MaxIDsp]=max(binsProb);
 MaxIDsp(:,isnan(MaxProbs))=NaN;
 idx=MaxProbsp>0 & ~isnan(MaxProbs);
 TrueTime1=TrueTime(idx);
 MaxIDsp=MaxIDsp(idx); 
 if ~isempty(BinsDecoded{nn,rr})
 ErrorShuff{nn}(rr,ss1)=nanmedian(abs(MaxIDsp-TrueTime1')./4);
 else
 ErrorShuff{nn}(rr,ss1)=NaN;
  end
 end
 end
 end
 MedErr{ii,ll}=ErrorMaxProb;
 Predictions{ii,ll}=BinsDecoded;
 TotPredicitions{ii,ll}=TotBinsDecoded;
 TotPredicitionsp{ii,ll}=TotBinsDecodedp;
 MedErrShuff{ii,ll}=ErrorShuff;
end
end
toc
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
save('InstantaneousTimeErrHPC','MedErr','MedErrShuff','Predictions','TotPredicitions','TotPredicitionsp','-v7.3')


%% PFC instantaneous decoding

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
MedErr=[];Predictions=[];TotPredicitions=[];TotPredicitionsp=[];
TrueTime=[];MedErrShuff=[];
for i=1:25
 TrueTime=[TrueTime; ones(12,1).*i];

end
for ii=[1:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load maze
dn1=pwd;

cluster=importdata('clusters');
elecreg=importdata('elecreg');
% load SelectivityIndicesSignificance.mat
% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
elecreg=importdata('elecreg');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
% if exist('FrameDecode.mat')
% delete('FrameDecode.mat');
% end
tic;

FiringRates1=FiringRates1(:,1);


GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
dn=pwd;
fstruct = dir('*FramePFCShuff*.mat');
if ii==6;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
fstruct1=array2table(fstruct);
LapInfo1=LapInfo;
if ~isempty(fstruct)
   
  
for ll=1:length(maze)
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
    plf1=plf11(LapInfo(:,5)==ll & LapInfo(:,6) ==arms(jj) & LapInfo(:,9)==1,ind);
    plf12=[plf12, nanmean(plf1)];
    plf_forplot{jj}=[plf_forplot{jj}; nanmean(plf1)];
    
end
plfConcat=[plfConcat; (plf12)];
end

ErrorMaxProb=[];BinsDecoded=[];ErrorShuff=[];
Decoded=[];                 
 for nn=1:length(arms)
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & LapInfo(:,9)==1);
    ind1=find(LapInfo(:,5)==ll & ~isnan(LapInfo(:,10)) & LapInfo(:,9)==1);
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(LapInfo(:,6)==arms(nn) & idx,:);  
    
   bins=[1:25:151];          
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

 [dprobp dprob, binedges]=BayesDecode_FramesWithPlot(tvec2,plfConcat,cltimevec,tbin);
 binsofInterest=bins(nn):bins(nn+1)-1;

 binsProb = dprob(binsofInterest,:); 
 [MaxProbs MaxIDs]=max(binsProb);
 binsProb = dprobp(binsofInterest,:);   
 [MaxProbsp MaxIDsp]=max(binsProb);
 MaxIDsp(:,isnan(MaxProbs))=NaN;
 idx=MaxProbsp>0 & ~isnan(MaxProbs);
 TrueTime1=TrueTime(idx);
 MaxIDsp=MaxIDsp(idx);
 
 TotBinsDecoded{nn,rr}=dprob;
 TotBinsDecodedp{nn,rr}=dprobp;
 BinsDecoded{nn,rr}=binsProb;
  if ~isempty(BinsDecoded{nn,rr})
 ErrorMaxProb{nn}(rr)=nanmedian(abs(MaxIDsp-TrueTime1')./4);
 else
 ErrorMaxProb{nn}(rr)=NaN;
  end
 for ss1=1:500
 dprobs=[];dprobp=[];ntimebin=size(dprob,1);   
 for j=1:size(dprob,2)
  dprobs(:,j)=circshift( dprob(:,j), round(rand(1)*ntimebin),1);
  dprobsp(:,j)=circshift( dprob(:,j), round(rand(1)*ntimebin),1);
 end    
 binsofInterest=bins(nn):bins(nn+1)-1;

 binsProb = dprobs(binsofInterest,:); 
 [MaxProbs MaxIDs]=max(binsProb);
 binsProb = dprobsp(binsofInterest,:);   
 [MaxProbsp MaxIDsp]=max(binsProb);
 MaxIDsp(:,isnan(MaxProbs))=NaN;
 idx=MaxProbsp>0 & ~isnan(MaxProbs);
 TrueTime1=TrueTime(idx);
 MaxIDsp=MaxIDsp(idx); 
 if ~isempty(BinsDecoded{nn,rr})
 ErrorShuff{nn}(rr,ss1)=nanmedian(abs(MaxIDsp-TrueTime1')./4);
 else
 ErrorShuff{nn}(rr,ss1)=NaN;
  end
 end
 end
 end
 MedErr{ii,ll}=ErrorMaxProb;
 Predictions{ii,ll}=BinsDecoded;
 TotPredicitions{ii,ll}=TotBinsDecoded;
 TotPredicitionsp{ii,ll}=TotBinsDecodedp;
 MedErrShuff{ii,ll}=ErrorShuff;
end
end
toc
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig3/Figures');
cd(GenPath)
save('InstantaneousTimeErrPFC','MedErr','MedErrShuff','Predictions','TotPredicitions','TotPredicitionsp','-v7.3')

