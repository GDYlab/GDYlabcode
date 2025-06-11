
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('TrackProbsSessionsAllTrials50','TrackProbsStacked','Trackbins','Trackbinsp','Trackdprobp') 

GoalProbs=[];sel=1:121;
MeanProbs1=[];MeanProbsShuff1=[];Errors1=[];pp=1;
 for ii=1:25
     tic;
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
dn1=pwd;
load LapInfo  
LapInfo1=LapInfo;
  tt=cell2mat(cellfun(@(x)(~isempty(x)),TrackProbsStacked(ii,:),'UniformOutput',false));
  tt1=tt(tt~=0);
  if isempty(tt1); continue; end   
   if any(ii==[1:5,21:25]) | ii==6 | ii==15 ; dd=0; else dd=1;end
    gg=0;
     for ll=[1,2,length(tt1)-dd,length(tt1)]   
     gg=gg+1;
         
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
    last_info= [NaN; LapInfo1(ind1(2:end)-1,6)];
      
    
NumberofErrors=LapInfo(:,11);
binns=[1,7,13,length(LapInfo)+1]; 
ProbOther1=[];  ProbOther2=[]; ProbOther3=[];MeanProbTarg=[]; MeanProbNext=[]; MeanProblast=[]; Err=[];      
 MeanProbTargShuff=[]; MeanProbNextShuff=[];MeanProblastShuff=[];
  arms=[1,2,0,3,4,0,5,6];    
  nn=1;
for rr=1:size(target_info,1)
 dprob=TrackProbsStacked{ii,ll}{rr};
%          dprob(isnan(dprob))=0;
 bins=[1:70:421;421:70:841];
 MeanProb=[];
 for ss=1:6
rr1=[];     
meanPP1=((dprob(bins(1,ss):bins(1,ss+1)-1,sel)));
% meanPP1=movmean(meanPP1,[2 4],2,'omitnan');
rr1(1,:,:)=meanPP1;
meanPP1=((dprob(bins(2,ss):bins(2,ss+1)-1,sel)));
% meanPP1=movmean(meanPP1,[2 4],2,'omitnan');
rr1(2,:,:)=meanPP1;
rr1=nanmean(rr1,1);rr1=squeeze(rr1(1,:,:));
MeanProb(rr,ss,:,:)=rr1;
 end
Err(rr)=NumberofErrors(rr)-1;
MeanProbTarg{rr}= squeeze(MeanProb(rr,arms(target_info(rr)),:,:))';
if arms(Next_info(rr))~=0;
MeanProbNext{rr}=squeeze(MeanProb(rr,arms(Next_info(rr)),:,:))';
else
MeanProbNext{rr}=ones(size(squeeze(MeanProb(rr,arms(target_info(rr)),:,:))')).*NaN;
end

if rr > 1 && arms(last_info(rr))~=0; 
MeanProblast{rr}=squeeze(MeanProb(rr,arms(last_info(rr)),:,:))';
else
MeanProblast{rr}=ones(size(squeeze(MeanProb(rr,arms(target_info(rr)),:,:))')).*NaN;
end


Oth=setdiff([1,2,3,4,5,6],[arms(target_info(rr))]);
ProbOther=squeeze(MeanProb(rr,Oth,:,:));
ProbOther=nanmean(ProbOther,1);
ProbOther1{rr}=squeeze(ProbOther(1,:,:))';

Oth=setdiff([1,2,3,4,5,6],[arms(Next_info(rr))]);
ProbOther=squeeze(MeanProb(rr,Oth,:,:));
ProbOther=nanmean(ProbOther,1);
ProbOther2{rr}=squeeze(ProbOther(1,:,:))';


if rr > 1 && arms(last_info(rr))~=0;
Oth=setdiff([1,2,3,4,5,6],[arms(last_info(rr))]);
ProbOther=squeeze(MeanProb(rr,Oth,:,:));
ProbOther=nanmean(ProbOther,1);
ProbOther3{rr}=squeeze(ProbOther(1,:,:))';
else
ProbOther3{rr}=ones(size(squeeze(MeanProb(rr,arms(target_info(rr)),:,:))')).*NaN;
end
%% for shuffle
 dprob=TrackProbsStacked{ii,ll}{rr};
 
%          dprob(isnan(dprob))=0;
 bins=[1:70:421;421:70:841];MeanProb=[];
 for dd1=1:300
 shuffids=circshift([1:840], round(rand(1).*150));    
 dprob1=dprob(shuffids,:);
 for ss=1:6
rr1=[];     
meanPP1=((dprob1(bins(1,ss):bins(1,ss+1)-1,sel)));
% meanPP1=movmean(meanPP1,[2 4],2,'omitnan');
rr1(1,:,:)=meanPP1;
meanPP1=((dprob1(bins(2,ss):bins(2,ss+1)-1,sel)));
% meanPP1=movmean(meanPP1,[2 4],2,'omitnan');
rr1(2,:,:)=meanPP1;
rr1=nanmean(rr1,1);rr1=squeeze(rr1(1,:,:));
MeanProb(dd1,rr,ss,:,:)=rr1;
 end
 end
 xx1=squeeze(MeanProb(:,rr,arms(target_info(rr)),:,:));
 xx2=nansum(xx1,2);
xx3=nanmean(squeeze(xx2(:,1,:)));
MeanProbTargShuff{rr}= xx3;
if arms(Next_info(rr))~=0;
xx1= squeeze(MeanProb(:,rr,arms(Next_info(rr)),:,:));   
xx2=nansum(xx1,2);
xx3=nanmean(squeeze(xx2(:,1,:)));    
MeanProbNextShuff{rr}=xx3;
else
MeanProbNextShuff{rr}=ones(size(xx3)).*NaN;
end

if rr > 1 && arms(last_info(rr))~=0; 
xx1=squeeze(MeanProb(:,rr,arms(last_info(rr)),:,:));
xx2=nansum(xx1,2);
xx3=nanmean(squeeze(xx2(:,1,:)));       
MeanProblastShuff{rr}=xx3;
else
MeanProblastShuff{rr}=ones(size(xx3)).*NaN;
end

nn=nn+1;
  

end
MeanProbs1{ii,gg,1}=MeanProbTarg;
MeanProbs1{ii,gg,2}=MeanProbNext;
MeanProbs1{ii,gg,3}=MeanProblast;
MeanProbs1{ii,gg,4}=ProbOther1;
MeanProbs1{ii,gg,5}=ProbOther2;
MeanProbs1{ii,gg,6}=ProbOther3;
MeanProbsShuff1{ii,gg,1}=MeanProbTargShuff;
MeanProbsShuff1{ii,gg,2}=MeanProbNextShuff;
MeanProbsShuff1{ii,gg,3}=MeanProblastShuff;
Errors1{ii,gg}=Err;
toc
 end
display(strcat('DoneForSess',num2str(ii))) 
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('StartArmDecodingSpatial50Shuff','MeanProbs1','MeanProbsShuff1',...
    'Errors1','-v7.3') 



 
 
