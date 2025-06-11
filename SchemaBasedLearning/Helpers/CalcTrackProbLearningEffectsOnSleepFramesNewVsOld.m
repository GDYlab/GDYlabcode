 
% Change to appropriate directory using files
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
gg=0;rr=0;


nn1=0;Indices1=[];Indices2=[];Probability1=[];
for ii=[6,7,9,10,11,12,14,15]
tic;    
ReplaySleep=[];
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LFPDataV2
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
load SleepReplaysConcat
Sleep1=ReplaySleep;
time1=Sleep1.t;
sig1=Sleep1.pdf;
clear Sleep1 
tic;

ReplaySleep=[]; 
load HPCSleepRipples5thresh
% save('framelags.mat','sig1','sig2','time1','time2');

    tt=cell2mat(cellfun(@(x)(~isempty(x)),sig1,'UniformOutput',false));
    tt1=tt(tt~=0);
    gg=0;
    if ii==6 | ii==15; dd=0; else dd=1;end
   nn1=nn1+1;
     for jj=[1,2,length(tt1)-dd,length(tt1)]   
     gg=gg+1; xx1=[];xx2=[]; rr=0; 
     for aa=[1,2,length(tt1)-dd,length(tt1)]   
     rr=rr+1;
 
     
    binsData=[(LFPData(aa*2-1).t)];
    binsize=round(length(binsData)/10);
    n = numel(binsData);
    bins = mat2cell(binsData,diff([0:binsize:n-1,n]));
    bins=bins(1:10);
     time1sleep=[];
     for xx=1:length(time1{1, aa}{1, jj})
         time1sleep=[time1sleep; time1{1, aa}{1, jj}{1, xx}(1,:)];
     end
     
     RipTimes=[];
     for xx=1:length(Ripples{1, aa})
         RipTimes=[RipTimes; Ripples{1, aa}{1, xx}.Rip1(:,1)];
     end 
     
     
     diffind1=[];
     for xx=1:length(RipTimes)
         difftimes=RipTimes(xx)-time1sleep;
         [mm indmm]=nanmin(abs(difftimes));
         if any(abs(difftimes) <= 0.5)
          diffind1(xx)= time1sleep(indmm); 
         end
     end   
      
        diffind=[];
     for xx=1:length(time1sleep)
         difftimes=time1sleep(xx)-diffind1;
         diffind(xx)= logical(any(difftimes == 0)); 
       
     end   
     
     
for cc=1:length(bins)  
[t ind]=restrict_time(time1sleep(:,1),bins{cc}(1),bins{cc}(end));    



DecodedFrames1=sig1{1, aa}{1, jj}(ind & diffind);
   
if ii < 11 ; nnn1=[2,4,7]; nnn2=[1,5,8];
else nnn1=[1,5,8]; nnn2=[2,4,7]; end
    
bins1=[1:70:561];
bins2=[561:70:1121];
Decoded1=[];Index1=[];
for ll=1:length(DecodedFrames1)
for nn=1:length(bins1)-1
%     DecodedFrames1{ll}=zscore(DecodedFrames1{ll},0,1);
    Probs1=DecodedFrames1{ll}(bins1(1):bins1(end)-1,:);   
%     Probs1=Probs1./sum(Probs1);
    Decoded1(1,nn,ll)=nanmean(nansum(Probs1(bins1(nn):bins1(nn+1)-1,:))');
    Probs1=DecodedFrames1{ll}(bins2(1):bins2(end)-1,:);
%     Probs1=Probs1./sum(Probs1);
    Decoded1(2,nn,ll)=nanmean(nansum(Probs1(bins1(nn):bins1(nn+1)-1,:))');
    
end
end   
x1=nanmean(nanmean(squeeze(Decoded1(1,nnn1,:))'));
x2=nanmean(nanmean(squeeze(Decoded1(2,nnn1,:))'));
x3=nanmean(nanmean(squeeze(Decoded1(1,nnn2,:))'));
x4=nanmean(nanmean(squeeze(Decoded1(2,nnn2,:))'));
Probability1{1,nn1}{gg,rr}(cc,:)=[x1,x2,x3,x4];
end
end
  toc
 end
 clear sig1 sig2 time1 time2
    
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('TrackProbLearnEffectOnSleepFramesHPCNewVsOldRipp.mat','Probability1')

