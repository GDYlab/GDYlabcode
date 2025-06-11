%% Fig 3a and b

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)

gg=0;rr=0;
Probability1=[];
Sess={[1,2,4,5],[21,22,24,25],[6,7,9,10,11,12,14,15],[16,17,19,20]};
for sess=1:length(Sess)
nn1=0;Indices1=[];Indices2=[];
for ii=Sess{sess}
tic;    
ReplaySleep=[];
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LFPDataV2
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
load SleepReplaysConcatCombo
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
    if any(ii==[1:5,21:25]) | ii==6 | ii==15; dd=0; else dd=1;end
    nn1=nn1+1;
     for jj=1:2  
     gg=gg+1; xx1=[];xx2=[]; rr=0; 
     for aa=[1,2,length(tt1)-dd,length(tt1)]   
    rr=rr+1;
    binsData=[(LFPData(aa*2-1).t)];
    binsize=round(length(binsData)/10);
    n = numel(binsData);
    bins = mat2cell(binsData,diff([0:binsize:n-1,n]));
    bins=bins(1:10);
    time1sleep=[];
     for xx=1:length(time1{1, aa})
         time1sleep=[time1sleep; time1{1, aa}{1, xx}(:,1)];
     end
        RipTimes=[];
     for xx=1:length(Ripples{1, aa})
         RipTimes=[RipTimes; Ripples{1, aa}{1, xx}.Rip1(:,1)];
     end 
     
     
     diffind1=[];
     for xx=1:length(RipTimes)
         difftimes=RipTimes(xx)-time1sleep;
         [mm indmm]=nanmin(abs(difftimes));
         if any(abs(difftimes) < 0.50)
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
  
DecodedFrames1=sig1{1, aa}(ind & diffind);
   
nnn1=[1,2,4,5,7,8];

bins11=[1:70:561];
bins22=[561:70:1121];
bins33=[1121:70:1681];
bins44=[1681:70:2241];
Decoded1=[];Index1=[];
for ll=1:length(DecodedFrames1)
if jj==1; bins1=bins11; bins2=bins22; 
elseif jj==2; bins1=bins33; bins2=bins44;end

for nn=1:length(bins1)-1
%     DecodedFrames1{ll}=zscore(DecodedFrames1{ll},0,1);
    Probs1=DecodedFrames1{ll};   
%     Probs1=Probs1./sum(Probs1);
    Decoded1(1,nn,ll)=nanmean(nansum(Probs1(bins1(nn):bins1(nn+1)-1,:))');
    Probs1=DecodedFrames1{ll};
%     Probs1=Probs1./sum(Probs1);
    Decoded1(2,nn,ll)=nanmean(nansum(Probs1(bins2(nn):bins2(nn+1)-1,:))');
    
end
end   
x1=nanmean(nanmean(squeeze(Decoded1(1,nnn1,:))'));
x2=nanmean(nanmean(squeeze(Decoded1(2,nnn1,:))'));
Probability1{sess,nn1}{gg,rr}(cc,:)=[x1,x2];
end
end
 toc
 end
 clear sig1 sig2 time1 time2    
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('TrackProbLearnEffectOnSleepFramesHPCRipp.mat','Probability1')
