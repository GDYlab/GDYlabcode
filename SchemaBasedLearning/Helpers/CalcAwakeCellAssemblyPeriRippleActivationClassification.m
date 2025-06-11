
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
gg=0;rr=0;

load('AwakeCellAssemblyActivationInSleepPFC25.mat', ...
        'React');
for i=1
    for j=1:4
        ddt{i,j}=[];
        ddt1{i,j}=[];
    end
end
thresh=3;Activations=[];SActivations=[];
for ii=[1,3,4,5,6,8,9,10,11,13,14,15,16,18,19,20,21,23,24,25]
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})

% save('framelags.mat','sig1','sig2','time1','time2');
load HPCSleepRipples3thresh

tic;
    tt=cell2mat(cellfun(@(x)(~isempty(x)),Ripples,'UniformOutput',false));
    tt1=tt(tt~=0);
    gg=0;
    if any(ii==1:5) | any(ii==21:25)| ii==8 | ii==13 | ii==18 | ii==6 | ii==15; dd=0; else dd=1;end
    rr=0;
    for aa=[1,2,length(tt1)-dd,length(tt1)]   
     rr=rr+1;    
     RipTimes=[];
     for xx=1:length(Ripples{1, aa})
         RipTimes=[RipTimes; Ripples{1, aa}{1, xx}.Rip1(:,1)];
     end
     for jj=1:4
     Reactivations=React{ii,jj,rr};
     CellCounts=Reactivations.cc;
     ReactTimes=Reactivations.t;
     ReactStrength=Reactivations.reactstrength;ZReactStrength=[];
     for cc=1:size(ReactStrength,1)
     ZReactStrength(cc,:)=zscore(Reactivations.reactstrength(cc,:));
     end
     ZReactStrength(ReactStrength<thresh)=0;   
  
     time_frame22=RipTimes;
     RippAct=[];ff=0; SurrAct=[];
     for ll=1:size(time_frame22,1)
     randInd=0;    
     while (randInd-20) <= 0 | (randInd+20) >= size(ZReactStrength,2)    
     randInd=(randi(numel(ZReactStrength)));
     end
     [minValue,closestIndex] = min(abs(ReactTimes-time_frame22(ll)));
     if (closestIndex-20) <= 0 | (closestIndex+20) >= size(ZReactStrength,2); continue; end
     ff=ff+1;
     RippAct(ff,:,:)=ZReactStrength(:,closestIndex-20:closestIndex+20); 
     SurrAct(ff,:,:)=ZReactStrength(:,randInd-20:randInd+20);
     end
     Ripples1=((RippAct>0));
     Activations{ii,jj}{rr}=Ripples1;
     
     SurrAct=[];
     for ww=1:1000
     ff=0;    
     for ll=1:size(time_frame22,1)
     randInd=0;    
     while (randInd-20) <= 0 | (randInd+20) >= size(ZReactStrength,2)    
     randInd=(randi(numel(ZReactStrength)));
     end
     [minValue,closestIndex] = min(abs(ReactTimes-time_frame22(ll)));
     if (closestIndex-20) <= 0 | (closestIndex+20) >= size(ZReactStrength,2); continue; end
     ff=ff+1;     
     SurrAct(ww,ff,:,:)=ZReactStrength(:,randInd-20:randInd+20);
     end
     end
     Ripples1=((SurrAct>0));
     Ripples2=nanmean(Ripples1,1);
     Ripples1=squeeze(Ripples2(1,:,:,:));
     SActivations{ii,jj}{rr}=Ripples1; 
 
end
toc
end
display(strcat('DoneForFolder',num2str(ii)));
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakeCAActivationsAround3Ripples3ZThPFC25msV21000','Activations','SActivations','-v7.3');

