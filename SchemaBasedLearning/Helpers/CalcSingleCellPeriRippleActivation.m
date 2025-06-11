%% 25 ms

%% With Ripples
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
gg=0;rr=0;
for i=1
    for j=1:4
        ddt{i,j}=[];
        ddt1{i,j}=[];
    end
end
thresh=3;
% for ii=[1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,19,20,21,22,24,25]
for ii=[13,18,23]
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
dn=pwd;
load SleepRates
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})

% save('framelags.mat','sig1','sig2','time1','time2');
load HPCSleepRipples3thresh
Activations=[];SActivations=[];
tic;
    tt=cell2mat(cellfun(@(x)(~isempty(x)),Ripples,'UniformOutput',false));
    tt1=tt(tt~=0);
    
    if any(ii==1:5) | any(ii==[21:25,3,8,13,18,23]) | ii==6 | ii==15; dd=0; else dd=1;end
    rr=0;
    for aa=[1,2,length(tt1)-dd,length(tt1)]   
     rr=rr+1;    
     RipTimes=[];
     for xx=1:length(Ripples{1, aa})
         RipTimes=[RipTimes; Ripples{1, aa}{1, xx}.Rip1(:,1)];
     end
     Reactivations=FiringRates1{aa,1};
     ReactTimes=FiringRates1{aa,2}(1,:);
     ReactStrength=Reactivations;ZReactStrength=[];
     for cc=1:size(ReactStrength,1)
     ZReactStrength(cc,:)=zscore(ReactStrength(cc,:));
     end
     
     time_frame22=RipTimes;IndAround=100;
     RippAct=[];ff=0; 
     for ll=1:size(time_frame22,1)
     [minValue,closestIndex] = min(abs(ReactTimes-time_frame22(ll)));
     if (closestIndex-IndAround) <= 0 | (closestIndex+IndAround) >= size(ZReactStrength,2); continue; end
     ff=ff+1;
     RippAct(ff,:,:)=ZReactStrength(:,closestIndex-IndAround:closestIndex+IndAround); 
     end
     Activations{1,rr}=RippAct;
 
end
toc
display(strcat('DoneForFolder',num2str(ii)));
cd(dn);
save('SingleCellPeriRippActivationZsc100','Activations','-v7.3');
end

