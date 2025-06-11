clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('AwakePFCCAEachGoalActivationGeneralizeCPDShuff1000V4.mat','CPD1','CPDS')
load('PFCCARippleActivationsV2','RippleActivations','RippleActivationsIndex','CPDInd','CPDVal','CPDShuffles','CPD');
load('IndicesOdorPFC')
load('IndicesFlavorPFC')
load('IndicesPlaceArmsPFC')
load('AwakeCellAssemblyActivationInSleepPFC25.mat', ...
        'React');
MembersIndices3=[];sigpcs=[];CPDInd1=[];
DaysInc={[1,3,4,5],[21,23,24,25],[6,8,9,10,11,13,14,15],[16,18,19,20]};
for kk=1:4
rr=0;tot_sum=[];
for ll=1:4;
    sigpcs{kk}{ll}=[];
    MembersIndices3{kk}{ll}=[];
end
for ii=DaysInc{kk}
    tt=cell2mat(cellfun(@(x)(~isempty(x)),Indices_Flavor(ii,:),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    gg=0;
    if any(ii==[1:5,21:25,6,15,8,13,18]); dd=0; else dd=1;end
    ll=0;
    tot_sum1=[];
    for jj=[1,2,length(tt1)-dd,length(tt1)]   
      ll=ll+1;
     x1=[(squeeze(Indices_Flavor{ii,jj}.Pvalues(:,:,1)))<0.05 ] ;
%      x1=[(squeeze(Indices_Flavor{ii,jj}.Prop(:,:,1)))];
     for ww=1:size(React{ii,ll}.sigpcs,2)
     rr=rr+1;    
     MembersIndices3{kk}{ll}=[ MembersIndices3{kk}{ll}; {any(x1')}];
     sigpcs{kk}{ll}=[sigpcs{kk}{ll}; {React{ii,ll}.sigpcs(:,ww)}];
     end
     end
    
end
end
%% sig results

ripphigh=[];ripplow=[];
DaysInc={[1,3,4,5],[21,23,24,25],[6,8,9,10,11,13,14,15],[16,18,19,20]};
sessinc={[1:3],[1:3],[1:4],[1:4]};
Prop1=[];Prop2=[];
for Days=1:4
    for ll=1:2;
    Prop1{Days,ll}=[];
    Prop2{Days,ll}=[];
    Prop1S{Days,ll}=[];
    Prop2S{Days,ll}=[];
    end   
end
for Days=1:4
    rr=0; 
    
for jj=[1,3]
    rr=rr+1;
        perind=CPDInd{Days,1}{jj,1};
        perind=perind((perind~=0));

        Phigh = prctile(perind,50);
        Plow = prctile(perind,50);
        pind1=perind>=0.95;
        pind2=perind<0.95; 
           
        
        flavInd1=MembersIndices3{Days}{jj}(pind1,:);
        flavInd2=MembersIndices3{Days}{jj}(pind2,:);
        
        WightInd1=sigpcs{Days}{jj}(pind1,:);
        WightInd2=sigpcs{Days}{jj}(pind2,:);
        
        ZWeightInd1=cellfun(@(x) zscore(x)>1.5,WightInd1,'UniformOutput',false);
        ZWeightInd2=cellfun(@(x) zscore(x)>1.5,WightInd2,'UniformOutput',false);
       
        for ww=1:size(ZWeightInd1,1)
            Prop1{Days,rr}=[Prop1{Days,rr}; [sum(ZWeightInd1{ww}==1 & flavInd1{ww}'==1)./sum(ZWeightInd1{ww}==1)]];
%               Prop1{Days,rr}=[Prop1{Days,rr}; flavInd1{ww}(ZWeightInd1{ww}==1)'];
        end
        
        for ww=1:size(ZWeightInd2,1)
            Prop2{Days,rr}=[Prop2{Days,rr}; [sum(ZWeightInd2{ww}==1 & flavInd2{ww}'==1)./sum(ZWeightInd2{ww}==1)]];
%             Prop2{Days,rr}=[Prop2{Days,rr}; flavInd2{ww}(ZWeightInd2{ww}==1)'];
        end
     
        
       for cc=1:300  
%         pind1=pind1(randperm(size(pind1,1)),:); 
%         pind2=pind2(randperm(size(pind2,1)),:);
        WightInd1=sigpcs{Days}{jj}(pind1,:);
        WightInd2=sigpcs{Days}{jj}(pind2,:);
        
        ZWeightInd1=cellfun(@(x) zscore(x)>1.5,WightInd1,'UniformOutput',false);
        ZWeightInd2=cellfun(@(x) zscore(x)>1.5,WightInd2,'UniformOutput',false);
       
        for ww=1:size(ZWeightInd1,1)
            ZWeightInd1{ww}=ZWeightInd1{ww}(randperm(size(ZWeightInd1{ww},1)),:); 
            Prop1S{Days,rr}=[Prop1S{Days,rr}; [sum(ZWeightInd1{ww}==1 & flavInd1{ww}'==1)./sum(ZWeightInd1{ww}==1)]];
%               Prop1{Days,rr}=[Prop1{Days,rr}; flavInd1{ww}(ZWeightInd1{ww}==1)'];
        end
        
        for ww=1:size(ZWeightInd2,1)
            ZWeightInd2{ww}=ZWeightInd2{ww}(randperm(size(ZWeightInd2{ww},1)),:); 
            Prop2S{Days,rr}=[Prop2S{Days,rr}; [sum(ZWeightInd2{ww}==1 & flavInd2{ww}'==1)./sum(ZWeightInd2{ww}==1)]];
%             Prop2{Days,rr}=[Prop2{Days,rr}; flavInd2{ww}(ZWeightInd2{ww}==1)'];
        end
       end
  
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('CueCellParticipationPFCCAsRawData','Prop1','Prop2','Prop1S','Prop2S');





