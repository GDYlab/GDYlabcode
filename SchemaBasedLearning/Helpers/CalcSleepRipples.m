
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
n=1;zreact=[];weight=[];sigpcs=[];
thresh=1;rr=0;indices=[];
g1=0;Contribution=[];zContribution=[];
for ii=[1:25]
    
% if ~any(ii== ([3,8,13])); continue; end
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load maze
load LFPDataV2
Filtered=[];Times=[];UnFilt=[];j=0;
for nn=1:2:length(LFPData)
    j=j+1;
Filtered{j}=LFPData(nn).firLFP{3};
Times{j}=LFPData(nn).t;
UnFilt{j}=LFPData(nn).LFP;
end
clear LFPData
dn=pwd; Ripples=[];
rr=rr+1;gg=0;
for kk=1:length(Times)
    tic;
    binsData=[(Times{kk})];
    binsize=round(length(binsData)/10);
    n = numel(binsData);
    bins = mat2cell(binsData,diff([0:binsize:n-1,n]));
    bins=bins(1:10);  
    gg=gg+1;   
     temptbins=[];zContribution1=[];Contribution1=[];indc=[];Filtered1=[];
     for cc=1:length(bins)       
     temptbins{cc}=Times{kk}(Times{kk}>=bins{cc}(1) & Times{kk}<=bins{cc}(end));
     Filtered1{cc}=Filtered{kk}(Times{kk}>=bins{cc}(1) & Times{kk}<=bins{cc}(end));
     Ripples1= RippleDetectRatGDlab_BRB(temptbins{cc},Filtered1{cc},5);
     zContribution1{cc}=Ripples1;
     end
Ripples{gg}=zContribution1; 
toc
end
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
save('HPCSleepRipples5thresh.mat','Ripples','-v7.3');
clear Ripples
end

