%% Fig 2e
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SessionWiseRDMMatHPCV2','spk3');
spk{1,1}=spk3;
load('SessionWiseRDMMatPFCV2','spk3');
spk{1,2}=spk3;
corrMat=[];
for i=1:4;
    for j=1:2
        for k=1:36
        corrMat{i,j}{k}=[];
        end
    end
end


Sess={[1:5],[21:25],[6:15],[16:20]};
for ee=2
rr=0;
for ll=1:3       
for kk=1:2 
for jj=1:6 
 rr= rr+1;
for Days=1:length(Sess)
for ii=Sess{Days}
for Reg=1:2
corrMat{Days,Reg}{rr}=[ corrMat{Days,Reg}{rr}; spk{1,Reg}{ii,ll}{jj,kk}];
end
end
end
end
end
end

for Reg=1:2
for Days=1:4

for ii=1:size(corrMat{Days,Reg},2)
    for jj=1:size(corrMat{Days,Reg},2)
         corrMat1(ee,Reg,Days,ii,jj)=nanmean(diag(corr(corrMat{Days,Reg}{1,ii},corrMat{Days,Reg}{1,jj},'rows','pairwise')));   
     end
end

end
end

end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SessionWiseRDMMatHPCSpatialDir1V2','spk3');
spk{1,1}=spk3;
load('SessionWiseRDMMatPFCSpatialDir1V2','spk3');
spk{1,2}=spk3;
load('SessionWiseRDMMatHPCSpatialDir2V2','spk3');
spk{3,1}=spk3;
load('SessionWiseRDMMatPFCSpatialDir2V2','spk3');
spk{3,2}=spk3;
for ee=[1,3]
corrMat=[];
for i=1:4;
    for j=1:2
        for k=1:36
        corrMat{i,j}{k}=[];
        end
    end
end


Sess={[1:5],[21:25],[6:15],[16:20]};

rr=0;
for ll=1:3       
for kk=1:2 
for jj=1:6 
 rr= rr+1;
for Days=1:length(Sess)
for ii=Sess{Days}
for Reg=1:2
corrMat{Days,Reg}{rr}=[ corrMat{Days,Reg}{rr}; spk{ee,Reg}{ii,ll}{jj,kk}];
end
end
end
end
end
end

for Reg=1:2
for Days=1:4

for ii=1:size(corrMat{Days,Reg},2)
    for jj=1:size(corrMat{Days,Reg},2)
       
        corrMat1(ee,Reg,Days,ii,jj)=nanmean(diag(corr(corrMat{Days,Reg}{1,ii},corrMat{Days,Reg}{1,jj},'rows','pairwise')));
             
     end
end

end
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('PlotExampleRSMDataAllDays','corrMat1');