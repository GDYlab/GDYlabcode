clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath);n=1;PlFieldTemp=struct; 

Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LFPDataV2
load mazel
dn=pwd;
cluster=importdata('clusters');
elecreg=importdata('elecreg');
SleepID=[1,3,5,7,9];
FiringRates1=[];
for i=1:length(cluster)
   c1G=cl2mat(cluster{i,1}); 
   t_spike = c1G.featuredata (:, 8);
   clear c1G
   for j=1:size(mazel,2)
   t_range=[LFPData(SleepID(j)).t(1) LFPData(SleepID(j)).t(end)];
   t_bin=0.01;
   Psth = CalcFiringRatesInSleep(t_spike, t_range, t_bin);
   FiringRates1{j,1}(i,:)=Psth(:,2)';  
   FiringRates1{j,2}(i,:)=Psth(:,1)';  
   end
   
end
save('SleepRates.mat', ...
        'FiringRates1','-v7.3');
display(strcat('Folder',num2str(ii),'completed'));    
     
    
end


