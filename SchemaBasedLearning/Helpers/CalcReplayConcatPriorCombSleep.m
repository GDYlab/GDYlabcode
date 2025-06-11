
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
gg=0;
for ii=[1:25]
ReplaySleep=[];    
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load maze
load LapInfo
dn=pwd;
LapInfo1=LapInfo;
ss=ii-5;

GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})

fstruct = dir('*FrameHPCDecodeConcatComb*.mat');
fstruct1=array2table(fstruct);

if ii==6;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end

if ~isempty(fstruct1)
 gg=gg+1;   
for i=1:2:(length(maze)*2)
     filename=fstruct(i).name;
     load(filename);
     dotLocations = find(filename == '.')
     filename1 = filename(1:dotLocations(1)-1)
     sesslocs=find(filename == 'b');
     sess1 = str2double(filename(sesslocs(1)+1:dotLocations(1)-1));
     assignin ( 'caller', 'framefile',eval(filename1));
     clear((filename1));
     sess1=i/2+0.5;
    
     
     

       tt=0;start_pdf=[];time=[];
           
         for a=1:length(framefile.HPCt) 
         tempt=(framefile.HPCt{a});        
         if size(tempt,2) <= 5; continue; end
         tt=tt+1;
         start_pdf{tt}=(framefile.HPCProb{a});
         time{tt}=tempt;
         
         end
 



ReplaySleep.t{sess1}=time;
ReplaySleep.pdf{sess1}=start_pdf;

end
clear framefile
save('SleepReplaysConcatCombo','ReplaySleep', '-v7.3');
end
end


    



%% For PFC

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
gg=0;
for ii=[1:25]
ReplaySleep=[];    
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load maze
load LapInfo
dn=pwd;
LapInfo1=LapInfo;
ss=ii-5;
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})

fstruct = dir('*FramePFCDecodeConcatComb*.mat');
fstruct1=array2table(fstruct);

if ii==6;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end

if ~isempty(fstruct1)
 gg=gg+1;   
for i=1:2:(length(maze)*2)
     filename=fstruct(i).name;
     load(filename);
     dotLocations = find(filename == '.')
     filename1 = filename(1:dotLocations(1)-1)
     sesslocs=find(filename == 'b');
     sess1 = str2double(filename(sesslocs(1)+1:dotLocations(1)-1));
     assignin ( 'caller', 'framefile',eval(filename1));
     clear((filename1));
     sess1=i/2+0.5;
    
     
   

       tt=0;start_pdf=[];time=[];
           
         for a=1:length(framefile.PFCt) 
         tempt=(framefile.PFCt{a});        
         if size(tempt,2) <= 5; continue; end
         tt=tt+1;
         start_pdf{tt}=(framefile.PFCProb{a});
         time{tt}=tempt;
         
         end
 



ReplaySleep.t{sess1}=time;
ReplaySleep.pdf{sess1}=start_pdf;

end
clear framefile
save('SleepReplaysConcatPFCComb','ReplaySleep', '-v7.3');
end
end
    



