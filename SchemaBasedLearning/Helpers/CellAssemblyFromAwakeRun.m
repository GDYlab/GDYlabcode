
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)

ReactStrength=[];
for ii=[1:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load mazel
dn=pwd;
LapInfo1=LapInfo;
 tic;

% load SelectivityIndicesSignificance.mat
% reading and outputting celltype
cd(dn);
[c1 c2] = textread('celltype','%s %s');
elecreg=importdata('elecreg');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % food is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_PFC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;

meanreactstrength=[];

if sum(ind_HPC) > 10
b=textread('clusters','%s');
b=b(ind_HPC);clg=[];
for c=1:size(b,1)
       clg{c}=cl2mat(b{c});              
       clear clG
end       
    
for kk=1:length(mazel)
   
    Laps=LapInfo(LapInfo(:,5)==kk,:);    
    time_frames1=LapInfo(LapInfo(:,5)==kk,[1,4]);
    if isnan(time_frames1(1))
    time_frames1(1)=Laps(1,3);
    end
    if isnan(time_frames1(end))
    time_frames1(end)=Laps(end,2);
    end
    

     
temptbins=time_frames1;
ReactStrength1= CellAssembliesSigPCsFromRunPFC(temptbins,clg,0.025);      
ReactStrength{ii,kk}=ReactStrength1;

end
end
toc
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakeCellAssembliesPFC25','ReactStrength');


