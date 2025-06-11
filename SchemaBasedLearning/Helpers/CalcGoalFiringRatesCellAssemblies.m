clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('AwakeCellAssembliesPFC25','ReactStrength');
AssemblyTuning=[];
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,5,5,5,5,5,4,4,4,4,4]);
for ii=1:25
   
GenPath=strcat(paths,files{ii});
cd(GenPath{1})

load maze
load LapInfo

cluster=importdata('clusters');
elecreg=importdata('elecreg');
% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % good is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;


b=textread('clusters','%s');
b=b(ind_HPC);clg=[];
for c=1:size(b,1)
       clg{c}=cl2mat(b{c});              
       clear clG
end 

tt=cell2mat(cellfun(@(x)(~isempty(x)),ReactStrength(ii,:),'UniformOutput',false));
tt1=tt(tt~=0);
if isempty(tt1); continue; end  
if any(ii==[1:5,21:25,8,13,18]) | ii==6 | ii==15 ; dd=0; else dd=1;end
gg=0;
 for jj=[1,2,length(tt1)-dd,length(tt1)] 
     gg=gg+1;
    tic;
    ReactStrength1=ReactStrength{ii, jj};
    PlFields=[];
    for kk =1:size(LapInfo,1)        
    lap_time=[LapInfo(kk,[7])-1 LapInfo(kk,[7])+5];
    
    if ~any(isnan(lap_time))
    Psth=CellAssembliesToFiringRates(lap_time,ReactStrength1,clg,0.025);
    PlFields(kk,:,:)=Psth.reactstrength;
    else      
    PlFields(kk,:,:)=ones(size(ReactStrength1.sigpcs,2),241).*NaN; 
    end
    end
AssemblyTuning{ii,gg}=PlFields;
    toc
end
display(strcat('DoneForFolder',num2str(ii)));
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakeCAAssemblyGoalTuningPFC','AssemblyTuning','-v7.3');