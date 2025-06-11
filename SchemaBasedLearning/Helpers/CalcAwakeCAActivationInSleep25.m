clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
cd('/media/baburam/DiskStorage1/Database/CellAssemblies/Data');
load('AwakeCellAssembliesPFC25');
sd=4;
Fss = 1;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel2 = kernel./sum(kernel);
React=[];
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load mazel
load LapInfo
load LFPDataV2

dn=pwd;
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
tic;

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
tic;
for ll=[1,2,length(tt1)-dd,length(tt1)]; 
ReactStrength1=ReactStrength{ii,ll};     
gg=gg+1; rr=0;
kk1=1:2:(length(mazel)*2)-1;
for kk=[1,2,length(tt1)-dd,length(tt1)];
rr=rr+1;
binsData=[nanmin(LFPData(kk1(kk)).t) nanmax(LFPData(kk1(kk)).t)];
   
temptbins=binsData;React1=[];
React1 = CellAssembliesToFiringRates(temptbins,ReactStrength1,clg,0.025);
React{ii,gg,rr}=React1;   
end
toc  
end
display(strcat('Folder',num2str(ii),'completed'));  
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakeCellAssemblyActivationInSleepPFC25.mat', ...
        'React','-v7.3');
    
