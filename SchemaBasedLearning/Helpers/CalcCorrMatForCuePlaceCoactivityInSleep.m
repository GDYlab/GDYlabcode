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
thresh=1;rr=0;indices=[];tbins=[];corrcounts=[];
g1=0;CorrMat=[];
for ii=[1:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load mazel
load LFPDataV2
dn=pwd;
LapInfo1=LapInfo;

ss=ii-5;
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
dn1=pwd;
fstruct = dir('*FrameHPCShuff*.mat');
fstruct1=array2table(fstruct);
if ii==6 ;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
aa1=length(fstruct);
if ii==6; 
    aa=8
elseif ii==15; 
    aa=4; 
else aa=6; 
end
if ~isempty(fstruct1);
g1=g1+1;    
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
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
rates=FiringRates1(a12(:,1)==1 & a12(:,2)==1,2);
tt=[1,2,4,5,7,8];
meanreactstrength=[];
gg=0;
if sum(ind_HPC)~0
    rr=rr+1;    
for kk=1:2:length(fstruct)-1;
    tic;
    binsData=[(LFPData(kk).t)];
    binsize=round(length(binsData)/10);
    n = numel(binsData);
    bins = mat2cell(binsData,diff([0:binsize:n-1,n]));
    bins=bins(1:10);
    gg=gg+1;
    cd(dn1);
     filename=fstruct(kk).name;
     load(filename);
     dotLocations = find(filename == '.');
     filename1 = filename(1:dotLocations(1)-1);
     sesslocs=find(filename == 'f');
     sess1 = str2double(filename(sesslocs(end)+1:dotLocations(1)-1));
     assignin ( 'caller', 'framefile',eval(filename1));
     clear((filename1));
     cd(dn);
     framelength=length(framefile);tempt=[];
     
     for a=1:framelength
     tempt(a,:)=[nanmin(squeeze(framefile(a).t(1,1,1,:))), max(squeeze(framefile(a).t(1,1,1,:)))];
     end
     clear framefile 
     temptbins=[];zContribution1=[];Contribution1=[];indc=[];tbins1=[];corrcounts1=[];
     for cc=1:length(bins)
     [t ind]=restrict_time(tempt(:,1),bins{cc}(1),bins{cc}(end));    
     temptbins{cc}=tempt(ind,:);
     ReactStrength= CellAssembliesCorrMatSleepFrames(temptbins{cc},ind_HPC);
%      zContribution1{cc}=zscore(ReactStrength.sigpcs,[],1);
     Contribution1{cc}=ReactStrength.corrmat;
     tbins1{cc}=ReactStrength.tbins;
     corrcounts1{cc}=ReactStrength.corrcounts;
%      indc{cc}=ReactStrength.inds;
     end
CorrMat{ii,gg}=Contribution1;
tbins{ii,gg}=tbins1;
corrcounts{ii,gg}=corrcounts1;
% zContribution{rr,gg}=zContribution1; 
% indices{rr,gg}=indc;
toc
end
end
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SleepFramesCorrMat10binsHPC','CorrMat','corrcounts','tbins');

%% For PFC
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
thresh=1;rr=0;indices=[];tbins=[];corrcounts=[];
g1=0;CorrMat=[];
for ii=[1:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load mazel
load LFPDataV2
dn=pwd;
LapInfo1=LapInfo;

ss=ii-5;
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
dn1=pwd;
fstruct = dir('*FramePFCShuff*.mat');
fstruct1=array2table(fstruct);
if ii==6 ;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
aa1=length(fstruct);
if ii==6; 
    aa=8
elseif ii==15; 
    aa=4; 
else aa=6; 
end
if ~isempty(fstruct1);
g1=g1+1;    
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
rates=FiringRates1(a12(:,1)==1 & a12(:,2)==1,2);
tt=[1,2,4,5,7,8];
meanreactstrength=[];
gg=0;
if sum(ind_HPC)>10
    rr=rr+1;    
for kk=1:2:length(fstruct)-1;
    tic;
    binsData=[(LFPData(kk).t)];
    binsize=round(length(binsData)/10);
    n = numel(binsData);
    bins = mat2cell(binsData,diff([0:binsize:n-1,n]));
    bins=bins(1:10);
    gg=gg+1;
    cd(dn1);
     filename=fstruct(kk).name;
     load(filename);
     dotLocations = find(filename == '.');
     filename1 = filename(1:dotLocations(1)-1);
     sesslocs=find(filename == 'f');
     sess1 = str2double(filename(sesslocs(end)+1:dotLocations(1)-1));
     assignin ( 'caller', 'framefile',eval(filename1));
     clear((filename1));
     cd(dn);
     framelength=length(framefile);tempt=[];
     
     for a=1:framelength
     tempt(a,:)=[nanmin(squeeze(framefile(a).t(1,1,1,:))), max(squeeze(framefile(a).t(1,1,1,:)))];
     end
     clear framefile 
     temptbins=[];zContribution1=[];Contribution1=[];indc=[];tbins1=[];corrcounts1=[];
     for cc=1:length(bins)
     [t ind]=restrict_time(tempt(:,1),bins{cc}(1),bins{cc}(end));    
     temptbins{cc}=tempt(ind,:);
     ReactStrength= CellAssembliesCorrMatSleepFrames(temptbins{cc},ind_HPC);
%      zContribution1{cc}=zscore(ReactStrength.sigpcs,[],1);
     Contribution1{cc}=ReactStrength.corrmat;
     tbins1{cc}=ReactStrength.tbins;
     corrcounts1{cc}=ReactStrength.corrcounts;
%      indc{cc}=ReactStrength.inds;
     end
CorrMat{ii,gg}=Contribution1;
tbins{ii,gg}=tbins1;
corrcounts{ii,gg}=corrcounts1;
% zContribution{rr,gg}=zContribution1; 
% indices{rr,gg}=indc;
toc
end
end
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath);
save('SleepFramesCorrMat10binsPFC','CorrMat','corrcounts','tbins');


%% All sleep
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath);
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
n=1;zreact=[];weight=[];sigpcs=[];
thresh=1;rr=0;indices=[];
g1=0;Contribution=[];zContribution=[];tbins=[];CorrMat=[];corrcounts=[];
for ii=[1:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load mazel
load LFPDataV2
dn=pwd;
LapInfo1=LapInfo;

ss=ii-5;
GenPath=strcat(paths,Replayfile{ii});
cd(GenPath{1})
dn1=pwd;
fstruct = dir('*FrameHPCShuff*.mat');
fstruct1=array2table(fstruct);
if ii==6 ;
fstruct([1:10], :) = fstruct([1,[3:10],2], :);
end
aa1=length(fstruct);
if ii==6; 
    aa=8
elseif ii==15; 
    aa=4; 
else aa=6; 
end
if ~isempty(fstruct1);
g1=g1+1;    
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
ind_HPC=a12(:,1)==1 & a12(:,2)==1;
rates=FiringRates1(a12(:,1)==1 & a12(:,2)==1,2);
tt=[1,2,4,5,7,8];
meanreactstrength=[];
gg=0;
if sum(ind_HPC)~0
    rr=rr+1;    
for kk=1:2:length(fstruct)-1;
    tic;
    binsData=[(LFPData(kk).t)];
    binsize=round(length(binsData)/10);
    n = numel(binsData);
    bins = mat2cell(binsData,diff([0:binsize:n-1,n]));
    bins=bins(1:10); 
    gg=gg+1;   
     temptbins=[];zContribution1=[];Contribution1=[];indc=[];tbins1=[];corrcounts1=[];
     for cc=1:length(bins)     
     temptbins{cc}=[bins{cc}(1), bins{cc}(end)];
     ReactStrength= CellAssembliesCorrMatSleepFrames20P(temptbins{cc},ind_HPC);
%      zContribution1{cc}=zscore(ReactStrength.sigpcs,[],1);
     Contribution1{cc}=ReactStrength.corrmat;
     tbins1{cc}=ReactStrength.tbins;
     corrcounts1{cc}=ReactStrength.corrcounts;
%      indc{cc}=ReactStrength.inds;
     end
CorrMat{ii,gg}=Contribution1;
tbins{ii,gg}=tbins1;
corrcounts{ii,gg}=corrcounts1;

toc
end
end
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath);
save('SleepFramesCorrMat10binsAllSleep20MultiHPCPFC','CorrMat','tbins','corrcounts');

