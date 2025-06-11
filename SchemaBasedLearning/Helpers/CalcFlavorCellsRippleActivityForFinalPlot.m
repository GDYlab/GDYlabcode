%% 100 ms calculations
clear all; close all;
Sess={[1,3,4,5],[21,23,24,25],[6,8,9,10,11,13,14,15],[16,18,19,20]};
for Days=1:length(Sess)
    
        for kk=1:4
            for nn=1:2
            RippleActivations{Days,nn}{1,kk}=[];
            RippleActivationsIndex{Days,nn}{1,kk}=[];
            FID{Days,kk}=[];
            end
        end
    
end
sd=5;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);

cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('IndicesFlavorPFC','Indices_Flavor')
% 
for Days=1:length(Sess)
tic;    
for ii=Sess{Days}
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
dn=pwd;  
load mazel
cluster=importdata('clusters');
elecreg=importdata('elecreg');
tic;
ss=ii-5;
% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % good is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;

load ('SingleCellPeriRippActivationZsc100');

        for kk=1:4
            xx1=squeeze(Activations{1,kk}(:,ind_HPC,:));
            xx1=nanmean(xx1,1);
            xx1=squeeze(xx1(1,:,:));
            xx1=movmean(xx1,[2 2],2);
           
            RippleActivations{Days,1}{1,kk}=[RippleActivations{Days,1}{1,kk}; xx1];
      
            
        end
 if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end         
MembersIndices33=[];
gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
gg=gg+1;       
MembersIndices33=(squeeze(Indices_Flavor{ii,kk}.Pvalues(:,:,1)))<0.05 ;
FID{Days,gg}=[FID{Days,gg};  any(MembersIndices33')'];
end  
   
toc 
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('PFCFlavorID&RippleActivations100','RippleActivations','FID');

