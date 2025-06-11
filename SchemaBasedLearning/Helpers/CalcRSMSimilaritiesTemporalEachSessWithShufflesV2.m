
% Change to appropriate directory using files

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
glmstats=[];spk3=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
for ii=1:25
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ii,ll}{i,j}=[];
    end
end
end
end
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
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
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates1(ind_HPC,1);
sel=16:40;
if sum(ind_HPC) < 10; continue; end
clear PTables firing11 spk11 sig_cells spk spk1 LapInfo
PlDir=[];plFields=[];Errors=[];
for kk=1:length(PlFieldTemp)
    PlDir(kk)=PlFieldTemp(kk).Direction;
    plFields{kk}=PlFieldTemp(kk).PlaceFields;
    Errors(kk,1)=PlFieldTemp(kk).TrialToReward;
end

PlFields1=plFields(PlDir==1);
PlFields2=plFields(PlDir==2);
error1=Errors(PlDir==1);
if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end 
    gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
    gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
spk2=[]; rr=0;
for j=1:length(rates) 
     rr=rr+1;
  spp1=rates1{j,1}(idx,sel);spp=[];
  for yy=1:size(spp1,1)
      spp(yy,:)=conv2(spp1(yy,:),kernel1,'same');
  end
  spk2(rr,:,:)= spp;
  
end 
    for bb=1:length(nn1)
    rw_spikes=[]; rw_spikes_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;    
    rw_spikes=(squeeze(spk2(:,id1 & id2,:)));
    rw_spikes_mean=nanmean(rw_spikes,2);
    spk3{ii,gg}{bb,1}=[spk3{ii,gg}{bb,1}; (squeeze(rw_spikes_mean(:,1,:)))];
    size1=size( (squeeze(rw_spikes_mean(:,1,:))));
    rw_spikes=[];rw_spikes_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_spikes=(squeeze(spk2(:,id1 & id2,:)));
    rw_spikes_mean=nanmean(rw_spikes,2);
    elseif sum(id1 & id2)==1
%     rw_spikes(:,1,:)=(squeeze(spk2(:,id1 & id2,:)));
%     rw_spikes_mean=rw_spikes;
    rw_spikes_mean(:,1,:)=ones(size1).*NaN; 
    else
    rw_spikes_mean(:,1,:)=ones(size1).*NaN;    
    end
    spk3{ii,gg}{bb,2}=[spk3{ii,gg}{bb,2}; (squeeze(rw_spikes_mean(:,1,:)))];
    end
end
toc
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SessionWiseRDMMatHPCV2','spk3');




cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
glmstats=[];spk3=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
for ii=1:25
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ii,ll}{i,j}=[];
    end
end
end
end
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
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
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates1(ind_HPC,1);
sel=16:40;
if sum(ind_HPC) < 10; continue; end
clear PTables firing11 spk11 sig_cells spk spk1 LapInfo
PlDir=[];plFields=[];Errors=[];
for kk=1:length(PlFieldTemp)
    PlDir(kk)=PlFieldTemp(kk).Direction;
    plFields{kk}=PlFieldTemp(kk).PlaceFields;
    Errors(kk,1)=PlFieldTemp(kk).TrialToReward;
end

PlFields1=plFields(PlDir==1);
PlFields2=plFields(PlDir==2);
error1=Errors(PlDir==1);
if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end 
    gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
    gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
spk2=[]; rr=0;
for j=1:length(rates) 
     rr=rr+1;
  spp1=rates1{j,1}(idx,sel);spp=[];
  for yy=1:size(spp1,1)
      spp(yy,:)=conv2(spp1(yy,:),kernel1,'same');
  end
  spk2(rr,:,:)= spp;
  
end 
    for bb=1:length(nn1)
    rw_spikes=[]; rw_spikes_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;    
    rw_spikes=(squeeze(spk2(:,id1 & id2,:)));
    rw_spikes_mean=nanmean(rw_spikes,2);
    spk3{ii,gg}{bb,1}=[spk3{ii,gg}{bb,1}; (squeeze(rw_spikes_mean(:,1,:)))];
    size1=size( (squeeze(rw_spikes_mean(:,1,:))));
    rw_spikes=[];rw_spikes_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_spikes=(squeeze(spk2(:,id1 & id2,:)));
    rw_spikes_mean=nanmean(rw_spikes,2);
    elseif sum(id1 & id2)==1
%     rw_spikes(:,1,:)=(squeeze(spk2(:,id1 & id2,:)));
%     rw_spikes_mean=rw_spikes;
    rw_spikes_mean(:,1,:)=ones(size1).*NaN; 
    else
    rw_spikes_mean(:,1,:)=ones(size1).*NaN;    
    end
    spk3{ii,gg}{bb,2}=[spk3{ii,gg}{bb,2}; (squeeze(rw_spikes_mean(:,1,:)))];
    end
end
toc
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SessionWiseRDMMatPFCV2','spk3');


GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SessionWiseRDMMatHPCV2','spk3');
spk{1,1}=spk3;
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
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
corrMat1=[];
for Reg=1:2
for Days=1:4

for ii=1:size(corrMat{Days,Reg},2)
    for jj=1:size(corrMat{Days,Reg},2)
        for kk=1:size(corrMat{Days,Reg}{1,ii},2)
        corrMat1(Reg,Days,ii,jj,kk)=(corr(corrMat{Days,Reg}{1,ii}(:,kk),corrMat{Days,Reg}{1,jj}(:,kk),'rows','pairwise'));
        end        
     end
end

end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('BigRSMTemporalCorrMatV2.mat','corrMat1')

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BigRSMTemporalCorrMatV2.mat','corrMat1')
SpecificInfo=[];SessInfo=[];Outcome=[];NovelInfo=[];SessID=[];
outcomei=[ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1];
outcomej=[ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1];

Noveli=[ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,...
    ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1];
Novelj=[ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,...
    ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1];
Across=zeros(36,36);
SessIDi=[ones(1,12),ones(1,12).*3,ones(1,12).*5];
SessID=[];
for ii=1:36
    for jj=1:36
       
        if ii==jj | jj==ii+6 | jj==ii+12 | jj==ii+18 | jj==ii+24 | jj==ii+30 | jj==ii+36 ...
         | jj==ii | jj==ii-6 | jj==ii-12 | jj==ii-18 | jj==ii-24 | jj==ii-30 | jj==ii-36          
        SpecificInfo(ii,jj)=1;
       else
         SpecificInfo(ii,jj)=0;  
       end
       
       if (ii <=12 & jj <= 12) | ((ii > 12 & ii <= 24) &  ((jj > 12 & jj <= 24)) ) | (ii > 24 & jj > 24) 
       SessInfo(ii,jj)=1;
       else
       SessInfo(ii,jj)=-1; 
       end
       Outcome(ii,jj)=outcomei(ii)*outcomej(jj);     
       NovelInfo(ii,jj)=Noveli(ii)*Novelj(jj); 
       if jj>ii
           Across(ii,jj)=logical(1);
       end
       
       SessID(ii,jj)=SessIDi(ii)*SessIDi(jj);
     end
end

SpecificInfo1=SpecificInfo;

SpecificInfo11=SpecificInfo1;
SpecificInfo11(SessInfo==1)=0;

Outcome11=Outcome;
Outcome11(SessInfo==1)=0;

SpecificInfo1(SpecificInfo1==-1)=0; 
Interact=SpecificInfo1.*Outcome;
Interact11=Interact;
Interact11(SessInfo==1)=0;


Pos=[];Neg=[];CPD1=[];CPDS=[];SessIDIDX=[3,5,15];
for Days=1:4
    
    for Reg=1:2
    CorrMatInc=squeeze(corrMat1(Reg,Days,:,:,:));
    for nn=1:size(CorrMatInc,3)
    CorrMatIncNN=squeeze(CorrMatInc(:,:,nn));
    PopVec=CorrMatIncNN(logical(Across & SessInfo==-1 ));
    Group=[];
    Group(:,1)=Outcome(Across & SessInfo==-1 );
    Group(:,2)=SpecificInfo(Across & SessInfo==-1 );     
    Group(:,3)=SessID(Across & SessInfo==-1 );  
       
    
    Response1=PopVec;
    models=[0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1];     
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE= cell2mat(tbl(6,2)); 
          
     models=[0 0 0; 0 1 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPD1(Days,Reg,nn,1)=(SSE1-SSE)./SSE1; 
         
    
    
     models=[1 0 0; 0 0 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;        
    CPD1(Days,Reg,nn,2)=(SSE1-SSE)./SSE1; 
    
     models=[1 0 0; 0 1 0; 0 0 0; 0 0 1];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPD1(Days,Reg,nn,3)=(SSE1-SSE)./SSE1; 

    
    for ww=1:1000
    x=1:size(Group,1);
    y=Group(randperm(numel(x)),:);    
        
    Response1=PopVec;
    models=[0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1];     
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE= cell2mat(tbl(6,2)); 
          
     models=[0 0 0; 0 1 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPDS(Days,Reg,nn,ww,1)=(SSE1-SSE)./SSE1;         
    
    
     models=[1 0 0; 0 0 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;        
    CPDS(Days,Reg,nn,ww,2)=(SSE1-SSE)./SSE1; 
    
     models=[1 0 0; 0 1 0; 0 0 0; 0 0 1];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPDS(Days,Reg,nn,ww,3)=(SSE1-SSE)./SSE1; 
    
        
    end    
           
    end  
    toc
     display(strcat('DoneForDay',num2str(Days),'Reg',num2str(Reg)))
    end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('CPDRSMTemporalV2.mat','CPD1','CPDS','-v7.3')


%% Spatial

% Change to appropriate directory using files
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
glmstats=[];spk3=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
for ii=1:25
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ii,ll}{i,j}=[];
    end
end
end
end
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
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
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates1(ind_HPC,1);
sel=16:40;
if sum(ind_HPC) < 10; continue; end
clear PTables firing11 spk11 sig_cells spk spk1 LapInfo
PlDir=[];plFields=[];Errors=[];
for kk=1:length(PlFieldTemp)
    PlDir(kk)=PlFieldTemp(kk).Direction;
    plFields{kk}=PlFieldTemp(kk).PlaceFields;
    Errors(kk,1)=PlFieldTemp(kk).TrialToReward;
end

PlFields1=plFields(PlDir==1);
PlFields2=plFields(PlDir==2);
plFields11=[];plFields22=[];
for rr=1:size(PlFields1{1,1},1)
    Pl1=(cellfun(@(x)x(rr,:),PlFields1,'UniformOutput',false)) ;
    Pl2=(cellfun(@(x)x(rr,:),PlFields2,'UniformOutput',false)) ;
    for tt=1:length(Pl1)
    plFields11(rr,tt,:)=Pl1{1,tt};
    plFields22(rr,tt,:)=Pl2{1,tt};
    end
    
end
error1=Errors(PlDir==1);
if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end 
    gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
    gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;plFields1=plFields11;plFields2=plFields22;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    plFields1=squeeze(plFields1(:,idx,:));
    plFields2=squeeze(plFields2(:,idx,:));
    for bb=1:length(nn1)
    rw_Place=[]; rw_Place_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;   
      
    rw_Place=squeeze(plFields2(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    spk3{ii,gg}{bb,1}=[spk3{ii,gg}{bb,1}; pl_means];
    size1=size( (squeeze(rw_Place_mean(:,1,:))));
    rw_Place=[];rw_Place_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_Place=squeeze(plFields2(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    elseif sum(id1 & id2)==1
%     rw_Place(:,1,:)=movmean(squeeze(plFields2(ind_HPC,id1 & id2,:)),[1 1],2);
%     rw_Place_mean=rw_Place;
    rw_Place_mean(:,1,:)=ones(size1).*NaN; 
    pl_means=squeeze(rw_Place_mean(:,1,:));
    else
    rw_Place_mean(:,1,:)=ones(size1).*NaN;    
     pl_means=squeeze(rw_Place_mean(:,1,:));
    end
    
    spk3{ii,gg}{bb,2}=[spk3{ii,gg}{bb,2}; pl_means];
    end
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SessionWiseRDMMatHPCSpatialDir2V2','spk3');


% Dir 1
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
glmstats=[];spk3=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
for ii=1:25
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ii,ll}{i,j}=[];
    end
end
end
end
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
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
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==1;
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates1(ind_HPC,1);
sel=16:40;
if sum(ind_HPC) < 10; continue; end
clear PTables firing11 spk11 sig_cells spk spk1 LapInfo
PlDir=[];plFields=[];Errors=[];
for kk=1:length(PlFieldTemp)
    PlDir(kk)=PlFieldTemp(kk).Direction;
    plFields{kk}=PlFieldTemp(kk).PlaceFields;
    Errors(kk,1)=PlFieldTemp(kk).TrialToReward;
end

PlFields1=plFields(PlDir==1);
PlFields2=plFields(PlDir==2);
plFields11=[];plFields22=[];
for rr=1:size(PlFields1{1,1},1)
    Pl1=(cellfun(@(x)x(rr,:),PlFields1,'UniformOutput',false)) ;
    Pl2=(cellfun(@(x)x(rr,:),PlFields2,'UniformOutput',false)) ;
    for tt=1:length(Pl1)
    plFields11(rr,tt,:)=Pl1{1,tt};
    plFields22(rr,tt,:)=Pl2{1,tt};
    end
    
end
error1=Errors(PlDir==1);
if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end 
    gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
    gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;plFields1=plFields11;plFields2=plFields22;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    plFields1=squeeze(plFields1(:,idx,:));
    plFields2=squeeze(plFields2(:,idx,:));
    for bb=1:length(nn1)
    rw_Place=[]; rw_Place_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;   
      
    rw_Place=squeeze(plFields1(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    spk3{ii,gg}{bb,1}=[spk3{ii,gg}{bb,1}; pl_means];
    size1=size( (squeeze(rw_Place_mean(:,1,:))));
    rw_Place=[];rw_Place_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_Place=squeeze(plFields1(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    elseif sum(id1 & id2)==1
%     rw_Place(:,1,:)=movmean(squeeze(plFields1(ind_HPC,id1 & id2,:)),[1 1],2);
%     rw_Place_mean=rw_Place;
    rw_Place_mean(:,1,:)=ones(size1).*NaN; 
    pl_means=squeeze(rw_Place_mean(:,1,:));
    else
    rw_Place_mean(:,1,:)=ones(size1).*NaN;    
     pl_means=squeeze(rw_Place_mean(:,1,:));
    end
    
    spk3{ii,gg}{bb,2}=[spk3{ii,gg}{bb,2}; pl_means];
    end
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SessionWiseRDMMatHPCSpatialDir1V2','spk3');


% PFC Dir 2
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
glmstats=[];spk3=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
for ii=1:25
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ii,ll}{i,j}=[];
    end
end
end
end
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
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
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates1(ind_HPC,1);
sel=16:40;
if sum(ind_HPC) < 10; continue; end
clear PTables firing11 spk11 sig_cells spk spk1 LapInfo
PlDir=[];plFields=[];Errors=[];
for kk=1:length(PlFieldTemp)
    PlDir(kk)=PlFieldTemp(kk).Direction;
    plFields{kk}=PlFieldTemp(kk).PlaceFields;
    Errors(kk,1)=PlFieldTemp(kk).TrialToReward;
end

PlFields1=plFields(PlDir==1);
PlFields2=plFields(PlDir==2);
plFields11=[];plFields22=[];
for rr=1:size(PlFields1{1,1},1)
    Pl1=(cellfun(@(x)x(rr,:),PlFields1,'UniformOutput',false)) ;
    Pl2=(cellfun(@(x)x(rr,:),PlFields2,'UniformOutput',false)) ;
    for tt=1:length(Pl1)
    plFields11(rr,tt,:)=Pl1{1,tt};
    plFields22(rr,tt,:)=Pl2{1,tt};
    end
    
end
error1=Errors(PlDir==1);
if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end 
    gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
    gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;plFields1=plFields11;plFields2=plFields22;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    plFields1=squeeze(plFields1(:,idx,:));
    plFields2=squeeze(plFields2(:,idx,:));
    for bb=1:length(nn1)
    rw_Place=[]; rw_Place_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;   
      
    rw_Place=squeeze(plFields2(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    spk3{ii,gg}{bb,1}=[spk3{ii,gg}{bb,1}; pl_means];
    size1=size( (squeeze(rw_Place_mean(:,1,:))));
    rw_Place=[];rw_Place_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_Place=squeeze(plFields2(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    elseif sum(id1 & id2)==1
%     rw_Place(:,1,:)=movmean(squeeze(plFields2(ind_HPC,id1 & id2,:)),[1 1],2);
%     rw_Place_mean=rw_Place;
    rw_Place_mean(:,1,:)=ones(size1).*NaN; 
    pl_means=squeeze(rw_Place_mean(:,1,:));
    else
    rw_Place_mean(:,1,:)=ones(size1).*NaN;    
     pl_means=squeeze(rw_Place_mean(:,1,:));
    end
    
    spk3{ii,gg}{bb,2}=[spk3{ii,gg}{bb,2}; pl_means];
    end
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SessionWiseRDMMatPFCSpatialDir2V2','spk3');


% Place Fields Dir 1 PFC
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
Animals=([1,2,3,4,5,1,2,3,4,5,1,2,3,4,5]);
Days=([1,1,1,1,1,2,2,2,2,2,3,3,3,3,3]);
glmstats=[];spk3=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
for ii=1:25
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ii,ll}{i,j}=[];
    end
end
end
end
for ii=1:25
GenPath=strcat(paths,files{ii});
cd(GenPath{1});
load LapInfo
load PlFieldTempNLaps
dn=pwd;
LapInfo1=LapInfo;
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
rates=FiringRates1(ind_HPC,1);
rates1=FiringRates1(ind_HPC,1);
sel=16:40;
if sum(ind_HPC) < 10; continue; end
clear PTables firing11 spk11 sig_cells spk spk1 LapInfo
PlDir=[];plFields=[];Errors=[];
for kk=1:length(PlFieldTemp)
    PlDir(kk)=PlFieldTemp(kk).Direction;
    plFields{kk}=PlFieldTemp(kk).PlaceFields;
    Errors(kk,1)=PlFieldTemp(kk).TrialToReward;
end

PlFields1=plFields(PlDir==1);
PlFields2=plFields(PlDir==2);
plFields11=[];plFields22=[];
for rr=1:size(PlFields1{1,1},1)
    Pl1=(cellfun(@(x)x(rr,:),PlFields1,'UniformOutput',false)) ;
    Pl2=(cellfun(@(x)x(rr,:),PlFields2,'UniformOutput',false)) ;
    for tt=1:length(Pl1)
    plFields11(rr,tt,:)=Pl1{1,tt};
    plFields22(rr,tt,:)=Pl2{1,tt};
    end
    
end
error1=Errors(PlDir==1);
if any(ii==[1:5,21:25,18,13,8]) | ii==6 | ii==15; dd=0; else dd=1;end 
    gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
    gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;plFields1=plFields11;plFields2=plFields22;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    plFields1=squeeze(plFields1(:,idx,:));
    plFields2=squeeze(plFields2(:,idx,:));
    for bb=1:length(nn1)
    rw_Place=[]; rw_Place_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;   
      
    rw_Place=squeeze(plFields1(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    spk3{ii,gg}{bb,1}=[spk3{ii,gg}{bb,1}; pl_means];
    size1=size( (squeeze(rw_Place_mean(:,1,:))));
    rw_Place=[];rw_Place_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_Place=squeeze(plFields1(ind_HPC,id1 & id2,:));
    rw_Place_mean=nanmean(rw_Place,2);
    pl_means=movmean((squeeze(rw_Place_mean(:,1,:))),[1 1],2);
    elseif sum(id1 & id2)==1
    rw_Place(:,1,:)=movmean(squeeze(plFields1(ind_HPC,id1 & id2,:)),[1 1],2);
    rw_Place_mean=rw_Place;
%     rw_Place_mean(:,1,:)=ones(size1).*NaN; 
    pl_means=squeeze(rw_Place_mean(:,1,:));
    else
    rw_Place_mean(:,1,:)=ones(size1).*NaN;    
     pl_means=squeeze(rw_Place_mean(:,1,:));
    end
    
    spk3{ii,gg}{bb,2}=[spk3{ii,gg}{bb,2}; pl_means];
    end
end
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('SessionWiseRDMMatPFCSpatialDir1V2','spk3');

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SessionWiseRDMMatHPCSpatialDir1V2','spk3');
spk{1,1}=spk3;

load('SessionWiseRDMMatPFCSpatialDir1V2','spk3');
spk{1,2}=spk3;

load('SessionWiseRDMMatHPCSpatialDir2V2','spk3');
spk{2,1}=spk3;

load('SessionWiseRDMMatPFCSpatialDir2V2','spk3');
spk{2,2}=spk3;
corrMat1=[];
for ee=1:2
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
        for kk=1:size(corrMat{Days,Reg}{1,ii},2)
        corrMat1(ee,Reg,Days,ii,jj,kk)=(corr(corrMat{Days,Reg}{1,ii}(:,kk),corrMat{Days,Reg}{1,jj}(:,kk),'rows','pairwise'));
        end        
     end
end

end
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('BigRSMSpatialCorrMatV2.mat','corrMat1')
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('BigRSMSpatialCorrMatV2.mat','corrMat1')
SpecificInfo=[];SessInfo=[];Outcome=[];NovelInfo=[];SessID=[];
outcomei=[ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1];
outcomej=[ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1,ones(1,6),ones(1,6).*-1];

Noveli=[ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,...
    ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1];
Novelj=[ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,...
    ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1,ones(1,3),ones(1,3).*-1];
Across=zeros(36,36);
SessIDi=[ones(1,12),ones(1,12).*3,ones(1,12).*5];
SessID=[];
for ii=1:36
    for jj=1:36
       
        if ii==jj | jj==ii+6 | jj==ii+12 | jj==ii+18 | jj==ii+24 | jj==ii+30 | jj==ii+36 ...
         | jj==ii | jj==ii-6 | jj==ii-12 | jj==ii-18 | jj==ii-24 | jj==ii-30 | jj==ii-36          
        SpecificInfo(ii,jj)=1;
       else
         SpecificInfo(ii,jj)=0;  
       end
       
       if (ii <=12 & jj <= 12) | ((ii > 12 & ii <= 24) &  ((jj > 12 & jj <= 24)) ) | (ii > 24 & jj > 24) 
       SessInfo(ii,jj)=1;
       else
       SessInfo(ii,jj)=-1; 
       end
       Outcome(ii,jj)=outcomei(ii)*outcomej(jj);     
       NovelInfo(ii,jj)=Noveli(ii)*Novelj(jj); 
       if jj>ii
           Across(ii,jj)=logical(1);
       end
       
       SessID(ii,jj)=SessIDi(ii)*SessIDi(jj);
     end
end

SpecificInfo1=SpecificInfo;

SpecificInfo11=SpecificInfo1;
SpecificInfo11(SessInfo==1)=0;

Outcome11=Outcome;
Outcome11(SessInfo==1)=0;

SpecificInfo1(SpecificInfo1==-1)=0; 
Interact=SpecificInfo1.*Outcome;
Interact11=Interact;
Interact11(SessInfo==1)=0;


Pos=[];Neg=[];CPD1=[];CPDS=[];SessIDIDX=[3,5,15];
for ee=1:2
for Days=1:4
    
    for Reg=1:2
    CorrMatInc=squeeze(corrMat1(ee,Reg,Days,:,:,:));
    for nn=1:size(CorrMatInc,3)
    CorrMatIncNN=squeeze(CorrMatInc(:,:,nn));
    PopVec=CorrMatIncNN(logical(Across & SessInfo==-1 ));
    Group=[];
    Group(:,1)=Outcome(Across & SessInfo==-1 );
    Group(:,2)=SpecificInfo(Across & SessInfo==-1 );     
    Group(:,3)=SessID(Across & SessInfo==-1 );  
       
    
    Response1=PopVec;
    models=[0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1];     
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE= cell2mat(tbl(6,2)); 
          
     models=[0 0 0; 0 1 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPD1(ee,Days,Reg,nn,1)=(SSE1-SSE)./SSE1; 
         
    
    
     models=[1 0 0; 0 0 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;        
    CPD1(ee,Days,Reg,nn,2)=(SSE1-SSE)./SSE1; 
    
     models=[1 0 0; 0 1 0; 0 0 0; 0 0 1];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPD1(ee,Days,Reg,nn,3)=(SSE1-SSE)./SSE1; 

    
    for ww=1:1000
    x=1:size(Group,1);
    y=Group(randperm(numel(x)),:);    
        
    Response1=PopVec;
    models=[0 0 0; 1 0 0; 0 1 0; 1 1 0; 0 0 1];     
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE= cell2mat(tbl(6,2)); 
          
     models=[0 0 0; 0 1 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPDS(ee,Days,Reg,nn,ww,1)=(SSE1-SSE)./SSE1;         
    
    
     models=[1 0 0; 0 0 0; 1 1 0; 0 0 1];   
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;        
    CPDS(ee,Days,Reg,nn,ww,2)=(SSE1-SSE)./SSE1; 
    
     models=[1 0 0; 0 1 0; 0 0 0; 0 0 1];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(5,2)) ;      
    CPDS(ee,Days,Reg,nn,ww,3)=(SSE1-SSE)./SSE1; 
    
        
    end    
           
    end  
    toc
     display(strcat('DoneForDay',num2str(Days),'Reg',num2str(Reg)))
    end
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('CPDRSMSpatialV2.mat','CPD1','CPDS','-v7.3')




