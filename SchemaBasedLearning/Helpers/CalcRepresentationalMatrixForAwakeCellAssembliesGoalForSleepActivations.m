
%%PFC
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
glmstats=[];spk3=[];
load AwakeCAAssemblyGoalTuningPFC
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
sess={[1,3,4,5],[21,23,24,25]};
fig1=figure;
oo=0;
for nn=1:length(sess)
    oo=oo+1;
  sessinc=[1:3];
for pp=sessinc  
for ll=1:4
for i=1:6
    for j=1:2
        spk3{oo,pp}{ll}{i,j}=[];
    end
end
end
for ii=sess{nn}
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load PlFieldTempNLaps

dn=pwd;
LapInfo1=LapInfo;
load mazel
cluster=importdata('clusters');
elecreg=importdata('elecreg');
tic;

% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % good is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;


if sum(ind_HPC) < 10; continue; end

clear PTables firing11 spk11 sig_cells spk spk1 LapInfo

if any(ii==[1:5,21:25,8,13,18]) | ii==6 | ii==15 ; dd=0; else dd=1;end

gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
rates1=AssemblyTuning{ii,pp}; rates=[];
bins=[1:10:241];rates1(rates1<3)=0;  
for ww=1:length(bins)-1
 r=squeeze(rates1(:,:,bins(ww):bins(ww+1)-1)); 
 r= nanmean(r,3);
 rates(:,:,ww)=r;
end
    
gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    spk2=squeeze(rates(idx,:,:));
    for bb=1:length(nn1)
    rw_spikes=[]; rw_spikes_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;    
    rw_spikes=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=nanmean(rw_spikes,1);
    spk3{oo,pp}{gg}{bb,1}=[spk3{oo,pp}{gg}{bb,1}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    size1=size( (squeeze(rw_spikes_mean(1,:,:))));
    rw_spikes=[];rw_spikes_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_spikes=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=nanmean(rw_spikes,1);
    elseif sum(id1 & id2)==1
    rw_spikes(1,:,:)=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=rw_spikes;
    else
    rw_spikes_mean(1,:,:)=ones(size1).*NaN;    
    end
    spk3{oo,pp}{gg}{bb,2}=[spk3{oo,pp}{gg}{bb,2}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    end
end
end


end
end





glmstats=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
sess={[6,8,9,10,11,13,14,15]}
fig1=figure;
for nn=1:length(sess)
    oo=oo+1;
   
  sessinc=[1:4];
for pp=sessinc  
for ll=1:4
for i=1:6
    for j=1:2
        spk3{oo,pp}{ll}{i,j}=[];
    end
end
end
for ii=sess{nn}
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load PlFieldTempNLaps

dn=pwd;
LapInfo1=LapInfo;
load mazel
cluster=importdata('clusters');
elecreg=importdata('elecreg');
tic;

% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % good is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;


if sum(ind_HPC) < 10; continue; end

clear PTables firing11 spk11 sig_cells spk spk1 LapInfo

if any(ii==[1:5,21:25,8,13,18]) | ii==6 | ii==15 ; dd=0; else dd=1;end

gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
rates1=AssemblyTuning{ii,pp}; rates=[];
bins=[1:10:241];rates1(rates1<3)=0;  
for ww=1:length(bins)-1
 r=squeeze(rates1(:,:,bins(ww):bins(ww+1)-1)); 
 r= nanmean(r,3);
 rates(:,:,ww)=r;
end
    
gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    spk2=squeeze(rates(idx,:,:));
    for bb=1:length(nn1)
    rw_spikes=[]; rw_spikes_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;    
    rw_spikes=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=nanmean(rw_spikes,1);
    spk3{oo,pp}{gg}{bb,1}=[spk3{oo,pp}{gg}{bb,1}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    size1=size( (squeeze(rw_spikes_mean(1,:,:))));
    rw_spikes=[];rw_spikes_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_spikes=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=nanmean(rw_spikes,1);
    elseif sum(id1 & id2)==1
    rw_spikes(1,:,:)=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=rw_spikes;
    else
    rw_spikes_mean(1,:,:)=ones(size1).*NaN;    
    end
    spk3{oo,pp}{gg}{bb,2}=[spk3{oo,pp}{gg}{bb,2}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    end
end
end

end
end



glmstats=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
sess={[16,18,19,20]}
fig1=figure;
for nn=1:length(sess)
    oo=oo+1;
  sessinc=[1:4];
for pp=sessinc  
for ll=1:4
for i=1:6
    for j=1:2
        spk3{oo,pp}{ll}{i,j}=[];
    end
end
end
for ii=sess{nn}
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
load LapInfo
load PlFieldTempNLaps

dn=pwd;
LapInfo1=LapInfo;
load mazel
cluster=importdata('clusters');
elecreg=importdata('elecreg');
tic;

% reading and outputting celltype
[c1 c2] = textread('celltype','%s %s');
a1 = strrep(c1,'p','1'); % pyrs is 1
a1 = str2num(cell2mat(strrep(a1,'i','0'))); % int is 0, output is number matrix
a2 = strrep(c2,'b','0'); % bad is 0
a2 = strrep(a2,'g','1'); % good is 1
a2 = str2num(cell2mat(strrep(a2,'f','2'))); % good is 2, output is number matrix
a12 = [a1 a2]; % concatenating the 2 vectors
ind_HPC=a12(:,1)==1 & a12(:,2)==1 & elecreg(:,2)==2;


if sum(ind_HPC) < 10; continue; end

clear PTables firing11 spk11 sig_cells spk spk1 LapInfo

if any(ii==[1:5,21:25,8,13,18]) | ii==6 | ii==15 ; dd=0; else dd=1;end

gg=0;
for kk=[1,2,length(mazel)-dd,length(mazel)] 
rates1=AssemblyTuning{ii,pp}; rates=[];
bins=[1:10:241];rates1(rates1<3)=0;  
for ww=1:length(bins)-1
 r=squeeze(rates1(:,:,bins(ww):bins(ww+1)-1)); 
 r= nanmean(r,3);
 rates(:,:,ww)=r;
end
    
gg=gg+1;
if ii< 11; nn1=[2,4,7,1,5,8]; else nn1=[1,5,8,2,4,7];  end  
    LapInfo = LapInfo1;
    sess_ind=(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    ind1=find(LapInfo(:,5)==kk & ~isnan(LapInfo(:,10)) & (LapInfo(:,9)~=0));
    idx=(~isnan(LapInfo(:,10)) & sess_ind); 
    LapInfo=LapInfo(idx,:); 
    spk2=squeeze(rates(idx,:,:));
    for bb=1:length(nn1)
    rw_spikes=[]; rw_spikes_mean=[];   
    id1=LapInfo(:,6)==nn1(bb);
    id2=LapInfo(:,9)==1;    
    rw_spikes=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=nanmean(rw_spikes,1);
    spk3{oo,pp}{gg}{bb,1}=[spk3{oo,pp}{gg}{bb,1}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    size1=size( (squeeze(rw_spikes_mean(1,:,:))));
    rw_spikes=[];rw_spikes_mean=[];
    id2=LapInfo(:,9)==2; 
    if sum(id1 & id2)>1      
    rw_spikes=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=nanmean(rw_spikes,1);
    elseif sum(id1 & id2)==1
    rw_spikes(1,:,:)=(squeeze(spk2(id1 & id2,:,:)));
    rw_spikes_mean=rw_spikes;
    else
    rw_spikes_mean(1,:,:)=ones(size1).*NaN;    
    end
    spk3{oo,pp}{gg}{bb,2}=[spk3{oo,pp}{gg}{bb,2}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    end
end
end

end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakeCellAssembliesRDMMatPFC25GoalForSleepCorrV2','spk3','-v7.3')

