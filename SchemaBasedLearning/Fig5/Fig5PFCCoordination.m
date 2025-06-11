%% Fig 5a
% Change to appropriate directory using files
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
spk3=[];
load AwakeCAAssemblyGoalTuningPFC
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
sess={[6:10,11:15]}
fig1=figure;
for nn=1:length(sess)
     
  sessinc=[1:4];
for pp=sessinc  
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ll}{i,j}=[];
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
    spk3{gg}{bb,1}=[spk3{gg}{bb,1}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
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
    spk3{gg}{bb,2}=[spk3{gg}{bb,2}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    end
end
end
labels={'BR',...
       'BUr',...
       'AnR',...
       'AnUr',...
       'AdR',...
       'AdUr'};
rr=0;corrMat=[];
 for ii=1:3       
      for kk=1:2 
       for jj=1:6 
         rr= rr+1;
            corrMat{rr}=spk3{ii}{jj,kk};
        end
    end
end
corrMat1=[];
for ii=1:length(corrMat)
    for jj=1:length(corrMat)
        corrMat1(ii,jj)=nanmean(diag(corr(corrMat{1,ii},corrMat{1,jj},'rows','pairwise')));
    end
end
left=0.1+(pp-1).*0.20; bot=0.65-(nn-1).*0.3; width=0.17; height=0.15;
ax(nn,pp)=axes('position',[left bot width height])
mm=corrMat1;
mm(logical(eye(size(mm))))=NaN;
corrMat2=corrMat1;
corrMat2(logical(eye(size(mm))))=NaN;
rang=nanmax(nanmax(corrMat2))-nanmin(nanmin(corrMat2));
minx=nanmin(nanmin(corrMat2))+rang/5;
maxx=nanmax(nanmax(corrMat2))-rang/5;
imagesc(corrMat2,[0.1 0.3])
% colormap hot
cmap=getNCLColmap('GMT_no_green.rgb',21)
colormap(gca,cmap);
if pp==4
cb=colorbar;
set(cb,'Position',[left+width+0.01 bot 0.01 0.05])
cb.FontSize=3;
end
hold on;

hold on;
plot([0 36],[12.5 12.5],'linestyle','--','Marker','none','linewidth',1,'color','b');
plot([0 36],[24.5 24.5],'linestyle','--','Marker','none','linewidth',1,'color','b');
plot([12.5 12.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','b');
plot([24.5 24.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','b');

plot([0 36],[6.5 6.5],'linestyle','--','Marker','none','linewidth',1,'color','k');
plot([0 36],[18.5 18.5],'linestyle','--','Marker','none','linewidth',1,'color','k');
plot([0 36],[30.5 30.5],'linestyle','--','Marker','none','linewidth',1,'color','k');

plot([6.5 6.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','k');
plot([18.5 18.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','k');
plot([30.5 30.5],[0 36],'linestyle','--','Marker','none','linewidth',1,'color','k');
set(gca,'box','off','tickdir','out','fontsize',5)
set(gca,'xtick',[3.5:6:36],'xticklabel',labels,'fontsize',5,'xticklabelrotation',90);
if pp==1
set(gca,'ytick',[3.5:6:36],'yticklabel',labels,'fontsize',5);
else
set(gca,'ytick',[],'yticklabel','none','fontsize',5);    
end
end
end
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,3.5,3.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300','AwakeCellAssembliesReactivationSimilarityLrn25msGoalPFC');
print(fig1,'-painters','-depsc','-r300','AwakeCellAssembliesReactivationSimilarityLrn25msGoalPFC'); 
close all;

%% Fig 5b, h and i
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('PFCCARippleActivationsV2','RippleActivations','RippleActivationsIndex','CPDInd','CPDVal','CPDShuffles','CPD','sessIDs');


col1={rgb('Red'),rgb('Blue'),rgb('Black'),rgb('Grey')};
%% comaprison figures
ripphigh=[];ripplow=[];
sessinc={[1,3],[1,3],[1,3],[1,3]};
fig1=figure
for Days=1:4
    rr=0;
for jj=[3]
    rr=rr+1;
%     left=0.1+(Days-1).*0.15; bot=0.80-(jj-1).*0.15; width=0.12; height=0.12;
    ss1=0
    for kk=sessinc{Days}
        ss1=ss1+1;
        left=0.1+(rr-1).*0.2+(ss1-1).*0.12; bot=0.60-(Days-1).*0.11; width=0.06; height=0.08;        
        perind=CPDInd{Days,1}{jj,1};
        perind=perind(~(perind==0));
%         pos=CPDVal{Days,1}{jj,1}(~isnan(perind));neg=CPDVal{Days,1}{jj,2}(~isnan(perind));
%         ripp1=(pos);ripp2=neg;
%         perind1=(ripp1-ripp2)./ripp2; 
        Phigh = prctile(perind,50);
        Plow = prctile(perind,50);
        pind1=perind>=0.95;
        pind2=perind<0.95;  
        Ind1{Days,rr}=perind(perind>Phigh);Ind2{Days,rr}=perind(perind<Plow);
        ripp=RippleActivations{Days,1}{jj,kk};sripp=RippleActivations{Days,2}{jj,kk};
        ripphigh{Days,rr}=[(ripp(pind1,:)-sripp(pind1,:))]./0.025;
        ripplow{Days,rr}=[(ripp(pind2,:)-sripp(pind2,:))]./0.025;
    ax1(Days,rr,ss1)=axes('position',[left bot width height]);
    if size(ripphigh{Days,rr},1)==1; 
    plot(1:size(ripphigh{Days,rr},2),(ripphigh{Days,rr}),'color','r','linestyle','-','marker','none','linewidth',1);hold on; 
    else
    plot(1:size(ripphigh{Days,rr},2),nanmean(ripphigh{Days,rr}),'color','r','linestyle','-','marker','none','linewidth',1);hold on;
    pp=patch([[1:size(ripphigh{Days,rr},2)] fliplr([1:size(ripphigh{Days,rr},2)])], ...
    [nanmean(ripphigh{Days,rr})+nansem(ripphigh{Days,rr}) fliplr(nanmean(ripphigh{Days,rr})-nansem(ripphigh{Days,rr}))], 'r')
    pp.EdgeColor = 'none';
    alpha(0.1); hold on
    end
    plot(1:size(ripplow{Days,rr},2),nanmean(ripplow{Days,rr}),'color','b','linestyle','-','marker','none','linewidth',1);
    pp=patch([[1:size(ripplow{Days,rr},2)] fliplr([1:size(ripplow{Days,rr},2)])], ...
    [nanmean(ripplow{Days,rr})+nansem(ripplow{Days,rr}) fliplr(nanmean(ripplow{Days,rr})-nansem(ripplow{Days,rr}))], 'b')
    pp.EdgeColor = 'none';
    alpha(0.1); hold on
    plot([21 21],[0 nanmax(nanmax(nanmean(ripphigh{Days,rr})))],'--k'); 
  
    set(gca,'fontsize',5);
    if ss1==1
    ylabel('activation','fontsize',4);
    end
    yl=get(gca,'ylim');
    yl1=[nanmin(yl) nanmax(yl)];
    set(gca,'ytick',yl1,'yticklabel',[round(yl1(1),1), round(yl1(2),1)],'fontsize',5);
   
    set(gca,'xtick',[1 21 41],'xticklabel',[-0.5 0 0.5],'fontsize',5);
    if Days==4
    xlabel('peri ripple time (s)','fontsize',5);
    end
    xlim([0 42]);
    set(gca,'box','off','tickdir','out')
%     ylim([-0.001 0.007]);  
    end
  
end
end


%% Comparison with baseline increases
col1={rgb('Black'),rgb('Grey'),rgb('Red'),rgb('Blue')};
sessinc={[1,3],[1,3],[1,3],[1,3]};
 left=left+0.2;ax1=[];
rr=0;Ind1=[];Ind2=[];ripphigh=[];ripplow=[];
for jj=[3]
    rr=rr+1;
%    
    ss=0;  
    for kk=[1,3] 
         ss=ss+1;
        
        left=left+(ss-1).*0.25; bot=0.20-(rr-1).*0.10; width=0.15; height=0.15;        
        ax1(ss)=axes('position',[left bot width height]);
              
        for Days=1:4
        
        perind=CPDInd{Days,1}{jj,1};
        perind=perind(~(perind==0));

        Phigh = prctile(perind,50);
        Plow = prctile(perind,50);
        pind1=perind>=0.95;
        pind2=perind<0.95;  
        ripp=RippleActivations{Days,1}{jj,kk};sripp=RippleActivations{Days,2}{jj,kk};
        ripphigh{Days,ss}=[(ripp(pind1,:)-sripp(pind1,:))]./0.025;
        ripplow{Days,ss}=[(ripp(pind2,:)-sripp(pind2,:))]./0.025;
        
    
    if size(ripphigh{Days,ss},1)==1; 
    x1=(ripphigh{Days,ss}-nanmean(ripplow{Days,ss}));
    
    indices1{Days,ss}=x1; 
    plot(1:size(x1,2),(x1),'color',col1{Days},'linestyle','-','marker','none','linewidth',1);hold on; 
    else
    x1=((ripphigh{Days,ss})-nanmean(ripplow{Days,ss}));  
    indices1{Days,ss}=x1;
    plot(1:size(x1,2),nanmean(x1),'color',col1{Days},'linestyle','-','marker','none','linewidth',1);hold on;
    pp=patch([[1:size(x1,2)] fliplr([1:size(x1,2)])], ...
    [nanmean(x1)+nansem(x1) fliplr(nanmean(x1)-nansem(x1))],col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on ;

    end
 
    
    set(gca,'fontsize',5);
    ylabel('\Delta activation','fontsize',5);
          
    set(gca,'xtick',[1 21 41],'xticklabel',[-0.5 0 0.5],'fontsize',4);
    
    xlabel('peri ripple time (s)','fontsize',4);
   
    xlim([0 41]);
    set(gca,'box','off','tickdir','out')
        end
   
    end
end
siglabels={'*','**','***'};
col1={rgb('Black'),rgb('Grey'),rgb('Red'),rgb('Blue')};
P=[];
    for ee=1:2
    for nn=1:4    
    yl(ee,:) = (get(ax1(ee), 'Ylim')) ;
    

     Pval=[];sig=[];
     for ss=1:size(ripphigh{nn,ee},2)   
         [h Pval(ss) ]=ttest2(ripphigh{nn,ee}(:,ss),ripplow{nn,ee}(:,ss));
         if Pval(ss) <0.05 & Pval(ss) >= 0.01 ;
         sig{ss}=siglabels{1};
         elseif Pval(ss) <0.01 & Pval(ss) >= 0.001; 
         sig{ss}=siglabels{2};  
         elseif Pval(ss) <0.001 & Pval(ss) >= 0;
         sig{ss}=siglabels{3}; 
         else
         sig{ss}=NaN;     
         end
     if ~isnan(sig{ss})  
     t=text('Parent', ax1(ee), 'Position', [ss yl(ee,end)+(nn-1).*0.005], 'String', sig{ss});
     t.FontSize=4; t.Color=col1{nn};
     t.Rotation=270;
     end
     end
     if ee==1; offset=0.003; else offset=0.001; end
     [ m mmax]=nanmax(nanmean(indices1{nn,ee}));
     plot(ax1(ee),[round(mmax,1) round(mmax,1)],[yl(ee,1) yl(ee,end)],'linestyle','--','color',col1{nn});
     t=text('Parent', ax1(ee), 'Position', [22 yl(ee,end)-(nn-1).*offset], 'String', strcat('Peak =',num2str(round(mmax,1).*0.025-0.525), ' s'));
     t.FontSize=4; t.Color=col1{nn};
     
    end 
if ee==1; offset=0.001; else offset=0.001; end    
x = [21 31 31 21];
y = [yl(ee,1) yl(ee,1) yl(ee,end)-offset yl(ee,end)-offset];
pp=patch(x,y,'m');  
pp.EdgeColor = 'none';
pp.Parent=ax1(ee);
alpha(pp,0.1); hold on

x = [11 21 21 11];
y = [yl(ee,1) yl(ee,1) yl(ee,end)-offset yl(ee,end)-offset];
pp=patch(x,y,'k');  
pp.EdgeColor = 'none';
pp.Parent=ax1(ee);
alpha(pp,0.1); hold on
    
end
P=[];
for nn=1:2
set(ax1(nn),'ylim',[yl(nn,1) yl(nn,end)]);   
plot(ax1(nn),[21 21],[yl(nn,1) yl(nn,end)],'-k'); 
sig=[];
Groups=[1,2;1,3;1,4;2,3;2,4;3,4];Pval1=[];
for ee=1:size(Groups,1)
    x1=(nanmean(indices1{Groups(ee,1),nn}(:,21:31))-nanmean(indices1{Groups(ee,1),nn}(:,10:20)));
    x2=(nanmean(indices1{Groups(ee,2),nn}(:,21:31))-nanmean(indices1{Groups(ee,2),nn}(:,10:20)));
    x11=((indices1{Groups(ee,1),nn}(:,21:31))-(indices1{Groups(ee,1),nn}(:,10:20)));
    x22=((indices1{Groups(ee,2),nn}(:,21:31))-(indices1{Groups(ee,2),nn}(:,10:20)));
    
    
    
    [Pval1(ee) h stats]=ranksum(x1(:),x2(:),'method','approximate');
    P{1,nn}{ee,1}=stats.ranksum; P{1,nn}{ee,2}=Pval1(ee);P{1,nn}{ee,3}=stats.zval;
    P{1,nn}{ee,4}=[length(x11(:)), length(x22(:))];  P{1,nn}{ee,5}=[nanmean(x1(:)), nanmean(x2(:))];
    P{1,nn}{ee,6}=[nansem(x1(:)), nansem(x2(:))]; 
       if nn==1; offset=0.001; else offset=0.001; end
       if Pval1(ee) <0.05 & Pval1(ee) >= 0.01 ;
         sig{ee}=siglabels{1};
         elseif Pval1(ee) <0.01 & Pval1(ee) >= 0.001; 
         sig{ee}=siglabels{2};  
         elseif Pval1(ee) <0.001 & Pval1(ee) >= 0;
         sig{ee}=siglabels{3}; 
         else
         sig{ee}='+';     
         end
    t=text('Parent', ax1(nn), 'Position', [2 yl(nn,end)-(ee-1).*offset], 'String', strcat(sig{ee}));
    t.FontSize=4; 

    
end
end


%% calculate percentile of CPDInd cdf

perind1{1,1}=[CPDInd{1,1}{1,1}; CPDInd{2,1}{1,1};CPDInd{3,1}{1,1};CPDInd{4,1}{1,1}];
perind1{1,2}=[CPDInd{1,1}{2,1};CPDInd{2,1}{2,1}];
perind1{1,3}=[CPDInd{1,1}{3,1};CPDInd{2,1}{3,1}];
perind1{1,4}=[CPDInd{3,1}{2,1};CPDInd{4,1}{2,1}];
perind1{1,5}=[CPDInd{3,1}{3,1};CPDInd{4,1}{3,1}];
% perind1{1,4}=[CPDInd{3,1}{3,1}];
% perind1{1,5}=[CPDInd{4,1}{3,1}];
left=left-0.2; bot=bot+0.25; width=0.1; height=0.1; 
ax1(5)=axes('position',[left bot width height]);
col1={rgb('Cyan'),rgb('Black'),rgb('Grey'),rgb('Red'),rgb('Blue')};
Ind1=[];Tot=[];
for ii=1:5
    perind=perind1{1,ii};
    perind=perind(perind~=0);
        
            
        Ind1{ii}=sum(perind>=0.95);
        if ii==5
        Tot{ii}=82;
        else
        Tot{ii}=length(perind);
        end
h2(ii)=bar(ii,Ind1{ii}./Tot{ii}.*100);hold on;
h2(ii).FaceColor='none';
h2(ii).EdgeColor=col1{ii};
h2(ii).BarWidth=0.5;
end
set(gca,'xtick',[1 2 3 4 5],'xticklabel',{'Ba-1','Ctrl','Ba-2','Asn','Asd'},'fontsize',5,'xticklabelrotation',45);
set(gca,'fontsize',4,'tickdir','out','box','off');
xlabel('Sessions','fontsize',5);
ylabel('proportion of assemblies','fontsize',5);
Ind11{1}=Ind1{1};
Ind11{2}=Ind1{2};
Ind11{3}=Ind1{3};
Tot11{1}=Tot{1};
Tot11{2}=Tot{2};
Tot11{3}=Tot{3};

Groups=combnk([1:5],2)
Groups=[1,2;1,3;1,4;1,5;2,3;2,4;2,5;3,4;3,5;4,5];zstat1=[];
xpos=[1,2.5;1,4.5;2.5,4.5];
P=[];
% Groups=combnk([1:5],2)
for ii=1:size(Groups,1)
     [zstat1(ii)] = myztest(Ind1{Groups(ii,1)},Ind1{Groups(ii,2)},...
     Tot{Groups(ii,1)},Tot{Groups(ii,2)});
     p_one(ii) = normcdf(zstat1(ii));
     if p_one(ii)>0.5; p_one(ii)=1-p_one(ii); end
     if p_one(ii)<0.05
         hh=sigstar(Groups(ii,:),p_one(ii),0)
     end
     P{ii,1}{1}=Groups(ii,:); P{ii,2}=zstat1(ii); P{ii,3}=[(Ind1{Groups(ii,1)}), (Ind1{Groups(ii,2)}), ...
         Tot{Groups(ii,1)}, Tot{Groups(ii,2)}, (Ind1{Groups(ii,1)})./Tot{Groups(ii,1)}, (Ind1{Groups(ii,2)})./Tot{Groups(ii,2)}]; 
     P{ii,4}=p_one(ii);
end


 
set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,3.5,3.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',('PFCCAActivationRipplesNormAssdCompared')); 
print(fig1,'-painters','-depsc','-r300',('PFCCAActivationRipplesNormAssdCompared')); 
close(fig1)

%% Fig 4c-g
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('AwakePFCCAEachGoalActivationGeneralizeCPDShuff1000V4.mat','CPD1','CPDS')
CPDInt=squeeze(CPD1(3,3,:,1));

load('PFCCARippleActivationsV2');
Act=(RippleActivations{3,1}{3,3});
SAct=(RippleActivations{3,2}{3,3});
NAct=(Act-SAct)./0.025;
ActSEM=(RippleActivationsSEM{3,1}{3,3});
IndCPD1=CPDInd{3,1}{3,1};
CPDShuff1=CPDShuffles{3,1}{3,1};
CPDVal1=CPD{3,1}{3,1};

load AwakeCAAssemblyGoalTuningPFC

load('AwakeCellAssemblyActivationInSleepPFC25.mat', ...
        'React');
glmstats=[];
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
sess={[6,8,9,10,11,13,14,15]};spk3=[];
ind_CA=[];sigpcs=[];
for nn=1:length(sess)
   
    
for pp=3
for ll=1:4
for i=1:6
    for j=1:2
        spk3{ll}{i,j}=[];
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
if any(ii==[14])
ind_CA1=ones(size(AssemblyTuning{ii,pp},2),1); 
else
ind_CA1=zeros(size(AssemblyTuning{ii,pp},2),1);
end

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
    spk3{gg}{bb,1}=[spk3{gg}{bb,1}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
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
    spk3{gg}{bb,2}=[spk3{gg}{bb,2}; conv2((squeeze(rw_spikes_mean(1,:,:))),kernel1,'same')];
    end
end
if any(ii==[14])
sigpcs=[sigpcs; React{ii,pp}.sigpcs];
end
ind_CA=[ind_CA; ind_CA1];
end

end
end

col1={rgb('Red'),rgb('Blue'),rgb('Green'),rgb('Black')};
col2={rgb('Salmon'),rgb('LightBlue'),rgb('LightGreen'),rgb('Black')};
fig1=figure;
CPDInt=CPDInt(~isnan(CPDInt));
IndCPD1=IndCPD1(~isnan(CPDInt));
CPDVal1=CPDVal1(~isnan(CPDInt));
CPDShuff1=CPDShuff1(~isnan(CPDInt),:);
CPDInt1=IndCPD1(ind_CA==1);
IndCPD=IndCPD1(ind_CA==1,:);
CPDVal=CPDVal1(ind_CA==1,:);
CPDShuff=(CPDShuff1(ind_CA==1,:));
[mn, ids]=sort(CPDInt1,'descend');
ids=ids([1,2,3,4,6,11]);
sigpcs1=sigpcs(:,ids);
IndCPD=IndCPD(ids,:);
CPDVal=CPDVal(ids,:);
CPDShuff=(CPDShuff(ids,:));
weightids=zeros(size(sigpcs1,1),1);
rr=0;
for cc=1:size(sigpcs1,2)
    rr=rr+100;
    zz=zscore(sigpcs1(:,cc));
    idds=find(zz>1.5);
    [mm iddsort]=sort(zz(idds),'descend');
    sortedid=idds(iddsort);
    for ss=1:length(idds)
    weightids(sortedid(ss))=1000-rr-ss;
    end
end

CPDInt1=CPDInt1(ids);

[mm idweights]=sort(weightids,'descend')
sigpcs1=sigpcs1(idweights,:);
ss=0;
for nn=1:length(CPDInt1)
ss=ss+1;
left=0.05+(ss-1).*0.06; bot=0.8; width=0.04; height=0.05; 
ax1(1,ss)=axes('position',[left bot width height]);



for cc=1:(size(sigpcs1,1))
zz=zscore((sigpcs1(:,nn)'));    
h(cc)=bar(cc,zz(cc)'); hold on

if zz(cc)>1.5 & sigpcs1(cc,nn)>0
    h(cc).EdgeColor='r'; 
else
    h(cc).EdgeColor=rgb('Gray'); 
end
h(cc).FaceColor='none'; 
end
view(90,90);
set(gca,'xtick',[0 (size(sigpcs1,1))],'xticklabel',[0 (size(sigpcs1,1))],...
    'fontsize',5);
if nn==1
xlabel('Cell ID','fontsize',5);
end

ylabel('Weights','fontsize',5)
%% RSM plots
rr=0;corrMat=[];
     for ii=1:3       
      for kk=1:2 
       for jj=1:6 
         rr= rr+1;
            x1=spk3{ii}{jj,kk}(ind_CA==1,:);
            corrMat{rr}=x1(ids(nn),:);
        end
    end
end
corrMat1=[];
for ii=1:length(corrMat)
    for jj=1:length(corrMat)
       corrMat1(ii,jj)=((corr(corrMat{1,ii}',corrMat{1,jj}','rows','pairwise')));
    end
end
left=0.05+(ss-1).*0.06; bot=bot-0.07; width=0.04; height=0.04; 
ax1(2,ss)=axes('position',[left bot width height])
mm=corrMat1;
mm(logical(eye(size(mm))))=NaN;
corrMat2=corrMat1;
corrMat2(logical(eye(size(mm))))=NaN;
% corrMat2=(corrMat2)-nanmean(nanmean(corrMat2));
imagesc(corrMat2,[nanmin(nanmin(corrMat2))+0.4 nanmax(nanmax(corrMat2))-0.6])
% colormap jet
cmap=getNCLColmap('GMT_no_green.rgb',21);
colormap(gca,cmap);

hold on;
plot([0 36],[12.5 12.5],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([0 36],[24.5 24.5],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([12.5 12.5],[0 36],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([24.5 24.5],[0 36],'linestyle','-','Marker','none','linewidth',0.5,'color','k');

plot([0 36],[6.5 6.5],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([0 36],[18.5 18.5],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([0 36],[30.5 30.5],'linestyle','-','Marker','none','linewidth',0.5,'color','k');

plot([6.5 6.5],[0 36],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([18.5 18.5],[0 36],'linestyle','-','Marker','none','linewidth',0.5,'color','k');
plot([30.5 30.5],[0 36],'linestyle','-','Marker','none','linewidth',0.5,'color','k');

set(gca,'box','on','tickdir','out','fontsize',4)
if nn==1
set(gca,'xtick',[1:6:31]+2.5,'xticklabel',{'BR','BUr','AnR','AnUr','AdR','AdUr'},...
    'ytick',[1:6:31]+2.5,'yticklabel',{'BR','BUr','AnR','AnUr','AdR','AdUr'},...
    'fontsize',4,'xticklabelrotation',45);
else
set(gca,'xtick',[1:6:31]+2.5,'xticklabel',{'BR','BUr','AnR','AnUr','AdR','AdUr'},...
    'ytick',[],...
    'fontsize',4,'xticklabelrotation',45);
end
% ylabel('Cell ID','fontsize',5)
% xlabel('Weights','fontsize',5)
% set(gca,'xtick',[],'fontsize',3,'xticklabelrotation',90);
% set(gca,'ytick',[],'fontsize',3);

left=0.05+(ss-1).*0.06; bot=bot-0.07; width=0.04; height=0.04; 
ax1(3,ss)=axes('position',[left bot width height])
x1=CPDShuff(nn,:);
bins=min(x1):0.001:max(x1);
binedges=(bins(1:end-1)+bins(2:end))/2;
[hh, bb]=histc(x1,bins);
h(1)=bar(bins,hh); hold on
h(1).FaceColor='k';
h(1).BarWidth=0.8;
h(2)=bar(CPDVal(nn),nanmax(hh)+5); hold on
h(2).BarWidth=0.001;
h(2).FaceColor='r';
h(2).EdgeColor='none';
xlim([-0.001 max([CPDVal(nn)+0.01 nanmax(bins)+0.01])]); hold on
annotation('textbox', [left-0.001 bot+height+0.01 0, 0], 'string', ...
    strcat('CPDInd= ',num2str(round(IndCPD(nn),2)))...
    ,'fontsize',5);
set(gca,'fontsize',5);
if nn==1
ylabel('CPD','fontsize',5)
end

xlabel('Counts','fontsize',5)


sp2=[];sp1=[];
for i=1:3
    for j=1:6
       x1=spk3{i}{j,1}(ind_CA==1,:);
       x2=spk3{i}{j,2}(ind_CA==1,:);
            sp1=[sp1; x1(ids(nn),:)];
            sp2=[sp2; x2(ids(nn),:)];
     
    end
end
left=0.05+(ss-1).*0.06; bot=bot-0.07; width=width; height=0.04; 
ax1(4,ss)=axes('position',[left bot width height]);
imagesc(sp1,[0 0.3]); hold on

cmap=getNCLColmap('MPL_Greys.rgb',128);
colormap(gca,cmap);
plot([0 25],[6.5 6.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 25],[12.5 12.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([4.5 4.5],[0 18],'linestyle','--','Marker','none','linewidth',1,'color','k');
set(gca,'ydir','normal');
set(gca,'xtick',[0 4.5 25],'xticklabel',[-1 0 5],'fontsize',5);
set(gca,'ytick',[2.5 8.5 14.5],'yticklabel',{'Ba','Assn','Assd'},'fontsize',4);
if nn==1
ylabel('Trials','fontsize',5)
else
set(gca,'ytick',[]);
end
xlabel('Time at goal (s)','fontsize',5)

% ylim([0.5 12.5]);
left=0.05+(ss-1).*0.06; bot=bot-0.07; width=width; height=0.04;
ax1(5,ss)=axes('position',[left bot width height]);
    
imagesc(sp2,[0 0.3]);  hold on
cmap=getNCLColmap('MPL_Greys.rgb',128);
colormap(gca,cmap);
plot([0 25],[6.5 6.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([0 25],[12.5 12.5],'linestyle','-','Marker','none','linewidth',1,'color','r');
plot([4.5 4.5],[0 18],'linestyle','--','Marker','none','linewidth',1,'color','k');
set(gca,'ydir','normal');
set(gca,'xtick',[0 4.5 25],'xticklabel',[-1 0 5],'fontsize',5);
set(gca,'ytick',[2.5 8.5 14.5],'yticklabel',{'Ba','Assn','Assd'},'fontsize',4);
if nn==1
ylabel('Trials','fontsize',5)
else
set(gca,'ytick',[]);
end
xlabel('Time at goal (s)','fontsize',5)

% ylim([0.5 12.5]);

% ylim([0.5 12.5]);
left=0.05+(ss-1).*0.06; bot=bot-0.07; width=width; height=0.04;
ax1(6,ss)=axes('position',[left bot width height]);
sp2=[];sp1=[];
for i=1:3
    sp2=[];sp1=[];
    for j=1:6
       x1=spk3{i}{j,1}(ind_CA==1,:);
       x2=spk3{i}{j,2}(ind_CA==1,:);
            sp1=[sp1; x1(ids(nn),:)];
            sp2=[sp2; x2(ids(nn),:)];
     
    end
ddt_mean=movmean(nanmean(sp1),[2 2],2);
ddt_sem=movmean(nansem(sp1),[2 2],2);
hb1(i)=plot(1:size(sp1,2),movmean(nanmean(sp1),[2 2],2),'linestyle','-','marker','none','color',col1{i}); hold on
pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{i})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on 
yl=get(gca,'ylim');
% plot(1:size(sp2,2),movmean(nanmean(sp2),[2 2],2),'linestyle','-','marker','none','color',col2{i});
end
plot([4.5 4.5],[nanmin(yl)-0.1 nanmax(yl)+0.1],'linestyle','--','Marker','none','linewidth',1,'color','k');
set(gca,'xtick',[0 4.5 25],'xticklabel',[-1 0 5],'fontsize',5);
if nn==1
ylabel('Activation','fontsize',5)
end
xlabel('Time at goal (s)','fontsize',5)
h8=legend ([hb1(1) hb1(2) hb1(3)], 'Ba', 'Assn','Assd','Location','Northwest')
legend('boxoff')
set(h8,'FontSize',4);
h8.ItemTokenSize = [5,5];


left=0.05+(ss-1).*0.06; bot=bot-0.07; width=width; height=0.04;
ax1(7,ss)=axes('position',[left bot width height]);
sp2=[];sp1=[];
for i=1:3
    sp2=[];sp1=[];
    for j=1:6
       x1=spk3{i}{j,1}(ind_CA==1,:);
       x2=spk3{i}{j,2}(ind_CA==1,:);
            sp1=[sp1; x1(ids(nn),:)];
            sp2=[sp2; x2(ids(nn),:)];
     
    end
ddt_mean=movmean(nanmean(sp2),[2 2],2);
ddt_sem=movmean(nansem(sp2),[2 2],2);    
% plot(1:size(sp1,2),movmean(nanmean(sp1),[2 2],2),'linestyle','-','marker','none','color',col1{i}); 
hb1(i)=plot(1:size(sp2,2),movmean(nanmean(sp2),[2 2],2),'linestyle','-','marker','none','color',col1{i}); hold on
pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{i})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on 
yl=get(gca,'ylim');
% plot(1:size(sp2,2),movmean(nanmean(sp2),[2 2],2),'linestyle','-','marker','none','color',col2{i});
end
plot([4.5 4.5],[nanmin(yl)-0.1 nanmax(yl)+0.1],'linestyle','--','Marker','none','linewidth',1,'color','k');
set(gca,'xtick',[0 4.5 25],'xticklabel',[-1 0 5],'fontsize',5);
if nn==1
ylabel('Activation','fontsize',5)
end
xlabel('Time at goal (s)','fontsize',5)

left=0.05+(ss-1).*0.06; bot=bot-0.07; width=width; height=0.04;
ax1(8,ss)=axes('position',[left bot width height]);
x1=NAct(ind_CA==1,:); 
x2=SAct(ind_CA==1,:);
x3=ActSEM(ind_CA==1,:);
ddt_mean=conv2(x1(ids(nn),:),kernel1,'same');
ddt_sem=conv2(x3(ids(nn),:),kernel1,'same');
plot(1:size(Act,2),ddt_mean,'linestyle','-','marker','none','color','r'); hold on
pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], 'r')
    pp.EdgeColor = 'none';
    alpha(0.1); hold on 
yl=get(gca,'ylim');
% plot(1:size(Act,2),movmean(x2(ids(nn),:),[1 3],2),'linestyle','-','marker','none','color','b');
plot([20.5 20.5],[nanmin(yl)-0.1 nanmax(yl)+0.1],'linestyle','--','Marker','none','linewidth',1,'color','k');
set(gca,'xtick',[1 20.5 41],'xticklabel',[-0.5 0 0 5],'fontsize',5);
if nn==1
ylabel('Norm Activation','fontsize',5)
else
set(gca,'ytick',[]);
end
xlim([0 42])
xlabel('Per-ripple time (s)','fontsize',5)

end

linkaxes([ax1(8,1) ax1(8,2) ax1(8,3) ax1(8,4) ax1(8,5) ...
    ax1(8,6) ...
    ],'xy');

linkaxes([ax1([6:7],1)],'xy');
linkaxes([ax1([6:7],2)],'xy');
linkaxes([ax1([6:7],3)],'xy');
linkaxes([ax1([6:7],4)],'xy');
linkaxes([ax1([6:7],5)],'xy');
linkaxes([ax1([6:7],6)],'xy');

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,8,8])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',('ExamplePFCCAActivationCPDRipples14')); 
print(fig1,'-painters','-depsc','-r300',('ExamplePFCCAActivationCPDRipples14')); 
close(fig1)


%% Fig 5j and k
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CueCellParticipationPFCCAsRawData.mat','Prop1','Prop2','Prop1S','Prop2S');

col1={rgb('Blue'),rgb('Grey'),rgb('Red'),rgb('Grey'),rgb('Cyan')};
labels={'Day 1','Day 4','Day 2&3','Day 5'};
P=[];
fig1=figure
for i=1:2
    left=0.1; bot=0.80-(i-1).*0.15; width=0.12; height=0.12;
    ax1(i,1)=axes('position',[left bot width height]);
    Means1=cellfun(@nanmean,(Prop1(:,i)));
    Sem1=cellfun(@nansem,(Prop1(:,i)));
    Means2=cellfun(@nanmean,(Prop2(:,i)));
    Sem2=cellfun(@nansem,(Prop2(:,i))); 
    
    Means1S=cellfun(@nanmean,(Prop1S(:,i)));
    Sem1S=cellfun(@nansem,(Prop1S(:,i)));
    Means2S=cellfun(@nanmean,(Prop2S(:,i)));
    Sem2S=cellfun(@nansem,(Prop2S(:,i))); 
P{i}{1,1}=[(Means1), Means1S, Means2, Means2S]; P{i}{1,2}=[(Sem1), Sem1S, Sem2, Sem2S];
    
    
h1(i)=bar([1:4],[Means1]);hold on
h2(i)=bar([1:4]+0.5,[Means1S]);hold on
h1(i).BarWidth=0.3;
h2(i).BarWidth=0.3;
h1(i).FaceColor='none';
h2(i).FaceColor='none';
h1(i).EdgeColor=col1{1};
h2(i).EdgeColor=col1{4};
errorbar([1:4],Means1,Sem1,'linestyle','none','color',col1{1},'capsize',3)
errorbar([1:4]+0.5,Means1S,Sem1S,'linestyle','none','color',col1{4},'capsize',3)
set(gca,'fontsize',5);
set(gca,'xtick',[1:4],'xticklabel',labels,'fontsize',4','xticklabelrotation',45,'box','off','tickdir','out');
ylabel('proportion','fontsize',5)


Groups=[1,1.5;2,2.5;3,3.5;4,4.5];
pval=[];hv1=[];
for tt=1:size(Groups,1)
    [hv1(tt) pval(tt) ci stats]=ttest2(Prop1{Groups(tt,1),i},Prop1S{Groups(tt,1),i});
    P{i}{1,3}(tt,:)=Groups(tt,:); P{i}{1,4}(tt)= pval(tt); P{i}{1,5}(tt)= stats.df; P{i}{1,6}(tt)= stats.tstat;
    if pval(tt) < 0.05
    hh1=sigstar10(Groups(tt,:),pval(tt),0);
    end
end


left=0.3; bot=0.80-(i-1).*0.15; width=0.12; height=0.12;
ax1(i,2)=axes('position',[left bot width height]);
h3(i)=bar([1:4],[Means2]);hold on
h4(i)=bar([1:4]+0.5,[Means2S]);hold on
h3(i).BarWidth=0.3;
h4(i).BarWidth=0.3;
h3(i).FaceColor='none';
h4(i).FaceColor='none';
h3(i).EdgeColor=col1{1};
h4(i).EdgeColor=col1{4};
errorbar([1:4],Means2,Sem2,'linestyle','none','color',col1{1},'capsize',2)
errorbar([1:4]+0.5,Means2S,Sem2S,'linestyle','none','color',col1{4},'capsize',2)
set(gca,'fontsize',5);
set(gca,'xtick',[1:4],'xticklabel',labels,'fontsize',4','xticklabelrotation',45,'box','off','tickdir','out');
ylabel('proportion','fontsize',5)
Groups=[1,1.5;2,2.5;3,3.5;4,4.5];
pval=[];hv1=[];
for tt=1:size(Groups,1)
  [hv1(tt) pval(tt) ci stats]=ttest2(Prop2{Groups(tt,1),i},Prop2S{Groups(tt,1),i});
   P{i}{1,7}(tt,:)=Groups(tt,:); P{i}{1,8}(tt)= pval(tt); P{i}{1,9}(tt)= stats.df; P{i}{1,10}(tt)= stats.tstat;
  if pval(tt) < 0.05
    hh2=sigstar10(Groups(tt,:),pval(tt),0);
    end
end

ylim([0 1])
end

linkaxes(ax1([1,2],[1,2]));

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,3.5,3.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',('CAPosGoalSimilarityCodingPropFlavorCellsSig')); 
print(fig1,'-painters','-depsc','-r300',('CAPosGoalSimilarityCodingPropFlavorCellsSig')); 
close(fig1)


%% Fig 5l and m
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('PFCFlavorID&RippleActivations100','RippleActivations','FID');

sd=7;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);

col1={rgb('Red'),rgb('Blue'),rgb('Black'),rgb('Grey')};
%% comaprison figures
ripphigh=[];ripplow=[];
sessinc={[1,3],[1,3],[1,3],[1,3]};
fig1=figure
for Days=1:4
    rr=0;
for jj=[3]
    rr=rr+1;
%     left=0.1+(Days-1).*0.15; bot=0.80-(jj-1).*0.15; width=0.12; height=0.12;
    ss1=0
    for kk=sessinc{Days}
        ss1=ss1+1;
        left=0.1+(rr-1).*0.2+(ss1-1).*0.12; bot=0.60-(Days-1).*0.11; width=0.06; height=0.08;        
        perind=FID{Days,jj};
        perind=perind(~isnan(perind));
        Phigh = prctile(perind,50);
        Plow = prctile(perind,50);
        pind1=perind>=1;
        pind2=perind<1;  
        Ind1{Days,rr}=perind(perind>Phigh);Ind2{Days,rr}=perind(perind<Plow);
        ripp=RippleActivations{Days,1}{1,kk};
        ripphigh{Days,rr}=[(ripp(pind1,:))];
        ripplow{Days,rr}=[(ripp(pind2,:))];
    ax1(Days,rr,ss1)=axes('position',[left bot width height]);
    if size(ripphigh{Days,rr},1)==1; 
    plot(1:size(ripphigh{Days,rr},2),(ripphigh{Days,rr}),'color','r','linestyle','-','marker','none','linewidth',1);hold on; 
    else
    plot(1:size(ripphigh{Days,rr},2),nanmean(ripphigh{Days,rr}),'color','r','linestyle','-','marker','none','linewidth',1);hold on;
    pp=patch([[1:size(ripphigh{Days,rr},2)] fliplr([1:size(ripphigh{Days,rr},2)])], ...
    [nanmean(ripphigh{Days,rr})+nansem(ripphigh{Days,rr}) fliplr(nanmean(ripphigh{Days,rr})-nansem(ripphigh{Days,rr}))], 'r')
    pp.EdgeColor = 'none';
    alpha(0.1); hold on
    end
    plot(1:size(ripplow{Days,rr},2),nanmean(ripplow{Days,rr}),'color','b','linestyle','-','marker','none','linewidth',1);
    pp=patch([[1:size(ripplow{Days,rr},2)] fliplr([1:size(ripplow{Days,rr},2)])], ...
    [nanmean(ripplow{Days,rr})+nansem(ripplow{Days,rr}) fliplr(nanmean(ripplow{Days,rr})-nansem(ripplow{Days,rr}))], 'b')
    pp.EdgeColor = 'none';
    alpha(0.1); hold on
    plot([101 101],[0 nanmax(nanmax(nanmean(ripphigh{Days,rr})))],'--k'); 

    set(gca,'fontsize',5);
    ylabel('activation','fontsize',5);
    yl=get(gca,'ylim');
    yl1=[nanmin(yl) nanmax(yl)];
    set(gca,'ytick',yl1,'yticklabel',[round(yl1(1),3), round(yl1(2),3)],'fontsize',5);

    set(gca,'xtick',[1 101 201],'xticklabel',[-1 0 1],'fontsize',5);
    if Days==4
    xlabel('peri ripple time (s)','fontsize',5);
    end
    xlim([0 201]);
    set(gca,'box','off','tickdir','out')

    end
  
end
end


%% Comparison with baseline increases
col1={rgb('Black'),rgb('Grey'),rgb('Red'),rgb('Blue')};
sessinc={[1,3],[1,3],[1,3],[1,3]};
 left=left+0.2;ax1=[];
rr=0;Ind1=[];Ind2=[];ripphigh=[];ripplow=[];
for jj=[3]
    rr=rr+1;
%    
    ss=0;  
    for kk=[1,3] 
         ss=ss+1;
        
        left=left+(ss-1).*0.2; bot=0.20-(rr-1).*0.10; width=0.15; height=0.15;        
        ax1(ss)=axes('position',[left bot width height]);
              
        for Days=1:4
        
        perind=FID{Days,jj};
        perind=perind(~isnan(perind));
        pind1=perind>=1;
        pind2=perind<1;  
        ripp=RippleActivations{Days,1}{1,kk};
        ripphigh{Days,ss}=[ripp(pind1,:)];
        ripplow{Days,ss}=[ripp(pind2,:)];
        
    
    if size(ripphigh{Days,ss},1)==1; 
    x1=conv2(ripphigh{Days,ss}-nanmean(ripplow{Days,ss}),kernel1,'same'); 
    
    indices1{Days,ss}=x1; 
    plot(1:size(x1,2),(x1),'color',col1{Days},'linestyle','-','marker','none','linewidth',1);hold on; 
    else
    x1=conv2(ripphigh{Days,ss}-nanmean(ripplow{Days,ss}),kernel1,'same');  
    x1=(x1);
    indices1{Days,ss}=x1;
    plot(1:size(x1,2),nanmean(x1),'color',col1{Days},'linestyle','-','marker','none','linewidth',1);hold on;
    pp=patch([[1:size(x1,2)] fliplr([1:size(x1,2)])], ...
    [nanmean(x1)+nansem(x1) fliplr(nanmean(x1)-nansem(x1))],col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on ;

    end
%     plot(1:size(ripplow,2),nanmean(ripplow),'color','b','linestyle','-','marker','none','linewidth',2);
    plot([1 201],[0 0],'--k');
    
    set(gca,'fontsize',5);
    ylabel('\Delta activation','fontsize',5);
          
    set(gca,'xtick',[1 101 201],'xticklabel',[-1 0 1],'fontsize',4);
    
    xlabel('peri ripple time (s)','fontsize',4);
   
    xlim([0 201]);
    set(gca,'box','off','tickdir','out')
    end
    
    end
end
siglabels={'*','**','***'};
col1={rgb('Black'),rgb('Grey'),rgb('Red'),rgb('Blue')};

    for ee=1:2
    for nn=1:4    
    yl(ee,:) = (get(ax1(ee), 'Ylim')) ;
    

     Pval=[];sig=[];
     for ss=1:size(ripphigh{nn,ee},2)   
         [h Pval(ss) ]=ttest2(ripphigh{nn,ee}(:,ss),ripplow{nn,ee}(:,ss),'Tail','Right');
         if Pval(ss) <0.05 & Pval(ss) >= 0.01 ;
         sig{ss}=siglabels{1};
         elseif Pval(ss) <0.01 & Pval(ss) >= 0.001; 
         sig{ss}=siglabels{2};  
         elseif Pval(ss) <0.001 & Pval(ss) >= 0;
         sig{ss}=siglabels{3}; 
         else
         sig{ss}=NaN;     
         end
     if ~isnan(sig{ss})  
     t=text('Parent', ax1(ee), 'Position', [ss yl(ee,end)+(nn-1).*0.001], 'String', sig{ss});
     t.FontSize=4; t.Color=col1{nn};
     t.Rotation=270;
     end
     end
     if ee==1; offset=0.003; else offset=0.003; end
     [ m mmax]=nanmax(nanmean(indices1{nn,ee}));
     plot(ax1(ee),[round(mmax,2) round(mmax,2)],[yl(ee,1) yl(ee,end)],'linestyle','--','color',col1{nn});
     t=text('Parent', ax1(ee), 'Position', [105 yl(ee,end)-(nn-1).*offset], 'String', strcat('Peak =',num2str(round(mmax,2).*0.01-1.01), ' s'));
     t.FontSize=4; t.Color=col1{nn};
     t.Rotation=270;
     
    end 
if ee==1; offset=0.001; else offset=0.001; end    
x = [101 126 126 101];
y = [yl(ee,1) yl(ee,1) yl(ee,end)-offset yl(ee,end)-offset];
pp=patch(x,y,'m');  
pp.EdgeColor = 'none';
pp.Parent=ax1(ee);
alpha(pp,0.1); hold on

x = [76 101 101 76];
y = [yl(ee,1) yl(ee,1) yl(ee,end)-offset yl(ee,end)-offset];
pp=patch(x,y,'k');  
pp.EdgeColor = 'none';
pp.Parent=ax1(ee);
alpha(pp,0.1); hold on
    
    end
P=[];
for nn=1:2
set(ax1(nn),'ylim',[yl(nn,1) yl(nn,end)]);   
plot(ax1(nn),[101 101],[yl(nn,1) yl(nn,end)],'-k'); 
sig=[];
Groups=[1,2;1,3;1,4;2,3;2,4;3,4];Pval1=[];
for ee=1:size(Groups,1)
    x1=(nanmean(indices1{Groups(ee,1),nn}(:,101:126))-nanmean(indices1{Groups(ee,1),nn}(:,75:100)));
    x2=(nanmean(indices1{Groups(ee,2),nn}(:,101:126))-nanmean(indices1{Groups(ee,2),nn}(:,75:100)));
    x11=((indices1{Groups(ee,1),nn}(:,101:126))-(indices1{Groups(ee,1),nn}(:,75:100)));
    x22=((indices1{Groups(ee,2),nn}(:,101:126))-(indices1{Groups(ee,2),nn}(:,75:100)));
  [Pval1(ee) h stats]=ranksum(x1(:),x2(:),'method','approximate');
    P{1,nn}{ee,1}=stats.ranksum; P{1,nn}{ee,2}=Pval1(ee);P{1,nn}{ee,3}=stats.zval;
    P{1,nn}{ee,4}=[length(x11(:)), length(x22(:))];  P{1,nn}{ee,5}=[nanmean(x1(:)), nanmean(x2(:))];
    P{1,nn}{ee,6}=[nansem(x11(:)), nansem(x22(:))]; 
       if nn==1; offset=0.001; else offset=0.003; end
       if Pval1(ee) <0.05 & Pval1(ee) >= 0.01 ;
         sig{ee}=siglabels{1};
         elseif Pval1(ee) <0.01 & Pval1(ee) >= 0.001; 
         sig{ee}=siglabels{2};  
         elseif Pval1(ee) <0.001 & Pval1(ee) >= 0;
         sig{ee}=siglabels{3}; 
         else
         sig{ee}='+';     
         end
    t=text('Parent', ax1(nn), 'Position', [51 yl(nn,end)-(ee-1).*offset], 'String', strcat(sig{ee}));
    t.FontSize=4;  
    
end
end

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.2,0.2,3.5,3.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',('FlavorCellsActivationRipplesNormAssdCompared')); 
print(fig1,'-painters','-depsc','-r300',('FlavorCellsActivationRipplesNormAssdCompared')); 
close(fig1)


%% Fig 5n
clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('GoalProbConsolidationEffectOnSleepFramesPFCNewVsOldRipp.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:6
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,8,9,10]};
for Days=1
for sess=1:4
    for jj=1:5
        if jj==1
     x=(([(Probability1{Days,sess}{1,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,2))]'))';
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
        else
      x=(([(Probability1{Days,sess}{3,jj-1}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{3,jj-1}(:,2))]'))';
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
        end
end
end

end

Sess={[11,13,14,15]};
for Days=1
for sess=5:8
    for jj=1
     if sess < 16   
     x=(([(Probability1{Days,sess}{1,jj}(:,2))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,1))]'))';
     indices{1,1}{6}=[indices{1,1}{6}; ((x-y))'];
     else
     x=(([(Probability1{Days,sess}{1,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,2))]'))';
     indices{1,1}{6}=[indices{1,1}{6}; ((x-y))'];  
     end
end
end

end
Sess={[11,13,14,15]};
for Days=1
for sess=5:8
    for jj=1:4
     x=(([(Probability1{Days,sess}{3,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{3,jj}(:,2))]'))';
     indices1{1,1}{jj}=[indices1{1,1}{jj}; ((x-y))'];
    
     
end
end

end

Sess={[21,23,24,25]};
for Days=1
for sess=9:12
    for jj=1
     if sess < 9
     x=(([(Probability1{Days,sess}{1,jj}(:,2))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,1))]'))';
     indices1{1,1}{5}=[indices1{1,1}{5}; ((x-y))'];
     else
     x=(([(Probability1{Days,sess}{1,jj}(:,1))]'))'; 
     y=(([(Probability1{Days,sess}{1,jj}(:,2))]'))';
     indices1{1,1}{5}=[indices1{1,1}{5}; ((x-y))'];  
     end
end
end

end


%% Combined consolidation Merged Day2 and 3

col1={rgb('Red'),...
    rgb('Blue'),rgb('Green'),rgb('LightBlue')};
fig1=figure(1)
left=0.05; bottom=0.75; width=0.20; height=0.2;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'Pre-Ba1','Pre-Nocue','Pre-Ba2'};
labels2={'Pre-Ba','Pre-Assn','Pre-Assd'};
for Days=1
    
    ddt1=[];ddt2=[];ddt3=[];
    for jj=1:6
        if jj==1
    ddt1=[ddt1, movmean([indices{1,Days}{jj}; ones(size(indices1{1,Days}{jj})).*NaN],[1,1],2)];
        else
    ddt1=[ddt1, movmean([indices{1,Days}{jj}; indices1{1,Days}{jj-1}],[1,1],2)];
        end
    end
   
    ddt_mean=nanmean(ddt1);
    ddt_sem=nansem(ddt1);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
    if Days==1
    plot([0 100],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([30.5 30.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([40.5 40.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([50.5 50.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([60.5 60.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([70.5 70.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([80.5 80.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([90.5 90.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([100.5 100.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',5,'tickdir','out',...
        'xtick',[5:10:64],'xticklabel',{'Ctrl','Pre-Ba','Pre-Assn','Pre-Assd','Post-Assd','Two Days Later'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
    ylabel('Frame prob (Changed-Unchanged)','fontsize',5)
    xlim([0, 65])
    ylim([-0.08 0.08]);
end


pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};

col1={rgb('Black'),rgb('Cyan'),rgb('Gray'),rgb('Blue'),rgb('Green'),rgb('Red')};
pval=[]; P=[]; 
left=0.05; bottom=0.55; width=0.15; height=0.10;
ax1=axes('Position',[left, bottom, width, height])  
  for jj=1:6
    ddt1=[];ddt2=[];ddt3=[];    

    
       for Days=1 
        if jj==1
    ddt_mean{Days}{jj}=[((indices{1,Days}{jj}(:)')')];
        else
        ddt_mean{Days}{jj}=[((indices{1,Days}{jj}(:)')'); ((indices1{1,Days}{jj-1}(:)')')];
        end
    [h1 stats]=ecdf(ddt_mean{Days}{jj});hold on
    plot(stats,h1,'color',col1{jj});
    plot([nanmedian(ddt_mean{Days}{jj}) nanmedian(ddt_mean{Days}{jj})],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{jj});
    end
%     [pval(jj) h ]=signrank(ddt_mean{1}{jj},ddt_mean{2}{jj});
%         if pval(jj) < 0.05 & pval(jj) > 0.01
%           text(-0.01,0.8,siglabels{1},'fontsize',4,'color','k');
%         elseif  pval(jj) < 0.01 & pval(jj) > 0.001
%           text(-0.01,0.8,siglabels{2},'fontsize',4,'color','k');
%         elseif  pval(jj) < 0.001 
%           text(-0.01,0.8,siglabels{3},'fontsize',4,'color','k');
%         end
    
    
    
    
       [pval1(jj) h stats]=signrank(ddt_mean{1}{jj});  
       P{1,1}{jj,1}=1; P{1,1}{jj,2}=pval1(jj);P{1,1}{jj,3}=stats.zval;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj});  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj});
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}); 
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.1,0.95-(jj-1).*0.1,siglabels{1},'fontsize',6,'color',col1{jj});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.1,0.95-(jj-1).*0.1,siglabels{2},'fontsize',6,'color',col1{jj});
        elseif  pval1(jj) < 0.001 
          text(-0.1,0.95-(jj-1).*0.1,siglabels{3},'fontsize',6,'color',col1{jj});
        end
    
%          [pval2(jj) h ]=signrank(ddt_mean{2}{jj});
%          if pval2(jj) < 0.05 & pval2(jj) > 0.01
%           text(-0.01,0.9,siglabels{1},'fontsize',6,'color',col1{2});
%         elseif  pval2(jj) < 0.01 & pval2(jj) > 0.001
%           text(-0.01,0.9,siglabels{2},'fontsize',6,'color',col1{2});
%         elseif  pval2(jj) < 0.001 
%           text(-0.01,0.9,siglabels{3},'fontsize',6,'color',col1{2});
%         end    
      
        
    set(gca,'box','off','tickdir','out','fontsize',5);
    if jj==3
    xlabel('Frame prob (Changed-Unchanged)','fontsize',5)
    end
    if Days==1
    ylabel('cdf','fontsize',5);    
    end
end
  


set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,5,3.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('GoalProbConsolidationEffectsOnSleepFramesNewVsOldPFCCombRippMerged'));
print(fig1,'-painters','-depsc','-r300',strcat('GoalProbConsolidationEffectsOnSleepFramesNewVsOldPFCCombRippMerged')); 
close all;



%% Fig 5 o

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('SleepFramesCorrMat10binsAllSleep20MultiHPCPFC');
load('IndicesPlaceArmsHPCPFC')
load('IndicesFlavorHPCPFC')
load('IndicesOdorHPCPFC')

sd=2;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);gg=0;
Coact=[];MembersZcontOld=[];
rr=0;
for ii=[6:15,21:25]
GenPath=strcat(paths,files{ii});
cd(GenPath{1})
 
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
ind=a12(:,1)==1 & a12(:,2)==1 & (elecreg(:,2)==1 | elecreg(:,2)==2);
ind_all=elecreg;
ind_all=elecreg(ind_HPC | ind_PFC,2);
if sum(ind_HPC) < 10 | sum(ind_PFC)<10; continue; end 
  
    tt=cell2mat(cellfun(@(x)(~isempty(x)),Indices_Place(ii,:),'UniformOutput',false));
    tt1=tt(tt~=0);
    if isempty(tt1); continue; end
    gg=0;rr=rr+1;
    if any(ii==[1:5,21:25,13,18,8]) ii==6 | ii==15; dd=0; else dd=1;end
    
   
    MembersIndices1=[];comm1=[];ll=0;
    MembersIndices2=[];comm2=[];
    MembersIndices3=[];comm3=[];
    MembersIndices4=[];comm4=[];mm1=0;
     for jj=[1,2,length(tt1)-dd,length(tt1)]   
      ll=ll+1; 
     MembersIndices1{ll}=(squeeze(Indices_Place{ii,jj}.Pvalues1(:,:,1,1)))<0.05;
     comm1{ll}=~(any(MembersIndices1{ll}(:,[2,3,5])')' & any(MembersIndices1{ll}(:,[1,4,6])')');
     
     MembersIndices2{ll}=(squeeze(Indices_Place{ii,jj}.Pvalues1(:,:,2,1)))<0.05;
     comm2{ll}=~(any(MembersIndices2{ll}(:,[2,3,5])')' & any(MembersIndices2{ll}(:,[1,4,6])')');
     
     MembersIndices3{ll}=(squeeze(Indices_Flavor{ii,jj}.Pvalues(:,:,1)))<0.05;
     comm3{ll}=~(any(MembersIndices3{ll}(:,[2,3,5])')' & any(MembersIndices3{ll}(:,[1,4,6])')');
     
     MembersIndices4{ll}=(squeeze(Indices_Odor{ii,jj}.Pvalues(:,:,1)))<0.05;
     comm4{ll}=~(any(MembersIndices4{ll}(:,[2,3,5])')' & any(MembersIndices4{ll}(:,[1,4,6])')');
     end
     
     comb=[1,3];
     if ii < 11; nn=[2,3,5];
     elseif ii > 10; nn=[1,4,6];
     end
     id11=[];id22=[];id1=[];id2=[];id33=[];id44=[];
     id3=[];id4=[];id5=[];id6=[];id7=[];id8=[];gg=0;
     for jj=1:size(comb,2);
     gg=gg+1; mm1=0;     
     for aa=[1,2,length(tt1)-dd,length(tt1)]  
     mm1=mm1+1;    
     zassemblycont=CorrMat{ii,aa};     
     for bb=1:length(nn)
         id11{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==1);
         id22{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==2);
 
         id1{1,jj}{bb}=(id11{jj}{bb});
         id1{2,jj}{bb}=id22{jj}{bb};
                 
      end
     comb1=[1 2];
     for ww=1:size(comb1,1)
      
     for cc=1:length(zassemblycont)   
     xx=[];    
     for bb=1:length(nn)
     indices=[];    
     for s=1:length(id1{comb1(ww,1),jj}{bb})
     for s1=1:length(id1{comb1(ww,2),jj}{bb})
     if id1{comb1(ww,1),jj}{bb}(s)== id1{comb1(ww,2),jj}{bb}(s1) ; continue; end   
     indices=[indices; (zassemblycont{cc}(id1{comb1(ww,1),jj}{bb}(s),id1{comb1(ww,2),jj}{bb}(s1)))];
     end
     end
%      if size(zassemblycont{cc},2) > 10; zz=1:10; else zz=1:size(zassemblycont{cc},2); end
%      indices=nanmean(indices')';
     indices=indices(:);
     xx=[xx; nanmean(indices)];      
     end  
     Coact{1,rr}{gg,mm1}(cc,:)=nanmean(xx);
     end
     end
    
     end
     end
         
     if ii < 11; nn=[1,4,6];
     elseif ii > 10; nn=[2,3,5];
     end
     id11=[];id22=[];id1=[];id2=[];id33=[];id44=[];
     id3=[];id4=[];id5=[];id6=[];id7=[];id8=[];gg=0;
     for jj=1:size(comb,2);
     gg=gg+1; mm1=0;     
     for aa=[1,2,length(tt1)-dd,length(tt1)]  
     mm1=mm1+1;    
     zassemblycont=CorrMat{ii,aa};     
     for bb=1:length(nn)
         id11{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==1);
         id22{jj}{bb}=find((MembersIndices3{comb(jj)}(:,nn(bb)) ) & comm3{comb(jj)} & ind_all==2);
 
         id1{1,jj}{bb}=id11{jj}{bb};
         id1{2,jj}{bb}=id22{jj}{bb};
                 
        end
     comb1=[1 2];
     for ww=1:size(comb1,1)
      
     for cc=1:length(zassemblycont)   
     xx=[];    
     for bb=1:length(nn)
     indices=[];    
     for s=1:length(id1{comb1(ww,1),jj}{bb})
     for s1=1:length(id1{comb1(ww,2),jj}{bb})
     if id1{comb1(ww,1),jj}{bb}(s)== id1{comb1(ww,2),jj}{bb}(s1) ; continue; end   
     indices=[indices; (zassemblycont{cc}(id1{comb1(ww,1),jj}{bb}(s),id1{comb1(ww,2),jj}{bb}(s1)))];
     end
     end
%      if size(zassemblycont{cc},2) > 10; zz=1:10; else zz=1:size(zassemblycont{cc},2); end
%      indices=nanmean(indices')';
     indices=indices(:);
     xx=[xx; nanmean(indices)];      
     end  
     Coact{2,rr}{gg,mm1}(cc,:)=nanmean(xx);
     end
     end
    
     end
     end
        
end


GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('CueCueCoactivityConsolidationEffectOnSleepFramesHPCPFCNewVsOld.mat','Coact')



clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('CueCueCoactivityConsolidationEffectOnSleepFramesHPCPFCNewVsOld.mat');

%% Calculating the indices
indices=[];indices1=[];indices2=[];
for i=1:2
    for k=1:2
    for j=1:6
        indices{i,k}{j}=[];
        indices1{i,k}{j}=[];
     
    end
    end
end
Sess={[6,9,10]};
for Days=1
for sess=1:3
    for jj=1:5
    if jj==1    
    x=(([(Coact{2,sess}{1,jj})]'))'; 
    y=(([(Coact{1,sess}{1,jj})]'))';
    x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
    else
    x=(([(Coact{1,sess}{2,jj-1})]'))'; 
    y=(([(Coact{2,sess}{2,jj-1})]'))';
    x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,1}{jj}=[indices{1,1}{jj}; ((x-y))'];
    end
     
end
end

end

Sess={[11,14,15]};
for Days=1
for sess=4:6
    for jj=1
     if sess < 7
     x=(([(Coact{2,sess}{1,jj})]'))'; 
     y=(([(Coact{1,sess}{1,jj})]'))';
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,1}{6}=[indices{1,1}{6}; ((x-y))'];
     else
     x=(([(Coact{1,sess}{1,jj})]'))'; 
     y=(([(Coact{2,sess}{1,jj})]'))';
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices{1,1}{6}=[indices{1,1}{6}; ((x-y))'];
     end
end
end

end


Sess={[11,14,15]};
for Days=1
for sess=[4:6]
    for jj=1:4
    x=(([(Coact{1,sess}{2,jj})]'))'; 
    y=(([(Coact{2,sess}{2,jj})]'))';
    x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
    indices1{1,1}{jj}=[indices1{1,1}{jj}; ((x-y))'];
     
end
end

end

Sess={[21,24,25]};
for Days=1
for sess=[6:8]
    for jj=1
     if sess < 7
     x=(([(Coact{2,sess}{1,jj})]'))'; 
     y=(([(Coact{1,sess}{1,jj})]'))';
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices1{1,1}{5}=[indices1{1,1}{5}; ((x-y))'];
     else
     x=(([(Coact{1,sess}{1,jj})]'))'; 
     y=(([(Coact{2,sess}{1,jj})]'))';
     x(:,all(isnan(x)))=0; y(:,all(isnan(y)))=0;
     indices1{1,1}{5}=[indices1{1,1}{5}; ((x-y))'];
     end
end
end

end



%% Merged

col1={rgb('Red'),...
    rgb('Blue'),rgb('Green'),rgb('LightBlue')};
fig1=figure(1)
left=0.05; bottom=0.75; width=0.20; height=0.2;
ax1(Days)=axes('Position',[left, bottom, width, height])
pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
labels1={'Pre-Ba1','Pre-Nocue','Pre-Ba2'};
labels2={'Pre-Ba','Pre-Assn','Pre-Assd'};
for Days=1
    
     ddt1=[];ddt2=[];ddt3=[];
    for jj=1:6
        if jj==1
    ddt1=[ddt1, movmean([indices{1,Days}{jj}; ones(size(indices1{1,Days}{jj})).*NaN],[2,2],2)];
        else
    ddt1=[ddt1, movmean([indices{1,Days}{jj}; indices1{1,Days}{jj-1}],[2,2],2)];
        end
    end
   
    ddt_mean=movmean(nanmean(ddt1),[1,2],2);
    ddt_sem=movmean(nansem(ddt1),[1,2],2);
    hh1(Days)=plot([1:length(ddt_mean)],ddt_mean,'linestyle','-','marker','none','linewidth',0.5,'color',col1{Days}); hold on;    
    pp=patch([[1:length(ddt_sem)] fliplr([1:length(ddt_sem)])], ...
       [ddt_sem+ddt_mean fliplr(ddt_mean-ddt_sem)], col1{Days})
    pp.EdgeColor = 'none';
    alpha(0.1); hold on  
    yl(Days,:)=get(gca,'ylim');
     if Days==1
    plot([0 100],[0 0],'--k');
    plot([10.5 10.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([20.5 20.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([30.5 30.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([40.5 40.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    plot([50.5 50.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([60.5 60.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([70.5 70.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([80.5 80.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([90.5 90.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
%     plot([100.5 100.5],[nanmin(yl(:,1)) nanmax(yl(:,2))],'--k');
    end
    set(gca,'box','off','fontsize',5,'tickdir','out',...
        'xtick',[5:10:64],'xticklabel',{'Ctrl','Pre-Ba','Pre-Assn','Pre-Assd','Post-Assd','Two Days Later'},'xticklabelrotation',45)
    xlabel('sleeps','fontsize',5);
      xlim([0, 65])
    ylabel('Coactivity (Changed-Unchanged)','fontsize',5)
end


pval=[]; cons=[];nn=0; pvalInd=[]; consInd=[];
xpos=[0,4,8,12];ddt_mean=[];siglabels={'*','**','***'};

col1={rgb('Black'),rgb('Cyan'),rgb('Gray'),rgb('Blue'),rgb('Green'),rgb('Red')};
pval=[]; P=[]; 
left=0.05; bottom=0.55; width=0.15; height=0.10;
ax1=axes('Position',[left, bottom, width, height])  
  for jj=1:6
    ddt1=[];ddt2=[];ddt3=[];    

    
       for Days=1 
        if jj==1
    ddt_mean{Days}{jj}=[((indices{1,Days}{jj}(:)')')];
        else
        ddt_mean{Days}{jj}=[((indices{1,Days}{jj}(:)')'); ((indices1{1,Days}{jj-1}(:)')')];
        end
    [h1 stats]=ecdf(ddt_mean{Days}{jj});hold on
    plot(stats,h1,'color',col1{jj});
    plot([nanmedian(ddt_mean{Days}{jj}) nanmedian(ddt_mean{Days}{jj})],[0 0.5],'linestyle','--','marker','none','linewidth',1,'color',col1{jj});
    end
%     [pval(jj) h ]=signrank(ddt_mean{1}{jj},ddt_mean{2}{jj});
%         if pval(jj) < 0.05 & pval(jj) > 0.01
%           text(-0.01,0.8,siglabels{1},'fontsize',4,'color','k');
%         elseif  pval(jj) < 0.01 & pval(jj) > 0.001
%           text(-0.01,0.8,siglabels{2},'fontsize',4,'color','k');
%         elseif  pval(jj) < 0.001 
%           text(-0.01,0.8,siglabels{3},'fontsize',4,'color','k');
%         end
    
    
    
    
       [pval1(jj) h stats]=signrank(ddt_mean{1}{jj});  
       P{1,1}{jj,1}=1; P{1,1}{jj,2}=pval1(jj);P{1,1}{jj,3}=stats.zval;
    P{1,1}{jj,4}=length(ddt_mean{1}{jj});  P{1,1}{jj,5}=nanmean(ddt_mean{1}{jj});
    P{1,1}{jj,6}=nansem(ddt_mean{1}{jj}); 
        if pval1(jj) < 0.05 & pval1(jj) > 0.01
          text(-0.1,0.95-(jj-1).*0.1,siglabels{1},'fontsize',6,'color',col1{jj});
        elseif  pval1(jj) < 0.01 & pval1(jj) > 0.001
          text(-0.1,0.95-(jj-1).*0.1,siglabels{2},'fontsize',6,'color',col1{jj});
        elseif  pval1(jj) < 0.001 
          text(-0.1,0.95-(jj-1).*0.1,siglabels{3},'fontsize',6,'color',col1{jj});
        end
    
%          [pval2(jj) h ]=signrank(ddt_mean{2}{jj});
%          if pval2(jj) < 0.05 & pval2(jj) > 0.01
%           text(-0.01,0.9,siglabels{1},'fontsize',6,'color',col1{2});
%         elseif  pval2(jj) < 0.01 & pval2(jj) > 0.001
%           text(-0.01,0.9,siglabels{2},'fontsize',6,'color',col1{2});
%         elseif  pval2(jj) < 0.001 
%           text(-0.01,0.9,siglabels{3},'fontsize',6,'color',col1{2});
%         end    
      
        
    set(gca,'box','off','tickdir','out','fontsize',5);
    if jj==3
    xlabel('Frame prob (Changed-Unchanged)','fontsize',5)
    end
    if Days==1
    ylabel('cdf','fontsize',5);    
    end
end
  

set(fig1,'paperunits','inches')
set(fig1,'papertype','usletter')
set(gcf,'paperposition',[0.1,0.1,5,3.5])
GenPath=strcat(paths1{2},'Results/PlotFigures/Fig5/Figures');
cd(GenPath);
print(fig1,'-djpeg','-r300',strcat('CueCueCoactivityConsolidationEffectsOnSleepCombMerged'));
print(fig1,'-painters','-depsc','-r300',strcat('CueCueCoactivityConsolidationEffectsOnSleepCombMerged')); 
close all;