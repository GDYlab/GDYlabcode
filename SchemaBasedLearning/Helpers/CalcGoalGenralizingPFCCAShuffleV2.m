%% calculation

clear all; close all;
cd('/media/baburam/DataStorage2/Results/PlotFigures')
paths=importdata('HomePath');
paths1=importdata('Paths');
files=importdata('FilesPath');
Replayfile=importdata('ReplayFilePath');
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('AwakeCellAssembliesRDMMatPFC25GoalForSleepCorrV2','spk3')

for i=1:4
spk{i,1}=spk3{i,1};
end

for i=1:4
spk{i,2}=spk3{i,2};
end

for i=1:4
spk{i,3}=spk3{i,3};
end

sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
corrMat=[];
bins=[1:5:71];

for Days=1:4
    for Dir=1:3
rr=0;
     for ii=1:3       
      for kk=1:2 
       for jj=1:6 
         rr= rr+1;
            xx1=spk{Days,Dir}{ii}{jj,kk};
            corrMat{Days,Dir}{rr}=xx1;
        end
        end
end
end
end


corrMat1=[];

for Days=1:4
for Dir=1:3 
for ii=1:size(corrMat{Days,Dir},2)
    for jj=1:size(corrMat{Days,Dir},2)
        for kk=1:size(corrMat{Days,Dir}{1,ii},1)
        corrMat1(Days,Dir,ii,jj,kk)=(corr(corrMat{Days,Dir}{1,ii}(kk,:)',corrMat{Days,Dir}{1,jj}(kk,:)'));
        end
        
     end
end
end
end


GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakeCAIdentificationGoalV2.mat','corrMat1','-v7.3')


GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('AwakeCAIdentificationGoalV2.mat','corrMat1')

SpecificInfo=[];SessInfo=[];Outcome=[];NovelInfo=[];
outcomei=[ones(1,6),ones(1,6).*0,ones(1,6),ones(1,6).*0,ones(1,6),ones(1,6).*0];
outcomej=[ones(1,6),ones(1,6).*0,ones(1,6),ones(1,6).*0,ones(1,6),ones(1,6).*0];

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
         SpecificInfo(ii,jj)=-1;  
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

Pos=[];Neg=[];CPD1=[];SessIDIDX=[3,5,15];CPDS=[];PosS=[];NegS=[];
for Days=1:4
    for Dir=1:3   
    tic;       
    CorrMatInc=squeeze(corrMat1(Days,Dir,:,:,:));
    for nn=1:size(CorrMatInc,3)
    CorrMatIncNN=squeeze(CorrMatInc(:,:,nn));
   
  
    PopVec=CorrMatIncNN(logical(Across & SessInfo==-1 ));
    
    Group=[];
    Group(:,1)=Outcome(Across & SessInfo==-1 );
    Group(:,2)=SpecificInfo(Across & SessInfo==-1 );
    Group(:,3)=SessID(Across & SessInfo==-1 );  
   
    Response1=PopVec;
    models=[0 0 0; 1 0 0; 0 1 0; 0 0 1];     
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE= cell2mat(tbl(5,2)); 
          
    models=[0 0 0; 0 1 0; 0 0 1]; 
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(4,2)) ;      
    CPD1(Days,Dir,nn,1)=(SSE1-SSE)./SSE1; 
         
        
    models=[0 0 0; 1 0 0; 0 0 1];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(4,2)) ;        
    CPD1(Days,Dir,nn,2)=(SSE1-SSE)./SSE1; 
    
    models=[0 0 0; 1 0 0; 0 1 0];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,Group(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(4,2)) ;        
    CPD1(Days,Dir,nn,3)=(SSE1-SSE)./SSE1; 
    
      if all(all(CorrMatIncNN==0)); 
       continue; end
    CPDS1=[];PosS1=[];NegS1=[];
    for ww=1:1000
        
    x=1:size(Group,1);
    y=Group(randperm(numel(x)),:); 
 
    
    models=[0 0 0; 1 0 0; 0 1 0; 0 0 1];     
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE= cell2mat(tbl(5,2)); 
          
    models=[0 0 0; 0 1 0; 0 0 1]; 
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(4,2)) ;      
    CPDS1(ww,1)=(SSE1-SSE)./SSE1; 
         
        
    models=[0 0 0; 1 0 0; 0 0 1];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(4,2)) ;        
    CPDS1(ww,2)=(SSE1-SSE)./SSE1; 
    
    models=[0 0 0; 1 0 0; 0 1 0];  
    [pANCOVAN,tbl,stats,terms] = anovan(Response1,y(:,[1,2,3]),'model',models,...
               'varnames',{'Outcome','Specific','Sess'},'display','off'); 
    SSE1= cell2mat(tbl(4,2)) ;        
    CPDS1(ww,3)=(SSE1-SSE)./SSE1;  
             
    end
    CPDS(Days,Dir,nn,:,:)=(CPDS1);
   
    
    
    
     end
     display(strcat('DoneForDay',num2str(Days)))
     toc      
         
    end    
    
end

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('AwakePFCCAEachGoalActivationGeneralizeCPDShuff1000V4.mat','CPD1','CPDS','-v7.3')

GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load('AwakePFCCAEachGoalActivationGeneralizeCPDShuff1000V4.mat','CPD1','CPDS')


GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
load AwakeCAActivationsAround3Ripples3ZThPFC25msV21000

Sess={[1,3,4,5],[21,23,24,25],[6,8,9,10,11,13,14,15],[16,18,19,20]};
for Days=1:length(Sess)
    for jj=1:4
        for kk=1:4
            for nn=1:2
            RippleActivations{Days,nn}{jj,kk}=[];
            RippleActivationsSEM{Days,nn}{jj,kk}=[];
            RippleActivationsIndex{Days,nn}{jj,kk}=[];
            CPDInd{Days,1}{jj,1}=[];
            CPDVal{Days,1}{jj,nn}=[];
            sessIDs{Days,1}{jj,kk}=[];
            end
        end
    end
end
sd=3;
Fss = 1/2;
npoints = round(2.*sd.*Fss);
kernel = normpdf(linspace(-4*sd, 4*sd, 2*npoints+1), 0, sd);
kernel1 = kernel./sum(kernel);
% 
for Days=1:length(Sess)
for ii=Sess{Days}
      
    if all(cellfun(@isempty,Activations(ii,:))); continue; end
    for jj=1:3
        for kk=1:4
            xx1=(Activations{ii,jj}{kk});
            xx1=nanmean(xx1,1);
            xx1=squeeze(xx1(1,:,:));
            xx11=(SActivations{ii,jj}{kk});
            xx11=nanmean(xx11,1);
            xx11=squeeze(xx11(1,:,:));
            yy1=((Activations{ii,jj}{kk})-(SActivations{ii,jj}{kk}))./(SActivations{ii,jj}{kk});
            yy11=[];
            for ee=1:size(yy1,2)
                yy11(ee,:)=nansem(squeeze(yy1(:,ee,:)));
            end
            yy1=squeeze(yy1(1,:,:));
            xx2=[];xx22=[];xx=[];
            for cc=1:size(xx1,1)
              
                xx(cc,:)=xx1(cc,:);
                xx1(cc,:)=movmean(xx1(cc,:),[2,2],2);
                xx2(cc,1)=(nanmean(xx(cc,19:23))-nanmean(xx(cc,:)))./nanmean(xx(cc,:));
                
                xx(cc,:)=xx1(cc,:);
                xx11(cc,:)=movmean(xx11(cc,:),[2,2],2);
                xx22(cc,1)=(nanmean(xx(cc,19:23))-nanmean(xx(cc,:)))./nanmean(xx(cc,:));
            end
            RippleActivations{Days,1}{jj,kk}=[RippleActivations{Days,1}{jj,kk}; xx1];
            RippleActivationsSEM{Days,1}{jj,kk}=[RippleActivationsSEM{Days,1}{jj,kk}; yy11];
%             xx2=
            RippleActivationsIndex{Days,1}{jj,kk}=[RippleActivationsIndex{Days,1}{jj,kk}; xx2];
            
            RippleActivations{Days,2}{jj,kk}=[RippleActivations{Days,2}{jj,kk}; xx11];
%             xx2=
            RippleActivationsIndex{Days,2}{jj,kk}=[RippleActivationsIndex{Days,2}{jj,kk}; xx22];
            
            sessIDs{Days,1}{jj,kk}=[sessIDs{Days,1}{jj,kk}; ones(size(xx1,1),1).*ii];
            
        end
        xx=squeeze(CPD1(Days,jj,:,1));
        xxS=squeeze(CPDS(Days,jj,:,:,1));
        xxInd=[];
        for ss=1:size(xx,1)
        xxInd(ss,1)=sum(xx(ss)>=xxS(ss,:))./size(xxS,2);
        end
        CPDShuffles{Days,1}{jj,1}=xxS;
        CPD{Days,1}{jj,1}=xx;
        CPDInd{Days,1}{jj,1}=[xxInd];
    
    end
end
end
GenPath=strcat(paths1{2},'Results/PlotFigures/Database');
cd(GenPath)
save('PFCCARippleActivationsV2','RippleActivations','RippleActivationsSEM','RippleActivationsIndex','CPDInd','CPDVal','CPDShuffles','CPD','sessIDs');

