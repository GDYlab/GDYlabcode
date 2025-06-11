function [ReactStrength]= CellAssembliesCorrMatSleepFrames(TimeInterval1,Ind_cellID)

%dir is directory for analyis
%TimeInterval1 and TimeInterval2 are when coactivations are analysed
%Ind_cellID have logical indices for which cells are included
b=textread('clusters','%s');
b=b(Ind_cellID);

trackcltime=[];

for c=1:size(b,1)
       clG=cl2mat(b{c});              
       trackcltime{c}= [restrict_time(clG.featuredata(:,8), TimeInterval1(:,1),TimeInterval1(:,2))];
       clear clG
end

t=0.02; %timescale
tracktimebins=[];
for a=1:size(TimeInterval1,1)
tracktimebins=[tracktimebins, TimeInterval1(a,1):t:TimeInterval1(a,2)];
end

tracktimebins=sort(tracktimebins,'ascend');

trackcount=[];
for c=1:size(trackcltime,2)
    x=[];
    [x ind]=hist(trackcltime{c},tracktimebins);
    trackcount{c}=x';
end

j=0;trackplcounts2=[];
for c=1:size(trackcount,2)
    if sum(trackcount{c})>=0
        j=j+1;
        trackplcounts2(j,:)=trackcount{c};
        inds2(j)=c;
    end
    
end
ztrackplcounts=[];
for c=1:size(trackplcounts2,1)
    ztrackplcounts(c,:)=zscore(trackplcounts2(c,:));
   
end

trackcorrmatnew=(ztrackplcounts*ztrackplcounts')/size(ztrackplcounts,2);
% eigenvals=eig(trackcorrmatnew);
% 
% [~,eigenind]=sort(eigenvals,'descend');
% eigenvals=eigenvals(eigenind);
% eigmax=((1+sqrt((size(ztrackplcounts,1)/size(ztrackplcounts,2))))^2);
% nsigpcas=sum((eigenvals>=eigmax));
% pctracknew=pca(trackcorrmatnew);

ReactStrength.corrmat=trackcorrmatnew;
ReactStrength.corrcounts=ztrackplcounts;
ReactStrength.tbins=tracktimebins;
% ReactStrength.inds=inds2;
% ReactStrength.weightpc=[];
% for c=1:size(ReactStrength.sigpcs,2)
% weights=[];
% weights=kron(ReactStrength.sigpcs(:,c)',ReactStrength.sigpcs(:,c));
% weights(logical(eye(size(weights))))=0;
% ReactStrength.weightpc(c,:,:)=weights;
% end
% 
% ReactStrength.reactstrength=[];
% for c=1:size(ReactStrength.sigpcs,2)
%    x1=[];x2=[];x3=[]; 
% x1=zpostplcounts'*squeeze(ReactStrength.weightpc(c,:,:));
% 
% x2=x1.*(zpostplcounts');
% 
% x3=sum(x2,2)';
% 
% ReactStrength.reactstrength(c,:)=x3;
% end

clearvars -except ReactStrength
fclose('all');











    
    




