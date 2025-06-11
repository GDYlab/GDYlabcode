function [ReactStrength]= CellAssembliesToFiringRates(TimeInterval2,ReactStrength1,clg,tbin)

%dir is directory for analyis
%TimeInterval1 and TimeInterval2 are when coactivations are analysed
%Ind_cellID have logical indices for which cells are included
% b=textread('clusters','%s');
% b=b(Ind_cellID);

postcltime=[];   

for c=1:size(clg,2)
       clG=clg{c};              
       postcltime{c}= [restrict_time(clG.featuredata(:,8), TimeInterval2(:,1),TimeInterval2(:,2))];
       clear clG
end

t=tbin; %timescale
posttimebins=[];
for a=1:size(TimeInterval2,1)
posttimebins=[posttimebins, TimeInterval2(a,1):t:TimeInterval2(a,2)];
end

postcount=[];
for c=1:size(postcltime,2)
    x=[];
    [x]=hist(postcltime{c},posttimebins);
    postcount{c}=x';
end
j=0;postplcounts2=[];
for c=1:size(postcount,2)    
        j=j+1;
        postplcounts2(j,:)=postcount{c};
        inds2(j)=c; 
    
end
zpostplcounts=[];
cellcounts=sum(logical(postplcounts2));
for c=1:size(postplcounts2,1)   
    zpostplcounts(c,:)=zscore(postplcounts2(c,:));    
end
ReactStrength=ReactStrength1;
weightpc=[];
for c=1:size(ReactStrength.sigpcs,2)
weights=[];
weights=kron(ReactStrength.sigpcs(:,c)',ReactStrength.sigpcs(:,c));
weights(logical(eye(size(weights))))=0;
weightpc(c,:,:)=weights;
end

ReactStrength.reactstrength=[];
for c=1:size(ReactStrength.sigpcs,2)
   x1=[];x2=[];x3=[]; 
x1=zpostplcounts'*squeeze(weightpc(c,:,:));

x2=x1.*(zpostplcounts');

x3=sum(x2,2)';

ReactStrength.reactstrength(c,:)=x3;

end
ReactStrength.t=posttimebins;
ReactStrength.cc=cellcounts;
clearvars -except ReactStrength
fclose('all');











    
    




