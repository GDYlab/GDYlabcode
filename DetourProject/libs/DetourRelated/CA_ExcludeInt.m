function [M1,M2,ClF,PlFields,nc,pyrid,clu2pyr,pyrintid] = CA_ExcludeInt(M1,M2,ClF,PlFields)
% this function remove interneurons from M1,M2,ClF, and PlFields
% nc is the totla number of pyr cells
% pyrid map pyr id to cluster id
% clu2pyr map cluster id to pyr id, if it's int, then it will give nan
% pyrintid map cluster id to pyr id, if it's int, the value is negative 
% for example 1,2,3,-1,-2,4 meaning cluster 1-6 are p1,p2,p3,i1,i2,p4
clu2pyr = nan(1,size(ClF,1)); 
pyrintid = nan(1,size(ClF,1)); 

%Exclude interneuron in analysis
m2temp = squeeze(M2(:,:,1));
pyrid = 1:size(m2temp,1);
intidall = 1:size(m2temp,1);

intid =   m2temp(:,21)==0;  % int neuron index
PlFields(intid,:,:,:) = [];
M2(intid,:,:) = [];
M1(intid,:,:,:) = [];
ClF(intid,:) = [];
nc = size(ClF,1);
pyrid(intid) = [];
intidall(~intid) = [];

for i = 1:length(pyrid)
    clu2pyr(pyrid(i)) = i;
    pyrintid(pyrid(i)) = i;
end

for i = 1:length(intidall)
    pyrintid(intidall(i)) = -i;
end

end

