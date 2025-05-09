% function [n,ia]=histcounts_row_kf(data,edge)
% edge (nx2)
% bin the hitogram according to each row: n(i)=sum( data>=edge(i,1) & data<=edge(i+1,2))
% ia: index of each data point in which bin
% ATTENTION: put all the equal value to the right side
function [n,ia]=histcounts_row_kf(data,edge)

n=zeros(1,size(edge,1));
ia=zeros(size(data));
for i=1:length(n)
   ind = data >= edge(i,1) & data < edge(i,2);
   n(i) = sum(ind); 
   ia(ind)=i;
end

