function PlF = addnoisetoplf(PlF)
% addnoisetoplf this function add small noise to the placemap in order to
% reduce non-decodeable moments
% PlF D1 cells, D2 spatial bins

for ic = 1:size(PlF,1)
    plftemp = PlF(ic,:);
    if max(plftemp) == 0
        PlF(ic,:) = ones(size(plftemp))*1e-8;
    else
        PlF(ic,:) = PlF(ic,:)+nanmean(PlF(ic,:))/1e5;
    end
end
end