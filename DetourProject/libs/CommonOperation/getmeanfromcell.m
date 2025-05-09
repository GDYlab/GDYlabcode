function mv = getmeanfromcell(d)
% this function get mean values from a cell array
mv = nan(size(d));
for iel = 1:numel(d)
   mv(iel) = nanmean(d{iel}); 
    
end

end