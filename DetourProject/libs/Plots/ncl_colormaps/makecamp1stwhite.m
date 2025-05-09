function cmap = makecamp1stwhite(cmap)
% make the first element in the colormap being completely white

if iscell(cmap)
    for ie = 1:length(cmap)
        cmap{ie}(1,:) = [1,1,1];
    end
    
else
    
    if ismatrix(cmap)
        cmap(1,:) = [1,1,1];
    end
end

end