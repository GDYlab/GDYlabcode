function ccat = CellCatinAllD(c1,c2,catd)
% this function concatenate all the elements in a cell array
% inputs : c1, c2 cell arrays to concatenate, their dimension must match
%          catd, dimension to concatenate, if set to 0, we will take (:)
%          and cat in 1
%          in each element in the cell array, we will do ccat{} = cat(catd,c1,{},c2{})


if ~isequal(size(c1), size(c2))
   error('input cell array must have the same size') 
end

if nargin < 3
   catd = 1; 
end

sizec = size(c1);
ccat = cell(sizec);
for ie = 1:numel(c1)
    if isempty(c1{ie})
        c1{ie} = [];
    end
    
    if isempty(c2{ie})
        c2{ie} = [];
    end
    
    if catd == 0
        ccat{ie} = cat(1,c1{ie}(:),c2{ie}(:));
    else
        ccat{ie} = cat(catd,c1{ie},c2{ie});
    end
end



end