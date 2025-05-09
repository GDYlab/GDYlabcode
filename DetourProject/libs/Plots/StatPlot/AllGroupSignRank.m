function [groups,pvals] = AllGroupSignRank(d)
% this function get p values across all the groups 
% d must be a matrix with rows being realizations and columns being groups
% to compare
% output are groups and pvals which can be used in sigstar

% check inputs
if size(d,2) < 2
    warning('Input Must be matrix with rows being realizations and columns (>=2) being groups to compare')
    groups = {};
    pvals = [];
    return
end

ng = size(d,2);

groups = {};
pvals = [];
for ig = 1:ng-1
    for jg = ig+1:ng
        try
            pvnow = signrank(d(:,ig)-d(:,jg));
            groups = cat(2,groups,[ig,jg]);
            pvals = cat(2,pvals,pvnow);
        end
    end
end


end

