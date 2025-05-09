function [groups,pvals] = AllGroupTtest2(d)
% this function get p values across all the groups 
% d must be a cell array with length over 1
% output are groups and pvals which can be used in sigstar

% check inputs
if ~iscell(d)
    warning('Input must be cell array with length over 1')
    groups = {};
    pvals = [];
    return
end
if length(d) < 2
    warning('Input must be cell array with length over 1')
    groups = {};
    pvals = [];
    return
end

ng = length(d);

groups = {};
pvals = [];
for ig = 1:ng-1
    for jg = ig+1:ng
         groups = cat(2,groups,[ig,jg]);
         [~,pvnow] = ttest2(d{ig},d{jg});
         pvals = cat(2,pvals,pvnow);
    end
end
    

end

