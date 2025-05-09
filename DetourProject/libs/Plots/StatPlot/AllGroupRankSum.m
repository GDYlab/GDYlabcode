function [groups,pvals] = AllGroupRankSum(d,xval)
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

if nargin < 2
    xval = [];
end

ng = length(d);

groups = {};
pvals = [];
for ig = 1:ng-1
    for jg = ig+1:ng
        if length(d{ig}(:)) <= 4 || length(d{jg}(:)) <= 4
            continue
        end
        
        if isempty(xval)
            groups = cat(2,groups,[ig,jg]);
        else
            groups = cat(2,groups,[xval(ig),xval(jg)]);
        end
        
        pvnow = ranksum(d{ig}(:),d{jg}(:));
        
        pvals = cat(2,pvals,pvnow);
    end
end


end

