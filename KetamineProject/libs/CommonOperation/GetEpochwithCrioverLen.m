function [tse,maskori] = GetEpochwithCrioverLen(t,idx,lencri)
% GetEpochwithCrioverLen
% Return the start/end times of the epochs in which a logical criterion is
% continuously satisfied for at least lencri seconds.
%
% inputs:  t      - vector of time stamps (s)
%          idx    - logical vector, whether the criterion is met at each sample
%          lencri - minimum epoch duration (s)
% outputs: tse     - nEpoch x 2, [start end] time of each qualifying epoch
%          maskori - logical vector (same size as idx) flagging the samples
%                    that belong to a qualifying epoch
%
% Pure-MATLAB implementation (no compiled dependencies).

idx = logical(idx(:))';
lencriidx = lencri / nanmean(diff(t));   % minimum length in samples

% locate runs of consecutive true values
d      = diff([false, idx, false]);
starts = find(d == 1);
ends   = find(d == -1) - 1;

% keep only runs that are long enough
runlen = ends - starts + 1;
keep   = runlen >= lencriidx;
starts = starts(keep);
ends   = ends(keep);

% build outputs
maskori = false(size(idx));
for k = 1:numel(starts)
    maskori(starts(k):ends(k)) = true;
end
maskori = reshape(maskori,size(idx));

t = t(:);
if isempty(starts)
    tse = [];
else
    tse = [t(starts(:)), t(ends(:))];
end
end
