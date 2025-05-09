function [tse,maskori] = GetEpochwithCrioverLen(t,idx,lencri)
% this function return the start and end time of epoches with a criteria
% and length over length criteria
% inputs  t is a vector of time stamps
%         idx is a vector of index showing whether or not the criteria is met or not
%         lencri is a length criteria for epoches
% output  tse, start and end time for epoches, rows are qualified epoches,
%         columns are start and end time
%         maskori as a index of qualified time stamps
lencriidx = lencri./nanmean(diff(t));

[B, N]   = RunLength(idx);
B(N < lencriidx) = false;
maskori = RunLength(B,N);
masklp = [false,berow(maskori),false];
timetemp = [nan,berow(t),nan];

ii1=strfind(masklp,[0 1])+1; % start index for all the epoches
ii2=strfind(masklp,[1 0]); % end index

if isempty(ii1)
    tse = [];
else
    tse = [becolumn(timetemp(ii1)),becolumn(timetemp(ii2))];
end

end