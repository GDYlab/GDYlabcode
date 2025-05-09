
function y = nanflip(x)
%NANFLIP flip nan to the end of the vector x
%   for example, if x = [nan 5 nan 7 9 nan nan 0]
%      y = nanflip(x) = [5 7 9 0 nan nan nan nan]

%% check input
if ~isvector(x)
    error('input must be a vector')
end


nanind = isnan(x);
if sum(nanind) == 0
   y = x; return 
end
    
xnan = x(nanind);
xnum = x(~nanind);

if isrow(x)
    y = [xnum,xnan];
else
    y = [xnum;xnan];
end


end

