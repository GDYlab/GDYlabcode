function idx = idxinrange_v2(x,range,epi)
% this function give the index of x in the range specified in range
% consider the numerical precison
if nargin < 3
   epi = 1e-12;  
end
idx1 = x-range(1) >= -epi;

idx2 = x-range(2) <= epi;

idx = idx1 & idx2;

end