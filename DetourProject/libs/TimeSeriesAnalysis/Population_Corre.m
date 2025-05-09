function pcos = Population_Corre(x1,x2)
%Population_Corre. The initial purpose of this function is to compute the population 
%                  correlation between two high dimensional random variables (random vectors).
%                  The current function actually compute the cosine
%                  value of angle between two vectors.

%                  
%Inputs            x1 & x2 two vectors with same dimension, the inner product is their vector dot product
%Output            pcorr  cosine value of the angle between these two vectors
%Note              if one of x1 or x2 is all Nans, it will return 0; if both of them are all
%                  Nans, it will return NaN
%   Yuchen Zhou 2021.Sep

% check Nan
x1 = squeeze(x1);
x2 = squeeze(x2);
if nansum(x1) == 0 && nansum(x2) == 0
    pcos = nan; return
end

if nansum(x1) == 0 || nansum(x2) == 0
    pcos = 0; return
end

pcos = nansum(x1.*x2)./(sqrt(nansum(x1.*x1))*sqrt(nansum(x2.*x2))); % inner product normalized by their length

end

