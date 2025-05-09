function v = nanweightmean(x,weight)
%nanweightmean compute the weighted mean of vector x based on weight.
%Ignore nan values in x, and define corresponding weight as zero

%check isvector and length

if ~isvector(x) || ~isvector(weight)
   error('input must be vector')
end

if length(x) ~= length(weight)
    error('input should have the same length')
end

weight(isnan(x)) = 0;
x(isnan(x)) = 0;
v = sum(x.*weight)./sum(weight);

end