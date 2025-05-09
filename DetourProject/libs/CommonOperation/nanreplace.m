function v = nanreplace(v)
%replace nans in vin by the nearest non-value by position
nonnan = ~isnan(v);
vvalues = v(nonnan);
%bsxfun compute the distance from every nan to all non-nan.
%min then find which of the distance is the smallest
[~, mindistcol] = min(abs(bsxfun(@minus, 1:numel(v), find(nonnan)')));
%the row at which the minimum was picked is then the index of the non-nan value that is closest
v = mindistcol(vvalues);
end