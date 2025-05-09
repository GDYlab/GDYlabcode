function v = TakeUpTriElements(m)
% take uptriangle elements from a matrix (diag not included), form a column vector
[r,c] = size(m);

idx = true(r,c);
upidx =  triu(idx,1);
upidx = upidx(:);
v = m(:);
v = v(upidx);

end