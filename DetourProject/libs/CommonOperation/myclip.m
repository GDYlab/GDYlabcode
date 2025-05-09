function v = myclip(v,range)
% clip the value in a range

if v >= range(2)
   v =  range(2);
end

if v <= range(1)
   v =  range(1);
end

end