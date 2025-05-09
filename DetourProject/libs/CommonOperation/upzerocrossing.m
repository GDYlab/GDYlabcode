function idx = upzerocrossing(x)
% this function find the idx of up 0 crossing for a given vector, the idx
% is the current value is below or equal 0, but the next value is larger than 0
x = berow(x);

negidx = x <= 0;

xtmp = x;
xtmp(xtmp == 0) = -1;
sgnx = sign(xtmp);

crossing = (-1).*sgnx(1:end-1) .*  sgnx(2:end);

crossing = cat(2,crossing,0);

negcrossing = crossing .* negidx;

idx = find(negcrossing == 1);


end