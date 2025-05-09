function idx = downzerocrossing(x)
% this function find the idx of down 0 crossing for a given vector, the idx
% is the current value is larger or equal 0, but the next value is smaller than 0
x = berow(x);

posidx = x >= 0;

xtmp = x;
xtmp(xtmp == 0) = 1;
sgnx = sign(xtmp);

crossing = (-1).*sgnx(1:end-1) .*  sgnx(2:end);

crossing = cat(2,crossing,0);

negcrossing = crossing .* posidx;

idx = find(negcrossing == 1);


end