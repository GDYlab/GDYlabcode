function pctdisplay(i,n,step)
% this function print the progress with i being current loop, n being total
% number of iterations, it will display every step percentages

if nargin < 3
    step = 10;
end

pctdisp = step:step:100;
ndisp = round(n.*pctdisp./100);

if ismember(i,ndisp)
    dispidx = ndisp == i;
    pctnum = pctdisp(dispidx);
    disp([num2str(pctnum),'% completed...'])
end

end