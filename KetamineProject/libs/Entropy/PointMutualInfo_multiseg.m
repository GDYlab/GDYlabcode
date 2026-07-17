function [MIXY,MIYX,hxp,hyp,hxn,hyn,hxpyn,hypxn,nMIXY,nMIYX] = PointMutualInfo_multiseg(xseg,yseg,tlen)
% PointMutualInfo_multiseg this function compute the mutual infomation
% between two time sereis  discrete processes x and y 
% with given the memory length tlen 
% in this verions both x and y have several segments, and the MI
% should be only computed within segments. For example, I don't
% want to concatenate frames as the end of one frame should not predict the
% start of next frame
% see     lyx file "InfoEntropy.lyx"
% inputs: xseg,yseg, should be cell array with the same length, each
%              element present the time series of one segment
%              In the element, if they are vectors, they should
%              have the same length. If they are matrix (multivariate),
%              they should match in row numbers.  In that case, we will
%              concatenate DOF of time series and delayed construction
%         tlen, memory of time series, shorter than length of x and y
% output: MIXY, MI of past of X and current Y
%         MIXY = I(Xpast;Ynow) = H(Xpast) + H(Ynow) - H(Xpast,Ynow)
%              = hxp + hyn - hxpyn
%         MIYX, MI of past of Y and current X
%         MIYX = I(Ypast;Xnow) = H(Ypast) + H(Xnow) - H(Ypast,Xnow)
%              = hyp + hxn - hypxn
%         nMIXY normalized MI = MIXY/H(Ynow)
%         nMIYX normalized MI = MIYX/H(Xnow)
% Yuchen Zhou  yuchenzhou93@gmail.com

%  see examples in the file tetptest.m
%% check inputs

if length(xseg) ~= length(yseg)
    error('Two time series segments should match in length')
end

nseg = length(xseg);
%% construct time series matrics with rows being observations and columns being DOF
ynowxpast = [];
xnowypast = [];
ypastall = [];
xpastall = [];
ynowall = [];
xnowall = [];
for iseg = 1:nseg
    x = xseg{iseg};
    y = yseg{iseg};
    alllen = size(x,1);
    % get number of observations
    nob = alllen-tlen;
    if nob <= 0
       continue 
    end
    % construct observation matrix Yi+1
    ynow = nan(nob,size(y,2));
    % construct observation matrix Xi+1
    xnow = nan(nob,size(x,2));
    % construct observation matrix YI
    ypast = nan(nob,tlen,size(y,2));
    % construct observation matrix XI
    xpast = nan(nob,tlen,size(x,2));
    for iob = 1:nob
        ynow(iob,:) = y(iob+tlen,:);
        xnow(iob,:) = x(iob+tlen,:);
        ypast(iob,:,:) = y(iob:iob+tlen-1,:);
        xpast(iob,:,:) = x(iob:iob+tlen-1,:);
    end
    
    % concatenate dimension 2&3 together
    ypast = reshape(ypast,nob,[]);
    xpast = reshape(xpast,nob,[]);
    
    ynowxpasttmp = cat(2,ynow,xpast);
    xnowypasttmp = cat(2,xnow,ypast);
    
    ypastall = cat(1,ypastall,ypast);
    xpastall = cat(1,xpastall,xpast);
    ynowxpast = cat(1,ynowxpast,ynowxpasttmp);
    xnowypast = cat(1,xnowypast,xnowypasttmp);
    ynowall = cat(1,ynowall,ynow);
    xnowall = cat(1,xnowall,xnow);
end
%% compute the joint entropy
% compute joint entropy

hxp = jointentropy(xpastall);
hyp = jointentropy(ypastall);

hxn = jointentropy(xnowall);
hyn = jointentropy(ynowall);

hxpyn = jointentropy(ynowxpast);
hypxn = jointentropy(xnowypast);

MIXY = hxp + hyn - hxpyn;
MIYX = hyp + hxn - hypxn;


nMIXY = MIXY/hyn;
nMIYX = MIYX/hxn;

nMIXY = myclip(nMIXY,[0,1]);
nMIYX = myclip(nMIYX,[0,1]);
end


function h = jointentropy(x)

% ndim = size(x,2);
% nobs = size(x,1);

[uniqueRows, ~, idx] = unique(x, 'rows');
occurrences = accumarray(idx, 1);

% Normalize to get probabilities
probabilities = occurrences / sum(occurrences);

% Compute entropy
h = -sum(probabilities .* log2(probabilities + eps)); 

end