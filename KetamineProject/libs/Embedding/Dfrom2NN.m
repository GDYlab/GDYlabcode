function [estd,mu] = Dfrom2NN(d)
% Dfrom2NN
% Estimate the intrinsic dimensionality of a point cloud from the ratio of
% the two nearest-neighbour distances ("TwoNN" estimator).
%   Facco, E., d'Errico, M., Rodriguez, A., & Laio, A. (2017). Estimating the
%   intrinsic dimension of datasets by a minimal neighborhood information.
%   Scientific Reports, 7, 12140.
%
% input:  d    - nobs x ndim matrix; rows are observations (time bins),
%                columns are the original dimensions (cells)
% outputs: estd - estimated intrinsic dimension (slope of the line fit)
%          mu   - 1 x nobs ratio r2/r1 of the 2nd- and 1st-nearest-neighbour
%                 distance for every observation (used for the CDF plot)
%
% The estimator uses the fact that, for points sampled from a locally uniform
% density of intrinsic dimension d, mu = r2/r1 follows a Pareto distribution
% with CDF F(mu) = 1 - mu^-d. Hence -log(1-F(mu)) is linear in log(mu) with
% slope d, recovered here by least squares (discarding the largest 10% of
% log(mu) for robustness, as in the reference).
%
% Self-contained R2017b port for the ketamine tutorial (no plotting / toolbox
% helpers). Requires pdist2 (Statistics Toolbox) and ecdf.
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

nobs = size(d,1);

% step 1: pairwise Euclidean distances
dis = pdist2(d,d);

% step 2: for each point find r1 (nearest) and r2 (2nd nearest); mu = r2/r1.
% Ignore self-distance and coincident points (distance 0).
dis(1:nobs+1:end) = inf;
dis(dis==0) = inf;
sorted_dis = sort(dis,1);
r1 = sorted_dis(1,:);
r2 = sorted_dis(2,:);
mu = (r2 + eps) ./ (r1 + eps);

% step 3: empirical CDF of mu
[F,musig] = ecdf(mu);
[unqmusig,idx] = unique(musig);
F = F(idx);

% step 4: line fit of log(mu) against -log(1-F(mu)) through the origin
xx = log(mu);
if numel(unqmusig) <= 2
    fmu = F(1) * ones(size(mu));
else
    fmu = interp1(unqmusig,F,mu);
end
yy = -log(1-fmu);

goodidx = ~isnan(xx) & ~isinf(xx) & ~isnan(yy) & ~isinf(yy);
xx = xx(goodidx);
yy = yy(goodidx);

% discard the largest 10% of log(mu) (heavy-tail trim, per the reference)
[~,newidx] = sort(xx,'ascend');
keplen = round(numel(newidx)*0.9);
kepidx = newidx(1:keplen);
xxkp = xx(kepidx);
yykp = yy(kepidx);

if numel(unique(xxkp)) == 1
    estd = xx(:) \ yy(:);   % fall back to all points if the trim left a single x
else
    estd = xxkp(:) \ yykp(:);
end
end
