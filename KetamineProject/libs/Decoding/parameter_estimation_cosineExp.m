function [pdf,normpdf] = parameter_estimation_cosineExp(d, tc, meanact, poslambda, alllambda)

% parameter_estimation_cosineExp - Computes pdf for position based on
% cosine angle difference and exponential distribution of angle difference

% [pdf,normpdf] = parameter_estimation_cosineExp(d, tc, meanact, posk, allk)
% Given a series of time bins of activities, this function will
% assign probability values to each of the possible states the encoded
% variable (position or head direction) could take in each time bin.
% Estimation is based on Bayes reconstruction method with chi square
% distribution of angle difference
% The inputs to this function are:
%     d, matrix with rows being cells and columns being time bins
%
%     tc, tuning curve or placemap,
%     rows are cells, columns are position bins
%
%     meanact, mean number of spikes from all the cells during the decoding
%     state
%
%     poslambda,  parameter of exponential distribution, a vector of position bins,
%
%     alllambda, parameter of exponential distribution, one scalar,


% The outputs of this function are:
%     pdf, decoding probability with Bayes reconstruction, rows are
%     position bins, columns are time bins

%     normpdf, normalized probability within each time bin


% based on d and tc, we will compute theta_pos, which is the angle
% difference with each averaged population vector, and theta_all, which is
% the angle difference with averaged population vector during the
% decoding state
% pdf = pb_prior * lambda_pos/lambda_all * ...
%        * exp(theta_all*lambda_all-theta_pos*lambda_pos)
% see lyxfile 'DecodingCosSpherical' for the derivation of equation
% Yuchen Zhou Jan 2024  yuchenzhou93@gmail.com

%% get basic parameters
dof = size(d,1);
tbins = size(d,2);
xybins = size(tc,2);

%% get angle difference based on cosine similarity

% optimize code to avoid for loops
meanact = reshape(meanact,[length(meanact),1]); % dof*1
meanact_rept = repmat(meanact,[1, tbins]); % dof*tbins

alladiff = AngleBetweenV(d,meanact_rept,1);% 1*tbins

dnow_repxy = reshape(d,[size(d,1),1,size(d,2)]); % dof*1*tbin
dnow_repxy = repmat(dnow_repxy,[1,xybins,1]); % dof*xybins*tbin

tcnow_rept = reshape(tc,[size(tc,1),size(tc,2),1]); %dof*xybins*1
tcnow_rept = repmat(tcnow_rept,[1,1,tbins]); % dof*xybins*tbin
posadiff = squeeze(AngleBetweenV(dnow_repxy,tcnow_rept,1));% xybins*tbins


%% compute components and reshape to desired size xybins * tbins
lambda_pos = reshape(poslambda,[length(poslambda),1]); % xybins*1

lambdacorr = lambda_pos./alllambda; % xybins*1
lambdacorr = repmat(lambdacorr,[1,tbins]); % xybins*tbins

lambda_pos_tbins = repmat(lambda_pos,[1, tbins]); % xybins*tbins
lambda_all_xytbins = repmat(alllambda,[xybins, tbins]); % xybins*tbins
alladiffrepxy = repmat(alladiff,[xybins, 1]);  % xybins*tbins


expcomp = exp(lambda_all_xytbins.*alladiffrepxy - lambda_pos_tbins.*posadiff);

%% get probability

pdf = lambdacorr .* expcomp;
sumpdf = nansum(pdf,1);
normpdf = pdf./repmat(sumpdf,[xybins 1]);
end


