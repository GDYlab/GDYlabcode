function [pdf, normpdf, tbins, nact, nspk] = reconstruct_CosineExp(ts, te, cells, tc, meanact, poslambda, alllambda, varargin)
% this function run the decoding based on cosine angle difference
% inputs:
%     ts: start time
%     te: stop time
%     cells : cells(i).time

%     tc, tuning curve or placemap,
%     rows are cells, columns are position bins
%
%     meanact, mean number of spikes from all the cells during the decoding
%     state
%
%     poslambda,  parameter of exponential distribution, a vector of position bins,
%
%     alllambda,  parameter of exponential distribution, one scalar,

%     varargin: tau, decoding time window, 20 ms by default
%               percent_overlap, overlaping window, 0 by default
%               max_pdf_comp_size, pdf segment size, 5000 by default

% outputs:
%     pdf, decoded probability, D1 spatial bins, D2 temporal bins
%     normpdf, normalized pdf with in time bins
%     tbins, n*2 time window range, D1 temporal bins, D2 start and end time
%     nact, get number of active units in each time bin

% based on d and tc, we will compute theta_pos, which is the angle
% difference with each averaged population vector, and theta_all, which is
% the angle difference with averaged population vector during the
% decoding state
% pdf = pb_prior * lambda_pos/lambda_all * ...
%        * exp(theta_all*lambda_all-theta_pos*lambda_pos)
% see lyxfile 'DecodingCosSpherical' for the derivation of equation
% Yuchen Zhou  yuchenzhou93@gmail.com

% % %   example to modify parameter:
% % %   reconstruct_CosineExp(ts, te, cells, tc, meanact, posk, allk, 'tau', 0.25)
% % %   now args.tau is overwrited to 0.02...

%% setting paramaters
args.tau = 0.02;
args.actidx = [];
args.percent_overlap = 0; % originally was 0
args.max_pdf_comp_size = 5000; % originally was 5000, then 2000

args = parseArgs(varargin, args);

tbins = get_tbins(ts, te, args.tau, args.percent_overlap);
spike_counts= zeros(length(cells), size(tbins,1));
for i = 1:length(cells)
    spike_counts(i,:) = histcounts_row_kf(cells(i).time, tbins);
end

% get number of active units in each time bin
acttbin = spike_counts >= 1;
nact = sum(acttbin,1);
nspk = sum(spike_counts,1);
    
%% compute probability
if size(spike_counts,2)>args.max_pdf_comp_size
    n_section = ceil(size(spike_counts,2)/args.max_pdf_comp_size);
    % (internal chunking for large decode ranges; message suppressed for tutorial)
    for i=0:n_section-1
        istart = i*args.max_pdf_comp_size+1;
        iend = min([(i+1)*args.max_pdf_comp_size, size(spike_counts,2)]);
        [pdf_short,n_pdf_s] = parameter_estimation_cosineExp(spike_counts(:,istart:iend), tc, meanact, poslambda, alllambda);

        switch i
            case 0
                pdf_temp = pdf_short;
                n_pdf = n_pdf_s;
            otherwise
                pdf_temp = cat(2,pdf_temp,pdf_short); 
                n_pdf = cat(2,n_pdf,n_pdf_s);
        end
    end
else
    [pdf_temp,n_pdf] =  parameter_estimation_cosineExp(spike_counts, tc, meanact, poslambda, alllambda);
end
pdf = pdf_temp;
normpdf = n_pdf;

end


function tbins = get_tbins(ts, te, tau, percent_overlap)
    tbins = ts:tau*(1-percent_overlap):te;
    if numel(tbins)==1
        tbins = [tbins, tbins+tau];
    end
    tbins = [tbins', tbins'+tau];
end


