
% Inputs:
% spike_counts (n by m matrix with n number of pyramidal cells and m time bins in sleep)
% tbins (vector containing time bins in sleep)
% frametimes (cells containing beginning and end of sleep frames)

% Outputs:
% frame pair correlations for the population vector across sleep frame pairs


rate_vec=[]; %population vector for a frame
for c=1:length(frametimes)
        rate_vec(:,c)=sum(spike_counts(:,tbins>=frametimes{c}(1,3) & tbins<=frametimes{c}(end,3)),2); 
end


corrs=[]; %ensemble frame pair correlations, resultant matrix is symmetric
for c=1:size(rate_vec,2)
      for d=1:size(rate_vec,2)
        corrs(c,d)=corr(rate_vec(:,c),rate_vec(:,d));
      end
end


