function tpm = FrameTPM_CMass(Seq,clu2pyr)
% FrameTPM_CMass find the transition probability matrix of cell activity in
% frames, in each frame, the cell activity sequence is obtained by sorting
% based on center of mass of spikes

% inputs:   Seq, cell array of frames of interets, within each cell is a
%           matrix with row being spk events, column 1 is spk times, column
%           2 is cluster id
%           clu2pyr, mapping cluster id to pyr id
% output:   tpm, n*n transition probability matrix with n being the number 
%           of pyr units (max value in clu2pyr). tpm(i,j) means the
%           possibility of observing pyr j as the next active cell
%           given pyr i fires. in tpm, sum over D2 gives value 1

% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

nc = nanmax(clu2pyr);
% create a vector containing all the spk events from all the frames
S = [0];
% we will use 0 to separate different frames
for il = 1:length(Seq)
    % make sure in Seq, spk events are ordered in time, also map the clu id to pyr id
    temp = Seq{il};
    [~,sind] = sort(temp(:,1),'ascend');
    temp = temp(sind,:);
    temp(:,2) = clu2pyr(temp(:,2));
    % now for each active pyr units in this frame, get it center of mass
    % spk time
    actp = unique(temp(:,2));
    actcm = nan(1,length(actp));
    for ic = 1:length(actp)
        cind = temp(:,2) == actp(ic);
        spkct = temp(cind,1);
        actcm(ic) = nanmean(spkct);
    end
    [~,cmsort] = sort(actcm,'ascend');
    actpst = actp(cmsort);
    S = cat(2,S,berow(actpst));
    S = cat(2,S,0);
end

n = numel(S);
y = zeros(nc,1);
p = zeros(nc,nc);
for k=1:n-1
    if S(k) ~= 0 &&  S(k+1) ~= 0
        y(S(k)) = y(S(k)) + 1;
        p(S(k),S(k+1)) = p(S(k),S(k+1)) + 1;
    end
end
tpm = bsxfun(@rdivide,p,y); tpm(isnan(tpm)) = 0;
%% we will further process 0 and 1 values in this matrix, to aviod 0 probability in shuffle
% make the diagonal nan
tpm = setdiag(tpm,nan);

% setting 0 and 1 values to the min or max value in (0,1)
tpmall = tpm(:);
tpm(tpm == 0) = nanmin(tpmall(tpmall>0));
tpm(tpm == 1) = nanmax(tpmall(tpmall<1));
% normalize again to make sure sum over D2 give value 1
for il = 1:size(tpm,1)
    tpm(il,:) = tpm(il,:)./nansum(squeeze(tpm(il,:)));
end

end