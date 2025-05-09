function pworatio = FramePWOrderPct(Seq,clu2pyr)
% FramePWOrderPct get the pairwise order based on given frames, the result
% is a matrix indicate probability of row fire before column


% inputs:   Seq, cell array of frames of interets, within each cell is a
%           matrix with row being spk events, column 1 is spk times, column
%           2 is cluster id
%           clu2pyr, mapping cluster id to pyr id
% output:   pworatio, n*n  probability matrix with n being the number 
%           of pyr units (max value in clu2pyr). tpm(i,j) means the
%           possibility of pyr i fires before pyr j given both of them are
%           active in a frame
% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

nc = nanmax(clu2pyr);

% get order matrix for all the frames
pwosnow = nan(nc,nc,length(Seq));
for il = 1:length(Seq)
    % make sure in Seq, spk events are ordered in time, also map the clu id to pyr id
    tempspk = Seq{il};
    [~,sind] = sort(tempspk(:,1),'ascend');
    tempspk = tempspk(sind,:);
    tempspk(:,2) = clu2pyr(tempspk(:,2));
    % now for each active pyr units in this frame, get it center of mass
    % spk time
    actp = unique(tempspk(:,2));
    actcm = nan(1,length(actp));
    for ic = 1:length(actp)
        cind = tempspk(:,2) == actp(ic);
        spkct = tempspk(cind,1);
        actcm(ic) = nanmean(spkct);
    end
    [~,cmsort] = sort(actcm,'ascend');
    actpst = actp(cmsort);
    pwosnow(:,:,il) = CA_PairWiseOrder(nc,actpst);
    % output:    pwo, pair wise order matrix, nc*nc. If the value is nan,
    %            meaning one of the pyr is not active in this sequence; on the
    %            diagonal, if the pyr is active, then the value is 0; If the
    %            value is 1, meaning the pyr in row fires before the pyr in
    %            column; if the value is -1, meaning the pyr in row fires after
    %            the pyr in column
end


% make a copy counting for each pair, the instance of being active in
% the same frame
pwocount = pwosnow;
pwocount(~isnan(pwocount)) = 1;
pwocount = nansum(pwocount,3);
% make a copy counting for each pair, the instance of row proceding
% column
pworc = pwosnow;
pworc(pworc ~= 1) = nan;
pworc = nansum(pworc,3);
pworatio = pworc./pwocount;

end


