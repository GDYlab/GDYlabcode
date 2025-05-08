function [SeqC,SeqSig] = function_DetourThetaCycle_PreSlpFrame_RankCorrPct(Seqses,ThetaCycOrder,clu2pyr,nc,nshf)
% function_DetourThetaCycle_PreSlpFrame_RankCorrPct
% this function run Bayesian decoding on detour tracks during post-detour run,
% we add detour segment and control segment to the decoding to
% detect flickering of detour segment during the reversal run
%
% Inputs:     Seqses, a cell array with length equals to number of
%                     sleep sessions, in each cell element, it's again a
%                     cell array with length equals to number of frames in
%                     that sleep. For each frame, it's a n*2 matrix, with
%                     row being all the spikes in this frame, and column
%                     being spike time and clusterID of spikes
%             for example, Seqses{1}{2} = [1.1,3;1.3,4;2.6,5;3.4,3] means in
%                     sleep1 frame 2, cluster 3,4,5 are active, where
%                     cluster 3 fires at 1.1, 3.4, cluster 4 fires at 1.3,
%                     and cluster 5 fires at 2.6;
%             ThetaCycOrder, spike sequence in theta cycles, it's a cell
%                     arrary with D1 being sessions and D2 being
%                     directions, in each cell element, it a cell array
%                     with length equals to number of theta cycle, in each
%                     element, it's a vector of pyr spike order 
%             for example, ThetaCycOrder{2,1}{7} = [5,2,1,10] means the 7th
%                       theta cycle in session 2 direction 1 has the pyr
%                       spike order of pyr5->pyr2->pyr1->pyr10
%             clu2pyr, is a vector mapping cluster ID to pyr ID, if the 
%                      cluster is interneuron, it gives nan,
%             for example, clu2pyr = [1,nan,2,3,nan,4] meaning cluster 2
%                      and 5 are interneurons, while cluster 1,3,4,6
%                      corresponding to Pyr1,2,3,4
%             clu2pyr can be obtained from command:
%             [M1,M2,ClF,PlFields,nc,pyrid,clu2pyr,pyrintid] = CA_ExcludeInt(M1,M2,ClF,PlFields)
%             nc,     total number of pyr units
%             nshf,   number of shuffle to compute rank order correlation
%             
% outputs:  SeqC is a cell array with dimension (session,direction,r and p value from rank order correlation)
%                     in this code, we only have value in detour session in
%                     the first dimension.
%                     in each cell element, it's a n*m matrix, with n being
%                     number of sleep frames, and m being number of theta
%                     cycles. The value is rank order correlation r or p
%                     value. If there are less than 5 common active cells,
%                     the value is NAN because we don't have enough shuffle
%                     to get a percentile
%
%           SeqSig  is an array with dimension (session,direction,number of
%                     signficant pairs and number of valid pairs)
%                     in this code, we only have value in detour session in
%                     the first dimension.
%                     to get ratio of significant correlation for session 3
%                     and direction 1, we can do
%                     ratio = SeqSig(3,1,1)/SeqSig(3,1,2)
%
% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

%% define parameters
detses = [2,3]; % detour sessions
presleep = [1,2]; % sleep session before detour experience
allpyr = 1:nc;

%% concatenate pre-detour sleep frames, get spike sequence based on center of mass
disp('Computing sleep frame sequence based on center of mass')
allslp = [];
nframe = [];
for islp = presleep
    slpseqnow = Seqses{islp};
    % convert cluid to pyrid
    pyrseq = cellfun(@(M) [M(:,1), arrayfun(@(x) clu2pyr(x), M(:,2))], slpseqnow, 'UniformOutput', false);
    pyrcomseq = cell(size(pyrseq));
    for iseq = 1:length(pyrseq)
        cacttime = nan(1,nc);
        for ic = 1:nc
            cidx =  pyrseq{iseq}(:,2) == ic;
            if sum(cidx) == 0
                continue
            end
            cacttime(ic) = nanmean(pyrseq{iseq}(cidx,1));
        end
        actcid = allpyr(~isnan(cacttime));
        actctime = cacttime(~isnan(cacttime));
        [~,sortidx] = sort(actctime,'ascend');
        pyrcomseq{iseq} = becolumn(actcid(sortidx));
    end
    allslp = cat(2,allslp,pyrcomseq);
    nframe = cat(2,nframe,length(pyrcomseq));
end
disp('Getting sleep frame sequence completed')

%% correlate spike order in theta cycles and pre-detour sleep frames
SeqC = cell(4,2,2); % D1 session, D2 directions, D3 r or p value
SeqSig = nan(4,2,2);  % D1 session, D2 directions, D3 sig and valid numbers,
for is = detses
    for idir = 1:2
        disp(['Correlate sleep frame with early detour lap theta cycles in Run',...
            num2str(is),'Direction',num2str(idir),'...'])
        % get the sequence
        seqnow = ThetaCycOrder{is,idir};
        seqcorrr = nan(length(allslp),length(seqnow));
        seqcorrpct = nan(length(allslp),length(seqnow));
        % run correlation with pre-detour sleep frames
        for irun = 1:length(seqnow)
            runseq = becolumn(seqnow{irun});
            for islp = 1:length(allslp)
                slpseq = becolumn(allslp{islp});
                % find overlap between theta cycle and pesleep frames
                ovset = intersect(runseq,slpseq);
                if length(ovset) < 5
                    continue
                end
                runidx = ismember(runseq,ovset);
                slpidx = ismember(slpseq,ovset);
                rnow = runseq(runidx);
                slpnow = slpseq(slpidx);
                r = corr(rnow,slpnow,'Type','Spearman');
                rshf = nan(1,nshf);
                for ishf = 1:nshf
                    rshf(ishf) = corr(rnow(randperm(length(rnow))),slpnow,'Type','Spearman');
                end
                seqcorrr(islp,irun) = r;
                seqcorrpct(islp,irun) = invprctile(rshf,r);
            end
        end
        SeqC{is,idir,1} = seqcorrr;
        SeqC{is,idir,2} = seqcorrpct;
        
        valididx = ~isnan(seqcorrr);
        validp = seqcorrpct(valididx);
        sigp = validp >= 95;
        SeqSig(is,idir,1) = sum(sigp);
        SeqSig(is,idir,2) = length(sigp);
    end
end

end
