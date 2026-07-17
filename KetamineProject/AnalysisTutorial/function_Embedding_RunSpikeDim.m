function emb = function_Embedding_RunSpikeDim(ses,mazel,times,M2,ClF,tbin,vcri)
% function_Embedding_RunSpikeDim
% Estimates the intrinsic dimensionality of the population spiking during fast
% run (manuscript Fig. 4). For each session and running direction it:
%   (1) bins pyramidal-cell spikes into tbin windows along every run lap that
%       crosses the track, keeping only bins above a velocity threshold,
%   (2) pools the bins across laps to form a (time bins) x (cells) activity
%       cloud,
%   (3) estimates that cloud's intrinsic dimension with the TwoNN estimator
%       (Facco et al. 2017; see Dfrom2NN).
% A higher intrinsic dimension means the population activity is less confined
% to a low-dimensional manifold. Ketamine (Run2Bi) is expected to raise it.
%
% inputs:
%   ses,mazel : track geometry / behavioral trajectory (see ../DataStructure.md)
%   times     : nsess x 2 cell; times{is,idir} is the nlap x 2 [start end] time
%               of each run lap for session is and direction idir (idir 1 = the
%               descending direction, idir 2 = the ascending direction)
%   M2        : per-cell parameters; column 21 is cell type (0 = interneuron)
%   ClF       : nc x nsess spikes per cell per session (pyr + int; see
%               DataStructure). Column 1 of each cell is the spike time.
%   tbin      : spike time bin, s (default 0.4)
%   vcri      : run velocity threshold, cm/s (default 10)
%
% output: emb struct with fields
%   emb.sesstr : 1 x nvalid session labels (effective sessions only)
%   emb.dim    : 2 x nvalid intrinsic dimensionality per direction per session
%   emb.mu     : 2 x nvalid cell; the r2/r1 ratios per direction per session
%                (their CDF is the manuscript example panel)
%   emb.exdir  : example direction index used by the manuscript panels (= 2)
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

if nargin < 6 || isempty(tbin); tbin = 0.4; end
if nargin < 7 || isempty(vcri); vcri = 10; end

%% preprocess
ses   = Detour_Ordertracks(ses);
mazel = uniquemazeltime(mazel);
ClF   = cellfun(@double,ClF,'uni',0);  % no-spike cells are stored as uint8
nsess = size(ClF,2);
dirs  = [-1,1];
lapcri = 0.1;   % a lap must cover >=10% of the track length to count
actcri = 1;     % keep only time bins with at least this many total spikes

% keep pyramidal cells only (cell type in M2 col 21, 0 = interneuron)
pyrid = find(M2(:,21,1) ~= 0);
ClF = ClF(pyrid,:);
nc  = numel(pyrid);

ntrack = numel(ses(1).tra_p);  % single linear track for the ketamine data

% adapt to the number of effective (non-empty) sessions this animal has
[validses,emb.sesstr] = KTM_ValidSessions(mazel);
slotnames = {'Run1','Run2Bi','Run3','Run2woBi'};
emb.dim = nan(2,nsess);
emb.mu  = cell(2,nsess);
emb.exdir = 2;

fprintf('function_Embedding_RunSpikeDim: %d sessions...\n',nsess);
for is = 1:nsess
    if isempty(mazel{is}); continue; end
    fprintf('  session %d/%d (%s)\n',is,nsess,slotnames{min(is,4)});
    smazel = mazel{is};
    mt   = smazel(:,6);          % time
    mlp  = smazel(:,4);          % linear position
    mvel = abs(smazel(:,7));     % running speed
    mdir = smazel(:,8);          % lap direction

    for idir = 1:2
        indld = mdir == dirs(idir);
        sesspk = [];             % nc x (pooled bins) spike counts
        sesvel = [];
        for it = 1:ntrack
            tralim = ses(is).tralim(it,:);
            indlpos = idxinrange(mlp,tralim);
            dirt = times{is,idir};
            for il = 1:size(dirt,1)
                indt = idxinrange(mt,dirt(il,:));
                allind = indt & indlpos & indld;   % this lap, track, direction
                if sum(allind) == 0; continue; end
                temppos = mlp(allind);
                % require the lap to actually traverse the track
                if (max(temppos)-min(temppos)) < lapcri*(max(tralim)-min(tralim))
                    continue;
                end
                [lt,uidx] = unique(mt(allind));
                lv = mvel(allind); lv = lv(uidx);
                if numel(lt) < 2; continue; end
                edges = min(lt):tbin:max(lt);
                if numel(edges) < 2; continue; end
                ctr = GetCenterfromEdge(edges);
                binvel = interp1(lt,lv,ctr);
                lapspk = zeros(nc,numel(ctr));
                for ic = 1:nc
                    st = ClF{ic,is}(:,1);
                    st = st(st>=min(lt) & st<=max(lt));
                    lapspk(ic,:) = histcounts(st,edges);
                end
                sesspk = cat(2,sesspk,lapspk);
                sesvel = cat(2,sesvel,binvel(:)');
            end
        end
        if isempty(sesspk); continue; end

        % drop low-velocity bins, then bins with no spikes at all
        sesvel(isnan(sesvel)) = 0;
        sesspk(:,sesvel < vcri) = [];
        sesspk(:,sum(sesspk,1) < actcri) = [];
        if size(sesspk,2) < 3; continue; end

        % intrinsic dimensionality of the (time bins) x (cells) cloud
        [estd,mu] = Dfrom2NN(sesspk');
        emb.dim(idir,is) = estd;
        emb.mu{idir,is}  = mu;
    end
end

% keep only the effective (non-empty) sessions, aligned with emb.sesstr
emb.dim = emb.dim(:,validses);
emb.mu  = emb.mu(:,validses);
end
