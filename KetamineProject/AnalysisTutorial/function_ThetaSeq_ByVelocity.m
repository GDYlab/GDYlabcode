function tseq = function_ThetaSeq_ByVelocity(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,thetapk,velbin,winlen)
% function_ThetaSeq_ByVelocity
% Computes the velocity-dependent theta sequence (manuscript Fig. 2) for one
% animal. For each run session and running direction it:
%   (1) Bayesian-decodes position during run in short (winlen) windows,
%   (2) finds each theta cycle (consecutive theta peaks) within run epochs,
%   (3) groups theta cycles by the mean running velocity within the cycle,
%   (4) re-centers every decoded time bin on the animal's actual position so
%       0 = current location, and
%   (5) averages the aligned posteriors across theta cycles within each
%       velocity bin.
% A theta sequence appears as a forward-sweeping diagonal (decoded position
% moving from behind to ahead of the animal within each theta cycle); the
% sweep grows with running speed.
%
% inputs:
%   ses,mazel : track geometry / behavioral trajectory (see ../DataStructure.md)
%   M1,M2,PlFields,PlFmesh,ClF : place fields and spikes (see ../DataStructure.md)
%   thetapk   : nsess x 1 cell; thetapk{is} is a vector of theta-peak times (s)
%   velbin    : edges of running-velocity bins, cm/s (default [5 15 30 120])
%   winlen    : decoding window length, s (default 0.04 = 40 ms for theta scale)
%
% output: tseq struct with fields
%   tseq.sesstr     : 1 x nsess session labels (data order)
%   tseq.velbin     : the velocity bin edges used
%   tseq.cttmesh    : 1 x nt time mesh relative to theta-cycle center (s)
%   tseq.relposmesh : 1 x nr relative-position mesh (cm); positive = decoded
%                     position AHEAD of the animal in its travel direction
%   tseq.meanpb     : {nsess, 2, nvbin} cell; each is nr x nt averaged posterior
%   tseq.ncyc       : nsess x 2 x nvbin count of theta cycles averaged
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

if nargin < 9 || isempty(velbin); velbin = [5 15 30 120]; end
if nargin < 10 || isempty(winlen); winlen = 0.04; end

%% preprocess
[~,~,ClF,PlFields,~] = CA_ExcludeInt(M1,M2,ClF,PlFields); % drop interneurons
ses   = Detour_Ordertracks(ses);
mazel = uniquemazeltime(mazel);
PlFmesh = double(PlFmesh);
nsess = length(mazel);
dirs  = [-1,1];
nvbin = length(velbin)-1;

%% parameters
overlap  = (winlen-0.004)/winlen; % 4 ms window step
thetaext = 0.16;                  % average over +/-160 ms around theta-cycle center
meshbin  = nanmean(diff(PlFmesh));
if meshbin >= 5
    pbext = 7;   % coarse (>=5 cm) place-field mesh
else
    pbext = 20;  % fine (~2 cm) place-field mesh -> wider relative-position window
end

% adapt to the number of effective (non-empty) sessions this animal has
[validses,tseq.sesstr] = KTM_ValidSessions(mazel);
slotnames = {'Run1','Run2Bi','Run3','Run2woBi'}; % full slot labels (for progress text)
tseq.velbin = velbin;
tseq.meanpb = cell(nsess,2,nvbin);
tseq.ncyc   = zeros(nsess,2,nvbin);
cttmesh = []; relposmesh = [];

%% main loop
fprintf('function_ThetaSeq_ByVelocity: decoding %d sessions x 2 directions...\n',nsess);
for is = 1:nsess
    smazel = mazel{is};
    if isempty(smazel); continue; end
    mt = smazel(:,6); mlp = smazel(:,4); mld = smazel(:,8);
    tl = ses(is).tralim(1,:);              % single linear track
    tpk = thetapk{is};
    for idir = 1:2
        fprintf('  session %d/%d (%s), direction %d/2 ...\n',is,nsess,slotnames{min(is,4)},idir);
        % run epochs: on track, correct direction, lasting >= 2 s
        tse = GetEpochwithCrioverLen(mt,(mlp>=tl(1)&mlp<=tl(2))&(mld==dirs(idir)),2);
        if isempty(tse); continue; end
        % decode ONLY within run epochs (with this direction's place fields).
        % Passing the epochs as the decoder time ranges skips immobile periods;
        % posteriors are window-local so results are identical to a full-session
        % decode, but far cheaper.
        plf = squeeze(PlFields(:,is,idir,:));
        sf  = BayePosDecode_LabCode_Lite(plf,PlFmesh,ClF(:,is),tse,tl,winlen,overlap);
        pdf = cat(2,sf.pdf);               % decoded posterior, space x time (all epochs)
        sb  = sf(1).spacebin(:);           % decoded position bins (cm)
        tcc = arrayfun(@(s) mean(s.tbin,2),sf,'uni',0);
        tc  = cat(1,tcc{:});               % window-center times
        dcddt   = winlen*(1-overlap);      % step between window centers
        meshsize = nanmean(diff(sb));
        npoints  = round(thetaext/dcddt);
        % accumulate aligned theta sequences, grouped by velocity bin
        acc = cell(1,nvbin);
        for iep = 1:size(tse,1)
            thi = tpk(tpk>=tse(iep,1) & tpk<=tse(iep,2));
            for k = 1:length(thi)-1
                ts = thi(k); te = thi(k+1); tcen = (ts+te)/2;
                % mean velocity within the theta cycle -> velocity bin
                sl = interp1(mt,mlp,ts); el = interp1(mt,mlp,te);
                vcyc = abs(sl-el)/(te-ts);
                vb = find(vcyc>=velbin(1:end-1) & vcyc<velbin(2:end),1);
                if isempty(vb); continue; end
                % regular time mesh spanning +/-npoints around theta center
                tnew = linspace(tcen-npoints*dcddt,tcen+npoints*dcddt,2*npoints+1);
                di = idxinrange(tc,[tnew(1) tnew(end)]);
                if sum(di) < 3; continue; end
                % interpolate posterior onto the regular time mesh
                pbnow = pdf(:,di); tbb = tc(di);
                tmp = zeros(numel(sb),numel(tnew));
                for ip = 1:numel(sb); tmp(ip,:) = interp1(tbb,pbnow(ip,:),tnew); end
                % re-center each time bin on the animal's actual position
                lpn = interp1(mt,mlp,tnew);
                seq = zeros(2*pbext+1,numel(tnew));
                for it = 1:numel(tnew)
                    sm = linspace(lpn(it)-meshsize*pbext,lpn(it)+meshsize*pbext,2*pbext+1);
                    seq(:,it) = interp1(sb,tmp(:,it),sm);
                end
                % orient the relative-position axis to the travel direction so
                % positive = AHEAD of the animal regardless of running direction
                % (otherwise the two directions are mirror images and cancel
                %  when pooled). seq rows increase with linear position; for the
                % descending direction, flip so ahead maps to positive.
                if dirs(idir) < 0; seq = flipud(seq); end
                if isempty(acc{vb}); acc{vb} = seq; else; acc{vb} = cat(3,acc{vb},seq); end
                if isempty(cttmesh)
                    cttmesh = tnew - tcen;
                    relposmesh = linspace(-meshsize*pbext,meshsize*pbext,2*pbext+1);
                end
            end
        end
        for vb = 1:nvbin
            tseq.ncyc(is,idir,vb) = size(acc{vb},3) * ~isempty(acc{vb});
            if ~isempty(acc{vb}); tseq.meanpb{is,idir,vb} = nanmean(acc{vb},3); end
        end
    end
end
tseq.cttmesh = cttmesh;
tseq.relposmesh = relposmesh;
% keep only the effective (non-empty) sessions, aligned with tseq.sesstr
tseq.meanpb = tseq.meanpb(validses,:,:);
tseq.ncyc   = tseq.ncyc(validses,:,:);
end
