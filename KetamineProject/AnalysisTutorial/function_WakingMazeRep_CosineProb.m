function rep = function_WakingMazeRep_CosineProb(ses,mazel,times,M1,M2,PlFields,PlFmesh,ClF,fratime,tbin)
% function_WakingMazeRep_CosineProb
% Quantifies how strongly awake-rest (immobile) population activity represents
% the maze, using a probability measure based on population-vector cosine
% similarity (manuscript Fig. 3). For one animal, per session, it:
%   (1) builds a lap-to-lap reference: for each track position bin, the angular
%       variability of the lap population vectors around the mean place-field
%       vector (this sets the per-bin sharpness lambda),
%   (2) decodes each awake-rest frame against the session's maze with a
%       cosine-exponential decoder: each 20 ms population vector is scored by
%       its angle to every position's place-field vector, turned into a
%       probability via an exponential whose rate comes from (1) and from the
%       rest-state angular spread, and
%   (3) summarizes each frame by the maze-representation log-probability
%       (mean over time bins of log max-probability).
% Ketamine (Run2Bi) is expected to weaken maze representation (lower values).
%
% inputs:
%   ses,mazel : track geometry / behavioral trajectory (see ../DataStructure.md)
%   times     : nsess x 2 cell; times{is,idir} is an nLap x 2 matrix of
%               directional lap start/end times (s)
%   M1,M2,PlFields,PlFmesh,ClF : place fields and spikes (see ../DataStructure.md)
%   fratime   : nsess x 1 cell; fratime{is} is an nFrame x 2 matrix of
%               awake-rest frame start/end times (s)
%   tbin      : decoding time bin, s (default 0.02 = 20 ms)
%
% output: rep struct with fields
%   rep.sesstr : 1 x nsess session labels (data order)
%   rep.logpb  : 1 x nsess cell; per-frame maze-representation log-probability
%                (pooled over the two decoding directions)
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

if nargin < 10 || isempty(tbin); tbin = 0.02; end

%% preprocess
[M1,M2,ClF,PlFields,nc] = CA_ExcludeInt(M1,M2,ClF,PlFields);
ClF = cellfun(@double,ClF,'uni',0); % cells with no spikes in a session may be
                                    % stored as empty non-double; force double
ses   = Detour_Ordertracks(ses);
mazel = uniquemazeltime(mazel);
PlFmesh = double(PlFmesh);
nsess = size(PlFields,2);
ntrack = length(ses(1).tra_p);
dirs = [-1,1];

% adapt to the number of effective (non-empty) sessions this animal has
[validses,rep.sesstr] = KTM_ValidSessions(mazel);
slotnames = {'Run1','Run2Bi','Run3','Run2woBi'}; % full slot labels (for progress text)
rep.logpb  = cell(1,nsess);

%% (1) lap-to-lap angular variability reference (per direction/session/track/bin)
% sdtangdiffvar{idir,is,it} : variance across laps of the angle between each
% lap's per-bin population vector and the mean place-field vector at that bin.
[sdtangdiffvar,segmesh] = local_LapAngleVar(ses,mazel,times,ClF,PlFields,PlFmesh,nc,ntrack,nsess,dirs);

%% (2-3) decode awake-rest frames and score maze representation
fprintf('function_WakingMazeRep_CosineProb: %d sessions...\n',nsess);
for is = 1:nsess
    if isempty(mazel{is}) || isempty(fratime{is}); continue; end
    fprintf('  session %d/%d (%s): %d rest frames\n',is,nsess,slotnames{min(is,4)},size(fratime{is},1));
    ft = fratime{is};
    sclf = ClF(:,is);

    % rest-state mean population vector and angular spread (sets alllambda)
    spkcount = [];
    for il = 1:size(ft,1)
        tb = local_tbins(ft(il,1),ft(il,2),tbin);
        sc = zeros(nc,size(tb,1));
        for ic = 1:nc; sc(ic,:) = histcounts_row_kf(sclf{ic}(:,1),tb); end
        spkcount = cat(2,spkcount,sc);
    end
    meanact = nanmean(spkcount,2);
    allang = nan(1,size(spkcount,2));
    for ib = 1:size(spkcount,2); allang(ib) = AngleBetweenV(spkcount(:,ib),meanact); end
    alllambda = sqrt(1/var(allang,0,2,'omitnan'));

    % decode each frame against this session's maze, both directions
    framelogpb = [];
    for idir = 1:2
        % build the decode template (place fields + per-bin lambda) for the track
        tl = ses(is).tralim(1,:);
        bidx = PlFmesh>=tl(1) & PlFmesh<=tl(2);
        plf = squeeze(PlFields(:,is,idir,bidx));
        plfmeshnow = PlFmesh(bidx);
        % per-bin lambda from the lap reference, interpolated onto plfmeshnow
        varpos = sdtangdiffvar{idir,is,1};
        meshpos = segmesh{idir,is,1};
        if isempty(varpos); continue; end
        poslambda = sqrt(1./varpos);
        poslambda(~isfinite(poslambda)) = nanmedian(poslambda(isfinite(poslambda)));
        poslambda = interp1(meshpos,poslambda,plfmeshnow,'linear','extrap');
        poslambda = poslambda(:);

        for il = 1:size(ft,1)
            sf = BayePosDecode_LabCode_Lite_CosineExp(plf,poslambda,plfmeshnow,sclf,...
                ft(il,:),[min(plfmeshnow) max(plfmeshnow)],tbin,0,meanact,alllambda);
            pdf = sf.pdf;                    % space x time
            if isempty(pdf); continue; end
            maxtbin = log(max(pdf,[],1));    % log max-prob per time bin
            framelogpb = cat(1,framelogpb,mean(maxtbin(isfinite(maxtbin))));
        end
    end
    rep.logpb{is} = framelogpb;
end
% keep only the effective (non-empty) sessions, aligned with rep.sesstr
rep.logpb = rep.logpb(validses);
end

%% ----------------------------------------------------------------------- %%
function [sdtangdiffvar,segmesh] = local_LapAngleVar(ses,mazel,times,ClF,PlFields,PlFmesh,nc,ntrack,nsess,dirs)
% lap-to-lap angular variability of the population vector around the mean
% place field, per direction/session/track/position-bin.
lapcri = 0.75;  % a lap must cover >=75% of the track length
grid = 2;       % place-field spatial bin (matches PlFmesh)
sdtangdiffvar = cell(2,nsess,ntrack);
segmesh = cell(2,nsess,ntrack);
for idir = 1:2
    for is = 1:nsess
        smazel = mazel{is};
        if isempty(smazel); continue; end
        mt = smazel(:,6); mlp = smazel(:,4); mtra = smazel(:,3); mld = smazel(:,8);
        for it = 1:ntrack
            posl = ses(is).tralim(it,:);
            tran = ses(is).tra_p(it);
            bidx = PlFmesh>=posl(1) & PlFmesh<=posl(2);
            aveplf = squeeze(PlFields(:,is,idir,bidx)); % cells x bins (mean place field)
            mesh = PlFmesh(bidx);
            dirt = times{is,idir};
            lapvecs = {};   % per-lap population vector (cells x bins) on the track
            for il = 1:size(dirt,1)
                sel = (mt>=dirt(il,1)) & (mt<dirt(il,2)) & (mtra==tran) & (mld==dirs(idir));
                if sum(sel)==0; continue; end
                lp = mlp(sel);
                if (max(lp)-min(lp)) < lapcri*(max(posl)-min(posl)); continue; end
                lapplf = CA_Getplacefield5cms(ClF,smazel(sel,:),is,grid);
                lpf = squeeze(lapplf(:,idir,bidx)); % cells x bins for this lap
                lapvecs{end+1} = lpf; %#ok<AGROW>
            end
            nlap = numel(lapvecs);
            if nlap < 2; continue; end
            nbin = numel(mesh);
            ang = nan(nbin,nlap);
            for ilap = 1:nlap
                for ip = 1:nbin
                    ang(ip,ilap) = AngleBetweenV(lapvecs{ilap}(:,ip),aveplf(:,ip));
                end
            end
            sdtangdiffvar{idir,is,it} = var(ang,0,2,'omitnan');
            segmesh{idir,is,it} = mesh;
        end
    end
end
end

%% ----------------------------------------------------------------------- %%
function tb = local_tbins(ts,te,tau)
% non-overlapping time bins of width tau spanning [ts te]
edges = ts:tau:te;
if numel(edges)==1; edges = [edges, edges+tau]; end
tb = [edges(:), edges(:)+tau];
end
