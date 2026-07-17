function beh = function_Beh_RunVelocityMeasures(mazel,ses,info)
% function_Beh_RunVelocityMeasures
% Computes the four run-behavior measures used to characterize the
% ketamine impact on locomotion (manuscript Fig. 1), for one animal,
% across the effective run sessions:
%   (1) distance travelled per minute  (overall locomotor output)
%   (2) angle between instantaneous velocity and the track direction
%       (straightness of the run path)
%   (3) maximum velocity per continuous-run epoch
%   (4) velocity profile vs distance to the track ends, separately while
%       approaching an end (decelerating) and departing an end (accelerating)
%
% inputs:
%   mazel : nsess x 1 cell, behavioral trajectory per session. Each element
%           is n x 8: [x, y, track, linpos, rundir, time, vel, lapdir].
%           x,y are in camera pixels; linpos and vel are in cm and cm/s.
%           See ../DataStructure.md.
%   ses   : 1 x nsess struct, track geometry. ses(i).tra_p lists the active
%           track id(s); ses(i).tralim gives [start end] linear position.
%   info  : dataset info struct; info.pxpercm converts pixels to cm.
%
% NOTE on session order: sessions are kept in DATA order. For an animal with
%   4 effective sessions (e.g. the URat4 sample) that order is
%   {Run1, Run2Bi, Run3, Run2woBi}; with 3 sessions it is {Run1, Run2Bi, Run3}.
%   beh.sesstr gives the correct label for each session index. (The manuscript
%   figure re-orders Run2woBi before Run3 for display only.)
%
% output: beh struct with fields
%   beh.sesstr      : 1 x nsess session labels (data order)
%   beh.dispermin   : 1 x nsess, distance travelled per minute (cm/min)
%   beh.maxvel      : 1 x nsess cell, max velocity (cm/s) per run epoch
%   beh.trackangle  : 1 x nsess cell, angle (deg) between velocity and track
%                     direction, pooled over run-epoch samples and directions
%   beh.disbinc     : 1 x nbin, distance-to-end bin centers (cm)
%   beh.velapproach : nsess x nbin cell, velocity samples while approaching an
%                     end (binned by distance to the upcoming end)
%   beh.veldepart   : nsess x nbin cell, velocity samples while departing an
%                     end (binned by distance from the previous end)
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

%% preprocess
mazel = uniquemazeltime(mazel);   % remove duplicate timestamps
ses   = Detour_Ordertracks(ses);  % order tracks by linear position
nsess = length(mazel);            % stored session slots
dirs  = [-1,1];                    % the two running directions
% adapt to the number of effective (non-empty) sessions this animal has
[validses,beh.sesstr] = KTM_ValidSessions(mazel);

%% (1) distance travelled per minute
% sum of |delta linear position| (excluding tracking jumps) over the session
% duration, scaled to per minute.
jumpcri = 100; % cm; ignore single-sample jumps larger than this (tracking errors)
dispermin = nan(1,nsess);
for is = 1:nsess
    smazel = mazel{is};
    if isempty(smazel); continue; end
    lpos = smazel(:,4); t = smazel(:,6);
    dlpos = abs(diff(lpos));
    dlpos(dlpos > jumpcri) = nan;
    dispermin(is) = nansum(dlpos)/(max(t)-min(t))*60;
end
beh.dispermin = dispermin;

%% (2) angle between velocity and track direction
% vcri/tlencri/discri define the continuous-run epochs used (final-result values)
angvcri = 20; angtlen = 1; angdiscri = 30;
trackangle = cell(1,nsess);
for is = 1:nsess
    smazel = mazel{is};
    if isempty(smazel); continue; end
    for idir = 1:2
        ag = local_VTrackAngle(ses,is,idir,smazel,info,angvcri,angtlen,angdiscri);
        trackangle{is} = cat(1,trackangle{is},ag(:));
    end
end
beh.trackangle = trackangle;

%% (3) maximum velocity per continuous-run epoch
maxvcri = 10; maxtlen = 2; maxdis = 50;
maxvel = cell(1,nsess);
for is = 1:nsess
    smazel = mazel{is};
    if isempty(smazel); continue; end
    t = smazel(:,6); lpos = smazel(:,4); lapdir = smazel(:,8);
    lpvel = local_lpvel(abs(smazel(:,7)),0.1); % low-pass velocity
    for idir = 1:2
        runidx = (lapdir == dirs(idir)) & (lpvel >= maxvcri);
        tse = GetEpochwithCrioverLen(t,runidx,maxtlen);
        for iep = 1:size(tse,1)
            ix = idxinrange(t,tse(iep,:));
            if (max(lpos(ix)) - min(lpos(ix))) < maxdis; continue; end % require real traversal
            maxvel{is} = cat(1,maxvel{is},max(lpvel(ix)));
        end
    end
end
beh.maxvel = maxvel;

%% (4) velocity profile vs distance to track ends (approach / depart)
d2vcri = 2; d2tlen = 1;
discorbin = 0:2:60;
disbinc = GetCenterfromEdge(discorbin);
velapproach = cell(nsess,length(disbinc)); % distance to upcoming end (decelerating)
veldepart   = cell(nsess,length(disbinc)); % distance from previous end (accelerating)
for is = 1:nsess
    smazel = mazel{is};
    if isempty(smazel); continue; end
    t = smazel(:,6); lpos = smazel(:,4); veldir = smazel(:,5);
    lpvel = local_lpvel(abs(smazel(:,7)),0.1); % low-pass velocity
    % distance to each track end (corner) in this session
    corners = unique(ses(is).tralim(:))';
    dd = bsxfun(@minus,lpos,corners);          % signed distance to each end
    dprev = dd; dprev(dprev<0) = nan;          % distance from a passed end
    dnext = dd; dnext(dnext>0) = nan; dnext = abs(dnext); % distance to an upcoming end
    dprevmin = nanmin(dprev,[],2);
    dnextmin = nanmin(dnext,[],2);
    for idir = 1:2
        runidx = (veldir == dirs(idir)) & (lpvel >= d2vcri);
        [~,linidx] = GetEpochwithCrioverLen(t,runidx,d2tlen);
        v = lpvel(linidx);
        dap = dnextmin(linidx); % approaching the upcoming end
        ddp = dprevmin(linidx); % departing the previous end
        for ib = 1:length(disbinc)
            ia = idxinrange(dap,discorbin([ib,ib+1]));
            velapproach{is,ib} = cat(1,velapproach{is,ib},v(ia));
            id = idxinrange(ddp,discorbin([ib,ib+1]));
            veldepart{is,ib} = cat(1,veldepart{is,ib},v(id));
        end
    end
end
beh.disbinc = disbinc;
beh.velapproach = velapproach;
beh.veldepart = veldepart;

%% keep only the effective (non-empty) sessions, aligned with beh.sesstr
beh.dispermin   = beh.dispermin(validses);
beh.trackangle  = beh.trackangle(validses);
beh.maxvel      = beh.maxvel(validses);
beh.velapproach = beh.velapproach(validses,:);
beh.veldepart   = beh.veldepart(validses,:);

end

%% ----------------------------------------------------------------------- %%
function lpvel = local_lpvel(absvel,cutoffHz)
% zero-phase low-pass filter of the velocity trace (FIR, ~cutoffHz)
absvel(isnan(absvel)) = 0;
MyFilt = fir1(30,cutoffHz/(30/2));
lpvel = Filter0(MyFilt,absvel);
end

%% ----------------------------------------------------------------------- %%
function allag = local_VTrackAngle(ses,is,idir,smazel,info,vcri,tlencri,discri)
% angle (deg) between the instantaneous velocity vector and the track vector,
% for run epochs in session is, direction idir.
dirsign = [-1,1];
t = smazel(:,6); lpos = smazel(:,4);
x = smazel(:,1)./info.pxpercm; % pixels -> cm
y = smazel(:,2)./info.pxpercm;
lpvel = local_lpvel(abs(smazel(:,7)),0.5);
allag = [];
for it = 1:length(ses(is).tra_p)
    trange = ses(is).tralim(it,:);
    onidx = idxinrange(lpos,trange);
    % track vector: mean of the 100 extreme-linpos samples at each end
    [~,ord] = sort(lpos(onidx),'ascend');
    xt = x(onidx); yt = y(onidx); xt = xt(ord); yt = yt(ord);
    npt = min(100,floor(numel(xt)/2));
    if npt < 1; continue; end
    tstart = [nanmean(xt(1:npt)), nanmean(yt(1:npt))];
    tend   = [nanmean(xt(end-npt+1:end)), nanmean(yt(end-npt+1:end))];
    tvec = tend - tstart; tvec = tvec./norm(tvec);
    % continuous-run epochs on this track in this direction
    runidx = idxinrange(lpos,trange) & (lpvel >= vcri);
    tse = GetEpochwithCrioverLen(t,runidx,tlencri);
    for iep = 1:size(tse,1)
        ix = idxinrange(t,tse(iep,:));
        lp = lpos(ix);
        if (lp(end)-lp(1))*dirsign(idir) < discri; continue; end % require real traversal in idir
        ept = t(ix);
        ex = smoothdata(x(ix),'movmedian',4); % suppress fast head jitter
        ey = smoothdata(y(ix),'movmedian',4);
        vx = diff(ex)./diff(ept); vy = diff(ey)./diff(ept);
        speed = sqrt(vx.^2 + vy.^2);
        bad = speed > 100; vx(bad) = []; vy(bad) = []; speed(bad) = [];
        valong = [vx(:) vy(:)] * tvec';
        cosang = abs(valong)./abs(speed);
        allag = cat(1,allag,rad2deg(acos(cosang)));
    end
end
end
