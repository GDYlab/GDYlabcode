function [detdecode,mscotpbdir] = function_DetRun_FlickeringofPreMobile(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,T2T4Laps,tbin,overlap)
% function_DetRun_FlickeringofPreMobile
% this function run Bayesian decoding on detour tracks during detour run,
% we add pre-detour mobile segment and control segment to the decoding to
% detect flickering of pre-detour mobile segment
% inputs:     ses, mazel, M1, M2, PlFields, PlFmesh, ClF, see documents: 
%             "DataStructure"
%             T2T4Laps is a cell array with D1 being track, D2 being
%                      directions, and D3 being seesion, for each element,
%                      it's a n*2 matrix containing start and end  time of
%                      laps on that track, session, and direction
%             for example, T2T4Laps{2,1,4} = [0 5;17 22] define the lap time of
%                          2 laps on track 2, direction 1, and session 4
%                          lap1 from 0 to 5, lap2 from 17 to 22
%             thetapk is a cell array with length equals to number of run
%                     sessions, in each element, it's a vector containing
%                     time of theta peaks 
%             tbin    is the decoding window length (use 40ms for flickering detection)
%             overlap is the decoding window overlap (ranging from 0 to 1)
% outputs:  detdecode is a structure with dimension (direction,track,session)
%                     detdecode(idir,tnum,is).seg(seslap,iseg) containing
%                     decoding results on direction idir, track tnum, and
%                     session is. Where we have 5 segments, and several
%                     laps on that track. 
%                     Segments are : 1. 50cm start stationary segment 
%                                    2. 50cm pre-detour mobile segment 
%                                    3. 150cm U shape detour segment 
%                                    4. 50cm control segment 
%                                    5. 50cm end stationary segment 
%                     detdecode(idir,tnum,is).seg(seslap,iseg).pdf is the
%                     decoding probability with D1 being spatial bins, and
%                     D2 being time bins
%                     detdecode(idir,tnum,is).seg(seslap,iseg).tbin specify
%                     the time bins (n*2 matrix as start and end of time bins)
%                     detdecode(idir,tnum,is).seg(seslap,iseg).tbin specify
%                     the spatial mesh 
%                     detdecode(idir,tnum,is).seg(seslap,iseg).nact is the
%                     number of active pyr cells in this bin
%                     detdecode(idir,tnum,is).seg(seslap,iseg).vel is the
%                     averaged run velocity of this time bin
%         mscotpbdir  is an array of averaged decoding probabilities on
%                     pre-detour mobile and control segments,
%                     D1 sessions, only detour session has value, 
%                     D2 directions, D3 pre-detour mobile and control pb
%
% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com


%% preprocess
% exclude int neurons in this analysis
[~,~,ClF,PlFields,~,~,~] = CA_ExcludeInt(M1,M2,ClF,PlFields);
% order tracks by linear pos
ses = Detour_Ordertracks(ses);
% get segment length for each detour track
rplens = Detour_GetDetourSegLen(ses,[2,4]);
% remove redundant time stamps in mazel
mazel = uniquemazeltime(mazel);

%% define parameters
dettra = [2,4]; % detoured tracks
detsesnum = [2,3]; % detour sessions
dirs = [-1,1];  % direction 1 has descending linear position, direction 2 had ascending linear position
lapcri = 0.50; % each lap should cover 50% of the detour track
vcri = 5;

%% decoding detour run laps, we will add current detour, two stationary linear segments
% pre-detour mobile segment, and control segment to the decoder

detdecode = struct;
% D1 direction, D2 track number, D3 run session, .seg decoding results
% on segments
for idir = 1:2
    for it = 1:length(dettra)
        tnum = dettra(it);
        disp(['Decoding D',num2str(idir),' T',num2str(tnum),'...'])
        
        othert = setdiff(dettra,tnum);
        detses = Det_FindDetTSes(tnum,ses);
        % activity session loop
        is = detses;
        % find laps on track of interest
        smazel = mazel{is};
        mazelt = smazel(:,6); % time stamp
        mazellpos = smazel(:,4); % linear position
        mazeltra = smazel(:,3); % track number
        mazelvel = smazel(:,7); % track velocity
        mazelldir = smazel(:,8); % lap direction
        % find lap runs in this session, this track
        posl = ses(is).tralim(tnum,:); % track linear pos limit
        % run direction loop
        dirt = T2T4Laps{tnum,idir,is}; % dir lap time
        indld = mazelldir == dirs(idir);
        indvel = abs(mazelvel) >= vcri;
        seslap = 1;
        % lap loop, for each lap
        for il = 1:size(dirt,1)
            indtra = mazeltra == ses(is).tra_p(tnum);
            indt = mazelt >= dirt(il,1) & mazelt < dirt(il,2);
            allind = indt & indtra & indld & indvel; % index of running for given lap, track and direction
            if sum(allind) == 0
                % no run on this track during this lap,
                continue;
            end
            temppos = mazellpos(allind);  % linear position of running for given lap, track and direction
            temptime = mazelt(allind);
            % make sure the run should cover 15% of the track
            % length
            lran = max(temppos)-min(temppos);
            if lran<lapcri*(max(posl)-min(posl)) || isempty(lran)
                continue;
            end
            % add segments to decoder: two stationary 50cm segment at
            % the start and the end of track; 50 cm pre-detour mobile
            % segment, 150 cm detour segment, 50 cm control segment
            % as the middle segment of the opposite track
            addses = [is,is-1,is,is,is];
            addlpos{1} = [ses(is).tralim(tnum,1), ses(is).tralim(tnum,1) + rplens{tnum}(1)];
            addlpos{2} = [ses(is-1).tralim(tnum,1) + rplens{tnum}(1), ses(is-1).tralim(tnum,2) - rplens{tnum}(2)];
            addlpos{3} = [ses(is).tralim(tnum,1) + rplens{tnum}(1), ses(is).tralim(tnum,2) - rplens{tnum}(2)];
            addlpos{4} = [ses(is).tralim(othert,1) + rplens{othert}(1), ses(is).tralim(othert,2) - rplens{othert}(2)];
            addlpos{5} = [ses(is).tralim(tnum,2) - rplens{tnum}(2), ses(is).tralim(tnum,2)];

            % now we will construct the deocoding template
            allsesmesh = [];
            allsesplf = [];
            meshlencount = 0;
            lposlimit = [];
            for iseg = 1:length(addses)
                [allsesmesh,allsesplf,meshlencount,lposlimit] = AddLpostoDecodeTemplate(addses(iseg),addlpos{iseg},idir,...
                    PlFields,PlFmesh,meshlencount,allsesmesh,allsesplf,lposlimit,ses);
            end
            % do decoding
            trangenow = [min(temptime),max(temptime)];
            ClFnow = ClF(:,is);
            % run the decoding for the concatenate plfield
            sframes = BayePosDecode_LabCode_Lite(allsesplf,allsesmesh,ClFnow,trangenow,[min(allsesmesh) max(allsesmesh)],tbin,overlap);
            % now we need to split the decoding results
            
            for iseg = 1:size(lposlimit,1)
                jssparange = [lposlimit(iseg,1),lposlimit(iseg,2)];
                % spatial bin range for this session
                dprob = sframes.pdf; %D1 is spatial bins, D2 is temporal bins
                spamesh = sframes.spacebin;
                spajsind = spamesh>=jssparange(1) & spamesh<=jssparange(2);
                dprob = dprob(spajsind,:);
                spamesh = spamesh(spajsind);
                % we want to map back to the linear position of the track in the session
                spamesh = spamesh - min(spamesh);
                orilimit1 = addlpos{iseg}(1);
                spamesh = spamesh + orilimit1;
                
                detdecode(idir,tnum,is).seg(seslap,iseg).pdf = dprob;
                detdecode(idir,tnum,is).seg(seslap,iseg).tbin = sframes.tbin;
                detdecode(idir,tnum,is).seg(seslap,iseg).spamesh = spamesh;
                detdecode(idir,tnum,is).seg(seslap,iseg).nact = sframes.nact;
                % get velocity
                ctbin = nanmean(sframes.tbin,2);
                detdecode(idir,tnum,is).seg(seslap,iseg).vel = interp1(mazelt,abs(mazelvel),ctbin);
            end
            seslap = seslap + 1;
        end
    end
end


%% get averaged probability on pre-detour mobile and control segments
vcri = 5; % only consider time bins with velocity larger than this value
actcri = 2; % only consider time bins with at least 2 active cells
mscotpbdir = nan(4,2,2);
% D1 sessions, only detour session has value, D2 directions, D3 pre-mobile
% and control pb

for is = detsesnum
    tnum = setdiff([1,2,3,4],ses(is).tra_p);
    for idir = 1:2
        dcdtmp = detdecode(idir,tnum,is).seg;
        nlap = size(dcdtmp,1);
        pdata = [];
        for ilap = 1:nlap
            goodbin = berow(dcdtmp(ilap,2).nact) >= actcri & berow(dcdtmp(ilap,2).vel) >= vcri;
            scpb = dcdtmp(ilap,2).pdf(:,goodbin);
            otpb = dcdtmp(ilap,4).pdf(:,goodbin);
            scpb = sum(scpb,1);
            otpb = sum(otpb,1);
            dtemp = [scpb',otpb'];
            pdata = cat(1,pdata,dtemp);
        end
        mscotpbdir(is,idir,1) = nanmean(pdata(:,1));
        mscotpbdir(is,idir,2) = nanmean(pdata(:,2));
    end
end



end

function [allsesmesh,allsesplf,meshlencount,lposlimit] = AddLpostoDecodeTemplate(is,lposb,idir,...
    PlFields,PlFmesh,meshlencount,allsesmesh,allsesplf,lposlimit,ses)

trainplf = squeeze(PlFields(:,is,idir,:));
trainmesh = PlFmesh;

lposb1 = lposb(1);
lposb2 = lposb(2);
minlposjs = min(ses(is).tralim(:,1));
maxlposjs = max(ses(is).tralim(:,2));
distol = 0;
lposb1 = max(minlposjs, lposb1-distol);
lposb2 = min(maxlposjs+nanmean(diff(trainmesh)), lposb2+distol);
% the end of maze might not reach a grid size(for example, 598 in 6 grid size)
% add one gird size to include the end of maze
lposb = [lposb1, lposb2];

posind = trainmesh>=lposb(1) & trainmesh<=lposb(2);
meshtemp = trainmesh(posind);
meshdiff = nanmean(diff(trainmesh)); % spatial mesh size
meshtemp = meshtemp-min(meshtemp); % map it to 0 -- segment length
% need to mark the spatial range for each session
plftemp = trainplf(:,posind);
templimit = [meshlencount,meshlencount + max(meshtemp) + meshdiff];
lposlimit = cat(1,lposlimit,templimit);

allsesmesh = cat(2,allsesmesh,meshtemp+meshlencount);
allsesplf = cat(2,allsesplf,plftemp);
meshlencount = meshlencount + max(meshtemp) + meshdiff;

end


