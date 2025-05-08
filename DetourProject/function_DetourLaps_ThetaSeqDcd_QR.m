function [lapdcd,tsig] = function_DetourLaps_ThetaSeqDcd_QR(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,seslaps,thetapk,tbin,overlap)
% function_DetourLaps_ThetaSeqDcd_QR
% this function run Bayesian decoding on detour tracks during detour run,
% correct the linear position based on animal's actual position, and
% compute the quadrant ratio for all the theta cycles grouped by laps
% inputs:     ses, mazel, M1, M2, PlFields, PlFmesh, ClF, see documents: 
%             "DataStructure"
%             seslaps is a cell array with length equals to number of run
%                     sessions, in each element, it's a vector containing
%                     lap time 
%             for example, seslaps{2} = [0 5 12 30] define the lap time of
%                          3 laps in session 2, lap1 from 0 to 5, lap2 from
%                          5 to 12, and lap3 from 12 to 30.
%             thetapk is a cell array with length equals to number of run
%                     sessions, in each element, it's a vector containing
%                     time of theta peaks 
%             tbin    is the decoding window length (use 20ms for theta sequence detection)
%             overlap is the decoding window overlap (ranging from 0 to 1)
% outputs:    lapdcd  is a structure with length equals to number of run
%                     sessions. lapdcd(1).Dlap(ilap,idir) containing
%                     decoding results for theta cycles with movement (mean velocity > 10)
%                     during ilap (1 to total number of laps) and in run direction idir (1 or 2)
%                     lapdcd(1).Dlap(ilap,idir).thetatime is a vector of
%                     middle time of all the decoded theta cycle (defined
%                     as peak to peak) in the lap
%                     lapdcd(1).Dlap(ilap,idir).allthetapb is a
%                     spacebin*timebin*thetacycles array containing
%                     decoding results for theta cycles in this lap, the
%                     space is corrected, so 0 always means animals current
%                     location across time bins
%                     lapdcd(1).Dlap(ilap,idir).allthetapb is a
%                     spacebin*timebin array averaging decoding results
%                     over all the theta cycles
%              tsig   is a structure with D1 being run sesions, D2 being
%                     run directions
%                     tsig(is,idir).lap(ilap).dataqr contains the quadrant
%                     ratio of all the theta cycles in session is,
%                     direction idir, and lap ilap
%                     tsig(is,idir).lap(ilap).shfqr is a array with number 
%                     of elements = shufflenumber*thetacycles
%                     quadrant ratios from time bin shuffle.
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
detses = [2,3]; % detour sessions
dirs = [-1,1]; 
% direction 1 means descending linear position, direction 2 means ascending
% linear position
vcri = 10; 
% velocity criteria of theta cycle, in the theta cycle, the mean velocity should larger than this value
thetaext = 0.04; % we will average over around +-40 ms
% was +-50ms, but some theta are fast, and there are phase reset effect at
% boundaries
pbbinextend = 20; % we will extend +-20 position bins from center of position bin
shfnum = 10; % for each theta cycle, we do 10 time bin shuffle to get time bin shuffle quadrant ratio

%% decoding laps on detour, correct for animal's actual position
lapdcd = struct;
% D1 run sessions,only s2s3 have values
for is = detses
    % get spikes and detour segment linear position range
    ClFnow = ClF(:,is);
    dettra = setdiff([1,2,3,4],ses(is).tra_p);
    detseglpos = ses(is).tralim(dettra,:);
    detseglpos(1) = detseglpos(1) + rplens{dettra}(1);
    detseglpos(2) = detseglpos(2) - rplens{dettra}(2);
    
    % get behavior data
    smazel = mazel{is};
    mazelt = smazel(:,6);
    mazellpos = smazel(:,4); % track velocity
    % smooth lpos
    mazellpos = smoothdata(mazellpos,'movmean',10);
    
    % load laps and theta cycles for this session
    lapnow = seslaps{is};
    thetapknow = thetapk{is};
    for ilap = 1:length(lapnow)-1
        laptimenow = [lapnow(ilap),lapnow(ilap+1)];
        timeidx = idxinrange(mazelt,laptimenow);
        lposidx = idxinrange(mazellpos,detseglpos);
        for idir = 1:2
            % pre-define decoding probability and theta cycle time
            allthetapb = [];
            allthetatime = [];
            
            % find time when animal run on the track in this lap
            % we will enforce the direction criteria later when we look at
            % each theta cycle
            allidx = timeidx & lposidx;
            
            % get continuous run on the detour segment in this lap
            % longer than 0.5 second
            trange = GetEpochwithCrioverLen(mazelt,allidx,0.5);
            if isempty(trange)
                lapdcd(is).Dlap(ilap,idir).avethetapb = [];
                continue
            end
            plftmp = squeeze(PlFields(:,is,idir,:));
            tempmesh = PlFmesh;
            
            % add tracks to decoder
            addses = is;
            addlpos{1} = ses(is).tralim(dettra,:);
            
            % now we will construct the deocoding template
            allsesmesh = [];
            allsesplf = [];
            meshlencount = 0;
            lposlimit = [];
            for iseg = 1:length(addses)
                [allsesmesh,allsesplf,meshlencount,lposlimit] = AddLpostoDecodeTemplate_Early(addses(iseg),addlpos{iseg},...
                    plftmp,tempmesh,meshlencount,allsesmesh,allsesplf,lposlimit,ses);
            end
            % do decoding with detoured track during continous runs on the
            % track in this lap
            sframes = BayePosDecode_LabCode_Lite(allsesplf,allsesmesh,ClFnow,trange,[min(allsesmesh) max(allsesmesh)],tbin,overlap);
            
            % process each continuous run
            for iep = 1:length(sframes)
                % spatial bin range for this session
                dprob = sframes(iep).pdf; %D1 is spatial bins, D2 is temporal bins
                spamesh = sframes(iep).spacebin;
                % map back to original linear position
                spamesh = spamesh - min(spamesh);
                orilimit1 = addlpos{1}(1);
                spamesh = spamesh + orilimit1;
                meshsize = nanmean(diff(spamesh));
                
                dcdt =  sframes(iep).tbin;
                alltimec = nanmean(dcdt,2);
                dcddt = nanmean(diff(alltimec));
                
                trangenow = trange(iep,:);
                % get averaged that sequence
                % find theta cycles in this interval
                thetanowidx = idxinrange(thetapknow,trangenow);
                thetaitv = thetapknow(thetanowidx);
                % theta period cycle
                for ith = 1:length(thetaitv)-1
                    % get the theta time range
                    thetast = thetaitv(ith);
                    thetaend = thetaitv(ith+1);
                    thetamid = mean([thetast,thetaend]);
                    % find linear position of this theta, should be in the
                    % detour segment
                    midtlpos = interp1(mazelt,mazellpos,thetamid);
                    lposidx = idxinrange(midtlpos,detseglpos);
                    
                    % this theta period should travel enough distance based
                    % on vel cri, and also match the direction
                    stlpos = interp1(mazelt,mazellpos,thetast);
                    endlpos = interp1(mazelt,mazellpos,thetaend);
                    thetadis = (endlpos-stlpos)*dirs(idir);
                    % if D1, stlpos should larger than end
                    % if D2, stlpos should smaller than end
                    discri = (thetaend-thetast)*vcri;
                    disidx = thetadis>=discri;
                    
                    if lposidx && disidx
                        % in this case, we will get the decoded probability
                        % find decoded time bins in this theta cycle
                        thetacenter = (thetast+thetaend)/2;
                        npoints = round(2*thetaext./dcddt);
                        
                        dcdtidx = idxinrange(alltimec,[thetacenter-thetaext*2,...
                            thetacenter+thetaext*2]);
                        % extend more to ensure accurate interpolation
                        pbnow = dprob(:,dcdtidx);
                        tbintmp = alltimec(dcdtidx);
                        % redefine the 2d interplation time mesh
                        tnewmesh = linspace(thetacenter-thetaext,...
                            thetacenter+thetaext,npoints);
                        % first, interpolate probability to the new
                        % time mesh
                        tmeshpb = zeros(size(pbnow,1),length(tnewmesh));
                        for ipos = 1:size(tmeshpb,1)
                            pbtmp=  interp1(tbintmp,pbnow(ipos,:),tnewmesh);
                            tmeshpb(ipos,:) = becolumn(pbtmp);
                        end
                        
                        % second, interpolate to relative pos related
                        % to animal's true position
                        % get the animal's actual position at this new
                        % time mesh
                        lposnewtbin = interp1(mazelt,mazellpos,tnewmesh);
                        
                        % for each time step, interpolation the
                        % position to the new space mesh
                        thetapb = zeros(2*pbbinextend+1,length(tnewmesh));
                        for itbin = 1:length(tnewmesh)
                            lposnow = lposnewtbin(itbin);
                            % redefine the 2d interplation spatial mesh
                            % current lpos+-10 space bins
                            spanewmesh = linspace(lposnow-meshsize*pbbinextend,...
                                lposnow+meshsize*pbbinextend,2*pbbinextend+1);
                            pbtmp = interp1(spamesh,tmeshpb(:,itbin),spanewmesh);
                            thetapb(:,itbin) = becolumn(pbtmp);
                        end
                        
                        % concatenate array to record this theta cycle
                        if isempty(allthetapb)
                            allthetapb = thetapb;
                            allthetatime = thetamid;
                        else
                            allthetapb = cat(3,allthetapb,thetapb);
                            allthetatime = cat(1,allthetatime,thetamid);
                        end
                    end
                end
            end
            lapdcd(is).Dlap(ilap,idir).thetatime = allthetatime;
            lapdcd(is).Dlap(ilap,idir).thetapb = allthetapb;
            % in theta pb, D1 space mesh, D2 time mesh
            lapdcd(is).Dlap(ilap,idir).avethetapb = nanmean(allthetapb,3);
        end
    end
end

%% compute quadrant ratio
tsig = struct;
for is = detses
    for idir = 1:2
        lapdcdnow = lapdcd(is).Dlap(:,idir);
        nlap = length(lapdcdnow);
        % get quadrant ratio of theta cycles in each lap
        for ilap = 1:nlap
            pbnow = lapdcdnow(ilap).thetapb;
            [dataqr,shfqr] = GetThetaQRSig_DShf(pbnow,idir,shfnum);
            dataqr = dataqr(~isnan(dataqr));
            shfqr = shfqr(~isnan(shfqr));
            tsig(is,idir).lap(ilap).dataqr = dataqr;
            tsig(is,idir).lap(ilap).shfqr = shfqr;
        end
    end
end

end

function [allsesmesh,allsesplf,meshlencount,lposlimit] = AddLpostoDecodeTemplate_Early(is,lposb,...
    PlFields,PlFmesh,meshlencount,allsesmesh,allsesplf,lposlimit,ses)

trainplf = PlFields;
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

function [dataqr,shfqr] = GetThetaQRSig_DShf(allpb,idir,nshf)

dataqr = nan(1,size(allpb,3));
shfqr = nan(nshf,size(allpb,3));

for ipb = 1:size(allpb,3)
    pbnow = squeeze(allpb(:,:,ipb));
    % get probability in each quadrant
    [Q1pb,Q2pb,Q3pb,Q4pb] = GetQPb(pbnow);
    Q13 = cat(1,becolumn(Q1pb(:)),becolumn(Q3pb(:)));
    Q13 = nansum(Q13);
    Q24 = cat(1,becolumn(Q2pb(:)),becolumn(Q4pb(:)));
    Q24 = nansum(Q24);
    if idir == 1
        qrval = (Q24 - Q13)/(Q24 + Q13);
    else
        qrval = (-Q24 + Q13)/(Q24 + Q13);
    end
    dataqr(ipb) = qrval;
    
    % compute quadrant ratio from time bin shuffle
    for ishf = 1:nshf
        shfpn = pbnow(:,randperm(size(pbnow,2)));
        [Q1pb,Q2pb,Q3pb,Q4pb] = GetQPb(shfpn);
        Q13 = cat(1,becolumn(Q1pb(:)),becolumn(Q3pb(:)));
        Q13 = nansum(Q13);
        Q24 = cat(1,becolumn(Q2pb(:)),becolumn(Q4pb(:)));
        Q24 = nansum(Q24);
        if idir == 1
            qrval = (Q24 - Q13)/(Q24 + Q13);
        else
            qrval = (-Q24 + Q13)/(Q24 + Q13);
        end
        shfqr(ishf,ipb) = qrval;
    end
end
end


function [Q1pb,Q2pb,Q3pb,Q4pb] = GetQPb(pb)
% in theta pb, D1 space mesh, D2 time mesh
[r,c] = size(pb);
rowbin = floor(r/2)-2;
colbin = floor(c/2);

Q1pb = pb(end-rowbin+1:end,end-colbin+1:end);
Q2pb = pb(end-rowbin+1:end,1:colbin);
Q3pb = pb(1:rowbin,1:colbin); 
Q4pb = pb(1:rowbin,end-colbin+1:end);

end
