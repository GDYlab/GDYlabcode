function [PCCG,tbbias] = function_DetourEarlylap_CCG_CorrXTimeScale(ses,mazel,M1,M2,PlFields,ClF,seslaps)
% function_DetourEarlylap_CCG_CorrXTimeScale
% this function first compute CCG across pyr pairs in based on spikes during
% detour session on detour segments for the first several laps. Then it
% will filter the CCG into behavioral time scale and theta time scale, and
% compute the temporal bias (as center of mass in a given time range)
% at these time scales.
% inputs:     ses, mazel, M1, M2, PlFields, ClF, see documents: 
%             "DataStructure"
%             seslaps is a cell array with length of run sessions.
%                      for each element, it's a vector defines the edges of
%                      lap time
%             for example, seslaps{3} = [0 5 17 22 45] define the lap time of
%                          4 in session 3
%                          lap1 from 0 to 5, lap2 from 5 to 17, lap3 from
%                          17 to 22, lap4 from 22 to 45
% outputs:   PCCG,  it's a structure with dimension (direction, session, lap number)
%                   we will only have values in detour sessions
%            PCCG(idir,is,ilap).lapt  is a vector gives the start and
%                                     end time of this lap on the detour segment
%            PCCG(idir,is,ilap).clf   is a cell array with length equals to
%                                     total number of pyr units, in each
%                                     element, it has the similar
%                                     structure with ClF, but only contains
%                                     spike of this cell on detour segment
%                                     in give lap with given run direction
%                                     and exceed run velocity criteria
%            PCCG(idir,is,ilap).ccgt  is the time bin for CCG
%            PCCG(idir,is,ilap).ccg(ip,jp,:) is the CCG for pair ip and jp
%                                            if ip is leading jp, we will
%                                            get CCG biased towards positive, 
%                                            otherwise, we will get CCG biased towards negative                  
%            PCCG(idir,is,ilap).ccgbh(ip,jp,:)  is the behavioral time scale CCG 
%            PCCG(idir,is,ilap).ccgbht          is the time bin for behavioral time scale CCG 
%            PCCG(idir,is,ilap).ccgthe(ip,jp,:) is the theta time scale CCG 
%            PCCG(idir,is,ilap).ccgthet         is the time bin for theta time scale CCG 
%
%            tbbias,  is a cell array with length equals to number of laps,
%                     in each element, it's a n*2 matrix which pool 
%                     CCG bias from cell pairs, detour sessions, and run
%                     directions. Rows are all the samples, columns are CCG
%                     temporal bias at theta or behavioral time scale
%                     we can directly get correlation between bias across
%                     CCG time scale for a given lap by run the command
%                     [r,p] = corr(tbbias{ilap}(:,1),tbbias{ilap}(:,2))
%
% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

%% setting parameters
nsess = length(ses); % number of run sessions
spkncri = 5; % spk number criteria
dirs = [-1,1]; % direction 1 has descending linear position, direction 2 had ascending linear position
vcri = 5; % only run speed over this criteria will be considered
ccgfre = 1000; % we will compute CCG at frequency of 1000 Hz (1 ms bin)
bhlpf = 1.5; % behaviour time scale CCG  low pass filter under 1.5 Hz
thelpf = 30; % theta time scale CCG  low pass filter under 30 Hz
detses = [2,3]; % detour sessions in detour project
nlap = 3; % we will compute CCG for the first 3 laps

pkcri = 0.2; % in computing CCG, we only consider CCG from pair with peak in CCG exceeding this value (Hz)
comrange = 0.2; % +- as the time range to compute center of mass of CCG (as our temporal bias) for theta scale

%% preprocess
% exclude int neurons in this analysis
[~,~,ClF,~,nc,pyrid] = CA_ExcludeInt(M1,M2,ClF,PlFields);
% order tracks by linear pos
ses = Detour_Ordertracks(ses);
mazel = uniquemazeltime(mazel);
% get segment length for each detour track
rplens = Detour_GetDetourSegLen(ses,[2,4]);

%% compute CCG for cell pairs
PCCG = struct;
for is = detses
    % get tracks on this session
    tra_p = ses(is).tra_p;
    dettra = setdiff([1,2,3,4],tra_p);
    smazel = mazel{is};
    mazelt = smazel(:,6); % time from pos recording
    mazellpos = smazel(:,4); % linear position
    it = dettra;
    detlim = [ses(is).tralim(it,1) + rplens{it}(1),...
        ses(is).tralim(it,2) - rplens{it}(2)];
    lapt = seslaps{is}; % this is a vector contain edges for laps
    
    % loop over lap and directions
    for ilap = 1:nlap
        for idir = 1:2
            disp(['Processing session',num2str(is),' lap',num2str(ilap),' direction',num2str(idir),'...'])
            PCCG(idir,is,ilap).pyrid = pyrid;

            indtra = idxinrange(mazellpos,detlim);
            indt = idxinrange(mazelt,lapt([ilap,ilap+1]));
            allind = indt & indtra;% index of running for given lap, track and direction
            
            laptratime_ori = mazelt(allind);
            laptime = [min(laptratime_ori), max(laptratime_ori)];
            
            % predefine clf where we will record spikes in given lap and
            % track, direction during run 
            clfs = cell(nc,1);
            for ic = 1:nc
                clfnow = clfs{ic};
                clftemp = ClF{ic,is};
                clft = clftemp(:,1); % spike time
                clfpos = clftemp(:,2); % spike lpos
                clfveldir = clftemp(:,5);  % spike vel dir
                clfvel = abs(clftemp(:,3)); % run vel at spike
                % cluster timestamps within lap run and exceed velocity
                % criteria
                tidx = idxinrange(clft,laptime);
                velidx = clfvel > vcri;
                posidx = idxinrange(clfpos,detlim);
                diridx = clfveldir == dirs(idir);
                idx = tidx & velidx & posidx & diridx;
                cltime = berow(clft(idx));
                if ~isempty(cltime)
                    clfnow = cat(2,clfnow,cltime);
                    clfs{ic} = clfnow;
                end
            end
            
            PCCG(idir,is,ilap).lapt = laptime;
            PCCG(idir,is,ilap).clf = clfs;
            
            disp('Getting spikes completed, computing CCG across cell pairs...')
            % look at pyr pairs, get CCGs
            for ip = 1:nc
                pctdisplay(ip,nc,50)
                cltimei = PCCG(idir,is,ilap).clf{ip};
                numspki = length(cltimei);
                for jp = ip:nc
                    cltimej = PCCG(idir,is,ilap).clf{jp};
                    numspkj = length(cltimej);
                    minspk = min([numspki,numspkj]);
                    if minspk < spkncri
                        PCCG(idir,is,ilap).ccg(ip,jp,:) = zeros(2001,1);
                        ccgt = 1000*(-ccgfre:ccgfre)*1/ccgfre/1000;
                        PCCG(idir,is,ilap).ccgt = ccgt;
                        continue
                    end
                    Time=[berow(cltimei), berow(cltimej)];   %concatenate spike time of both cells
                    Group=2*ones(length(Time),1);    %create a vectore of correct size filled with two WARNING ones won t work
                    Group(length(cltimei)+1:end)=3;    %filled up the times corresponding with the cluster l with threes
                    [ccg,ccgt,~]=CCG_zyc(Time,Group,1/ccgfre,ccgfre,1000);  %call the CCG function
                    PCCG(idir,is,ilap).ccg(ip,jp,:) = ccg(:,1,2);
                    PCCG(idir,is,ilap).ccgt = ccgt;
                    % we have ip vs jp ccg, if ip is leading jp, we will
                    % get CCG biased towards positive, otherwise, we
                    % will get CCG biased towards negative
                end
            end
        end
    end
end

%% filter ccg in corresponding bands
disp('Filtering CCG in behavioral and theta scales ...')
% behavior filter
BehFilt=fir1(400,bhlpf/(ccgfre/2));
% theta filter
TheFilt=fir1(100,thelpf/(ccgfre/2));
% loop over cell pair
for ip = 1:nc
    progress = round(ip/nc*100);
    if rem(progress,20) == 0
        disp(['Progress ',num2str(progress),'%...'])
    end
    for jp = ip:nc
        ccgmax = zeros(nsess,2,nlap);
        for is = detses
            for ilap = 1:nlap
                for idir = 1:2
                    ccgt = PCCG(idir,is,ilap).ccgt;
                    tidx = ccgt >= -0.05 & ccgt <= 0.05;
                    ccgd = squeeze(PCCG(idir,is,ilap).ccg(ip,jp,:));
                    smccgl = Filter0(TheFilt,ccgd);
                    smccgl = smccgl(tidx);
                    ccgmax(is,idir,ilap) = max(smccgl(:));
                end
            end
        end
        % loop over session, direction, and laps
        for is = detses
            for idir = 1:2
                for ilap = 1:nlap 
                    % raw ccg and time
                    ccgd = squeeze(PCCG(idir,is,ilap).ccg(ip,jp,:));
                    ccgt = PCCG(idir,is,ilap).ccgt;
                    % filter into behavioral time scale
                    smccgl = Filter0(BehFilt,ccgd);
                    PCCG(idir,is,ilap).ccgbh(ip,jp,:) = smccgl;
                    PCCG(idir,is,ilap).ccgbht = ccgt;
                    % filter into theta time scale
                    smccgl = Filter0(TheFilt,ccgd);
                    tidx = ccgt >= -0.05 & ccgt <= 0.05;
                    smccgl = smccgl(tidx);
                    ccgt = ccgt(tidx);
                    PCCG(idir,is,ilap).ccgthe(ip,jp,:) = smccgl;
                    PCCG(idir,is,ilap).ccgthet = ccgt;
                end
            end
        end
    end
end

%% get CCG bias at different scales
disp('Computing CCG temporal bias in behavioral and theta scales ...')
tbbias = cell(nlap,1);
% loop over cell pair
for ip = 1:(nc-1)
    for jp = (ip+1):1:nc
        % loop over session, direction, and laps
        for is = detses
            ccgmax = zeros(2,nlap);
            for ilap = 1:nlap
                for idir = 1:2
                    smccgl = PCCG(idir,is,ilap).ccgthe(ip,jp,:);
                    ccgmax(idir,ilap) = max(smccgl(:));
                end
            end
            allmax = sort(ccgmax(:),'descend');
            if allmax(1) <= pkcri
                continue
            end
            
            for ilap = 1:nlap
                for idir = 1:2
                    % get theta scale pk
                    smccgl = squeeze(PCCG(idir,is,ilap).ccgthe(ip,jp,:));
                    smmax = max(smccgl);
                    if  smmax <= pkcri
                       continue 
                    end

                    % get raw ccg com at -1 to 1s and -200ms to 200ms com
                    rawt = PCCG(idir,is,ilap).ccgt;
                    rawccg = squeeze(PCCG(idir,is,ilap).ccg(ip,jp,:));
                    behidx = idxinrange(rawt,[-1 1]);
                    thetaidx = idxinrange(rawt,[-comrange comrange]);
                    % compute center of mass in different time range
                    behcom = berow(rawt(behidx))*becolumn(rawccg(behidx))/sum(rawccg(behidx));
                    thetacom = berow(rawt(thetaidx))*becolumn(rawccg(thetaidx))/sum(rawccg(thetaidx));
                    tbbias{ilap,1} = cat(1,tbbias{ilap,1},[thetacom,behcom]);
                end
            end
        end
    end
end

end
