function [pctvar,shfpctvar] = function_PlaceMap_Mvregress(ses,M1,M2,PlFields,ClF,PlFmesh)
% function_PlaceMap_Mvregress
% this function extract tuning curves on the first and the last 50 cm
% segments on detoured tracks, and use a multi-variable regressor to study
% how post-detour tuning curves can be predicted based on pre-detour,
% detour, parallel track, and T1&T3. The result is quantified as the
% percentage residual difference between the full model and the reduced
% model where one regressor is removed. The result is also compared against
% shuffle datasets where the sample ID is shuffled in the mvregress
%
% inputs:     ses, M1, M2, PlFields, ClF, PlFmesh, see documents: "DataStructure"
%                     
% outputs:   pctvar,  is a 1*5 vector contain the percentage residual
%                     difference in the mvregress to predict post-detour 
%                     tuning curve if we remove:
%                     1. pre-detour tuning curve; 2. detour tuning curve;
%                     3. parralel track tuning curve; 4. T1 tuning curve;
%                     5. T3 tuning curve
%                     from the model
%            shfpctvar,  is a n*5 matrix, with row being shuffles and
%                        column being regressors
%
% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

%% preprocess and setting parameters
% exclude int neurons in this analysis
[~,~,~,PlFields,~,~] = CA_ExcludeInt(M1,M2,ClF,PlFields);
% order tracks by linear pos
ses = Detour_Ordertracks(ses);
% get segment length for each detour track
rplens = Detour_GetDetourSegLen(ses,[2,4]);

dettra = [2,4]; % detoured tracks 
klint = [1,3]; % kept linear tracks
segdof = 50; % for each segment, we will interpolate to make sure the tuning curve
% has the length of 50, so we will have the same vector length
shft = 500; % number of sample ID shuffle in mvregress
%% get tuning curves on kept segments on tracks
allplf = cell(1,6); 
% in allplf, we will have tuning curves for
% 1.Pre-detour, 2.Det, 3.Post-detour, 4.Parallel, 5.T1, 6.T3 tracks
% we only track the first and the last 50cm segments on tracks, and
% concatenate those tuning curves. We do this because there is no direct
% correpsondence between the 150cm U shape detour segment and 50cm removed
% or reversal segments.
% in each element of allplf, it's a matrix with
%  D1 being observations (cell * detour tracks * direction)
%  D2 being spatial bins (dimension of the data)

% find tuning curve on the first and the last 50 cm segments on detoured tracks
% in pre-detour, detour, and post detour sessions
for icat = 1:3
    % cat 1,2,3 are pre-detour,detour,post-detour sessions
    for idir = 1:2
        for it = dettra
            detses = Det_FindDetTSes(it,ses); 
            is = detses + icat-2;
            % get track limit and range of linear segment, rescale
            % them to have segdof bins
            tralim = ses(is).tralim(it,:);
            cornertol = 0;
            seg1 = [tralim(1)+cornertol,tralim(1) + rplens{it}(1)];
            seg2 = [tralim(2) - rplens{it}(2),tralim(2)-cornertol];
            idx1 = idxinrange(PlFmesh,seg1);
            idx2 = idxinrange(PlFmesh,seg2);
            plfm1 = PlFmesh(idx1);
            plfm2 = PlFmesh(idx2);
            
            % interpolate to make sure they have the same length
            newplfm1 = linspace(plfm1(1),plfm1(end),segdof);
            newplfm2 = linspace(plfm2(1),plfm2(end),segdof);
            plf1 = squeeze(PlFields(:,is,idir,idx1));
            plf2 = squeeze(PlFields(:,is,idir,idx2));
            newplf1 = interp1(plfm1,plf1',newplfm1);
            newplf2 = interp1(plfm2,plf2',newplfm2);
            newplf1 = newplf1';
            newplf2 = newplf2';

            nowplf = cat(2,newplf1,newplf2);
            % concatenate all the obeservations, D1 are observations (cell* detour tracks * direction)
            % D2 are spatial bins (dimension of the data)
            allplf{icat} = cat(1,allplf{icat},nowplf);
        end
    end
end

% find tuning curve on the first and the last 50 cm segments on the parallel tracks
% during detour session
for idir = 1:2
    for it = dettra
        detses = Det_FindDetTSes(it,ses);
        ot = setdiff(dettra,it);
        is = detses;
        % get track limit and range of linear segment, rescale
        % them to have segdof bins
        tralim = ses(is).tralim(ot,:);
        seg1 = [tralim(1)+cornertol,tralim(1) + rplens{ot}(1)];
        seg2 = [tralim(2) - rplens{ot}(2),tralim(2)-cornertol];
        idx1 = idxinrange(PlFmesh,seg1);
        idx2 = idxinrange(PlFmesh,seg2);
        plfm1 = PlFmesh(idx1);
        plfm2 = PlFmesh(idx2);
        
        % interpolate to make sure they have the same length
        newplfm1 = linspace(plfm1(1),plfm1(end),segdof);
        newplfm2 = linspace(plfm2(1),plfm2(end),segdof);
        plf1 = squeeze(PlFields(:,is,idir,idx1));
        plf2 = squeeze(PlFields(:,is,idir,idx2));
        newplf1 = interp1(plfm1,plf1',newplfm1);
        newplf2 = interp1(plfm2,plf2',newplfm2);
        newplf1 = newplf1';
        newplf2 = newplf2';

        nowplf = cat(2,newplf1,newplf2);
        % concatenate all the obeservations, D1 are observations (cell* detour tracks * direction)
        % D2 are spatial bins (dimension of the data)
        allplf{4} = cat(1,allplf{4},nowplf);
    end
end

% find tuning curve on the first and the last 50 cm segments on the T1 and
% T3 during detour session
for idir = 1:2
    % we still need this detour track loop, as for each detoured segment, we
    % want to explain tuning curve from T1 & T3
    for it = dettra
        for jt = 1:length(klint)
            is = detses;
            % get track limit and range of linear segment, rescale
            % them to have segdof bins
            tralim = ses(is).tralim(klint(jt),:);
            seg1 = [tralim(1)+cornertol,tralim(1) + 50];
            seg2 = [tralim(2) - 50,tralim(2)-cornertol];
            idx1 = idxinrange(PlFmesh,seg1);
            idx2 = idxinrange(PlFmesh,seg2);
            plfm1 = PlFmesh(idx1);
            plfm2 = PlFmesh(idx2);
            
            % interpolate to make sure they have the same length
            newplfm1 = linspace(plfm1(1),plfm1(end),segdof);
            newplfm2 = linspace(plfm2(1),plfm2(end),segdof);
            plf1 = squeeze(PlFields(:,is,idir,idx1));
            plf2 = squeeze(PlFields(:,is,idir,idx2));
            newplf1 = interp1(plfm1,plf1',newplfm1);
            newplf2 = interp1(plfm2,plf2',newplfm2);
            newplf1 = newplf1';
            newplf2 = newplf2';

            nowplf = cat(2,newplf1,newplf2);
            % concatenate all the obeservations, D1 are observations (cell* detour tracks * direction)
            % D2 are spatial bins (dimension of the data)
            allplf{4+jt} = cat(1,allplf{4+jt},nowplf);
        end
    end
end

%% build Multivariate linear regression full model
allobs = size(allplf{1},1);
X = cell(1,allobs); % this are the regressors 
% in X we have the following terms:
% 1. intercept, 2. pre-detour tuning curve, 3. detour tuning curve,
% 4. parallel track tuning curve, 5. T1 tuning curve, 6. T3 tuning curve,
% 7. averaged tuning turve
% The averaged tuning turve illustrate the spatial preference rather than 
% single cell tuning properties, for example, higher rate at corners
averate = nanmean(allplf{3},1);
DOF = size(allplf{1},2);

for iob = 1:allobs
   X{iob} =  [ones(DOF,1),allplf{1}(iob,:)',allplf{2}(iob,:)',...
       allplf{4}(iob,:)',allplf{5}(iob,:)',allplf{6}(iob,:)',becolumn(averate)];
       % intercept,pre,detour,parallel T,T1,T3,averate
end

% we need to predict post-detour tuning curve
Ynow = allplf{3};
[~,~,resi,~,~] = mvregress(X,Ynow,'algorithm','cwls');
% get the model residual
fullvar = norm(resi(:));

% build reduced model, get model residual, compare with sample ID shuffle
effectele = [2,3,4,5,6]; 
% we will remove each regressor once at a time, and see how much it
% contribute to the model residual
% we will remove pre,detour,parallel,T1,T3
pctvar = nan(1,5);
% this is the pct residual related to each regressor
for iele = 1:length(effectele)
    Xtmp = cell(1,allobs);
    for iob = 1:allobs
        % remove that regressor
        nowele = setdiff(effectele,effectele(iele));
        Xtmp{iob} = X{iob}(:,[1,nowele,7]);
    end
    [~,~,resinow,~,~] = mvregress(Xtmp,Ynow,'algorithm','cwls');
    % get the model residual
    resivar = norm(resinow(:));
    % get improvement from this regressor
    pctvar(effectele(iele)-1) = (resivar-fullvar)./fullvar*100;
end

% we will do a sample id shuffle, where we shuffle the sample id in post
% session
shfpctvar = nan(shft,5);
for ish = 1:shft
    shfplf = Ynow(randperm(size(Ynow,1)),:);
    [~,~,resinow,~,~]= mvregress(X,shfplf,'algorithm','cwls');
    fullshfvar = norm(resinow(:));
    for iele = 1:length(effectele)
        Xtmp = cell(1,allobs);
        for iob = 1:allobs
            % remove that regressor
            nowele = setdiff(effectele,effectele(iele));
            Xtmp{iob} = X{iob}(:,[1,nowele,7]);
        end
        % get the model residual
        [~,~,resinow,~,~] = mvregress(Xtmp,shfplf,'algorithm','cwls');
        resivar = norm(resinow(:));
        shfpctvar(ish,effectele(iele)-1) = (resivar-fullshfvar)./fullshfvar*100;
    end
end

end

