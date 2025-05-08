function pcorr = function_PplCorrXSpaceBin_PrevsDet(ses,M1,M2,PlFields,PlFmesh,ClF)
% function_PplCorrXSpaceBin_PrevsDet this function compute the population
% vector cosine similarity across spatial bin pairs between sessions
% inputs:     ses, mazel, M1, M2, PlFields, ClF, see documents: 
%             "DataStructure"
% output:     pcorr, is a structure with dimension (direction, sesi, sesj)
%             pcorr(idir,is,js).plfmesh is the tuning curve mesh
%             pcorr(idir,is,js).tracks  is the common track in sesi and sesj
%             pcorr(idir,is,js).istralim  is a n*2 matrix which gives the
%                                         linear position range (start to end)
%                                         of common tracks in sesi
%             pcorr(idir,is,js).jstralim  is a n*2 matrix which gives the
%                                         linear position range (start to end)
%                                         of common tracks in sesj
%             pcorr(idir,is,js).crcef     is a n*m matrix gives the
%                                         population vector cosine
%                                         similarity for spaital bin pairs
%                                         in sesi (row) and sesj (column)
%
% Yuchen Zhou 2025 Apr, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

%% preprocess
% exclude int neurons in this analysis
[~,~,~,PlFields,~,~] = CA_ExcludeInt(M1,M2,ClF,PlFields);
% order tracks by linear pos
ses = Detour_Ordertracks(ses);
nsess = length(ses);

%% compute population correlation between detour and pre-detour sessions acoss space bins
pcorr = struct;
% direction loop
for idir = 1:2
    % get the cell activites in reference session
    for is = 1:nsess % reference session loop
        % track limit of session is
        tralimis = ses(is).tralim;
        islim = [min(tralimis(:,1)), max(tralimis(:,2))];
        isind = find(PlFmesh>=islim(1) & PlFmesh<=islim(2));
        trais = ses(is).tra_p;
        for js = 1:nsess
            trajs = ses(js).tra_p;
            [alltra,indi,indj] = intersect(trais,trajs);
            limi = nan(length(alltra),2);
            limj = nan(length(alltra),2);
            % for each session, find the linear pos limit of overlap tracks
            for it = 1:length(alltra)
                % make sure they have the same length
                limi(it,:) = ses(is).tralim(indi(it),:);
                limj(it,:) = ses(js).tralim(indj(it),:);
            end
            
            % track limit of session js
            tralimjs = ses(js).tralim;
            jslim = [min(tralimjs(:,1)), max(tralimjs(:,2))];
            jsind = find(PlFmesh>=jslim(1) & PlFmesh<=jslim(2));
            %% find correlation matrix
            % find plfields on these sessions
            plfi = squeeze(PlFields(:,is,idir,isind));
            plfj = squeeze(PlFields(:,js,idir,jsind));

            % compute ppl correlation for all the spatial pairs as cosine
            % similarity
            corrtemp = nan(length(isind),length(jsind));
            for ip = 1:size(plfi,2)
                for jp = 1:size(plfj,2)
                    corrtemp(ip,jp) = Population_Corre(plfi(:,ip),plfj(:,jp));
                end
            end
            pcorr(idir,is,js).plfmesh = PlFmesh;
            pcorr(idir,is,js).crcef = corrtemp;
            pcorr(idir,is,js).tracks = alltra;
            pcorr(idir,is,js).istralim = limi;
            pcorr(idir,is,js).jstralim = limj;
        end
    end
end

end
