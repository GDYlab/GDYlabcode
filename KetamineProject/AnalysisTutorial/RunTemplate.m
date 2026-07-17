% RunTemplate.m
% This script demonstrates the key analyses of the ketamine project using a
% single-animal sample dataset (post-Run2-split; see ../DataStructure.md).
%
% To run: add the entire "Ketamine" folder and its subfolders to the MATLAB
% path. Select which analyses to run in the "specify which code to run"
% section by setting the flags to 1 (run) or 0 (skip).
%
% Developed and tested on MATLAB R2017b.
%
% The sample data lives in ../SampleData/. The variable SampleAnimal selects
% which animal's sample to load, so the sample can be swapped without editing
% the analysis code (a SampleData_<Animal>/ folder is expected per animal).
%
% To get access to the sample data, contact the Dragoi lab:
% george.dragoi@yale.edu
%
% Yuchen Zhou, yuchen.zhou@yale.edu, yuchenzhou93@gmail.com

clear; clc; close all

%% locate sample data (each analysis can use a different animal)
% Sample data is packaged per animal under SampleData/<animal>/. There are 
% 5 animals (URat1, URat2, URat3, URat4, XRat1) in the project
% (URat1/URat3 have 3 effective run sessions; the others have 4 — the Run2woBi
%  session is empty for the 3-session animals and is skipped automatically.)

animal_beh      = 'URat1'; % 1. behavior
animal_thetaseq = 'URat1'; % 2. theta sequence by velocity
animal_wakedcd  = 'URat1'; % 3. waking maze representation 
animal_minfo    = 'URat1'; % 4. spike mutual information
animal_embed    = 'URat1'; % 5. spike dimensionality
animal_lfpspec  = 'URat1'; % 6. bispectrum / theta asymmetry
currentScriptPath = mfilename('fullpath');
[currentFolder, ~, ~] = fileparts(currentScriptPath);
% helper: full path to a sample file for a given animal
spath = @(animal,file) [currentFolder, filesep, '..', filesep, 'SampleData', filesep, animal, filesep, file];

%% specify which code to run
runbeh      = 1; % 1. behavior: mean/max velocity, run-vs-track angle, end-approach velocity profile (Fig1)
runthetaseq = 1; % 2. theta sequence by velocity (Fig2)
runwakedcd  = 1; % 3. waking decoding via population-vector cosine-similarity probability (Fig3)
runminfo    = 1; % 4. mutual information of spikes during run (Fig4)
runembed    = 1; % 5. spike dimensionality during run (Fig4)
runlfpspec  = 1; % 6. bispectrum and theta wave asymmetry (Fig5)
% 7. LDS model training and inference (Fig6) is a separate Python tutorial: see python/

%% 1. behavior analysis (Fig1)
% estimated run time ~10 s
if runbeh
    tic
    load(spath(animal_beh,'Behavior.mat'),'mazel','ses','info');
    beh = function_Beh_RunVelocityMeasures(mazel,ses,info);
    nsess = numel(beh.sesstr);

    figure('Name','Behavior: ketamine impact on run velocity','color','w')
    % (a) distance travelled per minute
    subplot(2,2,1)
    bar(1:nsess,beh.dispermin,0.5,'FaceColor',[0.4 0.4 0.4]);
    set(gca,'xtick',1:nsess,'xticklabel',beh.sesstr); xlim([0.4 nsess+0.6])
    ylabel('Distance per minute (cm)'); title('Locomotor output')

    % (b) angle between velocity and track direction (mean +/- sem)
    subplot(2,2,2)
    m = cellfun(@nanmean,beh.trackangle);
    s = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),beh.trackangle);
    hold on; bar(1:nsess,m,0.5,'FaceColor',[0.4 0.4 0.4]);
    errorbar(1:nsess,m,s,'k','linestyle','none','CapSize',6); hold off
    set(gca,'xtick',1:nsess,'xticklabel',beh.sesstr); xlim([0.4 nsess+0.6])
    ylabel('Angle <track,velocity> (deg)'); title('Path straightness')

    % (c) maximum velocity per run epoch (mean +/- sem)
    subplot(2,2,3)
    m = cellfun(@nanmean,beh.maxvel);
    s = cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),beh.maxvel);
    hold on; bar(1:nsess,m,0.5,'FaceColor',[0.4 0.4 0.4]);
    errorbar(1:nsess,m,s,'k','linestyle','none','CapSize',6); hold off
    set(gca,'xtick',1:nsess,'xticklabel',beh.sesstr); xlim([0.4 nsess+0.6])
    ylabel('Max velocity (cm/s)'); title('Peak speed')

    % (d) velocity profile vs distance to track end (approaching / departing)
    subplot(2,2,4); hold on
    colors = lines(nsess);
    hp = gobjects(1,nsess);
    for is = 1:nsess
        ma = cellfun(@nanmean,beh.velapproach(is,:));
        md = cellfun(@nanmean,beh.veldepart(is,:));
        hp(is) = plot(beh.disbinc,ma,'-','color',colors(is,:),'linewidth',1.5); % solid: approaching end
        plot(beh.disbinc,md,'--','color',colors(is,:),'linewidth',1);           % dashed: departing
    end
    xlim([0 25]); ylim([0 45])
    hold off; xlabel('Distance to track end (cm)'); ylabel('Velocity (cm/s)')
    title('Approach (solid) / depart (dashed)'); legend(hp,beh.sesstr,'location','best')
    toc
end

%% 2. theta sequence by velocity (Fig2)
% estimated run time ~3-5 min (Bayesian decoding of run epochs at 4 ms steps)
if runthetaseq
    tic
    load(spath(animal_thetaseq,'Behavior.mat'),'ses','mazel');
    load(spath(animal_thetaseq,'ThetaSeq.mat'),'M1','M2','PlFields','PlFmesh','ClF','thetapk');
    velbin = [5 15 30 120];
    tseq = function_ThetaSeq_ByVelocity(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,thetapk,velbin,0.04);
    nsess = numel(tseq.sesstr); nvb = numel(velbin)-1;

    % grid of averaged theta sequences: rows = velocity bins, cols = sessions
    % (directions pooled). A forward-sweeping diagonal = theta sequence.
    figure('Name','Theta sequence by velocity','color','w')
    clim = 0.03;
    for vb = 1:nvb
        for is = 1:nsess
            % pool the two running directions
            stack = cat(3,tseq.meanpb{is,1,vb},tseq.meanpb{is,2,vb});
            ax = subplot(nvb,nsess,(vb-1)*nsess+is);
            if isempty(stack); axis(ax,'off'); continue; end
            mp = nanmean(stack,3);
            imagesc(ax,tseq.cttmesh([1 end]),tseq.relposmesh([1 end]),mp);
            set(ax,'ydir','normal'); caxis(ax,[0 clim]); ylim(ax,[-40 40]); hold(ax,'on')
            plot(ax,[0 0],ylim(ax),'w--'); plot(ax,xlim(ax),[0 0],'w--'); hold(ax,'off')
            if vb==1; title(ax,tseq.sesstr{is},'fontsize',9); end
            if is==1; ylabel(ax,sprintf('V%d-%d cm/s\nahead(+) cm',velbin(vb),velbin(vb+1))); end
            if vb==nvb; xlabel(ax,'Time (s)'); end
        end
    end
    toc
end

%% 3. waking decoding via population-vector cosine similarity (Fig3)
% estimated run time ~30 s
if runwakedcd
    tic
    load(spath(animal_wakedcd,'Behavior.mat'),'ses','mazel');
    load(spath(animal_wakedcd,'ThetaSeq.mat'),'M1','M2','PlFields','PlFmesh','ClF');
    load(spath(animal_wakedcd,'WakingRep.mat'),'fratime','times');
    rep = function_WakingMazeRep_CosineProb(ses,mazel,times,M1,M2,PlFields,PlFmesh,ClF,fratime,0.02);
    nsess = numel(rep.sesstr);

    % CDF of the per-frame maze-representation log-probability across sessions.
    % Ketamine (Run2Bi) shifts left = weaker waking maze representation.
    figure('Name','Waking maze representation','color','w'); hold on
    colors = lines(nsess); lh = []; ll = {};
    for is = 1:nsess
        d = rep.logpb{is}; d = d(isfinite(d));
        if isempty(d); continue; end
        [y,x] = ecdf(d);
        lh(end+1) = plot(x,y,'-','color',colors(is,:),'linewidth',1.5);  %#ok<SAGROW>
        ll{end+1} = rep.sesstr{is};                                      %#ok<SAGROW>
    end
    hold off; xlabel('Maze representation logPb'); ylabel('CDF')
    title('Waking maze representation (cosine-similarity probability)')
    legend(lh,ll,'location','southeast')
    toc
end

%% 4. mutual information of spikes during run (Fig4)
% estimated run time ~10-30 s
if runminfo
    tic
    load(spath(animal_minfo,'Behavior.mat'),'ses','mazel');
    load(spath(animal_minfo,'ThetaSeq.mat'),'M2','ClF');
    mi = function_RunSpikeMutualInfo(ses,mazel,M2,ClF,0.4,5,1);
    nsess = numel(mi.sesstr);

    figure('Name','Run spike mutual information','color','w')
    % (a) CDF of log normalized MI across sessions; ketamine (Run2Bi) shifts left
    subplot(1,2,1); hold on
    colors = lines(nsess); lh = []; ll = {};
    for is = 1:nsess
        d = mi.lognMI{is};
        if isempty(d); continue; end
        [y,x] = ecdf(d);
        lh(end+1) = plot(x,y,'-','color',colors(is,:),'linewidth',1.5); %#ok<SAGROW>
        ll{end+1} = mi.sesstr{is};                                      %#ok<SAGROW>
    end
    hold off; xlabel('log normalized MI'); ylabel('CDF')
    legend(lh,ll,'location','southeast'); title('Pairwise MI across sessions')

    % (b) example cell-by-cell MI matrix (first non-empty session); the red
    %     lines separate interneurons (first mi.nint cells) from pyramidal cells
    subplot(1,2,2)
    is = find(~cellfun(@isempty,mi.MImat),1);
    imagesc(mi.MImat{is},'AlphaData',~isnan(mi.MImat{is})); axis square; set(gca,'ydir','normal')
    caxis([0 0.2]); colorbar; hold on
    plot(xlim,[mi.nint+0.5 mi.nint+0.5],'r--'); plot([mi.nint+0.5 mi.nint+0.5],ylim,'r--'); hold off
    xlabel('present cell'); ylabel('past cell'); title([mi.sesstr{is},': nMI matrix'])
    toc
end

%% 5. spike dimensionality during run (Fig4)
% estimated run time ~10-30 s
if runembed
    tic
    load(spath(animal_embed,'Behavior.mat'),'ses','mazel');
    load(spath(animal_embed,'ThetaSeq.mat'),'M2','ClF');
    load(spath(animal_embed,'WakingRep.mat'),'times'); % run-lap timing (per session/direction)
    emb = function_Embedding_RunSpikeDim(ses,mazel,times,M2,ClF,0.4,10);
    nsess = numel(emb.sesstr);
    exdir = emb.exdir; % example direction shown in the manuscript

    figure('Name','Run spike dimensionality','color','w')
    % (a) CDF of mu = r2/r1 across sessions for the example direction. Higher
    %     intrinsic dimensionality makes the two nearest-neighbour distances
    %     more similar, so mu concentrates near 1 and the CDF rises faster;
    %     ketamine (Run2Bi) shifts the curve up/left.
    subplot(1,2,1); hold on
    colors = lines(nsess); lh = []; ll = {};
    for is = 1:nsess
        m = emb.mu{exdir,is};
        if isempty(m); continue; end
        [y,x] = ecdf(m);
        lh(end+1) = plot(x,y,'-','color',colors(is,:),'linewidth',1.5); %#ok<SAGROW>
        ll{end+1} = emb.sesstr{is};                                     %#ok<SAGROW>
    end
    hold off; xlabel('\mu = r2/r1'); ylabel('CDF(\mu)')
    legend(lh,ll,'location','southeast'); title('Nearest-neighbour ratio')

    % (b) intrinsic dimensionality per session (mean over the two run
    %     directions); ketamine (Run2Bi) is expected to raise it.
    subplot(1,2,2)
    m = nanmean(emb.dim,1); s = nanstd(emb.dim,0,1)./sqrt(sum(~isnan(emb.dim),1));
    hold on; bar(1:nsess,m,0.5,'FaceColor',[0.4 0.4 0.4]);
    errorbar(1:nsess,m,s,'k','linestyle','none','CapSize',6); hold off
    set(gca,'xtick',1:nsess,'xticklabel',emb.sesstr); xlim([0.4 nsess+0.6])
    ylabel('Intrinsic dimensionality'); title('Spike dimensionality during run')
    toc
end

%% 6. bispectrum and theta wave asymmetry (Fig5)
% estimated run time ~20-60 s
if runlfpspec
    tic
    load(spath(animal_lfpspec,'RunLFP.mat'),'seglfp','fs');
    lfp = function_LFPSpec_RunBispectrum(seglfp,fs,[7 10],[0 120]);
    nsess = numel(lfp.sesstr);

    figure('Name','Run LFP bispectrum / theta asymmetry','color','w',...
        'Units','normalized','Position',[0.25 0.25 0.4 0.4])
    % (a) imaginary part of the theta-band bicoherence for each session: a
    %     non-zero value reflects an asymmetric (saw-tooth) theta wave.
    sesplot = [1,2,4,5];
    for is = 1:nsess
        ax = subplot(2,3,sesplot(is));
        ss = lfp.bic{is};
        if isempty(ss); axis(ax,'off'); continue; end
        im = imag(ss); im(ss==0) = NaN;        % blank the redundant domain
        pcolor(ax,lfp.ff,lfp.ff,im); shading(ax,'interp'); hold(ax,'on')
        plot(ax,[0 60 120],[0 60 0],'--k');    % principal-domain boundary
        caxis(ax,[-0.6 0.6]); xlim(ax,[0 15]); ylim(ax,[0 15]); axis(ax,'square')
        xlabel(ax,'f1 (Hz)'); if is==1; ylabel(ax,'f2 (Hz)'); end
        title(ax,sprintf('%s\nasym=%.3f',lfp.sesstr{is},lfp.asym(is)),'fontsize',9)
        if is == nsess
           cb = addcolorbar; 
           cb.Label.String = 'Asymmetry';
        end
    end
    colormap(jet)

    % (b) theta wave asymmetry across sessions; ketamine (Run2Bi) stands out.
    ax = subplot('position',[0.75 0.25 0.2 0.5]);
    bar(ax,1:nsess,lfp.asym,0.5,'FaceColor',[0.4 0.4 0.4]);
    set(ax,'xtick',1:nsess,'xticklabel',lfp.sesstr); xlim(ax,[0.4 nsess+0.6])
    ylabel(ax,'Theta asymmetry'); title(ax,'Imag bicoherence (\theta)')
    toc
end
