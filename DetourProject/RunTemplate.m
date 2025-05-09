% This script includes code for several analyses presented in the manuscript
% and demonstrates how to apply them using sample data from a single animal.
% To run the script add the entire "Detour" fold and its subfolders into
% Matlab Path. You can select which analyses to run in the section 
% "specify which code to run" by setting those values to 1 and
% setting other to 0. Run all the scripts should take less than 10 min


clear;clc;close all

% Get the path of the current script
currentScriptPath = mfilename('fullpath');
[currentFolder, ~, ~] = fileparts(currentScriptPath);
% get the path of sample data
Datapath = [currentFolder,filesep,'SampleData',filesep];
% to get access to the sample data, please contact the lab: george.dragoi@yale.edu

% get nc and clu2pyr which can be used in several codes
load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2');
[~,~,~,~,nc,~,clu2pyr,~] = CA_ExcludeInt(M1,M2,ClF,PlFields);

%% specify which code to run
runqr = 1; % decoding in earlt detour lap theta cycles, get quadrant ratio distribution
runccgcorr = 1; % CCG temporal bias correlation between time scales in early laps
runpplcorr = 1; % population vector cosine similarity between detour and predetour segments
runplfreg = 1; % place map regression, explain post-detour place map by pre-detour and detour
runrankcorr = 1; % rank order correlation between pre-detour sleep frames and theta cycles
rundetflc = 1; % during detour run, the flickering of pre-detour mobile segment
runrevflc = 1; % during reversal run, the flickering of detour segment

%% code for early lap theta cycle quadrant ratio
if runqr
    % estimated run time ~2min
    tic
    load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2','PlFmesh');
    load([Datapath,'seslaps.mat'],'seslaps');
    load([Datapath,'thetapk.mat'],'thetapk');
    load([Datapath,'ses.mat'],'ses');
    load([Datapath,'mazel.mat'],'mazel');
    tbin = 0.02; % 20 ms decoding window
    overlap = 0; % 0 window overlap
    [lapdcd,tsig] = function_DetourLaps_ThetaSeqDcd_QR(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,seslaps,thetapk,tbin,overlap);
    % plot distibution of quadrant ratio for the first 3 laps, concatenate across
    % sessions and directions
    qr3laps = cell(1,3);
    detses = [2,3]; % detour sessions
    figure_letter(0.35,1.6)
    for ilap = 1:3
        for is = detses
            for idir = 1:2
                qr3laps{ilap} = cat(1,qr3laps{ilap},tsig(is,idir).lap(ilap).dataqr(:));
            end
        end
        ax = subplot(1,3,ilap);
        histogram(ax,qr3laps{ilap})
        p = signrank(qr3laps{ilap},0);
        hold(ax,'on')
        plot(ax,[0 0],ylim(ax),'-k','linewidth',1.5)
        hold(ax,'off')
        title(ax,['Lap',num2str(ilap),' P=',num2str(p,'%.3f')])
        xlabel(ax,'Quadrant Ratio')
        ylabel(ax,'Counts')
    end
    toc
end

%% code for early lap CCG correlation between theta and behavioral time scale
if runccgcorr
    % estimated run time ~4min
    tic
    load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2','PlFmesh');
    load([Datapath,'seslaps.mat'],'seslaps');
    load([Datapath,'ses.mat'],'ses');
    load([Datapath,'mazel.mat'],'mazel');
    [PCCG,tbbias] = function_DetourEarlylap_CCG_CorrXTimeScale(ses,mazel,M1,M2,PlFields,ClF,seslaps);
    
    % plot CCG bias at theta and behavioral scale in the first 3 laps
    figure_letter(0.35,1.6)
    for ilap = 1:3
        ax = subplot(1,3,ilap);
        dnow = tbbias{ilap};
        xlim(ax,[-0.15 0.15])
        ylim(ax,[-0.5 0.5])
        plot2vcorr(ax,dnow);
        title(ax,['Lap',num2str(ilap)])
        xlabel(ax,'Theta Temporal Bias')
        ylabel(ax,'Behaviroal Temporal Bias')
    end
    toc
end

%% code for population vector similarity between pre-detour mobile and detour segment 
if runpplcorr
    % estimated run time ~1min
    tic
    load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2','PlFmesh');
    load([Datapath,'seslaps.mat'],'seslaps');
    load([Datapath,'ses.mat'],'ses');
    % compute population vector cosine similarity across sessions
    pcorr = function_PplCorrXSpaceBin_PrevsDet(ses,M1,M2,PlFields,PlFmesh,ClF);
    
    % plot population vector similarity example between pre-detour and
    % detour session for one detour and driection
    % get segment length for each detour track
    rplens = Detour_GetDetourSegLen(ses,[2,4]);
    dettra = [2,4]; % detour tracks
    % define which detour track and direction to plot
    it = 4;
    idir = 2;
    detses = Det_FindDetTSes(it,ses);
    
    % plot parameters
    cmap = getNCLColmap('WhiteBlueGreenYellowRed.rgb',256);
    len = size(cmap,1);
    cmap = stripedCMap(cmap,len,round(len/5),round(len/5),0);
    trackcolors = {[0.83 0.42 0.4],[0.52 0.45 0.83],[0.92 0.70 0.20],[0.5 0.73 0.93],...
    [0.75 0.42 0.92],[0.92 0.12 0.75],[0.74 0.74 0.22],[0.22 0.74 0.74],[0.2 0.32 0.62]}; 
    colors2 = Color_FrenchDisp_modify;
    colors2{1} = colors2{3};
    colors2{2} = colors2{4};
    colors2{3} = trackcolors{1};
    colors2{4} = trackcolors{3};
    colors = trackcolors;
    
    % plot similarity, we will also add the color-coded track on x/yaxis
    figure_letter(0.35,0.9)
    ax = gca;
    axis square
    hold(ax,'on')
    plotPplVCosExample(ax,ses,rplens,it,idir,pcorr,cmap,colors,colors2)
    hold(ax,'off')
    toc
end

%% code for post-detour tuning curve predicted by pre-detour detour and others
if runplfreg
    % estimated run time ~1min
    tic
    load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2','PlFmesh');
    load([Datapath,'ses.mat'],'ses');
    [pctvar,shfpctvar] = function_PlaceMap_Mvregress(ses,M1,M2,PlFields,ClF,PlFmesh);
    colors = Color_FrenchDisp_modify;
    % plot each regressors contribution
    figure_letter(0.5)
    ax = gca;
    hold(ax,'on')
    for plotx = 1:5
        ifact = plotx;
        pctvar(ifact) = max([pctvar(ifact),1e-3]);
        % plot data
        stem(ax,plotx-0.15,pctvar(ifact),'filled','color',colors{plotx},'linewidth',1.5)
        % plot shuffle distribution
        shfd = shfpctvar(:,ifact);
        xbin = logspace(-3,2,60);
        [Nc,xbin] = histcounts(shfd,xbin);
        xbinc = (xbin(1:end-1) + xbin(2:end))/2;
        xx = berow(Nc)./max(Nc)*0.55;
        yy = berow(xbinc);
        xx = cat(2,xx,zeros(size(xx)));
        yy = cat(2,yy,fliplr(yy));
        fill(ax,xx+plotx+0.05,yy,1-0.6*(1-colors{plotx}),'LineStyle','none','FaceAlpha',0.6);
        pct = invprctile(shfd,pctvar(ifact));
        text(ax,plotx,max([0 pctvar(ifact)])+0.05,[num2str(pct,'%.3f'),'%'],'HorizontalAlignment','center')
    end
    set(ax,'yscale','log')
    xticks(ax,1:5)
    xlim(ax,[0 6])
    ylim(ax,[0.001 100])
    xticklabels(ax,{'Pre (\beta_1)','Det (\beta_2)', 'OT (\beta_3)','T1 (\beta_4)','T3 (\beta_5)'})
    xtickangle(ax,45)
    ylabel(ax,'Residual Decrease Pct')
    toc
end

%% code for rank order correlation between theta cycle and pre-detour sleep
if runrankcorr
    % estimated run time ~4min
    tic
    load([Datapath,'SleepSeq.mat'],'Seqses');
    load([Datapath,'ThetaCycOrder.mat'],'ThetaCycOrder');
    nshf = 100; % number of shuffles
    [SeqC,SeqSig] = function_DetourThetaCycle_PreSlpFrame_RankCorrPct(Seqses,ThetaCycOrder,clu2pyr,nc,nshf);
    % plot tools
    colors = Color_FrenchDisp_modify;
    cmap=getNCLColmap('MPL_rainbow.rgb',200);
    detses = [2,3]; % detour session
    
    % plot rank order correlation example 
    figure_letter(0.45,1.5)
    ax = subplot('position',[0.06 0.14 0.38 0.72]);
    dnow = SeqC{3,2,1};
    plotalpha = ones(size(dnow));
    plotalpha(isnan(dnow)) = 0;
    % if less than 5 common active cell, set that pixel transparent
    h = imagesc(ax,[1,size(dnow,2)],[1,size(dnow,1)],dnow);
    h.AlphaData = plotalpha;
    ax.YDir = 'normal';
    shading flat
    colormap(ax,cmap)
    caxis(ax,[-1 1])
    set(ax,'layer','top')
    cb = addcolorbar(ax,0.02);
    cb.Label.String = 'Rankcorr';
    xlabel(ax,'Theta cycles')
    ylabel(ax,'Pre-detour sleep frames')
    % plot significant ratio
    ax = subplot('position',[0.6 0.14 0.33 0.72]);
    hold(ax,'on')
    for js = 1:length(detses)
        is = detses(js);
        for idir = 1:2
            bar(ax,idir + (js-1)*2,SeqSig(is,idir,1)/SeqSig(is,idir,2),...
                0.25,'FaceColor',colors{1})
        end
    end
    plot(ax,xlim(ax),[0.05 0.05],'--k');
    hold(ax,'off')
    xticks(ax,[1,2,3,4])
    xticklabels(ax,{'R2D1','R2D2','R3D1','R3D2'})
    ylabel(ax,'Significant Ratio')
    toc
end


%% code for get flickering of pre-detour mobile segment during detour run
if rundetflc
    % estimated run time ~5s
    tic
    load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2','PlFmesh');
    load([Datapath,'ses.mat'],'ses');
    load([Datapath,'mazel.mat'],'mazel');
    load([Datapath,'T2T4Laps.mat'],'T2T4Laps');
    tbin = 0.04;
    overlap = 0;
    [detdecode,mscotpbdir] = function_DetRun_FlickeringofPreMobile(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,T2T4Laps,tbin,overlap);
    
    % plot rate of strong representation moments of pre-detour mobile
    % or contol segment
    % refine to detour sessions
    plotd = mscotpbdir([2,3],:,:);
    plotd = reshape(plotd,[],2); % D1 samples, D2 pre-mobile and control
    figure_letter(0.32,0.7)
    ax = gca;
    MyPairSEMBarPlot(plotd,ax,1:2,1,'sigplot',0,'markersize',8,'CapSize',14)
    xticks(ax,[1,2])
    xticklabels(ax,{'Mobile','Control'})
    ylabel(ax,'Averaged probability')
    toc
end

%% code for get flickering of detour segment during reversal run
if runrevflc
    % estimated run time ~5s
    tic
    load([Datapath,'M1M2PlFieldsPlFmeshClF.mat'],'M1','PlFields','ClF','M2','PlFmesh');
    load([Datapath,'ses.mat'],'ses');
    load([Datapath,'mazel.mat'],'mazel');
    load([Datapath,'T2T4Laps.mat'],'T2T4Laps');
    tbin = 0.04;
    overlap = 0;
    [revdecode,mdetotpbdir] = function_PostDetRun_FlickeringofDetour(ses,mazel,M1,M2,PlFields,PlFmesh,ClF,T2T4Laps,tbin,overlap);
    
    % plot rate of strong representation moments of pre-detour mobile
    % or contol segment
    % refine to post-detour sessions
    plotd = mdetotpbdir([3,4],:,:);
    plotd = reshape(plotd,[],2); % D1 samples, D2 pre-mobile and control
    figure_letter(0.32,0.7)
    ax = gca;
    MyPairSEMBarPlot(plotd,ax,1:2,1,'sigplot',0,'markersize',8,'CapSize',14)
    xticks(ax,[1,2])
    xticklabels(ax,{'Detour','Control'})
    ylabel(ax,'Averaged probability')
    toc
end
