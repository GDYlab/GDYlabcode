clear;clc;close all

%% define subpanel positions
testpos = 1;

% plot place map
[~,pos11] = tight_subplot(2,2,[0.025 0.025],[0.75 0.04],[0.065 0.775],0);
% plot tuning curves
[~,pos12] = tight_subplot(1,4,[0.0 0.04],[0.77 0.05],[0.32 0.04],0);

% theta sequence vs run vel
[~,pos21] = tight_subplot(1,2,[0.05 0.01],[0.48 0.365],[0.04 0.72],0);
[~,pos22] = tight_subplot(2,4,[0.015 0.019],[0.47 0.34],[0.355 0.335],0);
[~,pos23] = tight_subplot(1,2,[0.05 0.01],[0.48 0.365],[0.74 0.02],0);

% % dissociation between lfp and behaivor
% lfp examples
[~,pos31] = tight_subplot(1,2,[0.0 0.035],[0.22 0.62],[0.05 0.67],0);
% var vs vel
[~,pos32] = tight_subplot(1,4,[0.0 0.04],[0.22 0.62],[0.38 0.04],0);

% test position
if testpos
f = figure_letter(0.6);
testsubpos(pos11)
testsubpos(pos12)
testsubpos(pos21)
testsubpos(pos22)
testsubpos(pos23)
testsubpos(pos31)
testsubpos(pos32)
end
close all
f = figure_letter(0.6);
%% define parameters
subpanelfont = 12;
labelfont = 7;
legendfont = 6;
smallfont = 5;


global txtfonts sescolors sesstr pdfcolors cmap
txtfonts = [subpanelfont,labelfont,legendfont,smallfont];
sescolors = Color_Ketamine_v2;
pdfcolors = {[0.68 0.4 0.1],[0.1 0.4 0.68]};
sesstr = {'Run1','Run2Bi','Run2woBi','Run3'};
cmap = getNCLColmap('WhiteYellowOrangeRed.rgb',100);
%% add subpanel annotations
addsubp = 1;
if addsubp
    panellbpos = {[0 1],[0.27 1],...
        [0 0.71],[0.28 0.71],[0.685 0.71]...
        [0 0.43],[0.32 0.43]};
    panellb = {'a','b','c','d','e','f','g'};
    for isub = 1:length(panellbpos)
        posstr = panellbpos{isub};
        % add subpanel label
        annotation(f,'textbox',...
            [posstr(1) posstr(2) 0 0],...
            'String',panellb{isub},...
            'LineStyle','none',...
            'FontSize',txtfonts(1),...
            'fontweight','bold',...
            'FitBoxToText','off');
    end
end

%% define data paths
Inpath = cell(10,12); % 10 subpanels, each could have up to 6 data source
rnum = "XRat1";
setnum = 1;

[Datasets] = getCADatasets(rnum,setnum); %
Dataset = Datasets{1};

Data_Path = [getDataRoot(),'CAResults',filesep,'Preprocess',filesep,'KTM_Run',filesep];
Data_Path2 = [getDataRoot(),'CAResults',filesep,'Preprocess',filesep,'KTM_RunLFP',filesep];
Data_Path3 = [getDataRoot(),'CAResults',filesep,'Preprocess',filesep,'KTM_AwakeRest',filesep];
Data_Path4 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_Sleep',filesep];
Data_Path5 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_AwakeRest',filesep];
Data_Path6 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_Run',filesep];
Data_Path7 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_onMazeAll',filesep];

Inpath{1,1} = [Data_Path, Dataset,'Step06',filesep];% path to ses
Inpath{1,2} = [Data_Path, Dataset,'Step04',filesep];% path to M1, PlFields 2cm resolution
Inpath{1,3} = [Data_Path, Dataset,'Step03',filesep];% path to mazel
Inpath{1,4} = [Data_Path2, Dataset,'Thetapeak',filesep];% path to run LFP theta peaks
Inpath{1,5} = [Data_Path2, Dataset,'BestTT',filesep];% path to run LFP
Inpath{1,6} = [Data_Path2, Dataset,'OneTTRips',filesep];% path to run LFP ripple events
Inpath{1,7} = [Data_Path3, Dataset,'BSMUAstd2',filesep];% path to frames
Inpath{1,8} = [Data_Path6, Dataset,'RunSeq',filesep,'CosAngDiff_bylap',filesep];% path to mu and sigma of spherical cooridnations
Inpath{1,9} = [Data_Path5, Dataset,'BayeDcode',filesep,'All442',filesep,'FrameSes_Alllambda','Exp_Var',filesep,'BS_MUAstd2',filesep];% path alllambda
Inpath{1,10} = [Data_Path3, Dataset,'PyrISI_runvsrest_byDir_v2_KetImpactRun2',filesep,'Twin300Step120',filesep,'runvel10restvel2duration0.5',filesep];% path to rest/run mua
Inpath{1,11} = [getDataRoot(),'CAResults/Sets_Ave/KTM_AwakeRest/KTMUActNew/PyrISI_runvsrest_byDir_v2_KetImpactRun2/Twin300Step120/runvel10restvel2duration0.5/'];% path to rest/run ISI kldiv

Inpath{2,1} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/KTMUActNew/Run/LFPvsRunVel_v2_KetImpactRun2/Twin300Step120/']; % run vel vs lfp var
Inpath{2,2} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/RemapDrift/RunSeq/MazeSeq_EachSes_KetImpactRun2/Twin300Step120Vel10/'];% maze sequence of juvenile

Inpath{3,1} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/SpatialTuning/Pyr_SigSpatialInfo_Track_Abs_KetImpactRun2/Twin300Step120/']; % tuning feautre 
Inpath{3,2} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/RunDecoding/BayeDcode/BayesianRun_InSesErrvsShuffle_KetImpactRun2/Twin300Step120/']; % decoding results
Inpath{3,3} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/FiringRate/CellRate_EntireSes_KetImpactRun2/Twin300Step120/']; % rate in sessions

% Inpath{3,1} = [Data_Path7, Dataset, 'Decoding',filesep,'MazeAllLaps_TAlignLposC_Lfit_2cm_byVel_KetImpactRun2',...
%     filesep,'Twin300Step120',filesep,'WinLen',num2str(winlen),'Vel',num2str(velbin),filesep];
Inpath{4,1} = [getDataRoot(),'CAResults/Sets_Ave/KTM_onMazeAll/KTMUActNew/Decoding/MazeAllLaps_TAlignLposC_Lfit_2cm_byVel_KetImpactRun2/Twin300Step120/WinLen0.04Vel5   15   30  120/'];
Inpath{4,2} = [getDataRoot(),'CAResults/Sets_Ave/KTM_onMazeAll/KTMUActNew/Decoding/MazeAllLaps_TAlignLposC_Lfit_2cm_byVel_KetImpactRun2/Twin300Step120/WinLen0.04Vel5   15   30  120/'];
Inpath{4,3} = [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/SpatialTuning/Pyr_SigSpatialInfo_Maze_AbsbyVRange_KetImpactRun2/Twin300Step120/5   15   30  120/'];

Outpath =   '/media/yuchen/data14/CAResults/Figures/Ketamine/';
figurename = 'DissociationBehHPC_v9';
mkdir(Outpath)
velbin = [5,15,30,120];
%% make subplots
% plot placemap sort based on each session
PlotXsesPlaceMap(pos11,Inpath{2,2});

% plot spatial info by vel
PlotInfoByVel(pos21,Inpath{4,3})

% plot theta sequence
PlotThetaSeqSlope(pos23,Inpath{4,2})

PlotThetaSeqExample(pos22,Inpath{4,1},velbin)

% plot spatial tuning across sessions
PlotXsesSpatialTuning(pos12,Inpath{3,1},Inpath{3,2},Inpath{3,3});



% plot run1 run2 lfp examples
PlotRunLFPvsVel(pos31,Inpath{1,3},Inpath{1,5})

% plot lfp var vs run vel
PlotRunLFPVarvsVel(pos32,Inpath{2,1})

% % plot run/rest raster 
% PlotRestRunRaster(pos31,Inpath{1,1},Inpath{1,2},Inpath{1,3})


%% save the result
File_Path = strcat(Outpath,figurename,'.fig');
saveas(gcf, File_Path);
File_Path = strcat(Outpath,figurename,'.png');
print(gcf,File_Path,'-dpng','-r300'); 
% export_fig([Outpath,figurename,'.pdf'])
print('-painters',gcf,'-dpdf',[Outpath,figurename,'_AI.pdf']) % this is for illustrator
close all
%% subpanel plot function
function  PlotInfoByVel(pos,InPath)
global txtfonts sescolors sesstr pdfcolors
load([InPath,'ctuningbySesVbin.mat'],'cstab','cinfo','cpk','flen')
nsess = size(cstab,2);
nvbin = size(cstab,3);
vbin = [5,15,30,120];
pkcri = 1;

infoxvel = cell(nsess,nvbin);
for is = 1:nsess
    for ivbin = 1:nvbin
        pkidx = cpk(:,is,ivbin) >= pkcri;
        infoxvel{is,ivbin} = cinfo(pkidx,is,ivbin);
    end
end

% plot spatial info vs ses vel;
ax = subplot('position',pos{1});
hold(ax,'on')
set(ax,'FontSize',txtfonts(4))
% plot slop vs ses vel;
hold(ax,'on')
vbinmean = GetCenterfromEdge(vbin);
% vbinmean = [10,20,40,80];
% vbinmean = log(vbinmean);

vxlabel = [20 40 60];
% vxlabel = log(vxlabel);
% vbinstr = cell(size(vxlabel));
% for iv = 1:length(vxlabel)
%     vbinstr{iv} = num2str(exp(vxlabel(iv)));
% end

lgd = cell(1,nsess);
alltable = []; 
for is = 1:nsess
   sidx = KetSesMap(is);
   hd(is) = plot(ax,nan,nan,'s','markerfacecolor',sescolors{is},'markeredgecolor',sescolors{is},'markersize',3.5); 
   dxx = [];
   dyy = [];
   for ivbin = 1:nvbin
       dtmp = infoxvel{sidx,ivbin};
       velbin = vbinmean(ivbin).*ones(size(dtmp));
       
       dxx = cat(1,dxx,velbin);
       dyy = cat(1,dyy,dtmp);
       sesnum = is.*ones(size(dtmp));
       
       alltable = cat(1,alltable,[dtmp, sesnum, velbin]);
   end

   h = PlotCorrFitCIShade(ax, [dxx,dyy], sescolors{is}, 0.2,'doScatter',false,'doText',false);
   [R,p] = nancorr(dxx,dyy);
   lgd{is} = [sesstr{is},'; R=',num2str(R,'%.2f'),'; P=',num2str(p,'%.2f')];
end

[lgd,icon]=legend(hd,lgd,'fontsize',txtfonts(4)-0.5,'location','northwest');
icon(1).Position(1) = icon(1).Position(1) - 0.14;
icon(2).Position(1) = icon(2).Position(1) - 0.14;
icon(3).Position(1) = icon(3).Position(1) - 0.14;
icon(4).Position(1) = icon(4).Position(1) - 0.14;

legend boxoff
set(lgd,'Position',[0.025 0.605 0.155 0.06],'FontSize',4.5);
xlabel(ax,'Velocity (cm/s)')
ylabel(ax,'Spatial Info (bits/spk)')
xticks(ax,vxlabel)
% xticklabels(ax,vbinstr)
ylim(ax,[0.5 2.5])
yticks(ax,0.5:0.5:2.5)
xlim(ax,[10,70])
ax.TickLength = [0.035, 0.035];

% tt = 'Theta sequence slope';
% title(tt)

[p,tbls,stats] = anovan(alltable(:,1),{alltable(:,2),alltable(:,3)},'continuous',[2],'model','interaction',...
    'varnames',{'session','vel'},'Display','off');
[results,c,m,gnames] = multcompare(stats,'dimension',[1],'Display','off');
%
% use anova to test the effect of ses diff and no, det, crossdet
% pos{2}(4) = pos{2}(4) - 0.02;
ax = subplot('position',pos{2});
hold(ax,'on')
MySEMPlot_v2(c,ax,sescolors,1:nsess,0.25,'linewidth',1.4,'markersize',0)
ylim(ax,[0.5 2.5])
myoverlapsigstar_v3(ax,{[1,2],[2,3],[2,4]},[results(1,6) results(4,6) results(5,6)],'starty',1.8)
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
xtickangle(ax,20)
title(ax,['ANOVAN Multicompare'],'FontWeight','normal')
% ylabel(ax,'Spatial Info')
xlabel(ax,'Sessions')
xlim(ax,[0.3 nsess+0.7])
yticks(ax,0.5:0.5:2.5)
yticklabels(ax,[])
set(ax,'FontSize',txtfonts(4))
ax.TickLength = [0.035, 0.035];

end


function PlotXsesPlaceMap(pos,Inpath)
global txtfonts sescolors sesstr cmap
% plot juvenile map
pkcri = 1;
infopctcri = 0;
parastr = ['pkcri',num2str(pkcri),'infopct',num2str(infopctcri)];
idir = 1;
load([Inpath,'dirmappk',parastr,'.mat'],'dirmappk')

nsess = size(dirmappk,2);

for is = 1:nsess

    sidx = KetSesMap(is);
    % based on cell seq in each session
    pkloc = dirmappk{idir,sidx,2};
    % sort cells based on pk loc
    [~,sid] = sort(pkloc,'ascend');

    ax = subplot('position',pos{is});
    hold(ax,'on')
    
    zz = dirmappk{idir,sidx,1};
    if is == 3
        imagesc(ax,[0 1],[1, 60],zz(sid,:))
    else
        imagesc(ax,[0 1],[1, length(pkloc)],zz(sid,:))
    end
    
    set(ax,'Ydir','normal')
    xticks(ax,[0 1])
    %         caxis(ax,[-1 3])
    colormap(ax,cmap)
    caxis(ax,[0 1])
    title(ax,[sesstr{is}],'FontWeight','bold','FontSize',txtfonts(3)+0.5)
    axis tight
    
    yticks(ax,[20 40 60 80])
    if is == 1
        ylb = ylabel(ax,['Cell ID from each session']);
        ylb.Position(2) = ylb.Position(2) -60;
%     else
%         yticklabels(ax,[])
    end
    if is == 3
       xlb = xlabel(ax,'Normalized position');
       xlb.Position(1) = xlb.Position(1) + 0.7; 
       xlb.Position(2) = xlb.Position(2) + 6; 
       ylim(ax,[1 60])
    end
    if is <= 2
       xticks(ax,[]) 
    end
   
    
    if is == nsess
        cb = addcolorbar;
        cb.Label.String = 'Normalized rate';
        cb.Label.Position(1) = cb.Label.Position(1) - 2.5;
        cb.Ticks = [0,1];
    end
    set(ax,'fontsize',txtfonts(4))
    % we will save placemap and their pk location
end

end

function PlotXsesSpatialTuning(pos,Inpath,Inpath2,Inpath3)
global txtfonts sescolors sesstr
% plot juvenile map

load([Inpath,'ratd.mat'],'ratd')
load([Inpath3,'rates.mat'],'rates')
ratesnew = permute(rates, [1, 2, 3]);
ictype = 1;
seslb = sesstr;
mstrs = {'Log10 peak rate (Hz)','Spatial info (bits/spk)','Log10 field length (cm)'};
nsess = 4;
pkcri = 1;
poolm = cell(nsess,4); % D1 sessions, D2 measures
allylim = {[-2 3],[0 4.5],[0.5 3]};
allyticks = {[-2:2],[0:3],[1:0.5:2]};

for is = 1:nsess
    sidx = KetSesMap(is);
    for id = 1:length(ratd)
        cinfo = ratd(id).cinfo(:,sidx,:,:);
        cpk = ratd(id).cpk(:,sidx,:,:);
        valididx = cpk(:) >= log10(pkcri);
        flen = log10(ratd(id).flen{sidx});
        
        poolm{is,1} = cat(1,poolm{is,1},cpk(:));
        poolm{is,2} = cat(1,poolm{is,2},cinfo(valididx));
        poolm{is,3} = cat(1,poolm{is,3},flen(:));
    end
end

for im = 1:3
    ax = subplot('position',pos{im});
    hold(ax,'on')
    MyBoxPlot_v4(poolm(:,im),ax,'sigplot',0,'color',sescolors,'wislinewidth',1,'meanlinewidth',1.5)
    xticks(ax,1:nsess)
    xticklabels(ax,sesstr)    
    xtickangle(ax,35)
    ylabel(ax,mstrs{im})
    yticks(ax,allyticks{im})
    xlim(ax,[0.2 nsess+0.8])
    [groups,pvals] = AllGroupRankSum(poolm(:,im));
    pidx = [1,4,5];
    ylim(ax,allylim{im})
%     myoverlapsigstar(ax,groups,pvals)
    ystart = allylim{im}(2) * 0.8 + allylim{im}(1) * 0.2;
    yrange = diff(allylim{im}) * 0.03;
    myoverlapsigstar_v3(ax,groups(pidx),pvals(pidx),'starty',ystart,'ygap',yrange,'starsize',6)
    set(ax,'fontsize',txtfonts(3)-0.5)
    set(ax, 'YGrid', 'on', 'XGrid', 'off')
    ax.TickLength = [0.03 0.03];
%     xlb = xlabel(ax,'Sessions');
end

%% plot decoding error

load([Inpath2,'decodeerror'],'shfdecer','decer')
nsess = 4;
alld = cell(2*nsess,1);
plotx = nan(2*nsess,1);
pltc = cell(2*nsess,1);
gs1 = cell(1,nsess);
ps1 = nan(1,nsess);

ax = subplot('position',pos{4});
hold(ax,'on')
for is = 1:nsess
   sidx = KetSesMap(is);
   alld{2*is-1} = log10(decer{sidx});
   alld{2*is} = log10(shfdecer{sidx});
   plotx(2*is-1) = is-0.15;
   plotx(2*is) = is+0.15;
   pltc{2*is-1} = sescolors{is};
   pltc{2*is} = [0.6 0.6 0.6];
   gs1{is} = [is-0.15,is+0.15];
   ps1(is) = ranksum(decer{sidx},shfdecer{sidx});
   myoverlapsigstar_v3(ax,gs1(is),ps1(is),'starty',2.25,'ygap',0.13,'starsize',6)
end
gs2{1} = [1,2] - 0.12;
gs2{2} = [2,3] - 0.12;
gs2{3} = [2,4] - 0.12;
ps2(1) = ranksum(decer{1},decer{2});
ps2(2) = ranksum(decer{2},decer{3});
ps2(3) = ranksum(decer{2},decer{4});

MyBoxPlot_v4(alld,ax,'volinpdf',0,'plotx',plotx,'colors',pltc,'boxlen',0.1,'sigplot',0)
ylim(ax,[-1.5 4])

myoverlapsigstar_v3(ax,gs2,ps2,'starty',2.75,'ygap',0.13,'starsize',6)
% myoverlapsigstar_v3(ax,gs,ps)
% hd(1) = plot(ax,nan,nan,'-','linewidth',5,'color',colors{6});
% hd(2) = plot(ax,nan,nan,'-','linewidth',5,'color',[0.6 0.6 0.6]);
% legend(hd,{'Data','Shuffle'})
% legend boxoff

yticks(ax,[-1:2])
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
ylb = ylabel(ax,'Log10 decoding error (cm)');
% ylb.Position(1) = ylb.Position(1) + 0.2;
% ylb.Position(2) = ylb.Position(2) - 3.2;
set(ax,'fontsize',txtfonts(3)-0.5)
set(ax, 'YGrid', 'on', 'XGrid', 'off')
xtickangle(ax,35)
% xlb = xlabel(ax,'Sessions');
% xlb.Position(1) = xlb.Position(1) - 2.65;

plot(ax,3.15,3.75,'s','markerfacecolor',sescolors{1},'markeredgecolor',sescolors{1},'markersize',3)
plot(ax,3.3,3.75,'s','markerfacecolor',sescolors{2},'markeredgecolor',sescolors{2},'markersize',3)
plot(ax,3.45,3.75,'s','markerfacecolor',sescolors{3},'markeredgecolor',sescolors{3},'markersize',3)
plot(ax,3.6,3.75,'s','markerfacecolor',sescolors{4},'markeredgecolor',sescolors{4},'markersize',3)
text(ax,3.75,3.75,'Data','HorizontalAlignment','left','FontSize',txtfonts(3)-1)

plot(ax,3.375,3.45,'s','markerfacecolor',[0.6 0.6 0.6],'markeredgecolor',[0.6 0.6 0.6],'markersize',3)
text(ax,3.7,3.4,'Shuffle','HorizontalAlignment','left','FontSize',txtfonts(3)-1)

end

function PlotThetaSeqExample(pos,Inpath,velbin)
global txtfonts sescolors sesstr pdfcolors cmap
vbin = velbin;
% cmap = getNCLColmap('MPL_PuRd.rgb',100); 
% cmap = getNCLColmap('MPL_Reds.rgb',100); 

load([Inpath,'thetapbDS.mat'],'thetapbDS','thetalfpDS','meshsize','actextend','cttmesh')

lfitrange = 0.05; % we will compute quadrant ratio +-50ms

pbcutoff = 0.04;
if meshsize >=5
    pbbinextend = 7; % we will extend +-5 position bins from center of position bin
else
    pbbinextend = 20; % we will extend +-20 position bins from center of position bin
end

ntrack = size(thetapbDS,3);
vbin2plot = [2,3];
nvbin = size(thetapbDS,4);
nsess = size(thetapbDS,1);
for is = 1:nsess
    sid = KetSesMap(is);
    for ivbinidx = 1:2
        ivbin = vbin2plot(ivbinidx);
        pdfnow = nan(41,81,5000);
        count = 1;
        for idir = 1:2
            for it = 1:ntrack
                dtemp = thetapbDS{sid,idir,it,ivbin};
                
                if length(dtemp) <= 2
                    continue
                end
                
                for iseg = 1:length(dtemp)
                    pdfsample = dtemp{iseg};
                    if idir == 1
                        pdfsample = flipud(pdfsample);
                    end
                    pdfnow(:,:,count) = pdfsample;
                    count = count + 1;
                end
            end
        end
        meanpb = nanmean(pdfnow,3);
        
        % get qr test p value
        lefttidx = cttmesh >= -lfitrange & cttmesh <0;
        rightidx = cttmesh <= lfitrange & cttmesh >= 0;
        halfposb = round(size(pdfnow,1)/2);
        v1 = pdfnow(1:halfposb,lefttidx,:);
        v2 = pdfnow(1:halfposb,rightidx,:);
        v3 = pdfnow(halfposb+1:end,rightidx,:);
        v4 = pdfnow(halfposb+1:end,lefttidx,:);
        v13 = cat(1,v1(:),v3(:));
        v24 = cat(1,v2(:),v4(:));
        if idir == 1
            [~,ptmp] = ttest2(v13,v24,'tail','left');
        else
            [~,ptmp] = ttest2(v13,v24,'tail','right');
        end
        
        tmppb = meanpb;
        % we need to set center position bins as 0, other wise the
        % linefit will just detect the center bins
        timeidx = idxinrange_v2(cttmesh,[-lfitrange lfitrange]);
        
        midremovebin = 0; % we will remove in total 2*midremovebin+1 bins in the middle
        tmppb(pbbinextend-midremovebin+1:pbbinextend+midremovebin+1,:) = [];
        % get max jump
        jumpmesh = linspace(-meshsize*pbbinextend,meshsize*pbbinextend,2*pbbinextend-2*midremovebin);
        mjump = Decode_Maxjump(tmppb,jumpmesh);
        
        tmppb(tmppb>=pbcutoff) = pbcutoff;
        
        tmppb_rd = tmppb(:,timeidx);
        tmppb_rd = tmppb_rd./repmat(sum(tmppb_rd,1),size(tmppb_rd,1),1);
        
%         [slope,interc,lfs0,~]=est_line_detect(1:size(tmppb_rd,2),1:size(tmppb_rd,1),tmppb_rd,'kernelwidth',[10 0]);
%         % slope is posbins/timebins of best linefit
%         % to get sptial covergae we should do slope*number of time bins* spatial bin
%         % length
%         spacover = slope*(size(tmppb_rd,2)-1)*meshsize;
%         % also get the lpos of linefit at start and end of frame
%         lfsta = (interc-1)*meshsize - meshsize*(pbbinextend-midremovebin) + slope*1*meshsize;
%         lfend = (interc-1)*meshsize - meshsize*(pbbinextend-midremovebin) + slope*size(tmppb_rd,2)*meshsize;
        
        ax = subplot('position',pos{is + (ivbinidx-1)*nsess});
        hold(ax,'on')
        if idir == 2
            pmoveleft(0.01,ax)
        end
        if is == 1 
            ylabel(ax,{['Vel',num2str(vbin(ivbin)),'-',num2str(vbin(ivbin+1))],['position (cm)']})
        else
            yticklabels(ax,[])
        end
        
        caxis(ax,[0 pbcutoff])
        
        if ivbinidx == 1
            title(ax,[sesstr{is}])
        end
        
        if ivbinidx < 2
            xticklabels(ax,[])
        else
            xticks(ax,[-0.145,0,0.145])
            xticklabels(ax,{'-0.16','0','0.16'})
        end
        imagesc(ax,[-actextend actextend],[-meshsize*(pbbinextend-midremovebin) meshsize*(pbbinextend-midremovebin)],tmppb)
        colormap(ax,cmap);
        xlim(ax,[-actextend actextend])
        %             xticks(ax,[-actextend 0 actextend])
        ylim(ax,[-40 40])
        
        % plot the range where we compute linefit
        plot(ax,[-lfitrange -lfitrange],[-meshsize*(pbbinextend-midremovebin) meshsize*(pbbinextend-midremovebin)],'--k','linewidth',1)
        plot(ax,[lfitrange lfitrange],[-meshsize*(pbbinextend-midremovebin) meshsize*(pbbinextend-midremovebin)],'--k','linewidth',1)
        % plot the linefit
        %             if slopeidx
        %plot(ax,[-lfitrange lfitrange],[lfsta lfend],'-','linewidth',1,'color',[0.23 0.36 0.85])
        %             end
        
        %         ylim(ax,[-meshsize*pbbinextend meshsize*pbbinextend])
        if ivbinidx == 2  && is == 2
            xlb = xlabel(ax,'Time (s)');
            xlb.Position(1) = xlb.Position(1) + 0.2;
            xlb.Position(2) = xlb.Position(2) + 10;
        end
        
        set(ax,'FontSize',txtfonts(4))                
        if  is == nsess && ivbinidx == 2
            cb = addcolorbar;
            cb.Label.String = 'Probability';
            cb.Ticks = [0 pbcutoff];
            cb.Label.Position(1) = cb.Label.Position(1) - 1.5;
        end
    end
end
end

function  PlotThetaSeqSlope(pos,InPath)
global txtfonts sescolors sesstr pdfcolors
load([InPath,'tslope.mat'],'tslope','meanvel')
load([InPath,'thetapbDS.mat'],'thetapbDS')
vbin = [5,15,30,120];

% plot slop vs ses vel;
ax = subplot('position',pos{1});
ax.TickLength = [0.035, 0.035];
hold(ax,'on')
set(ax,'FontSize',txtfonts(4))
vbinmean = GetCenterfromEdge(vbin);
% vbinmean = [10,20,40,80];
vbinmean = log(vbinmean);
vbinstr = cell(size(vbinmean));
nrat = size(tslope,3);
nsess = size(tslope,1);
nvbin = length(vbin)-1;

sesnum = becolumn(1:nsess);
vbinnum = berow(1:nvbin);

sesnum = repmat(sesnum,[1,nvbin,nrat]);
vbinnum = repmat(vbinnum,[nsess,1,nrat]);
vbinnum = meanvel;
for iv = 1:length(vbinmean)
    vbinstr{iv} = num2str(exp(vbinmean(iv)));
end


lgd = cell(1,nsess);
for is = 1:nsess
   sidx = KetSesMap(is);
   hd(is) = plot(ax,nan,nan,'s','markerfacecolor',sescolors{is},'markeredgecolor',sescolors{is},'markersize',3.5); 
   dtmp = squeeze(tslope(sidx,:,:));
   velbin = squeeze(vbinnum(sidx,:,:));
   dxx = dtmp(:);
   dyy = velbin(:);
   h = PlotCorrFitCIShade(ax, [dyy,dxx], sescolors{is}, 0.2,'doScatter',false,'doText',false);
   [R,p] = nancorr(dxx,dyy);
   lgd{is} = [sesstr{is},'; R=',num2str(R,'%.2f'),'; P=',num2str(p,'%.2f')];
end

[lgd,icon]=legend(hd,lgd,'fontsize',txtfonts(4)-0.5,'location','northwest');
icon(1).Position(1) = icon(1).Position(1) - 0.14;
icon(2).Position(1) = icon(2).Position(1) - 0.14;
icon(3).Position(1) = icon(3).Position(1) - 0.14;
icon(4).Position(1) = icon(4).Position(1) - 0.14;

legend boxoff
set(lgd,'Position',[0.723 0.605 0.155 0.06],'FontSize',4.5);

% tt = 'Theta sequence slope';
% title(tt)

xlabel(ax,'Velocity (cm/s)')
ylabel(ax,'Theta sequence slope')
ylim(ax,[0 0.4])
xlim(ax,[10 70])
xticks(ax,[20,40,60])
yticks(ax,[0:0.1:0.5])
% xlim(ax,[log(5),log(150)])
% xticks(ax,vbinmean)
% xticklabels(ax,vbinstr)
set(ax,'FontSize',txtfonts(4))

% build anova to check the contribution of session and vel

alltable = [];
for iel = 1:numel(tslope)
    if ~isnan(tslope(iel))
        alltable = cat(1,alltable,[tslope(iel), KetSesMap(sesnum(iel)), vbinnum(iel)]);
    end
end

[p,tbls,stats] = anovan(alltable(:,1),{alltable(:,2),alltable(:,3)},'continuous',[2],'model','interaction',...
    'varnames',{'session','vel'},'Display','off');
[results,c,m,gnames] = multcompare(stats,'dimension',[1],'Display','off');
%
% use anova to test the effect of ses diff and no, det, crossdet
% pos{2}(4) = pos{2}(4) - 0.02;
ax = subplot('position',pos{2});
hold(ax,'on')
MySEMPlot_v2(c,ax,sescolors,1:nsess,0.25,'linewidth',1.4,'markersize',0)
ylim(ax,[0 0.4])
myoverlapsigstar_v3(ax,{[1,2],[2,3],[2,4]},[results(1,6) results(4,6) results(5,6)],'starty',0.25)
ylim(ax,[0 0.3])
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
xtickangle(ax,20)
title(ax,['ANOVAN Multicompare'],'FontWeight','normal')
% ylabel(ax,'Theta sequence slope')
xlabel(ax,'Sessions')
xlim(ax,[0.3 nsess+0.7])
yticks(ax,[0:0.1:0.5])
yticklabels(ax,[])
ylim(ax,[0 0.4])
ax.TickLength = [0.035, 0.035];
set(ax,'FontSize',txtfonts(4))
end

function PlotRunLFPvsVel(pos,Inpath1,Inpath2)
global txtfonts sescolors sesstr pdfcolors
load ([Inpath1,'Maze.mat'],'mazel','times')
load ([Inpath2,'Ptet.mat'],'Ptet')


ses2run = [1,2];
trange = {[9411,9483.42],[18923,18514]};
tspan = 0.5;
for is = ses2run
    ax = subplot('position',pos{is});
    hold(ax,'on')
    smazel = mazel{is};
    mazelt = smazel(:,6);
    mazelvel = abs(smazel(:,7));
    lfptnow = Ptet(is).t;
    lfpnow = Ptet(is).LFP;
    
    lgdstr = {'lowVel','highVel'};
    for iv = 1:2
        trangenow = [trange{is}(iv),trange{is}(iv)+tspan];
        lfpidx = idxinrange(lfptnow,trangenow);
        lfpplot = lfpnow(lfpidx);
        lfpplot = lfpplot - nanmean(lfpplot);
        lfptplot = lfptnow(lfpidx);
        lfptplot = lfptplot - min(lfptplot);
        plot(ax,lfptplot,lfpplot,'-','color',pdfcolors{iv})
        
        mtidx = idxinrange(mazelt,trangenow);
        meanv = nanmean(mazelvel(mtidx));
        lfpvar = var(lfpplot);
%         lgdstr{iv} = ['meanVel',num2str(meanv,'%.1f'),'; LFPVar',num2str(lfpvar,'%.2f')];
    end
    xlabel(ax,'Time (s)')
    if is == 1
        ylabel(ax,'LFP (mV)')
    else
        yticklabels(ax,[])
        pmoveleft(0.015,ax)
    end
    xlim(ax,[0 0.5])
    ylim(ax,[-5e-4 7e-4])
    [lgd,icon]=legend(lgdstr,'fontsize',txtfonts(4),'Location','northeast');
    lgd.Position(2) = lgd.Position(2) + 0.0205;
    icon(3).XData(1) = icon(3).XData(1)+0.42;
    icon(5).XData(1) = icon(5).XData(1)+0.42;
    legend boxoff
    title(ax,sesstr{is})
    set(ax,'fontsize',txtfonts(4))
end



end

function PlotRunLFPVarvsVel(pos,Inpath1)
global txtfonts sescolors sesstr pdfcolors
load ([Inpath1,'lfpseg.mat'],'lfpseg')
nsess = length(lfpseg);
vcri = 10;
for is = 1:nsess
    sid = KetSesMap(is);
    dnow = lfpseg{sid};
    
    velnow = dnow(:,7);
    valididx = velnow >= vcri;
    velnow = velnow(valididx);
    im = 1;
    ax = subplot('position',pos{is});
    xnow = dnow(valididx,im);
    plot2vcorr(ax,cat(2,log(velnow),log(xnow)),'markersize',0.5,'marker','.','rpval',0,'color',sescolors{is})
    xx = becolumn(log(velnow));
    yy = becolumn(log(xnow));
    goodidx = ~isnan(xx) & ~isnan(yy);
    [r,p] = corr(xx(goodidx),yy(goodidx));
    xlabel(ax,'LogVel')
    if is == 1
    ylabel(ax,'LogVar')
    end
    title(ax,[' ',sesstr{is},' R=',num2str(r,'%.2f'),' P=',num2str(p,'%.3f')],'HorizontalAlignment','center')
    set(ax,'fontsize',txtfonts(4))
end


end


