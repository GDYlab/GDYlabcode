clear;clc;close all

dbstop if error
%% define subpanel positions
testpos = 1;

% % dissociation between run and waking rest representation
% replay examples
[~,pos11] = tight_subplot(1,3,[0.0 0.018],[0.84 0.05],[0.05 0.07],0);
[~,pos12] = tight_subplot(1,3,[0.0 0.018],[0.52 0.37],[0.05 0.07],0);

[~,pos21] = tight_subplot(2,1,[0.015 0.065],[0.05 0.70],[0.055 0.81],0);
[~,pos22] = tight_subplot(2,2,[0.025 0.01],[0.05 0.70],[0.26 0.57],0);
[~,pos23] = tight_subplot(1,2,[0.04 0.065],[0.065 0.71],[0.56 0.05],0);
% test position
if testpos
f = figure_letter(0.55);
testsubpos(pos11)
testsubpos(pos12)
testsubpos(pos21)
testsubpos(pos22)
testsubpos(pos23)
end
close all
f = figure_letter(0.55);
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
    panellbpos = {[0 1],[0 0.68],[0 0.34],[0.20 0.34],[0.50 0.34],[0.72 0.34]};
    panellb = {'a','b','c','d','e','f'};
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

rnum = "URat1";
setnum = 1;

[Datasets2] = getCADatasets(rnum,setnum); %
Dataset2 = Datasets2{1};

Data_Path = [getDataRoot(),'CAResults',filesep,'Preprocess',filesep,'KTM_Run',filesep];
Data_Path2 = [getDataRoot(),'CAResults',filesep,'Preprocess',filesep,'KTM_RunLFP',filesep];
Data_Path3 = [getDataRoot(),'CAResults',filesep,'Preprocess',filesep,'KTM_AwakeRest',filesep];
Data_Path4 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_Sleep',filesep];
Data_Path5 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_AwakeRest',filesep];
Data_Path6 = [getDataRoot(),'CAResults',filesep,'Solver',filesep,'KTM_Run',filesep];

Inpath{1,1} = [Data_Path, Dataset,'Step06',filesep];% path to ses
Inpath{1,2} = [Data_Path, Dataset,'Step04',filesep];% path to M1, PlFields 2cm resolution
Inpath{1,3} = [Data_Path, Dataset,'Step03',filesep];% path to mazel
Inpath{1,4} = [Data_Path2, Dataset,'Thetapeak',filesep];% path to run LFP theta peaks
Inpath{1,5} = [Data_Path2, Dataset,'BestTT',filesep];% path to run LFP
Inpath{1,6} = [Data_Path2, Dataset,'OneTTRips',filesep];% path to run LFP ripple events
Inpath{1,7} = [Data_Path3, Dataset,'BSMUAstd2',filesep];% path to frames
Inpath{1,8} = [Data_Path6, Dataset,'RunSeq',filesep,'CosAngDiff_bylap',filesep];% path to mu and sigma of spherical cooridnations
Inpath{1,9} = [Data_Path5, Dataset,'BayeDcode',filesep,'All442',filesep,'FrameSes_Alllambda','Exp_Var',filesep,'BS_MUAstd2',filesep];% path alllambda
Inpath{1,10} = [Data_Path3, Dataset,'MUARestRun_KetImpactRun2',filesep,'Twin300Step120',filesep,'runvel10restvel2duration2',filesep];% path to rest/run mua
Inpath{1,11} = [getDataRoot(),'CAResults/Sets_Ave/KTM_AwakeRest/KTMUActNew/MUARestRun_KetImpactRun2/Twin300Step120/runvel10restvel2duration2/'];% path to rest/run mua kldiv


Inpath{2,1} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/KTMUActNew/Run/LFPvsRunVel_v2_KetImpactRun2/Twin300Step120/']; % run vel vs lfp var

twin = 300; %  window length
tstep = 120; %  window step
vel = 10;
runorsplit = 1;
[outpstr,pm,pplf,pses,plfp,prip,ptpk,ptph,parf,ppispk] = KTMRun_SesorR2Split_Path(Dataset,runorsplit,twin,tstep,vel);
Inpath{3,1} = pses;% path to ses
Inpath{3,2} = pplf;% path to M1, PlFields 2cm resolution
Inpath{3,3} = pm;% path to mazel
Inpath{3,4} = ptpk;% path to run LFP theta peaks
Inpath{3,5} = plfp;% path to run LFP
Inpath{3,6} = prip;% path to run LFP ripple events
Inpath{3,7} = parf;% path to frames
Inpath{3,8} = [Data_Path6, Dataset,'RunSeq',filesep,'CosAngDiff_bylap',outpstr,filesep];% path to mu and sigma of spherical cooridnations
Inpath{3,9} = [Data_Path5, Dataset,'BayeDcode',filesep,'All442',filesep,'FrameSes_Alllambda',outpstr,'Exp_Var',filesep,'BS_MUAstd2',filesep];% path alllambda
AllInpath31 = {Inpath{3,1},Inpath{3,2},Inpath{3,3},Inpath{3,4},Inpath{3,5},Inpath{3,6},Inpath{3,7},Inpath{3,8},Inpath{3,9}};

[outpstr,pm,pplf,pses,plfp,prip,ptpk,ptph,parf,ppispk] = KTMRun_SesorR2Split_Path(Dataset2,runorsplit,twin,tstep,vel);
Inpath{3,1} = pses;% path to ses
Inpath{3,2} = pplf;% path to M1, PlFields 2cm resolution
Inpath{3,3} = pm;% path to mazel
Inpath{3,4} = ptpk;% path to run LFP theta peaks
Inpath{3,5} = plfp;% path to run LFP
Inpath{3,6} = prip;% path to run LFP ripple events
Inpath{3,7} = parf;% path to frames
Inpath{3,8} = [Data_Path6, Dataset2,'RunSeq',filesep,'CosAngDiff_bylap',outpstr,filesep];% path to mu and sigma of spherical cooridnations
Inpath{3,9} = [Data_Path5, Dataset2,'BayeDcode',filesep,'All442',filesep,'FrameSes_Alllambda',outpstr,'Exp_Var',filesep,'BS_MUAstd2',filesep];% path alllambda
AllInpath32 = {Inpath{3,1},Inpath{3,2},Inpath{3,3},Inpath{3,4},Inpath{3,5},Inpath{3,6},Inpath{3,7},Inpath{3,8},Inpath{3,9}};

Inpath{4,1} = [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/KTMUActNew/Run/RipFeatureXsession_KetImpactRun2/Twin300Step120/'];
Inpath{4,2} = [getDataRoot(),'CAResults/Sets_Ave/KTM_AwakeRest/KTMUActNew/FrameStatsXSes/BSMUAstd2Vcri5_KetImpactRun2/Twin300Step120/'];
Inpath{4,3} =  '/media/yuchen/data14/CAResults/Sets_Ave/KTM_AwakeRest/KTMUActNew/BayeDcode/MazeRep/NowSesMaze_CosAngDiff_KetImpactRun2/Twin300Step120/Exp_Var/BS_MUAstd2/'; % dcd pb
Inpath{4,4} =  '/media/yuchen/data14/CAResults/Sets_Ave/KTM_AwakeRest/KTMUActNew/BayeDcode/MazeRep/MazeRep_SigRatio_KetImpactRun2/Twin300Step120/BS_MUAstd2/'; % dcd sig ratio
Inpath{4,5} = [getDataRoot(),'CAResults/Sets_Ave/KTM_AwakeRest/KTMUActNew/BayeDcode/MazeRep/MazeRep_SigRatio_AbsM_KetImpactRun2/Twin300Step120/BS_MUAstd2/'];

Outpath =   '/media/yuchen/data14/CAResults/Figures/Ketamine/';
figurename = 'WakingRestDissociation_v6';
mkdir(Outpath)
%% make subplots
% plot decoding probability
PlotDcdProb(pos23{2},Inpath{4,3})

% plot decoding sig ratio
PlotDcdSigRatio(pos23{1},Inpath{4,4})

% plot ripple rate
PlotRipRate(pos21{1},Inpath{4,1})

% plot ripple zscore power in frames
PlotFrameRipZ(pos21{2},Inpath{4,2})


% plot two parameter replay
Plot2DReplay(pos22,Inpath{4,5})

% plot e-ket session decoding
PlotSes2Dcd1(pos12,AllInpath31)

% plot e-ket session decoding
PlotSes2Dcd2(pos12,AllInpath32)

% plot e-ket session decoding
PlotSes2Dcd3(pos12,AllInpath32)

% plot control session decoding
PlotSes1Dcd1(pos11,AllInpath31)

% plot control session decoding
PlotSes1Dcd2(pos11,AllInpath32)

% plot control session decoding
PlotSes1Dcd3(pos11,AllInpath31)


%% save the result
File_Path = strcat(Outpath,figurename,'.fig');
saveas(gcf, File_Path);
File_Path = strcat(Outpath,figurename,'.png');
print(gcf,File_Path,'-dpng','-r300'); 
% export_fig([Outpath,figurename,'.pdf'])
print('-painters',gcf,'-dpdf',[Outpath,figurename,'_AI.pdf']) % this is for illustrator
close all
%% subpanel plot function

function Plot2DReplay(pos,Inpath)
global txtfonts sescolors sesstr
colors = Color_FrenchDisp;
cellstr = 'NCell5';

load([Inpath,'bsframes',cellstr,'.mat'],'allf','allsfcorr');
spacri = 0;
nsess = 4;

%% make the plots
im = 2;
seqrange = 0:0.1:0.9;
jprange = 0.1:0.1:1;
jumpcol = [11,12];
for is = 1:nsess
    sesnow = KetSesMap(is);
    sesidx = allf(:,1) == sesnow & abs(allf(:,13)) >= spacri ;
    tempf = allf(sesidx,:);
    sesidx = allsfcorr(:,1) == sesnow & abs(allsfcorr(:,13)) >= spacri ;
    tempshff = allsfcorr(sesidx,:);
    
    ax1 = subplot('position',pos{is});
    TwoProportionZtest2D_v3(tempf(:,[jumpcol(im),9]),tempshff(:,[jumpcol(im),9]),...
        jprange,seqrange,0,1,ax1);
    
    set(ax1,'fontsize',txtfonts(4))
    if is == 1
        ylb = ylabel(ax1,'Absolute correlation','FontSize',txtfonts(3));
        ylb.Position(2) = ylb.Position(2)+0.6;
        ylb.Position(1) = ylb.Position(1)-0.05;
    end
    title(ax1,sesstr{is},'FontWeight','normal','FontSize',txtfonts(3))
    ax1.YDir = 'reverse';
    xticks(ax1,[0.2 0.5 0.8])
    yticks(ax1,[0 0.3 0.6 0.9])
    
    if is == 1 || is == 3
        yticklabels(ax1,{'>0','>0.3','>0.6','>0.9'})
    else
        yticklabels(ax1,[])
    end
    if is > 2
        xticklabels(ax1,{'<0.2','<0.5','<0.8'})
%         xtickangle(ax1,45)
    else
        xticklabels(ax1,[])
    end
    
    hold(ax1,'on')
    if is == 3
        xlb = xlabel(ax1,'Max jump','FontSize',txtfonts(3));
        xlb.Position(1) = xlb.Position(1)+0.6;
%         xlb.Position(1) = xlb.Position(1);
    end
    set(ax1,'TickLength',[0.05, 0.05])
    set(ax1,'layer','top')
    rectangle(ax1,'Position',[0.05 0.55 0.4 0.4],'EdgeColor',[0.05 0.05 0.05],'linewidth',1.5)
    if is == 4
        cb = colorbar;
        cb.Location = 'eastoutside';
        cb.Position = [0.434 0.051 0.0105 0.091];
        cb.Label.String = '   p data>shuffle';
        
        cb.FontSize = txtfonts(3);
        cb.Ticks = [0 0.01 0.02 0.03 0.04 0.05];
        hl = plot(ax1,nan,nan,'s','markeredgecolor','k','markersize',8);
        [lgd,icon] = legend(hl,'NoData');
        icon(1).Position(1) = 0.4;
        legend boxoff
        set(lgd,'Position',[0.4065 0.145 0.1 0.0196]);
    end
end
end

function PlotDcdProb(pos,Inpath)
global txtfonts sescolors sesstr cmap
framespk = 2;
tdrcri = 1000; % frames with theta/delta ratio < ...
ripzcri = -100; % frames with rip zscore > ...
nsess = 4;
flfpstr = ['TDcri',num2str(tdrcri),'RipZcri',num2str(ripzcri)];
% cmap=getNCLColmap('MPL_BuPu.rgb',100);
% juvenile group

load([Inpath,'allf.mat'],'allf')
mcols = [15,16,17,18,19]; % mean, median, sum, max, min logpb
im = 1;
ax = subplot('position',pos);
hold(ax,'on')
dcdf = cell(1,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    sesidx = allf(:,1) == sid &  allf(:,5) >= framespk & allf(:,20) < tdrcri & allf(:,21) > ripzcri;
    tempf = allf(sesidx,:);
    lgpb = tempf(:,mcols(im));
    dcdf{is} = lgpb;
end
[groups,pvals] = AllGroupRankSum(dcdf);
pidx = [1,4,5];
MyBoxPlot_v3(dcdf,ax,'sigplot',0,'color',sescolors)
ylim(ax,[2 5.6])
yl = ylim(ax);
xticks(ax,1:nsess)
% xlabel(ax,'Sessions')
xticklabels(ax,sesstr)
ylabel(ax,'Log sum probability ')
myoverlapsigstar_v3(ax,groups(pidx),pvals(pidx),'ygap',0.02*diff(yl))
set(ax,'fontsize',txtfonts(3))
set(ax, 'YGrid', 'on', 'XGrid', 'off')
xtickangle(ax,30)
% xlabel(ax,'Sessions')
end

function PlotDcdSigRatio(pos,Inpath)
global txtfonts sescolors sesstr
medjp = 100;
covercri = 0;
cellstr = 'NCell5';
tdrcri = 1000; % frames with theta/delta ratio < ...
ripzcri = -100; % frames with rip zscore > ...
flfpstr = ['TDcri',num2str(tdrcri),'RipZcri',num2str(ripzcri)];
parastr = ['_MedJPct',num2str(medjp),'SpaCover',num2str(covercri),cellstr,flfpstr];
cmap=getNCLColmap('MPL_BuPu.rgb',100);
nsess = 4;
% juvenile group

load([Inpath,'sigr',parastr,'.mat'],'sigr2d','shfsigr2d');
sigr2d = permute(sigr2d,[2,1,3]);
shfsigr2d = permute(shfsigr2d,[2,1,3]);
sigr2d = reshape(sigr2d,nsess,[]);
shfsigr2d = reshape(shfsigr2d,nsess,[]);

ax = subplot('position',pos);
hold(ax,'on')
p = nan(1,nsess);
g = cell(1,nsess);
for is = 1:nsess
   sid = KetSesMap(is);
   bar(ax,is+0.15,nanmean(sigr2d(sid,:)),0.3,'facecolor',sescolors{is},'edgecolor',sescolors{is},'facealpha',0.6)
   bar(ax,is-0.15,nanmean(shfsigr2d(sid,:)),0.3,'facecolor',[0.7 0.7 0.7],'edgecolor',[0.7 0.7 0.7],'facealpha',0.6)
   for isample = 1:size(sigr2d,2)
      plotx1 = is+0.15;
      plotx1 = plotx1 + (rand-0.5)/8;
      
      plotx2 = is-0.15;
      plotx2 = plotx2 + (rand-0.5)/8;
      plot(ax,plotx1,sigr2d(sid,isample),'o',...
          'markerfacecolor',0.8*sescolors{is},'markeredgecolor',0.8*sescolors{is},'markersize',2)
      plot(ax,plotx2,shfsigr2d(sid,isample),'o',...
          'markerfacecolor',[0.4 0.42 0.4],'markeredgecolor',[0.4 0.42 0.4],'markersize',2)
   end
   addSEMBartoPlot(ax,is+0.15,{sigr2d(sid,:)},'capsize',0.15,'linewidth',1)
   addSEMBartoPlot(ax,is-0.15,{shfsigr2d(sid,:)},'capsize',0.15,'linewidth',1)
   [~,ptmp] = ttest(sigr2d(sid,:)-shfsigr2d(sid,:),0);
   %ptmp = signrank(sigr2d(is,:)-shfsigr2d(is,:),0);
   p(is) = ptmp;
   g{is} = [is-0.15,is+0.15];
end
sigstar(g,p,0,7)

xticks(ax,1:nsess)
xlim(ax,[0.5,nsess+0.5])
xticklabels(ax,sesstr)
ylabel(ax,'Significant ratio')
%xlabel(ax,'Sessions')
xtickangle(ax,30)
% title(ax,'Juvenile','FontWeight','normal','fontsize',txtfonts(3))
set(ax,'fontsize',txtfonts(3))
set(ax, 'YGrid', 'on', 'XGrid', 'off')

% adult group
% load([Inpath2,'sigr',parastr,'.mat'],'sigr2d','shfsigr2d');
% sigr2d = permute(sigr2d,[2,1,3]);
% shfsigr2d = permute(shfsigr2d,[2,1,3]);
% sigr2d = reshape(sigr2d,nsess,[]);
% shfsigr2d = reshape(shfsigr2d,nsess,[]);
% 
% ax = subplot('position',pos{2});
% hold(ax,'on')
% p = nan(1,nsess);
% g = cell(1,nsess);
% for is = 1:nsess
%    bar(ax,is+0.15,nanmean(sigr2d(is,:)),0.3,'facecolor',sescolors{is},'edgecolor',sescolors{is})
%    bar(ax,is-0.15,nanmean(shfsigr2d(is,:)),0.3,'facecolor',[0.7 0.7 0.7],'edgecolor',[0.7 0.7 0.7])
%    for isample = 1:size(sigr2d,2)
%       plotx = [is+0.15,is-0.15];
%       plotx = plotx + (rand-0.5)/8;
%       plot(ax,plotx,[sigr2d(is,isample),shfsigr2d(is,isample)],'o-','color',[0.6 0.6 0.6],'markersize',3)
%    end
%    addSEMBartoPlot(ax,is+0.15,{sigr2d(is,:)},'capsize',0.15,'linewidth',1)
%    addSEMBartoPlot(ax,is-0.15,{shfsigr2d(is,:)},'capsize',0.15,'linewidth',1)
%    [~,ptmp] = ttest(sigr2d(is,:)-shfsigr2d(is,:),0);
%    %ptmp = signrank(sigr2d(is,:)-shfsigr2d(is,:),0);
%    p(is) = ptmp/2;
%    g{is} = [is-0.15,is+0.15];
% end
% sigstar(g,p,0,7)
% 
% xticks(ax,1:nsess)
% xticklabels(ax,{'Ctrl','E-Ket','L-Ket'})
% % ylabel(ax,'Significant Ratio')
% xlabel(ax,'Sessions')
% title(ax,'Adult','FontWeight','normal','fontsize',txtfonts(3))
% set(ax,'fontsize',txtfonts(3))

plot(ax,0.74,0.39,'s','markerfacecolor',sescolors{1},'markeredgecolor',sescolors{1},'markersize',3.5)
plot(ax,0.92,0.39,'s','markerfacecolor',sescolors{2},'markeredgecolor',sescolors{2},'markersize',3.5)
plot(ax,1.08,0.39,'s','markerfacecolor',sescolors{3},'markeredgecolor',sescolors{3},'markersize',3.5)
plot(ax,1.26,0.39,'s','markerfacecolor',sescolors{4},'markeredgecolor',sescolors{4},'markersize',3.5)
text(ax,1.4,0.39,'Data','HorizontalAlignment','left','FontSize',txtfonts(3)-0.5)

plot(ax,1,0.36,'s','markerfacecolor',[0.7 0.7 0.7],'markeredgecolor',[0.7 0.7 0.7],'markersize',3.5)
text(ax,1.4,0.358,'Shuffle','HorizontalAlignment','left','FontSize',txtfonts(3)-0.5)
end


function PlotSes1Dcd1(pos,AllInpath1)
global txtfonts sescolors sesstr
load ([AllInpath1{5},'Ptet.mat'],'Ptet')

winlen = 0.04;
allses = [1];

alldir = [1];
alltime = {[10510 10511]};

posnow = pos{1};

addcb = 0;
addylb = 1;
iexample = 1;
addxlb = 0;
ttstr = [];
timerangedecoding(posnow,AllInpath1,Ptet,allses(iexample),alltime{iexample},alldir(iexample),winlen,addcb,addylb,ttstr,addxlb)

end

function PlotSes1Dcd2(pos,AllInpath1)
global txtfonts sescolors sesstr
load ([AllInpath1{5},'Ptet.mat'],'Ptet')

winlen = 0.04;
allses = [1];

alldir = [2];
alltime = {[14494.3 14495.3]};
posnow = pos{2};

addcb = 0;
addylb = 0;
iexample = 1;
addxlb = 0;
ttstr = [sesstr{1},' Decoding Examples'];
timerangedecoding(posnow,AllInpath1,Ptet,allses(iexample),alltime{iexample},alldir(iexample),winlen,addcb,addylb,ttstr,addxlb)

end

function PlotSes1Dcd3(pos,AllInpath1)
global txtfonts sescolors sesstr
load ([AllInpath1{5},'Ptet.mat'],'Ptet')

winlen = 0.04;
allses = [1];

alldir = [2];
alltime = {[9734.2 9735.2]};
alltime = {[10434.8 10435.8]};
posnow = pos{3};

addcb = 1;
addylb = 0;
iexample = 1;
addxlb = 0;
ttstr = [];
timerangedecoding(posnow,AllInpath1,Ptet,allses(iexample),alltime{iexample},alldir(iexample),winlen,addcb,addylb,ttstr,addxlb)

end

function PlotSes2Dcd1(pos,AllInpath1)
global txtfonts sescolors sesstr
load ([AllInpath1{5},'Ptet.mat'],'Ptet')

winlen = 0.04;
allses = [2];

alldir = [2];
alltime = {[19235.3 19236.3]};

posnow = pos{1};

addcb = 0;
addylb = 1;
iexample = 1;
addxlb = 1;
ttstr = [];
timerangedecoding(posnow,AllInpath1,Ptet,allses(iexample),alltime{iexample},alldir(iexample),winlen,addcb,addylb,ttstr,addxlb)
end

function PlotSes2Dcd2(pos,AllInpath1)
global txtfonts sescolors sesstr
load ([AllInpath1{5},'Ptet.mat'],'Ptet')

winlen = 0.04;
allses = [2];

alldir = [2];
alltime = {[20414.2 20415.2]};
alltime = {[20395 20396]};
posnow = pos{2};

addcb = 0;
addylb = 0;
iexample = 1;
addxlb = 1;
ttstr = [sesstr{2},' Decoding Examples'];
timerangedecoding(posnow,AllInpath1,Ptet,allses(iexample),alltime{iexample},alldir(iexample),winlen,addcb,addylb,ttstr,addxlb)
end

function PlotSes2Dcd3(pos,AllInpath1)
global txtfonts sescolors sesstr
load ([AllInpath1{5},'Ptet.mat'],'Ptet')

winlen = 0.04;
allses = [2];

alldir = [2];
alltime = {[20563.5 20564.5]};
alltime = {[20589 20590]};
% alltime = {[19754.5 19755.5]};
posnow = pos{3};

addcb = 1;
addylb = 0;
iexample = 1;
addxlb = 1;
ttstr = [];
timerangedecoding(posnow,AllInpath1,Ptet,allses(iexample),alltime{iexample},alldir(iexample),winlen,addcb,addylb,ttstr,addxlb)
end


function PlotRipRate(pos,Inpath)
global txtfonts sescolors sesstr
ripstr = ['sup',num2str(5),'st',num2str(2),'_mcheck'];

% juvenile group
load([Inpath,'ripxses',ripstr,'.mat'],'riplenses','ripspec','ripwl','wlff','wltt','ripspkses','riprate','ripactcses','ripmlfp')
nsess = size(riprate,2);

ax = subplot('position',pos);
hold(ax,'on')
set(ax,'fontsize',txtfonts(4))
riprate_t = nan(5,nsess);
plotd = cell(1,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    riprate_t(:,is) = riprate(:,sid);
    plotd{is} = riprate(:,sid);
end
[groups,pvals] = AllGroupTtest(riprate_t);
pidx = [1,4,5];
MyBarSampleSEMPlot_v2(plotd,ax,'colors',sescolors,'barwidth',0.35)
% ylim(ax,[2 5.6])
% yl = ylim(ax);
xticks(ax,1:nsess)
xlim(ax,[0.3 nsess+0.7])
% xlabel(ax,'Sessions')
xticklabels(ax,[])
ylabel(ax,'Ripple (counts/sec)','fontsize',txtfonts(4))
myoverlapsigstar_v3(ax,groups(pidx),pvals(pidx),'starty',0.082,'ygap',0.003,'starsize',7,'nsgap',0.9)
set(ax, 'YGrid', 'on', 'XGrid', 'off')
% xlabel(ax,'Sessions')

end

function PlotFrameRipZ(pos,Inpath)
global txtfonts sescolors sesstr
cellstr = 'NCell5';
% juvenile group

load([Inpath,'fstats',cellstr,'.mat'],'fstats','frate','frate_state','muakurtsis','friprate'); 
nsess = size(fstats,2);

im = 4;
ax = subplot('position',pos);
hold(ax,'on')
set(ax,'fontsize',txtfonts(4))
plotd = cell(1,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    plotd{is} = fstats{im,sid};
end
[groups,pvals] = AllGroupRankSum(plotd);
pidx = [1,4,5];
MyBoxPlot_v3(plotd,ax,'sigplot',0,'color',sescolors)

xticks(ax,1:nsess)
xlim(ax,[0.3 nsess+0.7])
yl = ylim;
ylim(ax,[yl(1),6])
% xlabel(ax,'Sessions')
xticklabels(ax,sesstr)
ylb = ylabel(ax,'Zscore ripple power','fontsize',txtfonts(4));
ylb.Position(1) = ylb.Position(1) - 0.4; 
myoverlapsigstar_v3(ax,groups(pidx),pvals(pidx),'starty',4.8,'starsize',7)
set(ax, 'YGrid', 'on', 'XGrid', 'off')
% xlabel(ax,'Sessions','fontsize',txtfonts(3))
xtickangle(ax,30)
end




function timerangedecoding(pos,allinpath,Ptet,sesn,trange,idir,winlen,addcb,addylb,ttstr,addxlb)

global txtfonts sescolors cmap
method = 'Exp_Var';
ploth = 0.04;

ld = load ([allinpath{1},'ses.mat'],'ses');
ld2 = load ([allinpath{2},'M12ClFThrSpeedVarB_note.mat'],'M1','M1_note_param','PlFields','ClF','PlFmesh','M2');
load ([allinpath{3},'Maze.mat'],'mazel','times')
load ([allinpath{4},'thetapk.mat'],'thetapk')
% cmap{1} = getNCLColmap('MPL_BuGn.rgb',100);
% cmap{2} = getNCLColmap('MPL_PuRd.rgb',100); % this is for current detour
% cmap{3} = getNCLColmap('MPL_PuBu.rgb',100);
% cmap{4} = getNCLColmap('MPL_Purples.rgb',100);
% cmap{5} = getNCLColmap('MPL_Greys.rgb',100); % this is for all other tracks
cmap2=getNCLColmap('MPL_afmhot.rgb',100);
colors = Color_FrenchDisp;
nsess = size(ld2.M1,4);

nactcri = 1; % only show bins with active cells not less than this value
stdcri = 1; % log10(sumpb) mean + std criteria

% load ([allinpath{5},'Ptet.mat'],'Ptet')
try
    load ([allinpath{6},'sesrip_sup5st2_mcheck.mat'],'sesrip')
catch
    load ([allinpath{6},'sesrip_sup5st2.mat'],'sesrip')
end
load ([allinpath{7},'MazeNCell5.mat'],'fratime')
load ([allinpath{8},'sdtangdiff.mat'],'sdtangdiff','sdtangdiffmu','sdtangdiffvar','segmesh')
load ([allinpath{9},'seslambdaNCell5.mat'],'seslambda')
%% preprocess
% exclude int neurons in this analysis
[M1,M2,ClF,PlFields,nc,pyrid,clu2pyr] = CA_ExcludeInt(ld2.M1,ld2.M2,ld2.ClF,ld2.PlFields);
PlFmesh = ld2.PlFmesh;
% order tracks by linear pos
ses = Detour_Ordertracks(ld.ses);

ntrack = length(ses(1).tra_p);
lambdascale = 1; % this would adjust the sharpness of decoding result
nsample = 1e4;
%% define parameters
tbin = winlen; % decoding window length
% % dcdtinc = 0.004; % 4 ms window increasement
% % overlap = (tbin-dcdtinc)/tbin;
% 50% overlap
% if winlen > 0.1
%     overlap = 0.5;
% else
%     dcdtinc = 0.004; % 4 ms window increasement
%     overlap = (tbin-dcdtinc)/tbin;
% end
overlap = 0.9;
%% do the decoding

alldecode = struct;
% D1 direction, D2 run session,
is = sesn;
% ses(is).tralim(end,2) = ses(is).tralim(end,2) - 8;
% get behavior data
% find laps on track of interest
smazel = mazel{is};
mazelt = smazel(:,6);
% remove redundant timestamps due to track overlap
[mazelt,Iuni,~] = unique(mazelt);
smazel = smazel(Iuni,:);% time from pos recording
mazelt = smazel(:,6);
mazellpos = smazel(:,4); % linear position
mazeltra = smazel(:,3); % track number
mazelvel = smazel(:,7); % track velocity
mazelldir = smazel(:,8); % lap direction

% time bin spikes in this session
cells = struct;
for ic = 1:nc
    cells(ic).time = ClF{ic,is}(:,1);
end
tbins = get_tbins(mazelt(1), mazelt(end), tbin, overlap);
nbins = size(tbins,1);
% take 1e4 samples
subsample = randperm(nbins);
subsample = subsample(1:nsample);
subtbin = tbins(subsample,:);

spike_counts= zeros(length(cells), nsample);
for i=1:length(cells)
    spike_counts(i,:) = histcounts_row_kf(cells(i).time, subtbin);
end
meanact = nanmean(spike_counts,2);
% get angle with average activity
allbinang = nan(1,size(spike_counts,2));
for ibin = 1:size(spike_counts,2)
    allbinang(ibin) = AngleBetweenV(spike_counts(:,ibin),meanact);
end

alllambda = seslambda(is);

plftmp = squeeze(PlFields(:,is,idir,:));
tempmesh = PlFmesh;

% add tracks to decoder
addses = is*ones(1,ntrack);
addlpos = cell(1,ntrack);
addlambda = cell(1,ntrack);
lambdamesh = cell(1,ntrack);
for it = 1:ntrack
    mpos = sdtangdiffmu{idir,is,it};
    varpos = sdtangdiffvar{idir,is,it};
    meshnow = segmesh{idir,is,it};
    
    if strcmp(method,'Exp_Mean')
        poslambda = 1./mpos;
    end
    if strcmp(method,'Exp_Var')
        poslambda = sqrt(1./varpos);
    end
    if strcmp(method,'ChiSq')
        poslambda = mpos;
    end
    
    addlpos{it} = [ses(is).tralim(it,1), ses(is).tralim(it,2)];
    poslambda = poslambda * lambdascale;
    addlambda{it} = poslambda;
    lambdamesh{it} = meshnow;
end
% now we will construct the deocoding template
allsesmesh = [];
allsesplf = [];
allseslambda = [];
meshlencount = 0;
lposlimit = [];
for iseg = 1:length(addses)
    [allsesmesh,allsesplf,allseslambda,meshlencount,lposlimit] = AddLpostoDecodeTemplate(addses(iseg),addlpos{iseg},addlambda{iseg},...
        plftmp,tempmesh,lambdamesh{iseg},meshlencount,allsesmesh,allsesplf,allseslambda,lposlimit,ses);
end

% do decoding
trangenow = trange;
ClFnow = ClF(:,is);

% run the decoding for the concatenate plfield
if strcmp(method,'Exp_Mean') || strcmp(method,'Exp_Var')
    sframes = BayePosDecode_LabCode_Lite_CosineExp(allsesplf,allseslambda,allsesmesh,...
        ClFnow,trangenow,[min(allsesmesh) max(allsesmesh)],tbin,overlap,meanact,alllambda);
end
if strcmp(method,'ChiSq')
    sframes = BayePosDecode_LabCode_Lite_CosineChiSq(allsesplf,allseslambda,allsesmesh,...
        ClFnow,trangenow,[min(allsesmesh) max(allsesmesh)],tbin,overlap,meanact,alllambda);
end


normdcdpdf = BayePosDecode_LabCode_Lite(allsesplf,allsesmesh,ClFnow,trangenow,...
    [min(allsesmesh) max(allsesmesh)],tbin,overlap);


for iseg = 1:size(lposlimit,1)
    jssparange = [lposlimit(iseg,1),lposlimit(iseg,2)];
    % spatial bin range for this session
    dprob = sframes.pdf; %D1 is spatial bins, D2 is temporal bins
    dprobnorm = normdcdpdf.pdf; %D1 is spatial bins, D2 is temporal bins
    
    spamesh = sframes.spacebin;
    spajsind = spamesh>=jssparange(1) & spamesh<=jssparange(2);
    dprob = dprob(spajsind,:);
    dprobnorm = dprobnorm(spajsind,:);
    spamesh = spamesh(spajsind);
    % we want to map back to the linear position of the track in the session
    spamesh = spamesh - min(spamesh);
    orilimit1 = addlpos{iseg}(1);
    spamesh = spamesh + orilimit1;
    
    alldecode.seg(iseg).pdf = dprob;
    alldecode.seg(iseg).pdfnorm = dprobnorm;
    alldecode.seg(iseg).tbin = sframes.tbin;
    alldecode.seg(iseg).spamesh = spamesh;
    alldecode.seg(iseg).nact = sframes.nact;
    alldecode.seg(iseg).nspk = sframes.nspk;
    % get velocity
    ctbin = nanmean(sframes.tbin,2);
    alldecode.seg(iseg).vel = interp1(mazelt,abs(mazelvel),ctbin);
end


%% do the plot
ftimesnow = fratime{is};
ftimecenter = nanmean(ftimesnow,2);
smazel = mazel{is};
mazelt = smazel(:,6); % time from pos recording
mazellpos = smazel(:,4); % linear position
try
    goodrip = sesrip{is}(:,5) == 1;
    rippknow = sesrip{is}(goodrip,2); % peak time of ripples in this session
catch
    rippknow = sesrip{is}(:,2); % peak time of ripples in this session
end

tmpdcd = alldecode.seg;

nseg = length(tmpdcd);

dcdtbin = alldecode.seg(1).tbin;
alltc = nanmean(dcdtbin,2);
tcstep = nanmean(diff(alltc));
tplotrange = [alltc(1)-tcstep/2, alltc(end)+tcstep/2];

%% plot decoding results
js = 1;
pdf = tmpdcd(js).pdf;
nact = tmpdcd(js).nact;
less2idx = nact < nactcri;
pdf(:,less2idx) = nan;

binsum = sum(pdf,1);
if js == 1
    segsum = binsum;
else
    segsum = segsum + binsum;
end

lpos = tmpdcd(js).spamesh;
pdf = tmpdcd(js).pdfnorm;
nact = tmpdcd(js).nact;
less2idx = nact < nactcri;
pdf(:,less2idx) = nan;
postemp = pos;
postemp(2) = postemp(2) + 0.01;
postemp(4) = postemp(4) - 0.01;
ax = subplot('position',pos);
hold(ax,'on')

% normalize within time bins
% binsum2 = repmat(segsum,size(pdf,1),1);
% pdfnorm = pdf./binsum2;

imagesc(ax,[alltc(1) alltc(end)],lpos([1 end]),pdf)
caxis(ax,[0 0.1])
colormap(ax,cmap)
drawnow

tidx = mazelt >= tplotrange(1)-3 & mazelt <= tplotrange(end)+3;
posnow = mazellpos(tidx);
MyFilt=fir1(30,0.5/(30/2));
posnow = Filter0(MyFilt,posnow);
plottime = mazelt(tidx);
plot(ax,plottime,posnow,':','color',[0.22 0.35 0.52],'linewidth',2)
ylim(ax,[lpos(1)-5 lpos(end)+2])
if js == 1
    ylim(ax,[lpos(1)-6 lpos(end)])
end

if js == nseg
    ylim(ax,[lpos(1) lpos(end)+6])
end
xlim(ax,tplotrange)
yticks(ax,[0 50 100 150])
set(ax,'fontsize',txtfonts(4))
if addylb == 1
    ylabel(ax,'Position (cm)','fontsize',txtfonts(4))
else
    yticklabels(ax,[])
end
box(ax,'on')
set(ax,'layer','top')
xticks(ax,[])
set(ax,'TickLength',[0.002, 0.002])
if js == nseg && addcb
    cb = addcolorbar;
    cb.Label.String = ' Normalized Probability';
    cb.Ticks = [0 0.1];
end

if js == nseg && ~isempty(ttstr)
%     title(ax,ttstr,'FontWeight','normal','FontSize',txtfonts(2))
    title(ax,ttstr,'FontSize',txtfonts(2),'FontWeight','normal')
end
ylim(ax,[-5 155])
    
%% plot pb
postemp = pos;
postemp(2) = postemp(2) - ploth;
postemp(4) = ploth;
axpb = subplot('position',postemp);
hold(axpb,'on')
% plot the sum of pb
binsumlog = log10(segsum);
binsumlog = log(segsum);
plot(axpb,alltc,binsumlog,'linewidth',1.25,'color',cmap(70,:))
box(axpb,'on')
set(axpb,'layer','top')
set(axpb,'TickLength',[0.002, 0.002])
xlim(axpb,tplotrange)
xticklabels(axpb,[])
set(axpb,'fontsize',txtfonts(4))
if addylb == 1
    ylabel(axpb,'log\Sigma Probability','fontsize',txtfonts(4)-1.5)
else
    yticklabels(axpb,[])
end
%     ylim(axpb,[-2 2])
%     yticks(axpb,[-2 0 2])
%     plot(axpb,tplotrange,[0 0],'--k')
%     ylim(axpb,[1 2])
%     yticks(axpb,[1 1.5 2])
%     plot(axpb,tplotrange,[1.5 1.5],'--k')

ylim(axpb,[1.8 5.2])
yticks(axpb,2:4)
set(axpb, 'YGrid', 'on', 'XGrid', 'off')

%% plot raster
postemp = pos;
postemp(2) = postemp(2) - ploth*1.8;
postemp(4) = ploth*0.8;
ax = subplot('position',postemp);
hold(ax,'on')

sclf = ClF(:,is);
myrasterplot(ax,sclf,tplotrange,'markersize',1.5)
yl = ylim(ax);

box(ax,'on')
set(ax,'layer','top')
set(ax,'TickLength',[0.002, 0.002])
xlim(ax,tplotrange)
xticklabels(ax,[])
% high light frames in this time range
fidx = find(idxinrange(ftimecenter,tplotrange));
if ~isempty(fidx)
    for ift = 1:length(fidx)
        ftimethis =  ftimesnow(fidx(ift),:);
        xxfill = [ftimethis,fliplr(ftimethis)];
        yyfill = [yl(1),yl(1),yl(2),yl(2)];
        fill(ax,xxfill,yyfill,[0.42,0.32,0.68],'LineStyle','none','FaceAlpha',0.4)
    end
end
if  addylb == 1
    ylabel(ax,'Pyr','fontsize',txtfonts(4))
end

%% plot row LFP
postemp = pos;
postemp(2) = postemp(2) - ploth*2.6;
postemp(4) = ploth*0.8;
ax = subplot('position',postemp);
hold(ax,'on')

lfpt = Ptet(is).t;
lfpnow = Ptet(is).LFP;
fs = Ptet(is).fs;
% high pass filter 1 Hz
htband = [2 400];
MyFilt=fir1(20,htband/(fs/2));
lfpnow = Filter0(MyFilt,lfpnow);
                
                
tidx = idxinrange(lfpt,tplotrange);
plot(ax,lfpt(tidx),lfpnow(tidx),'-k','linewidth',1);
box(ax,'on')
set(ax,'layer','top')
set(ax,'TickLength',[0.002, 0.002])
yticks(ax,[])
set(ax,'fontsize',txtfonts(3))
if addylb == 1
    ylabel(ax,'LFP','fontsize',txtfonts(4))
else
    yticklabels(ax,[])
end
xlim(ax,tplotrange)
xticklabels(ax,[])
yl = ylim(ax);
% plot ripple
% tmprippktidx = idxinrange(rippknow,tplotrange);
% tmprippk = rippknow(tmprippktidx);
% if ~isempty(tmprippk)
%     yplot = ones(size(tmprippk)).* (yl(2)*0.9+yl(1)*0.1);
%     plot(ax,tmprippk,yplot,'r*','linewidth',1,'markersize',3)
% end


%% plot wavelet results
postemp = pos;
postemp(2) = postemp(2) - ploth*3.7;
postemp(4) = ploth*1.1;
ax = subplot('position',postemp);
hold(ax,'on')
wltange = [tplotrange(1) - 1,tplotrange(end) + 1];
% wltange = [tplotrange(1),tplotrange(end)];
wltidx =  idxinrange(lfpt,wltange);
[ss,ff] = cwt(lfpnow(wltidx),Ptet(is).fs);
ff = log10(ff);
% %             cwt(lfp,Ltemp.fs)
% %             imagesc([ts(1) ts(end)],[ff(1) ff(end)],abs(ss))
pcolor(ax,lfpt(wltidx),ff,abs(ss))
ytk = [8,32,128];
ytklabels = cell(1,length(ytk));
for ll = 1:length(ytk)
    ytklabels{ll} = num2str(ytk(ll));
end
shading flat
colormap(ax,cmap2)
% %             relat = [min(laptime)-min(ts) max(laptime)-min(ts)];
ylim(ax,log10([3 256]))
yticks(ax,log10(ytk))
yticklabels(ax,ytklabels)
caxis([1e-5 1e-4])
box(ax,'on')
set(ax,'layer','top')
set(ax,'TickLength',[0.002, 0.002])
xlim(ax,tplotrange)
set(ax,'fontsize',txtfonts(4))
if addylb == 1
    ylabel(ax,'Freq (Hz)','fontsize',txtfonts(4))
else
    yticklabels(ax,[])
end

if addcb
   cb = addcolorbar;
   cb.Label.String = 'Power';
end
xticks(ax,[tplotrange(1)+0.03 tplotrange(end)-0.03])
xticklabels(ax,{num2str(tplotrange(1),'%.1f'), num2str(tplotrange(end),'%.1f')})
if addxlb
    xlb = xlabel(ax,'Time (s)');
    xlb.Position(2) = xlb.Position(2) + 0.6;
end
end

function [allsesmesh,allsesplf,allseslambda,meshlencount,lposlimit] = AddLpostoDecodeTemplate(is,lposb,lambda,...
    PlFields,PlFmesh,lambdamesh,meshlencount,allsesmesh,allsesplf,allseslambda,lposlimit,ses)

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

lambdaidx = idxinrange(lambdamesh,lposb);
tmplambda = lambda(lambdaidx);

templimit = [meshlencount,meshlencount + max(meshtemp) + meshdiff];
lposlimit = cat(1,lposlimit,templimit);

allsesmesh = cat(2,allsesmesh,meshtemp+meshlencount);
allsesplf = cat(2,allsesplf,plftemp);
allseslambda = cat(1,allseslambda,tmplambda);
meshlencount = meshlencount + max(meshtemp) + meshdiff;

end

