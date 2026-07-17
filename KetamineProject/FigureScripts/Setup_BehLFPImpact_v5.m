clear;clc;close all


%% define subpanel positions
testpos = 1;

% rest run rest protocal for juvenile
pos11 = [0.02 0.79 0.32 0.15];
% showing behavior at different vel range
[~,pos12] = tight_subplot(3,1,[0.0085 0.01],[0.81 0.02],[0.40 0.435],0);
pos12{4} = [0.610 0.815 0.105 0.145];
[~,pos13] = tight_subplot(3,1,[0.0085 0.01],[0.81 0.02],[0.785 0.05],0);

% vel and angle of vel
[~,pos21] = tight_subplot(1,3,[0.035 0.05],[0.55 0.28],[0.052 0.42],0);
% deceleration and acceleration
[~,pos22] = tight_subplot(1,2,[0.035 0.05],[0.55 0.28],[0.67 0.03],0);

% LFP
[~,pos31] = tight_subplot(1,4,[0.035 0.01],[0.28 0.55],[0.052 0.27],0); 
pos32 = [0.80 0.28 0.17 0.19];
% test position
if testpos
f = figure_letter(0.55);
testsubpos({pos11})
testsubpos(pos12)
testsubpos(pos21)
testsubpos(pos22)
testsubpos(pos31)
testsubpos(pos32)
end
close all
f = figure_letter(0.55);
%% define parameters
subpanelfont = 13;
labelfont = 7;
legendfont = 5;
smallfont = 4;

global txtfonts sescolors sesstr cmap
txtfonts = [subpanelfont,labelfont,legendfont,smallfont];
sescolors = Color_Ketamine_v2;
sesstr = {'Run1','Run2Bi','Run2woBi','Run3'};
cmap = getNCLColmap('WhiteYellowOrangeRed.rgb',100);

%% add subpanel annotations
addsubp = 1;
if addsubp
    panellbpos = {[0 1],[0.32 1],[0 0.775],[0.605 0.77],[0 0.515],[0.75 0.515]};
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
Inpath = cell(10,6); % 10 subpanels, each could have up to 6 data source
rnum = "URat4";
setnum = 1;
[Datasets] = getCADatasets(rnum,setnum); %
Datasets = Datasets{1};
Inpath{1,1} =  [getDataRoot(),'CAResults/Preprocess/KTM_Run/',Datasets,'/Step06/']; % ses
Inpath{1,2} =  [getDataRoot(),'CAResults/Preprocess/KTM_Run/',Datasets,'/Step04/']; % plf
Inpath{1,3} =  [getDataRoot(),'CAResults/Preprocess/KTM_Run/',Datasets,'/Step03/']; % mazel
%Inpath{1,4} =  [getDataRoot(),'CAResults/Solver/KTM_Run/URat4/kd1_2024-05-23_16-33-37/InjectTAlign/VelCDFvsRun13/Twin300Step120/']; %run2split
Inpath{1,4} =  [getDataRoot(),'CAResults/Solver/KTM_Run/',Datasets,'/InjectTAlign/AllRun_VelCDFvsRun1/Twin120Step60/']; %run2split

Inpath{2,1} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_RunBeh/KTMUActNew/KTM_AbsDisperMinute_KetImpactRun2/Twin300Step120/']; % vel and angle of vel
Inpath{2,2} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_RunBeh/KTMUActNew/KTMVel_SigImpactRange_byDir_KetImpactRun2/Twin300Step120/']; % vel and angle of vel
Inpath{2,3} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_RunBeh/KTMUActNew/KTMVel_MaxVel_byDir_KetImpactRun2/Twin300Step120/']; % vel and angle of vel
Inpath{2,4} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/RemapDrift/RunSeq/MazeSeq_EachSes_KetImpactRun2/Twin300Step120Vel10/'];% maze sequence of juvenile
Inpath{2,5} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_RunBeh/KTMUXDNew/XSesVel_Dis2Corner_Decelerate_KetImpactRun2/Twin300Step120/']; % vel during deceleration
Inpath{2,6} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_RunBeh/KTMUXDNew/XSesVel_Dis2Corner_Accelerate_KetImpactRun2/Twin300Step120/']; % vel during deceleration

Inpath{3,1} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/SpatialTuning/Pyr_SigSpatialInfo_Track_Abs_KetImpactRun2/Twin300Step120/']; % tuning feautre 
Inpath{3,2} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/RunDecoding/BayeDcode/BayesianRun_InSesErrvsShuffle_KetImpactRun2/Twin300Step120/']; % decoding results
Inpath{3,3} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_Run/KTMUActNew/FiringRate/CellRate_EntireSes_KetImpactRun2/Twin300Step120/']; % rate in sessions

Inpath{4,1} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/KTMUActNew/Run/RunLFP_PsdXses_GoodSeg_KetImpactRun2/Twin300Step120/Vbins0    3   10   30  100/']; % load theta frequency by vel
Inpath{4,2} =  [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/KTMUActNew/Run/GammaCoh_LRXses_byDir_KetImpactRun2/Twin300Step120/Velcri10Gband50  100/']; % L/R hemisphere gamma coherence

Inpath{5,1} =  '/media/yuchen/data14/CAResults/cartoons/detour/'; % image of the rat cartoon
Outpath =   '/media/yuchen/data14/CAResults/Figures/Ketamine/';
figurename = 'Setup_BehLFPImpact_v5';
mkdir(Outpath)
%% make subplots
% plot rest run rest protocal
PlotProtocal(pos11,Inpath{5,1});
drawnow

% plot example showing run2 split
PlotRun2Split(pos12,pos13,Inpath{1,1},Inpath{1,2},Inpath{1,3},Inpath{1,4},rnum);
drawnow

% plot power spectrum across session by vel
PlotPSDbyVel(pos31,Inpath{4,1});

% plot run vel and angle distribution
PlotRunVelFXSes(pos21,pos22,Inpath{2,1},Inpath{2,2},Inpath{2,3},Inpath{2,5},Inpath{2,6});

% plot gamma coherence
PlotGammaCoh(pos32,Inpath{4,2});

% % plot theta frequency across session by vel
% PlotThetaFrqbyVel(pos32,Inpath{3,3},Inpath{3,4});

% % plot decoding result
% PlotXsesRunDecoding(pos23,Inpath{2,5},Inpath{2,6});

%% save the result
File_Path = strcat(Outpath,figurename,'.fig');
saveas(gcf, File_Path);
File_Path = strcat(Outpath,figurename,'.png');
print(gcf,File_Path,'-dpng','-r400'); 
% export_fig([Outpath,figurename,'.pdf'])
print('-painters',gcf,'-dpdf',[Outpath,figurename,'_AI.pdf']) % this is for illustrator
close all
%% subpanel plot function

function PlotProtocal(pos,inpath)
global txtfonts 
plotcolors = {[0.12 0.12 0.12],[0.12 0.12 0.12],[0.12 0.12 0.12]};
[ratimg,~,ratalpha] = imread([inpath,'rat1.png']);
% ratimg = flip(ratimg,1);
% ratalpha = flip(ratalpha,1);
ratimgt = nan(size(ratimg,2),size(ratimg,1),size(ratimg,3));
for ic = 1:3
    ratimgt(:,:,ic) = ratimg(:,:,1)';
end
ratimg = ratimgt;
ratimg(ratimg > 0.95) = 0.95;
ratimg = 1-0.8*(1-ratimg);
ratalpha = ratalpha';
ax = subplot('position',pos);
hold(ax,'on')

axis off
% plot juvenile group, day1
text(ax,-0.4,3.5,{'Day1'},'FontSize',txtfonts(3)+1,'HorizontalAlignment','center')
yoff = 2;
xgap = 3.2;
text(ax,-0.4,yoff,{'Juvenile', 'P30-31'},'FontSize',txtfonts(3)+1,'HorizontalAlignment','center')



for islp = 1:3
    rectangle(ax,'Position',[0.5+(islp-1)*xgap yoff-0.8 1.3 1.6],...
        'Curvature',0.04,'FaceColor',[0.9 0.9 0.9],'edgecolor',[0.15 0.15 0.15])
    text(ax,1+(islp-1)*xgap+0.15,yoff,['Sleep',num2str(islp)],'FontSize',txtfonts(4)+1,'HorizontalAlignment','center')
end

for irun = 1:3
   plot(ax,[2.1 3.5] +(irun-1)*xgap, [yoff yoff],'-','linewidth',4,'color',[0.12 0.12 0.12]) 
   if irun == 1
       image(ratimg,'alphadata',im2double(ratalpha),'XData',[2 3.1],'YData',[-0.62 0.62]+yoff)
   end
end

% plot injection time

plot(ax,[2,2],[0.3 3.6],'--','color',plotcolors{1}) 
plot(ax,2,3.8,'v','markersize',7,'markerfacecolor',plotcolors{1},'markeredgecolor',plotcolors{1}) 
plot(ax,[2,2]+xgap,[0.3 3.6],'--','color',plotcolors{2}) 
plot(ax,2+xgap,3.8,'v','markersize',7,'markerfacecolor',plotcolors{2},'markeredgecolor',plotcolors{2}) 

text(ax,2,5,{'IP: Saline','Injection'},'FontSize',txtfonts(2)-1.5,'HorizontalAlignment','center')
text(ax,2+xgap,5,{'IP: Ketamine','Injection (10mg/kg)'},'FontSize',txtfonts(2)-1.5,'HorizontalAlignment','center')

text(ax,2.9,-0.25,{'Run1'},'FontSize',txtfonts(2)+0.5,'HorizontalAlignment','center','color',0.9*plotcolors{1}) 
text(ax,2.9+xgap,-0.25,{'Run2'},'FontSize',txtfonts(2)+0.5,'HorizontalAlignment','center','color',0.9*plotcolors{2}) 
text(ax,2.9+xgap*2,-0.25,{'Run3'},'FontSize',txtfonts(2)+0.5,'HorizontalAlignment','center','color',0.9*plotcolors{3}) 

% plot juvenile group, day2
text(ax,-0.4,-1,{'Day2'},'FontSize',txtfonts(3)+1,'HorizontalAlignment','center')
yoff = -2.5;
xgap = 3.2;
text(ax,-0.4,yoff,{'Juvenile', 'P31-32'},'FontSize',txtfonts(3)+1,'HorizontalAlignment','center')


for islp = 1:2
    rectangle(ax,'Position',[0.5+(islp-1)*xgap yoff-0.8 1.3 1.6],...
        'Curvature',0.04,'FaceColor',[0.9 0.9 0.9],'edgecolor',[0.15 0.15 0.15])
    text(ax,1+(islp-1)*xgap+0.15,yoff,['Sleep',num2str(islp)],'FontSize',txtfonts(4)+1,'HorizontalAlignment','center')
end

for irun = 1:2
   plot(ax,[2.1 3.5] +(irun-1)*xgap, [yoff yoff],'-','linewidth',4,'color',[0.12 0.12 0.12]) 
   if irun == 1
       image(ratimg,'alphadata',im2double(ratalpha),'XData',[2 3.1],'YData',[-0.62 0.62]+yoff)
   end
end

xlim(ax,[-1 11.2])
ylim(ax,[-3.8 3.8])
end

function PlotRun2Split(pos,pos2,InPath,InPath2,InPath3,InPath4,rnum)
global txtfonts sescolors sesstr

load ([InPath,'ses.mat'],'ses')
load ([InPath2,'M12ClFThrSpeedVarB_note.mat'],'M1','M1_note_param','PlFields','M2','M2_note_param','PlFmesh','ClF')
load ([InPath3,'Maze.mat'],'mazel','times')
load ([InPath4,'tbinvel.mat'],'xsestbinvel','xsesksp','xsesvelbinc','run13vel','run1vel','run3vel','hyptime')

nsess = length(ses);
ktmtime = CA_KTM_KTMInjectionTime(rnum);
for is = 1:nsess
    ax = subplot('position',pos{is});
    hold(ax,'on')
    smazel = mazel{is};
    mazelt = smazel(:,6);
    mazelpos = smazel(:,4);

    plotx = mazelt - min(mazelt);
    
    sidx = KetSesMap(is);
    
    plot(ax,plotx,mazelpos,'-','linewidth',0.5,'color',sescolors{sidx})
    ylim(ax,[0 160])
    yticks(ax,[0 150])
    xlim(ax,[0 1600])
    
    % highlight sig impact time in session 2
    if is == 2
        sigidx = idxinrange(mazelt,hyptime);
        plot(ax,plotx(sigidx),mazelpos(sigidx),'-','color',sescolors{2},'linewidth',0.5)
        plot(ax,plotx(~sigidx),mazelpos(~sigidx),'-','color',sescolors{3},'linewidth',0.5)
        xlim(ax,[0 1600])
    end
    if is == 3
        xlabel(ax,'Session time (s)')
    else
        xticklabels(ax,[])
    end
    if is == 2
        ylabel(ax,{'Position (cm)',['Run',num2str(is)]})
    else
        ylabel(ax,['Run',num2str(is)])
    end
    set(ax,'fontsize',txtfonts(4)+0.5)
end

% plot CDFs
ax = subplot('position',pos{nsess+1});
hold(ax,'on')
% plot RUN13 CDF
for irun13 = 1:length(run1vel)
    dtmp = run1vel{irun13};
    dtmp(isnan(dtmp)) = [];
    [cdfy,cdfx] = ecdf(dtmp);
    plot(ax,cdfx,cdfy,'-','linewidth',0.5,'color',sescolors{1});
end
for irun13 = 1:length(run3vel)
    dtmp = run3vel{irun13};
    dtmp(isnan(dtmp)) = [];
    [cdfy,cdfx] = ecdf(dtmp);
    plot(ax,cdfx,cdfy,'-','linewidth',0.5,'color',sescolors{4});
end
% plot vel CDF from run2 in different time bins
psmallidx = 1;
tbinvel = xsestbinvel{2};
ksp = xsesksp{2};
relkbinc = xsesvelbinc{2};
sigidx = idxinrange(relkbinc+ktmtime,hyptime);
for ibin = 1:length(tbinvel)
    dtmp = tbinvel{ibin};
    dtmp(isnan(dtmp)) = [];
    [cdfy,cdfx] = ecdf(dtmp);
%     if ksp(ibin) > 0.05
%         psmallidx = 0;
%     end
    if sigidx(ibin)
        plot(ax,cdfx,cdfy,'-','color',sescolors{2},'linewidth',0.5);
%         plot(ax2,relkbinc(ibin),ksp(ibin),'o','markerfacecolor',sescolors{2},...
%             'markeredgecolor',sescolors{2},'linewidth',0.5,'markersize',3.5);
    else
        plot(ax,cdfx,cdfy,'-','color',sescolors{3},'linewidth',0.5);
%         plot(ax2,relkbinc(ibin),ksp(ibin),'o','markerfacecolor',sescolors{3},...
%             'markeredgecolor',sescolors{3},'linewidth',0.5,'markersize',3.5);
    end
end
xlabel(ax,'Velocity (cm/s)')
ylabel(ax,'CDF')
set(ax,'fontsize',txtfonts(4)+0.5)

hd(1) = plot(ax,nan,nan,'-','color',sescolors{1},'linewidth',0.5);
hd(2) = plot(ax,nan,nan,'-','color',sescolors{2},'linewidth',0.5);
hd(3) = plot(ax,nan,nan,'-','color',sescolors{3},'linewidth',0.5);
hd(4) = plot(ax,nan,nan,'-','color',sescolors{4},'linewidth',0.5);
[lgd,icon] = legend(ax,hd,sesstr,'fontsize',txtfonts(4),...
    'Position',[0.630 0.844 0.095 0.044]);
lgd.FontSize = txtfonts(4)-0.5;
icon(5).XData(1) = icon(5).XData(1) + 0.45;
icon(7).XData(1) = icon(7).XData(1) + 0.45;
icon(9).XData(1) = icon(9).XData(1) + 0.45;
icon(11).XData(1) = icon(11).XData(1) + 0.45;
legend boxoff

% plot value
sesmap = [1,0,4];
for is = 1:nsess

    smazel = mazel{is};
    mazelt = smazel(:,6);
    ksp = xsesksp{is};
    relkbinc = xsesvelbinc{is};
    ax = subplot('position',pos2{is});
    hold(ax,'on')
    if is == 2
        sigidx = idxinrange(relkbinc+ktmtime,hyptime);
        plot(ax,relkbinc(sigidx)+ktmtime-min(mazelt),ksp(sigidx),'o','markerfacecolor',sescolors{2},...
            'markeredgecolor',sescolors{2},'linewidth',0.5,'markersize',2.5);
        
        plot(ax,relkbinc(~sigidx)+ktmtime-min(mazelt),ksp(~sigidx),'o','markerfacecolor',sescolors{3},...
            'markeredgecolor',sescolors{3},'linewidth',0.5,'markersize',2.5);
    else
        plot(ax,relkbinc+ktmtime-min(mazelt),ksp,'o','markerfacecolor',sescolors{sesmap(is)},...
            'markeredgecolor',sescolors{sesmap(is)},'linewidth',0.5,'markersize',2.5);
    end
    xl = xlim(ax);
    yticks(ax,[0.05 0.5 1])
    plot(ax,xl,[0.05 0.05],'--k','linewidth',0.5);
    xlim(ax,[0 1600])
    set(ax,'fontsize',txtfonts(4)+0.5)
    ylim(ax,[-0.05 1.05])
    if is == 3
       xlabel(ax,'Session time (s)') 
    else
        xticklabels(ax,[])
    end
%     if is == 2
%         ylabel(ax,'Max pval vs. Run1')
%     end

    
    if is == 2
        ylabel(ax,{'Max pval vs. Run1',['Run',num2str(is)]})
    else
        ylabel(ax,['Run',num2str(is)])
    end
    
end

end

function PlotRunVelFXSes(pos,pos2,Inpath,Inpath2,Inpath3,Inpath4,Inpath5)
global txtfonts sescolors sesstr
% plot distance per min
load([Inpath,'dispermall.mat'],'dispermall')
nsess = size(dispermall,2);
ax = subplot('position',pos{1});
hold(ax,'on')
sesidx = 1:nsess;
plotsesidx = KetSesMap(sesidx);
hold(ax,'on')
plotd = dispermall(:,plotsesidx);

[g,p] = AllGroupTtest(plotd);
pidx = [1,4,5];
pidx = 1:6;
MyPairSEMBarPlot(plotd,ax,1:nsess,1,'colors',sescolors,'linewidth',1.5,'CapSize',6)
myoverlapsigstar_v3(ax,g(pidx),p(pidx),'starsize',6,'blockdis',3.5)
xlabel(ax,'Sessions')
ylabel(ax,'Distance per minute (cm)')
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
xlim(ax,[0.35 nsess + 0.65])
ylim(ax,[0 1800])
set(ax,'fontsize',txtfonts(4))

% plot angle 
vcri = 20;
tlencri = 1;
parastr = ['Vcri',num2str(vcri),'Tlen',num2str(tlencri)];
load([Inpath2,'vel',parastr,'.mat'],'vel','avecost')
ax = subplot('position',pos{3});
nsess = size(avecost,2);
sesidx = 1:nsess;
plotsesidx = KetSesMap(sesidx);
hold(ax,'on')
plotd = avecost(:,plotsesidx);

[g,p] = AllGroupSignRank(plotd);
pidx = [1,4,5];
pidx = 1:6;
MyPairSEMBarPlot(plotd,ax,1:nsess,1,'colors',sescolors,'linewidth',1.5,'CapSize',6)
ylim(ax,[0 20])
myoverlapsigstar_v3(ax,g(pidx),p(pidx),'starsize',6,'starty',14)
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
xlim(ax,[0.35 nsess + 0.65])
xlabel(ax,'Sessions')
ylabel(ax,'Angle <Track,Velocity> (deg)')
set(ax,'fontsize',txtfonts(4))
ylim(ax,[0 20])

% plot max vel
vcri = 10; % absolute velocity larger than...
tlencri = 2; % continuous run over period
rundis = 50; % continuous run over distance
parastr = ['Vcri',num2str(vcri),'Tlen',num2str(tlencri),'Dis',num2str(rundis)];
load([Inpath3,'lpvelmax',parastr,'.mat'],'avevel')
ax = subplot('position',pos{2});
nsess = size(avevel,2);
sesidx = 1:nsess;
plotsesidx = KetSesMap(sesidx);
hold(ax,'on')
plotd = avevel(:,plotsesidx);

[g,p] = AllGroupSignRank(plotd);
pidx = [1,4,5];
pidx = 1:6;
MyPairSEMBarPlot(plotd,ax,1:nsess,1,'colors',sescolors,'linewidth',1.5,'CapSize',6)
ylim(ax,[0 140])
myoverlapsigstar_v3(ax,g(pidx),p(pidx),'starsize',6,'blockdis',3.8)
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
xlim(ax,[0.35 nsess + 0.65])
xlabel(ax,'Sessions')
ylabel(ax,'Max velocity (cm/s)')
set(ax,'fontsize',txtfonts(4))
ylim(ax,[0 150])

% plot vel vs distance to end
vcri = 2;%0;
tlen = 1;%0.1;
parastr = ['Vcri',num2str(vcri),'Tlen',num2str(tlen)];
load([Inpath4,'vel2tend',parastr,'.mat'],'vel2tend','disbinc','ratmean')
nsess = 4;
ax = subplot('position',pos2{1});
set(ax,'fontsize',txtfonts(4))
pmoveleft(0.01,ax)
hold(ax,'on')
meanv = squeeze(nanmean(ratmean,3));
semv = squeeze(nansembydim(ratmean,3));
pvals = nan(size(disbinc));

for is = 1:nsess
    sesnow = KetSesMap(is);
    meannow = meanv(sesnow,:);
    semnow = semv(sesnow,:);
    PlotMeanSemShade(ax,disbinc,meannow,semnow,sescolors{is},0.2,'linewidth',1)
    hd(is) = plot(ax,nan,nan,'-','color',sescolors{is},'linewidth',1);
end
for ibin = 1:length(disbinc)
    d1 = squeeze(ratmean(1,ibin,:));
    d2 = squeeze(ratmean(2,ibin,:));
    d4 = squeeze(ratmean(3,ibin,:)); % 3 corresponds to run 3 in ratmean
    [~,p12] = ttest(d1-d2,0);
    [~,p24] = ttest(d4-d2,0);
    pvals(ibin) = max([p12,p24]);
    plotxsigstar(ax,pvals(ibin),disbinc(ibin),50,2.2)
end

myarrowplot(ax,34,58,-17,0,'linewidth',1,'markersize',3)
text(ax,25,66,{'Run toward corner','Deceleration'},'FontSize',txtfonts(4)+0.5,'HorizontalAlignment','center')

xlim(ax,[-0.5 50])
ylim(ax,[0 65])
[lgd,icon] = legend(hd,sesstr,'FontSize',txtfonts(4));
icon(5).XData(1) = icon(5).XData(1) + 0.45;
icon(7).XData(1) = icon(7).XData(1) + 0.45;
icon(9).XData(1) = icon(9).XData(1) + 0.45;
icon(11).XData(1) = icon(11).XData(1) + 0.45;
legend boxoff
set(lgd,'Position',[0.712 0.56 0.086 0.05]);
xlabel(ax,'  Distance to corner (cm)')
ylabel(ax,'Velocity (cm/s)')


% plot vel vs distance to end
vcri = 2;%0;
tlen = 1;%0.1;
parastr = ['Vcri',num2str(vcri),'Tlen',num2str(tlen)];
load([Inpath5,'vel2tend',parastr,'.mat'],'vel2tend','disbinc','ratmean')
nsess = 4;
ax = subplot('position',pos2{2});
set(ax,'fontsize',txtfonts(4))
pmoveleft(0.01,ax)
hold(ax,'on')
meanv = squeeze(nanmean(ratmean,3));
semv = squeeze(nansembydim(ratmean,3));
pvals = nan(size(disbinc));

for is = 1:nsess
    sesnow = KetSesMap(is);
    meannow = meanv(sesnow,:);
    semnow = semv(sesnow,:);
    PlotMeanSemShade(ax,disbinc,meannow,semnow,sescolors{is},0.2,'linewidth',1)
    hd(is) = plot(ax,nan,nan,'-','color',sescolors{is},'linewidth',1);
end
for ibin = 1:length(disbinc)
    d1 = squeeze(ratmean(1,ibin,:));
    d2 = squeeze(ratmean(2,ibin,:));
    d4 = squeeze(ratmean(3,ibin,:)); % 3 corresponds to run 3 in ratmean
    [~,p12] = ttest(d1-d2,0);
    [~,p24] = ttest(d4-d2,0);
    pvals(ibin) = max([p12,p24]);
    plotxsigstar(ax,pvals(ibin),disbinc(ibin),50,2.2)
end

myarrowplot(ax,17,58,17,0,'linewidth',1,'markersize',3)
text(ax,25,66,{'Run away from corner','Acceleration'},'FontSize',txtfonts(4)+0.5,'HorizontalAlignment','center')

xlim(ax,[-0.5 50])
ylim(ax,[0 65])
[lgd,icon] = legend(hd,sesstr,'FontSize',txtfonts(4));
icon(5).XData(1) = icon(5).XData(1) + 0.45;
icon(7).XData(1) = icon(7).XData(1) + 0.45;
icon(9).XData(1) = icon(9).XData(1) + 0.45;
icon(11).XData(1) = icon(11).XData(1) + 0.45;
legend boxoff
set(lgd,'Position',[0.884 0.56 0.086 0.05]);
xlabel(ax,'  Distance to corner (cm)')
ylabel(ax,'Velocity (cm/s)')


end

function PlotPSDbyVel(pos,Inpath)
global txtfonts sescolors sesstr
% plot juvenile map
cmap=getNCLColmap('MPL_BuPu.rgb',100);
vbin = [0,3,10,30,100];
load([Inpath,'psd.mat'],'psdses','ff')
nsess = 4;
allses = 1:nsess;
sesplotorder = [setdiff(allses,2),2];

ffplot = logspace(0,2.4,50);
for is = 1:size(psdses,1)
    for iv = 1:size(psdses,2)
        dnow = psdses{is,iv};
        ditp = nan(size(dnow,1),length(ffplot));
        for isample = 1:size(dnow,1)
            ditp(isample,:) = interp1(ff,dnow(isample,:),ffplot);
        end
        psdses{is,iv} = ditp;
    end
end


for iv = 1:size(psdses,2)
    ax = subplot('position',pos{iv});
    hold(ax,'on')
    for sesloop = 1:length(sesplotorder)
        ises = sesplotorder(sesloop);
        sidx = KetSesMap(ises);
        specnow = psdses{sidx,iv};
        if isempty(specnow)
           continue 
        end
        meanspec = nanmean(specnow,1);
        semspec = nansembydim(specnow,1);
        specrange(:,1) = meanspec - semspec;
        specrange(:,2) = meanspec + semspec;

        xx = [berow(ffplot),fliplr(berow(ffplot))];
        yy = [berow(specrange(:,1)),fliplr(berow(specrange(:,2)))];
        fill(xx,yy,1-0.8*(1-sescolors{ises}),'LineStyle','none','FaceAlpha',0.2)
        hd(ises) = plot(ffplot,meanspec,'linewidth',1,'color',sescolors{ises});
        set(ax,'xscale','log')
        set(ax,'yscale','log')
    end
    
    otherses = [1,3,4];
    pvals = nan(size(ffplot));
    
    for ifbin = 1:length(ffplot)
        dktm = psdses{2,iv}(:,ifbin);
        dother = [];
        for ises = otherses
            dother = cat(1,dother,psdses{ises,iv}(:,ifbin));
        end
        pvals(ifbin) = ranksum(dktm,dother);
        plotxsigstar(ax,pvals(ifbin),ffplot(ifbin),1e-8,5e-9)
    end
    
    set(ax,'fontsize',txtfonts(3)-0.5)
    xlim(ax,[2 250])
    ylim(ax,[1e-12 1e-7])
    
    if iv == 1
        ylb = ylabel(ax,'Power (mV^2/Hz)');
%         ylb.Position(2) = 5e-13;
    else
            yticks(ax,[])
    end
    
    if iv == 1
        [lgd,icon] = legend(hd,sesstr,'Location','southwest');
        lgd.FontSize = txtfonts(4);
        lgd.Position(1) = lgd.Position(1) - 0.05;
        lgd.Position(2) = lgd.Position(2) - 0.02;
        icon(5).XData(1) =  icon(5).XData(1) + 0.45;
        icon(7).XData(1) =  icon(7).XData(1) + 0.45;
        icon(9).XData(1) =  icon(9).XData(1) + 0.45;
        icon(11).XData(1) =  icon(11).XData(1) + 0.45;
        legend boxoff
    end
    title(ax,['Velocity:',num2str(vbin(iv)),'-',num2str(vbin(iv+1)),'cm/s'],'FontWeight','normal')
    set(ax,'Xscale','log')
    set(ax,'Yscale','log')
    xticks(ax,[4 8 16 32 64 128])
%     if iv <= 2
%        xticklabels(ax,[]) 
%     end
%     if iv == 3
    xlb = xlabel(ax, 'Frequency (Hz)');
%     xlb.Position(1) = xlb.Position(1) + 480;
%     end
    set(ax,'fontsize',txtfonts(3)-0.5)
    set(ax, 'XGrid', 'on', 'YGrid', 'off','XMinorGrid','off')
end

end

function PlotGammaCoh(pos,Inpath)
global txtfonts sescolors sesstr
load([Inpath,'cohave.mat'],'cohave')

nsess = size(cohave,2);
allses = 1:nsess;
plotses = KetSesMap(allses);
ax = subplot('position',pos);
set(ax,'fontsize',txtfonts(3)-0.5)
plotd = cohave(:,plotses);
[g,p] = AllGroupTtest(plotd);
pidx = [1,4,5];
pidx = [1:6];
MyPairSEMBarPlot(plotd,ax,1:nsess,1,'colors',sescolors,'linewidth',1.5,'CapSize',6)

xlabel(ax,'Session')
ylabel(ax,'Left-Right Gamma coherence')
xticks(ax,1:nsess)
xticklabels(ax,sesstr)

myoverlapsigstar_v3(ax,g(pidx),p(pidx)/2,'starty',0.65,'ygap',0.02,'starsize',6,'blockdis',3)
ylim(ax,[0 0.9])
end



function sigstr = sigstrfromp(p)

if p >= 0.05
    sigstr = 'n.s.';
end
if p < 0.05 && p >= 0.01
    sigstr = '*';
end
if p < 0.01 && p >= 0.001
    sigstr = '**';
end
if p < 0.001
    sigstr = '***';
end

end

function tmpcdfwithsigbar(ax,dplot)
global txtfonts sescolors
ydis = 0.12;
MyCDFsCompare_v3(dplot,ax,sescolors,{'Ctrl','E-Ket','L-Ket'},txtfonts(4),'unequal',[],ydis)

ldgstr = {'Ctrl','E-Ket','L-Ket'};
lgdfont = txtfonts(4);
ppairs = [1,3,1;2,2,3];
xl = xlim(ax);
xstart = 0.65;
xtend = 0.80;
ally = nan(1,length(dplot));
for id = 1:length(dplot)
    ally(id) = ydis*id;
end

xlgend = [xl(2)*xstart+xl(1)*(1-xstart),xl(2)*xtend+xl(1)*(1-xtend)];

xoffset = 0.012;
sigstr = cell(1,size(ppairs,2));
sigpval = nan(1,size(ppairs,2));
for ip = 1:size(ppairs,2)
    id = ppairs(1,ip);
    jd = ppairs(2,ip);
    [~,pt] = kstest2(dplot{id}(:),dplot{jd}(:));
%     [~,pt] = kstest2(dplot{id}(:),dplot{jd}(:),'tail','smaller');
    pt = ranksum(dplot{id}(:),dplot{jd}(:));
    sigpval(ip) = pt;
%     [~,pt] = ttest2(dplot{id}(:),dplot{jd}(:),'tail','right');
    sigstr{ip} = sigstrfromp(pt);
end
% plot sigbar
if strcmp(sigstr{1},sigstr{2})
    xdiff = 0.013;
    for ip = 1:size(ppairs,2)
        id = ppairs(1,ip);
        jd = ppairs(2,ip);
        xoffset = xoffset + xdiff*(ip-1);
        plot(ax,[xl(2)*(xstart-xoffset)+xl(1)*(1+xoffset-xstart) xl(2)*(xstart-xoffset)+xl(1)*(1+xoffset-xstart)],[ally(id) ally(jd)],'-k')
    end
    text(ax,xl(2)*(xstart-xoffset-0.01)+xl(1)*(1+xoffset+0.01-xstart),mean([ally(id),ally(jd)]),sigstr{1},'HorizontalAlignment','right','FontSize',lgdfont)
else
    xdiff = 0.013;
    for ip = 1:size(ppairs,2)
        id = ppairs(1,ip);
        jd = ppairs(2,ip);
        xoffset = xoffset + xdiff*(ip-1);
        plot(ax,[xl(2)*(xstart-xoffset)+xl(1)*(1+xoffset-xstart) xl(2)*(xstart-xoffset)+xl(1)*(1+xoffset-xstart)],[ally(id) ally(jd)],'-k')
        text(ax,xl(2)*(xstart-xoffset-0.01)+xl(1)*(1+xoffset+0.01-xstart),mean([ally(id),ally(jd)]),sigstr{ip},'HorizontalAlignment','right','FontSize',lgdfont)
    end
end


end


function plotxsigstar(ax,p,x,y,gap)
fs = 7;
hold(ax,'on')
% if p < 0.001
%     plot(ax,x,y+gap,'k*','markersize',size,'linewidth',lwidth)
%     plot(ax,x,y,'k*','markersize',size,'linewidth',lwidth)
%     plot(ax,x,y-gap,'k*','markersize',size,'linewidth',lwidth)
% elseif p < 0.01
%     plot(ax,x,y+gap,'k*','markersize',size,'linewidth',lwidth)
%     plot(ax,x,y,'k*','markersize',size,'linewidth',lwidth)
% elseif p < 0.05
%     plot(ax,x,y+gap,'k*','markersize',size,'linewidth',lwidth)
% end
% 

hold(ax,'on')
if p < 0.001
    text(ax,x,y+gap,'*','fontsize', fs,'HorizontalAlignment','Center');
    text(ax,x,y,'*','fontsize', fs,'HorizontalAlignment','Center');
    text(ax,x,y-gap,'*','fontsize', fs,'HorizontalAlignment','Center');
elseif p < 0.01
    text(ax,x,y+gap,'*','fontsize', fs,'HorizontalAlignment','Center');
    text(ax,x,y,'*','fontsize', fs,'HorizontalAlignment','Center');
elseif p < 0.05
    text(ax,x,y+gap,'*','fontsize', fs,'HorizontalAlignment','Center');
end



% if p < 0.001
%     plot(ax,x,y+gap,'ks','markersize',size,'linewidth',lwidth,'markerfacecolor',[0 0 0])
%     plot(ax,x,y,'ks','markersize',size,'linewidth',lwidth,'markerfacecolor',[0 0 0])
%     plot(ax,x,y-gap,'ks','markersize',size,'linewidth',lwidth,'markerfacecolor',[0 0 0])
% elseif p < 0.01
%     plot(ax,x,y+gap,'ks','markersize',size,'linewidth',lwidth,'markerfacecolor',[0 0 0])
%     plot(ax,x,y,'ks','markersize',size,'linewidth',lwidth,'markerfacecolor',[0 0 0])
% elseif p < 0.05
%     plot(ax,x,y+gap,'ks','markersize',size,'linewidth',lwidth,'markerfacecolor',[0 0 0])
% end

end