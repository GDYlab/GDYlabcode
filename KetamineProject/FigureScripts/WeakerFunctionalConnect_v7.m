clear;clc;close all
dbstop if error

%% define subpanel positions
testpos = 1;

% cofiring across sessions
[~,pos11] = tight_subplot(2,2,[0.038 0.018],[0.575 0.065],[0.055 0.695],0);
pos12 = [0.41 0.65 0.20 0.25];

% cartoon of transfer entropy
pos21{1} = [0.66 0.82 0.32 0.12];
pos21{2} = [0.68 0.61 0.09 0.13];
pos21{3} = [0.86 0.61 0.09 0.13];

% MUinfo examples
[~,pos31] = tight_subplot(2,2,[0.035 0.015],[0.07 0.545],[0.045 0.695],0);
% MUinfo across sessions
[~,pos32] = tight_subplot(2,2,[0.035 0.022],[0.07 0.545],[0.41 0.30],0); 

% dimensionality of spikes
pos41 = [0.76 0.07 0.21 0.185];
[~,pos42] = tight_subplot(1,2,[0.035 0.045],[0.325 0.545],[0.76 0.03],0); 

% test position
if testpos
f = figure_letter(0.5);
testsubpos(pos11)
testsubpos({pos12})
testsubpos(pos21)

testsubpos(pos31)
testsubpos(pos32)
testsubpos(pos41)
testsubpos(pos42)
% testsubpos(post44)
end
close all
f = figure_letter(0.5);
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
    panellbpos = {[0 1],[0.35 1],[0.62 1],[0 0.535],[0.36 0.535],[0.71 0.535]};
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

annotation(f,'textbox',...
    [0.07 0.945 0.2 0.04],...
    'String',{'Example of cofiring'},...
    'LineStyle','none','Fontsize',txtfonts(2),'HorizontalAlignment','center',...
    'FitBoxToText','off');


annotation(f,'textbox',...
    [0.03 0.475 0.3 0.04],...
    'String',{'Example of functional connectivity'},...
    'LineStyle','none','Fontsize',txtfonts(2),'HorizontalAlignment','center',...
    'FitBoxToText','off');


annotation(f,'textbox',...
    [0.405 0.475 0.3 0.04],...
    'String',{'Functional connectivity across sessions'},...
    'LineStyle','none','Fontsize',txtfonts(2),'HorizontalAlignment','center',...
    'FitBoxToText','off');

%% define data paths
Inpath = cell(10,6); % 10 subpanels, each could have up to 6 data source
% rnum = "Detour_R3";
% setnum = 2;
% 
% rnum = "URat1";
% setnum = 1;
% [Datasets] = getCADatasets(rnum,setnum); %
% Datasets = Datasets{1};

Inpath{1,1} =  '/media/yuchen/data14/CAResults/Solver/KTM_Run/URat2/kd1_2016-03-16_11-03-24/Cofiring/PyrInt_Cofiring_DirXses_KetImpactRun2/Twin300Step120/'; % cofiring example
Inpath{1,2} =  '/media/yuchen/data14/CAResults/Sets_Ave/KTM_Run/KTMUActNew/Cofiring/Cofiring/PyrInt_Cofiring_Xses_vsREMSlp_KetImpactRun2/Twin300Step120/BrainState_EMG_Vcri2_SlpLen10/'; % cofiring strength

Inpath{2,1} =  '/media/yuchen/data14/CAResults/Solver/KTM_Run/XRat1/kd1_2024-05-23_14-28-50/TEntropy/XCells/PyrInt_HighVelRunXsesMI_DirSeg_Point_v2_KetImpactRun2/Twin300Step120/Tbin0.2spkcri5/'; % MUinfo example
Inpath{2,2} =  '/media/yuchen/data14/CAResults/Sets_Ave/KTM_XRunSleep/KTMUActNew/TEntropy/XCells/PyrInt_RunREM_MIMag_PreCellVCorr_KetImpactRun2/Twin300Step120/RunTbin0.4spkcri5EMG_SlpLen10Tbin0.4Spkcri5/'; % MUinfo xses
Inpath{2,3} =  '/media/yuchen/data14/CAResults/Preprocess/KTM_Run/XRat1/kd1_2024-05-23_14-28-50/PyrIntSpk_KetImpactRun2/Twin300Step120/'; %path to pyr int number


Inpath{3,1} =  '/media/yuchen/data14/CAResults/Sets_Ave/KTM_XRunSleep/KTMUActNew/Embedding/XRunREMSlp_Dimfrom2NN_KetImpactRun2/Twin300Step120/RunTbin0.4Vcri10EMG_SlpLen10tbin0.4/'; % dimension of spikes 
Inpath{3,2} =  '/media/yuchen/data14/CAResults/Solver/KTM_Run/URat1/kd1_2016-01-30_12-23-12/Embedding/TDXSes_Dimfrom2NN_KetImpactRun2/Twin300Step120/tbin0.4vcri10/';

Outpath =   '/media/yuchen/data14/CAResults/Figures/Ketamine/';
figurename = 'WeakerFunctionalConnect_v7';
mkdir(Outpath)
%% make subplots
% plot transfer entropy cartoon
MUinfoCartoon(pos21);

% plot MUinfo example
PlotMUinfoExample(pos31,Inpath{2,1},Inpath{2,3});

% plot cofiring strength across sessions
PlotCofXSes(pos12,Inpath{1,2})

% plot dimension of data
PlotDfrom2NNXSes(pos41,Inpath{3,1})

% plot MUinfo across sessions
PlotMUinfoXses(pos32,Inpath{2,2});

% plot mu distribution example 
PlotMuExampleXSes(pos42,Inpath{3,2})

% plot cofiring example
PlotCofExample(pos11,Inpath{1,1})


%% save the result
File_Path = strcat(Outpath,figurename,'.fig');
saveas(gcf, File_Path);
File_Path = strcat(Outpath,figurename,'.png');
print(gcf,File_Path,'-dpng','-r300'); 
% export_fig([Outpath,figurename,'.pdf'])
print('-painters',gcf,'-dpdf',[Outpath,figurename,'_AI.pdf']) % this is for illustrator
close all
%% subpanel plot function

function MUinfoCartoon(pos)
global txtfonts sescolors cmap
% cmap = getNCLColmap('MPL_PuRd.rgb',100); % this is for detour
rng(83)

nc = 5;
nt = 20;
aid = 3;
bid = 2;

mrate = rand(1,nc);
allbspk = zeros(nc,nt);
mrate(aid) = 0.95;
for ic = 1:nc
    allbspk(ic,:) = poissrnd(mrate(ic),1,nt);
end

% use cell A to predict cell B

for tbin = 1:nt-1
    spknow =  allbspk(aid,tbin);
    allbspk(bid,tbin+1) = determineBfromA(spknow);
end

% plot all spikes
ax = subplot('position',pos{1});
hold(ax,'on')
% plot all spikes
for tbin = 1:nt
    for ic = 1:nc
        spknow = allbspk(ic,tbin);
        if spknow > 0
           xoff = 0.1+rand(1,spknow)/2;
           xoff = cumsum(xoff);
           if xoff(end) > 0.9
               xoff = xoff./xoff(end)*0.9;
           end
        else
            continue
        end
        switch ic
            case aid
                crgb = [0.8 0.35 0.45];
                markersize = 1.5;
                len = 0.15;
            case bid
                crgb = [0.45 0.35 0.8];
                markersize = 1.5;
                len = 0.15;
            otherwise
                crgb = [0.6 0.6 0.6];
                markersize = 0.5;
                len = 0.1;
        end
%         plot(ax,xoff+tbin-1,ic*ones(1,length(spknow)),'s',...
%             'markerfacecolor',crgb,'markeredgecolor',crgb,'markersize',markersize)
        xall = xoff+tbin-1;
        yall = ic*ones(1,length(xall));
        for ispk = 1:length(xall)
            plot(ax,[xall(ispk),xall(ispk)],[-len len]+yall(ispk),'-',...
                'color',crgb,'linewidth',markersize)
        end
    end
end

% plot time bins
for tbin = 1:nt+1
    plot(ax,[tbin-1 tbin-1],[0.5 nc+0.5],'-','color',[0.75 0.75 0.75]);    
end
axis tight
xtk = 0.5:4:19.5;
xticks(ax,xtk)
xtklbl = cell(1,length(xtk));
for itk = 1:length(xtk)
    xtklbl{itk} = num2str(ceil(xtk(itk)));
end
xticklabels(ax,xtklbl)
xlabel(ax,'Time bins')
% ylabel(ax,'Cells')
yticks(ax,[bid,aid])
yticklabels({'B','A'})
box on
set(ax,'fontsize',txtfonts(4))
ax.TickLength = [0 0];
title(ax,'Cartoon of functional connectivity','FontWeight','normal','FontSize',txtfonts(2),'HorizontalAlignment','center')

% plot joint distributon
ax2 = subplot('position',pos{2});
hold(ax2,'on')
% linex = 1.4;
% plot(ax2,[linex linex],[0 5],'-','linewidth',4,'color',[0.6 0.6 0.6 0.6])
aspike = normrnd(1.8,1.2,[1e4,1]);
aspike = round(aspike);
aspike(aspike<0) = 0;

bspike = -1/3*(aspike).^3+1/2*(aspike).^2+1.8333*aspike;
offsetb = normrnd(0,0.46,size(aspike));
bspike = bspike + offsetb;
bspike = round(bspike);
bspike(bspike<0) = 0;
% [n] = hist2d([aspike,bspike],1:5,1:5,[0 6 0 6]);
n = hist2d_v2(aspike(1:10:end),bspike(1:10:end),-0.5:1:5.5,-0.5:1:5.5);
% row of n corresponds to bspike, column of n corrspond to aspike
imagesc(ax2,[0 5],[0 5],n)
xlim(ax2,[-0.5 5.5])
ylim(ax2,[-0.5 5.5])
colormap(ax2,cmap)
set(ax2,'Ydir','normal')
xlabel(ax2,'Current spikes of A ')
ylabel(ax2,'Future spikes of B')
% [te,nte] = GaussianTransferEntropy(aspike(2:end),bspike(1:end-1),1);
% [MI,nMI] = MultualInfo([becolumn(aspike),becolumn(bspike)],10);
[MI,~,~,~,~,~,~,~,nMI,~] = PointMutualInfo_multiseg({becolumn(aspike(2:end))},{becolumn(bspike(1:end-1))},1);
% xlim(ax2,[0 4])
% ylim(ax2,[0 5])
title(ax2,['nMI(CurrentA;FutureB) ',num2str(nMI,'%.3f')],'FontWeight','normal')
set(ax2,'fontsize',txtfonts(4))
cb = addcolorbar;
cb.Ticks = [0 228];
cb.TickLabels = {'Low','High'};

% plot future spike distribution on yaxis
% tmppos = pos{2};
% tmppos(1) = tmppos(1) + tmppos(3) + 0.002;
% tmppos(3) = 0.05;
% ax4 = subplot('position',tmppos);
% hold(ax4,'on')
% linex = 1.4;
% [nb,bedg] = histcounts(bspike,-0.1:0.2:5.1,'Normalization','pdf');
% bc = (bedg(1:end-1)+bedg(2:end))/2;
% xx = [berow(nb),zeros(size(berow(nb)))];
% yy = [berow(bc),fliplr(berow(bc))];
% fill(ax4,xx,yy,[0.8 0.8 0.8],'LineWidth',1,'FaceAlpha',0.4)
% 
% aidx = idxinrange(aspike,[1.35 1.45]);
% tmpb = bspike(aidx);
% [nb,bedg] = histcounts(tmpb,-0.1:0.2:5.1,'Normalization','pdf');
% bc = (bedg(1:end-1)+bedg(2:end))/2;
% xx = [berow(nb),zeros(size(berow(nb)))];
% yy = [berow(bc),fliplr(berow(bc))];
% fill(ax4,xx,yy,[0.85 0.3 0.2],'LineStyle','none','FaceAlpha',0.8)
% ylim(ax4,[0 5])
% axis off

% plot joint distributon
ax3 = subplot('position',pos{3});
hold(ax3,'on')
aspike2 = normrnd(1.8,1.2,[1e4,1]);
aspike2 = round(aspike2);
aspike2(aspike2<0) = 0;

bspike2 = normrnd(2.0,1.2,[1e4,1]);
bspike2 = round(bspike2);
bspike2(bspike2<0) = 0;


[n] = hist2d_v2(bspike2(1:10:end),aspike2(1:10:end),-0.5:1:5.5,-0.5:1:5.5);
imagesc(ax3,[0 5],[0 5],n)
xlim(ax3,[-0.5 5.5])
ylim(ax3,[-0.5 5.5])
colormap(ax3,cmap)
set(ax3,'Ydir','normal')
xlabel(ax3,'Current spikes of B ')
ylabel(ax3,'Future spikes of A')
% [te,nte] = GaussianTransferEntropy(aspike(2:end),bspike(1:end-1),1);
% [MI,nMI] = MultualInfo([becolumn(aspike),becolumn(bspike)],10);
[MI,~,~,~,~,~,~,~,nMI,~] = PointMutualInfo_multiseg({becolumn(bspike2(2:end))},{becolumn(aspike2(1:end-1))},1);
% xlim(ax3,[0 4])
% ylim(ax3,[0 5])
title(ax3,['nMI(CurrentB;FutureA) ',num2str(nMI,'%.3f')],'FontWeight','normal')
set(ax3,'fontsize',txtfonts(4))
cb = addcolorbar;
cb.Label.String = 'Counts';
cb.Ticks = [0 110];
cb.TickLabels = {'Low','High'};
cb.Label.Position(1) = cb.Label.Position(1) - 1;

end

function PlotMUinfoExample(pos,Inpath,Inpath2)
global txtfonts sescolors sesstr cmap
tlag = 5;
% cmap=getNCLColmap('MPL_BuPu.rgb',100);
load([Inpath,'behspk_maxt1  2  3  4  5.mat'],'behspk')
load([Inpath2,'RunClF.mat'],'ClF','ClFPyr','ClFInt','pyridx')

npyr = size(ClFPyr,1);
nint = size(ClFInt,1);
celltype = cat(1,ones(nint,1),zeros(npyr,1));
% plot high vel
im = 1;% data
in = 2;% normalized
itl = 1;
nsess = 4;
idir = 1;
for is = 1:nsess
    sid = KetSesMap(is);
    dnow = behspk(sid,idir).xcellnte;
    ax = subplot('position',pos{is});
    hold(ax,'on')
    zz = squeeze(dnow(:,:,itl,im + (in-1)*2));
    zz = log(zz);
    % get nan values
    nanidx = isnan(zz);
    sumnan = sum(nanidx,1);
    goodidx = sumnan ~= size(zz,1);
%     goodidx = sumnan ~= 999;
    %zz = zz(goodidx,goodidx);
    %celltypenow = celltype(goodidx);
    celltypenow = celltype;
    nintnow = sum(celltypenow);
    npyrnow = length(celltypenow) - nintnow;
    
    imagesc([1 size(zz,2)],[1,size(zz,1)],zz')
    caxis([-4 -0])
    colormap(ax,cmap)
    
    title(ax,sesstr{is})
    axis tight
    axis square
    xl = xlim(ax);
    yl = ylim(ax);
    plot(ax,xl,[nintnow+0.5 nintnow+0.5],'-k')
    plot(ax,[nintnow+0.5 nintnow+0.5],yl,'-k')
    set(ax,'layer','top')
    set(ax,'Ydir','Normal')
    if is >= 3
        xticks(ax,[nintnow/2+0.5, nintnow + npyrnow/2])
        xticklabels(ax,{'Ints','Pyrs'})
        xlabel(ax,'Current cell')
    else
        xticks(ax,[])
    end
    set(ax,'Fontsize',txtfonts(4))
  
    if rem(is,2) == 1
    yticks(ax,[nintnow/2+0.5, nintnow + npyrnow/2])
    yticklabels(ax,{'Ints','Pyrs'})
    ylabel(ax,'Future cell')
    else
         yticks(ax,[nintnow/2+0.5, nintnow + npyrnow/2])
         yticklabels(ax,[])
    end

    if is == nsess
        cb = addcolorbar(ax,0.008);
        cb.Label.String = 'Log normalized MI';
%         cb.Position(4) = cb.Position(4) - 0.09;
%         cb.Position(2) = cb.Position(2) + 0.045;
    end
end


end

function PlotMUinfoXses(pos,Inpath)
global txtfonts sescolors sesstr
cmap=getNCLColmap('MPL_BuPu.rgb',100);


allylim = {[-5 1.5],[-5.5 0],[-5.5 0],[-5.2 0.5]};
load([Inpath,'xrunslpmi.mat'],'xrunslpmi','allmeanval');

allses = size(xrunslpmi,1);
nsess = allses - 1;

% plot different cell types, mean value from each rate
typestr = {'int→int','pyr→pyr','pyr→int','int→pyr'};
itl = 1;
in = 2 ;
tmpstr = sesstr;
tmpstr{allses} = 'REMS';

for ictype = 1:4
    ax = subplot('position',pos{ictype});
    hold(ax,'on')
    tmpd = squeeze(allmeanval(ictype,itl,in,:,1:nsess));
    tmpd2 = tmpd;
    sesids = 1:nsess;
    sid = KetSesMap(sesids);
    tmpd2(:,sesids) = tmpd2(:,sid);
    tmpd = tmpd2;
    
    MyPairSEMPlot_v2(tmpd,ax,1:nsess,1,'sigplot',0,'CapSize',5,'linewidth',1.5,'markersize',2.5,'colors',sescolors);
    
    plotdnow = squeeze(allmeanval(ictype,itl,in,:,allses));
    postketidx = [3:3:15];
    plotdnow(postketidx) = [];
    % only use preKTM sleep
    plotdnow(isnan(plotdnow)) = [];
    plotx = allses + (rand(size(plotdnow))-0.5)/10;
    plot(ax,plotx,plotdnow,'o','markersize',2.5,...
        'markerfacecolor',1-0.6*(1-sescolors{5}),'MarkerEdgeColor',1-0.6*(1-sescolors{5}))
    errorbar(ax,allses,nanmean(plotdnow),sem(plotdnow),'-','linewidth',1.5,'CapSize',5,'color',0.5*sescolors{5});
    xticks(ax,1:allses)
    if ictype >= 3
        xticklabels(ax,tmpstr)
    else
      xticklabels(ax,[])  
    end
    xtickangle(ax,45)
    set(ax, 'YGrid', 'on', 'XGrid', 'off')
    
    title(ax,typestr{ictype})
    if rem(ictype,2) == 1
        ylabel(ax,['log normalized MI'])
    end
    xlim(ax,[0.5,allses+0.5])
    ylim(ax,allylim{ictype})
%     [g,p] = AllGroupTtest(tmpd);
    [g,p] = AllGroupSignRank(tmpd);
    % [g,p] = AllGroupRankSum(plotd);
    gplot = {g{1},g{4},g{5}};
    pplot = [p(1),p(4),p(5)];
    yl = ylim(ax);
    myoverlapsigstar_v3(ax,gplot,pplot,'starty',0.25*yl(1)+0.75*yl(2),'starsize',7,'blockdis',2.5)

    set(ax,'fontsize',txtfonts(4))
end


end

function PlotCofExample(pos,Inpath)
global txtfonts sescolors sesstr
load([Inpath,'cofire.mat'],'cofire')
cmap3 = getNCLColmap('BlueWhiteOrangeRed.rgb',256);
idir = 2;
nsess = 4;
for is = 1:nsess
    sid = KetSesMap(is);
    ax = subplot('position',pos{is});
    hold(ax,'on')
    pmoveleft(0.01,ax)
    dnow = cofire(idir,sid).corr;
    dnow = setdiag(dnow,0);
    nc = size(dnow,1);
    dnow(isnan(dnow)) = 0;
    imagesc(ax,[1 nc],[1 nc],dnow)
    colormap(ax,cmap3)
    xlim(ax,[0.5 nc+0.5])
    ylim(ax,[0.5 nc+0.5])
    xticks(ax,[])
    yticks(ax,[])
    set(ax,'layer','top')
    box on
    if is >= 3
        xlabel(ax,'Cell ID')
    end
    if rem(is,2) == 1
        ylabel(ax,'Cell ID')
    end
    caxis(ax,[-0.15 0.15])
    title(ax,sesstr{is})
    set(ax,'fontsize',txtfonts(4))
    axis square
    
    if is == nsess
        cb = addcolorbar;
        cb.Label.String = 'Correlation';
    end
end


end

function PlotCofXSes(pos,Inpath)
global txtfonts sescolors sesstr
load([Inpath,'cfxses.mat'],'cfxses','ctypexses')
nsess = size(ctypexses,2)-1;
allrunses = nsess+1;

tmpstr = sesstr;
tmpstr{allrunses} = 'REMS';


tmpd = nan(10,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    cofnow = becolumn(cfxses{sid});
    tmpd(1:length(cofnow),is) = cofnow;
end

ax = subplot('position',pos);
set(ax,'fontsize',txtfonts(4))
MyPairSEMPlot_v2(tmpd,ax,1:nsess,1,'sigplot',0,'CapSize',5,'linewidth',1.5,'markersize',3,'colors',sescolors);

plotdnow = cfxses{allrunses};
postketidx = [3:3:15];
plotdnow(postketidx) = [];
% only use preKTM sleep

plotx = allrunses + (rand(size(plotdnow))-0.5)/10;
plot(ax,plotx,plotdnow,'o','markersize',3,...
        'markerfacecolor',1-0.6*(1-sescolors{5}),'MarkerEdgeColor',1-0.6*(1-sescolors{5}))
errorbar(ax,allrunses,nanmean(plotdnow),sem(plotdnow),'-','linewidth',1.5,'CapSize',5,'color',0.5*sescolors{5});
xticks(ax,1:allrunses)
xticklabels(ax,tmpstr)
ylabel(ax,'Log Abs(Correlation)')
set(ax, 'YGrid', 'on', 'XGrid', 'off')

% [g,p] = AllGroupTtest(tmpd);
[g,p] = AllGroupSignRank(tmpd);
    gplot = {g{1},g{4},g{5}};
    pplot = [p(1),p(4),p(5)];
myoverlapsigstar_v3(ax,gplot,pplot,'starty',-4.2,'ygap',0.06)

% nsess = 4;
% tmpd = squeeze(ratd(:,:,1));
% tmpd2 = tmpd;
% allses = 1:nsess;
% sid = KetSesMap(allses);
% tmpd2(:,allses) = tmpd2(:,sid);
% tmpd = tmpd2;

% ax = subplot('position',pos);
% MyPairSEMPlot(tmpd,ax,1:nsess,1,'sigplot',0,'markersize',3.5,'linewidth',1)
% xticks(ax,1:nsess)
% xticklabels(ax,sesstr)
% ylabel(ax,'Log Abs(Correlation)')

% [g,p] = AllGroupTtest(tmpd);
% % [g,p] = AllGroupRankSum(plotd);
% gplot = {g{1},g{4},g{5}};
% pplot = [p(1),p(4),p(5)];
% 
% myoverlapsigstar_v3(ax,gplot,pplot,'starty',-4.2,'ygap',0.06)
xtickangle(ax,45)
xlim(ax,[0.5,allrunses+0.5])

th = title(ax,'Cofiring strength across sessions','FontWeight','normal','FontSize',txtfonts(2),'HorizontalAlignment','center');
th.Position(2) = th.Position(2) + 0.2; 
end

function PlotDfrom2NNXSes(pos,Inpath)
global txtfonts sescolors sesstr
load([Inpath,'xrslpdSpkcri1.mat'],'xrslpd');
nsess = size(xrslpd,2)-1;
allrunses = nsess+1;

tmpstr = sesstr;
tmpstr{allrunses} = 'REMS';

tmpd = nan(10,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    cofnow = becolumn(xrslpd{sid});
    tmpd(1:length(cofnow),is) = cofnow;
end

ax = subplot('position',pos);
set(ax,'fontsize',txtfonts(4))
ylim(ax,[4 14])
MyPairSEMPlot_v2(tmpd,ax,1:nsess,1,'sigplot',0,'CapSize',5,'linewidth',1.5,'markersize',3,'colors',sescolors);

plotdnow = xrslpd{allrunses};
postketidx = [3:3:15];
plotdnow(postketidx) = [];
% only use preKTM sleep
plotx = allrunses + (rand(size(plotdnow))-0.5)/10;
plot(ax,plotx,plotdnow,'o','markersize',3,...
        'markerfacecolor',1-0.6*(1-sescolors{5}),'MarkerEdgeColor',1-0.6*(1-sescolors{5}))
errorbar(ax,allrunses,nanmean(plotdnow),sem(plotdnow),'-','linewidth',1.5,'CapSize',5,'color',0.5*sescolors{5});
xticks(ax,1:allrunses)
xticklabels(ax,tmpstr)
ylabel(ax,'Intrinsic dimensionality')
set(ax, 'YGrid', 'on', 'XGrid', 'off')

%[g,p] = AllGroupTtest(tmpd);
[g,p] = AllGroupSignRank(tmpd);
    gplot = {g{1},g{4},g{5}};
    pplot = [p(1),p(4),p(5)];
myoverlapsigstar_v3(ax,gplot,pplot,'starty',12.5,'ygap',0.2,'starsize',7,'blockdis',2.5)

xtickangle(ax,45)
xlim(ax,[0.5,allrunses+0.5])

end

function PlotMuExampleXSes(pos,Inpath)
global txtfonts sescolors sesstr
load([Inpath,'dspkSpkcri1.mat'],'dspk','dzspk','allmu')
nsess = 2;
ax = subplot('position',pos{1});
set(ax,'FontSize',txtfonts(4))
hold(ax,'on')
for is = 1:nsess
    munow = allmu{2,1,is};
    [F,musig] = ecdf(munow);
    plot(ax,musig,F,'-','color',sescolors{is})
end
axis tight
xlabel(ax,'\mu=r2/r1')
ylabel(ax,'CDF(\mu)')
[lgd,icon] = legend(ax,sesstr(1:nsess),'location','southeast');
icon(3).XData(1) = icon(3).XData(1) + 0.42;
icon(5).XData(1) = icon(5).XData(1) + 0.42;
set(lgd,...
    'Position',[0.768 0.335 0.085 0.036]);
legend boxoff
% line fit
tth = title(ax,'Intrinsic dimensionality of spikes','FontWeight','normal','FontSize',txtfonts(2),'HorizontalAlignment','center');
tth.Position(1) = tth.Position(1) + 1.18;
tth.Position(2) = tth.Position(2) + 0.25;

ax = subplot('position',pos{2});
set(ax,'FontSize',txtfonts(4))
hold(ax,'on')
lgdstr = cell(1,nsess);
for is = 1:nsess
    mu = allmu{2,1,is};
    
    [F,musig] = ecdf(mu);
    [unqmusig,idx,~] = unique(musig);
    F = F(idx);
    fmu = interp1(unqmusig,F,mu);
    
    xx = log(mu);
    yy = -log(1-fmu);
    % remove nan or inf values
    goodidx = ~isnan(xx) & ~isinf(xx) & ~isnan(yy) & ~isinf(yy);
    xx = xx(goodidx);
    yy = yy(goodidx);
    % discard 10% large xx as mentioned in ref
    [~,newidx] = sort(xx,'ascend');
    keplen = round(length(newidx)*0.9);
    kepidx = newidx(1:keplen);
    xxkp = xx(kepidx);
    yykp = yy(kepidx);
    
    % linefit from original
    estd = xxkp(:)\yykp(:);
    
%     plot(ax,xx, yy, 'o','color',[0.65 0.65 0.65],'markersize',1.5)
    plot(ax,xxkp, yykp, 'o','color',1-0.75*(1-sescolors{is}),'markersize',2)
    y_est = xxkp*estd;
    hd(is) = plot(ax,xxkp, y_est, '-','color',sescolors{is});
    lgdstr{is} = [sesstr{is},' Slope ',num2str(estd,'%.3f')];
end
xlabel(ax,'log\mu')
ylabel(ax,'-log(1-CDF(\mu))')
[lgd,icon] = legend(hd,lgdstr,'location','southeast');
icon(3).XData(1) = icon(3).XData(1) + 0.3;
icon(5).XData(1) = icon(5).XData(1) + 0.3;
set(lgd,...
    'Position',[0.862 0.435 0.11 0.0367]);
legend boxoff
    
end

function bspk = determineBfromA(aspk)
switch aspk
    case 0
        bspk = 0;
    case 1
        bspk = 2;
    case 2
        bspk = 3;
    case 3
        bspk = 1;
    otherwise
        bspk = 0;
end
end

function plotrandomspk(ax,allbspk,xoffset,yoffset,markersize,len)
nt = size(allbspk,2);
nc = size(allbspk,1);

% plot all spikes
for tbin = 1:nt
    for ic = 1:nc
        spknow = allbspk(ic,tbin);
        if spknow > 0
           xoff = 0.1+rand(1,spknow)/2;
           xoff = cumsum(xoff);
           if xoff(end) > 0.9
               xoff = xoff./xoff(end)*0.9;
           end
        else
            continue
        end

        crgb = [0.1 0.1 0.1];

        xall = xoff+tbin-1 + xoffset;
        yall = ic*ones(1,length(xall))+yoffset;
        for ispk = 1:length(xall)
            plot(ax,[xall(ispk),xall(ispk)],[-len len]+yall(ispk),'-',...
                'color',crgb,'linewidth',markersize)
        end
    end
end

end