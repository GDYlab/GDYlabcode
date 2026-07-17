clear;clc;close all
dbstop if error


% find REM examples in URat3
% find run examples in
%% define subpanel positions
testpos = 1;
% REM LDS mode and ELBO in fiting
[~,pos11] = tight_subplot(1,2,[0.03 0.03],[0.58 0.08],[0.015 0.62],0);
% REM fiting CA1 tsne examples
[~,pos12] = tight_subplot(3,4,[0.01 0.01],[0.585 0.085],[0.43 0.27],0);
% REM error vs shuffle
pos13 = [0.80 0.59 0.185 0.33];


% Run LDS cartoon
[~,pos21] = tight_subplot(1,2,[0.03 0.03],[0.08 0.58],[0.015 0.62],0);
% Run fiting CA1 tsne in run1 and ketamine
[~,pos22] = tight_subplot(3,4,[0.01 0.01],[0.085 0.585],[0.43 0.27],0);
% Run error across sessions
pos23 = [0.80 0.09 0.185 0.33];

% test position
if testpos
f = figure_letter(0.43);
testsubpos(pos11)
testsubpos(pos12)
testsubpos(pos13)
testsubpos(pos21)
testsubpos(pos22)
testsubpos(pos23)
end

close all
f = figure_letter(0.43);
%% define parameters
subpanelfont = 12;
labelfont = 7;
legendfont = 5;
smallfont = 4.5;

global txtfonts sescolors sesstr pdfcolors cmap
txtfonts = [subpanelfont,labelfont,legendfont,smallfont];
sescolors = Color_Ketamine_v2;
pdfcolors = {[0.68 0.4 0.1],[0.1 0.4 0.68]};
sesstr = {'Run1','Run2Bi','Run2woBi','Run3'};
cmap = getNCLColmap('WhiteYellowOrangeRed.rgb',100);
%% add subpanel annotations
addsubp = 1;
if addsubp
    panellbpos = {[0 1],[0.385 1],[0.755 1],[0 0.515],[0.385 0.515],[0.755 0.515]};
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

% REM LDS Predict example 1
Inpath{1,1} =  '/media/yuchen/data14/CAResults/Solver/KTM_Sleep/URat3/kd1_2019-12-18_11-19-51/LDS/TbinSmooth_FitEachREM/EMG_Vcri2_SlpLen10Tbin0.2Kstd2_latentpct0.99/'; 
% REM LDS Predict example 2
Inpath{1,2} =  '/media/yuchen/data14/CAResults/Solver/KTM_Sleep/XRat1/kd1_2024-05-23_14-28-50/LDS/TbinSmooth_FitEachREM/EMG_Vcri2_SlpLen10Tbin0.2Kstd2_latentpct0.99/'; 
% REM LDS Predict example 3
Inpath{1,3} =  '/media/yuchen/data14/CAResults/Solver/KTM_Sleep/URat4/kd1_2024-05-23_16-33-37/LDS/TbinSmooth_FitEachREM/EMG_Vcri2_SlpLen10Tbin0.2Kstd2_latentpct0.99/'; 
% REM one step error
Inpath{1,4} = [getDataRoot(),'CAResults/Sets_Ave/KTM_Sleep/KTMUActNew/LDS/TbinSmooth_FitEachREM_PredictStepErr/EMG_Vcri2_SlpLen10Tbin0.2Kstd2_latentpct0.99/'];

% REM LDS predict run example 1
Inpath{2,1} = [getDataRoot(),'CAResults/Solver/KTM_XRunSleep/URat1/kd1_2016-01-30_12-23-12/LDS/TbinSmooth_MazeLap_PredictTSNE_Catlap_KetImpactRun2/Twin300Step120/EMG_Vcri2_SlpLen10Tbin0.2Kstd5_latentpct0.99/Tbin0.2Kstd5Vcri10Tlen1Discri50/']; 
% REM LDS predict run example 2
Inpath{2,2} = [getDataRoot(),'CAResults/Solver/KTM_XRunSleep/XRat1/kd1_2024-05-23_14-28-50/LDS/TbinSmooth_MazeLap_PredictTSNE_Catlap_KetImpactRun2/Twin300Step120/EMG_Vcri2_SlpLen10Tbin0.2Kstd5_latentpct0.99/Tbin0.2Kstd5Vcri10Tlen1Discri50/']; 
% REM predict run one step error
Inpath{2,3} = [getDataRoot(),'CAResults/Sets_Ave/KTM_XRunSleep/KTMUActNew/LDS/TbinSmooth_MazeLap_PredictStepErr_KetImpactRun2/Twin300Step120/EMG_Vcri2_SlpLen10Tbin0.2Kstd2_latentpct0.99/'];

Outpath =   '/media/yuchen/data14/CAResults/Figures/Ketamine/';
figurename = 'LSDREMandRun_v5';

%% make subplots
% plot REM one step error
PlotREMPredictErr(pos13,Inpath{1,4})

% plot circuit model in REM
PlotCriCuitModelInREM(pos11)

% plot circuit model in run
PlotCriCuitModelInRun(pos21)

% plot LDS fitting and prediction
PlotLDSFitPredict(pos11)

% plot REM CA1 predict TSNE
PlotREMPredictTSNE(pos12,Inpath{1,1},Inpath{1,2},Inpath{1,3})

% plot REM CA1 predict Run TSNE
PlotREMPredictRunTSNE(pos22,{Inpath{2,1},Inpath{2,2}})

% plot REM to Run one step error
PlotREMPredictRunErr(pos23,Inpath{2,3})


%% save the result
File_Path = strcat(Outpath,figurename,'.fig');
saveas(gcf, File_Path);
File_Path = strcat(Outpath,figurename,'.png');
print(gcf,File_Path,'-dpng','-r300'); 
% export_fig([Outpath,figurename,'.pdf'])
print('-painters',gcf,'-dpdf',[Outpath,figurename,'_AI.pdf']) % this is for illustrator
close all
%% subpanel plot function
function PlotCriCuitModelInREM(pos)
global txtfonts sescolors sesstr cmap
% cmap=getNCLColmap('MPL_PuRd.rgb',100);
% load([Inpath,'preslpcorr',parastr,'.mat'],'preslpcorr','preslpcorrbinm','preslpcorrdir','oridforplot')
ax = subplot('position',pos{1});
hold(ax,'on')
xlim(ax,[0 3.2])
ylim(ax,[-0.25 1.75])
axis off
% plot intrinsic/external/read out
rectangle(ax,'Position',[0 1 1.2 0.5],'Curvature',0.2,'LineStyle','none',...
    'facecolor',[0.2 0.4 0.8 0.4])
rectangle(ax,'Position',[0 0 1.2 0.5],'Curvature',0.2,'LineStyle','none',...
    'facecolor',[0.2 0.4 0.8 0.4])
rectangle(ax,'Position',[2 1 1.2 0.5],'Curvature',0.2,'LineStyle','none',...
    'facecolor',[0.9 0.2 0.2 0.4])
text(ax,0.6,1.25,{'Intrinsic','Dynamics'},'HorizontalAlignment','center','FontSize',txtfonts(2),'FontWeight','normal')
text(ax,2.6,1.25,{'External','Inputs'},'HorizontalAlignment','center','FontSize',txtfonts(2),'FontWeight','normal')
text(ax,0.6,0.25,{'CA1','Activities'},'HorizontalAlignment','center','FontSize',txtfonts(2),'FontWeight','normal')

% % show CA3 recurrent connection
% circular_arrow(ax, 0.15, [0.5 1.65], 90, 300, 1, [0.2 0.4 0.8], 5,'cback3',1.5);
% text(ax,0.5,2,{'Recurrent','connection'},'HorizontalAlignment','center','FontSize',txtfonts(3)+1)

% plot CA3 to CA1 projection
% plot(ax,[0.5 0.5 0.8],[1 0.25 0.25],'-','linewidth',3,'color',[0.2 0.4 0.8])
myarrowplot(ax,0.6,1,0,-0.45,'color',[0.2 0.4 0.8],'linewidth',3,'markersize',4)
% text(ax,0.35,0.12,{'Projection'},'HorizontalAlignment','center','FontSize',txtfonts(3)+1)

% plot minimum EC inputs

myarrowplot(ax,0.8,0.75,0,-0.2,'color',[0.9 0.2 0.2],'linewidth',1,'linestyle',':','markersize',4)
plot(ax,[2.6 2.6 0.8 0.8],[1 0.75 0.75 0.55],':','linewidth',1,'color',[0.9 0.2 0.2])

text(ax,1.5,1.75,'REMS','HorizontalAlignment','center','FontSize',txtfonts(2)+1)

text(ax,1.6,0.75,'x','color',[0.9 0.1 0.1],'fontsize',12,'HorizontalAlignment','center','FontWeight','bold')

end

function PlotCriCuitModelInRun(pos)
global txtfonts sescolors sesstr cmap
% cmap=getNCLColmap('MPL_PuRd.rgb',100);
% load([Inpath,'preslpcorr',parastr,'.mat'],'preslpcorr','preslpcorrbinm','preslpcorrdir','oridforplot')
ecld = [3,1];
itld = [1,3];
ca1c = {[0.7 0.4 0.6 0.4],[0.5 0.5 0.7 0.4]};
actxoff = [1.85,1.15];
runstr = {'Normal Run','Ketamine'};

for im = 1:2
    ax = subplot('position',pos{im});
    hold(ax,'on')
    xlim(ax,[0 3.2])
    ylim(ax,[-0.25 1.75])
    axis off
    % plot three brain regions
    rectangle(ax,'Position',[0 1 1.2 0.5],'Curvature',0.2,'LineStyle','none',...
        'facecolor',[0.2 0.4 0.8 0.4])
    rectangle(ax,'Position',[0+actxoff(im) 0 1.2 0.5],'Curvature',0.2,'LineStyle','none',...
        'facecolor',ca1c{im})
    rectangle(ax,'Position',[2 1 1.2 0.5],'Curvature',0.2,'LineStyle','none',...
        'facecolor',[0.9 0.2 0.2 0.4])
    rectangle(ax,'Position',[0 0 1.1 0.5],'Curvature',0.2,'LineStyle','--',...
        'edgecolor',[0.2 0.4 0.8],'linewidth',1.5)
    
    text(ax,0.6,1.25,{'Intrinsic','Dynamics'},'HorizontalAlignment','center','FontSize',txtfonts(2),'FontWeight','normal')
    text(ax,2.6,1.25,{'External','Inputs'},'HorizontalAlignment','center','FontSize',txtfonts(2),'FontWeight','normal')
    text(ax,0.55,0.25,{'REMS LDS','Prediction'},'HorizontalAlignment','center','FontSize',txtfonts(2)-0.5,'FontWeight','normal')
    text(ax,0.6+actxoff(im),0.25,{'CA1','Activities'},'HorizontalAlignment','center','FontSize',txtfonts(2),'FontWeight','normal')
    
    myarrowplot(ax,0.5,1,0,-0.42,'color',[0.2 0.4 0.8],'linewidth',1,'markersize',3,'linestyle','--')
    % text(ax,0.35,0.12,{'Projection'},'HorizontalAlignment','center','FontSize',txtfonts(3)+1)
    
    % plot EC inputs
    
    myarrowplot(ax,0.65+actxoff(im),0.75,0,-0.2,'color',[0.9 0.2 0.2],'linewidth',ecld(im),'linestyle','-','markersize',4)
    plot(ax,[2.6 2.6 0.65+actxoff(im) 0.65+actxoff(im)],[1 0.75 0.75 0.55],'-','linewidth',ecld(im),'color',[0.9 0.2 0.2])
    myarrowplot(ax,0.55+actxoff(im),0.75,0,-0.2,'color',[0.2 0.4 0.8],'linewidth',itld(im),'linestyle','-','markersize',4)
    plot(ax,[0.6 0.6 0.55+actxoff(im) 0.55+actxoff(im)],[1 0.75 0.75 0.55],'-','linewidth',itld(im),'color',[0.2 0.4 0.8])
    
    text(ax,1.5,1.75,runstr{im},'HorizontalAlignment','center','FontSize',txtfonts(2)+1)
    
    if im == 2
        myarrowplot(ax,0.65+actxoff(im),-0.2,0,0.18,'color',[0.7 0.7 0.7],'linewidth',1,'linestyle','-','markersize',4)
        plot(ax,[0.5 0.5 0.65+actxoff(im) 0.65+actxoff(im)],[-0.05 -0.2 -0.2 -0.1],'-','linewidth',1,'color',[0.7 0.7 0.7])
        text(ax,1.2,-0.35,{'More similar','REMS ''intrusion'''},'HorizontalAlignment','center','FontSize',txtfonts(3))
    end
end
end

function PlotLDSFitPredict(pos)
global txtfonts sescolors sesstr cmap
% cmap=getNCLColmap('MPL_PuRd.rgb',100);
% load([Inpath,'preslpcorr',parastr,'.mat'],'preslpcorr','preslpcorrbinm','preslpcorrdir','oridforplot')
rng(83)

ax = subplot('position',pos{2});
hold(ax,'on')
xlim(ax,[-0.2 3.2])
ylim(ax,[-2.15 2.85])
axis off
% plot LDS
text(ax,1.5,1.8,{'Linear dynamical system','(LDS) Model of REMS'},'HorizontalAlignment','center','FontSize',txtfonts(2)-1,'FontWeight','normal')
rectangle(ax,'Position',[0 1 3 1.6],'Curvature',0.2,'LineStyle','-',...
    'edgecolor',[0 0 0],'facecolor',[0.2 0.4 0.8 0.15])

% plot CA1 activities
rectangle(ax,'Position',[0 -2.15 3 2],'Curvature',0.1,'LineStyle','none',...
    'facecolor',[0.4 0.6 0.8 0.1])
nc = 8;
rate = rand([nc,1])*10;
for ic = 1:nc
   ploty = ic*0.15-2.15;
   rnow = round(rate(ic));
   spkt = rand([rnow,1])*2.8+0.05;
   plot(ax,spkt,ploty*ones(size(spkt)),'ks','LineStyle','none','MarkerFaceColor','k','markersize',2)
end

text(ax,1.5,-0.4,{'Measured CA1','Activities'},'HorizontalAlignment','center','FontSize',txtfonts(2)-1,'FontWeight','normal')


myarrowplot(ax,0.9,0.05,0,0.6,'color',[0.2 0.2 0.2],'linewidth',3,'markersize',5)
myarrowplot(ax,2.1,0.8,0,-0.6,'color',[0.2 0.2 0.2],'linewidth',3,'markersize',5)
text(ax,0.2,0.45,'Inference','HorizontalAlignment','center','FontSize',txtfonts(2))
text(ax,2.8,0.45,'Prediction','HorizontalAlignment','center','FontSize',txtfonts(2))

end

function PlotREMPredictTSNE(pos,inpath1,inpath2,inpath3)
global txtfonts sescolors sesstr cmap
cmap1 = getNCLColmap('MPL_PuRd.rgb',100);
cmap1 = flipud(cmap1);
cmap1 = cmap1(1:60,:);

cmap2 = getNCLColmap('MPL_GnBu.rgb',100);
cmap2 = flipud(cmap2);
cmap2 = cmap2(1:60,:);
% othertc = [0.94 0.91 0.90];
% otherpc = [0.90 0.91 0.94];
othertc = [0.8 0.8 0.8];
otherpc = [0.8 0.8 0.8];

examses = [1,1,3,3];
examid = [7,17,19,20];
highmsize = 1.5;
for iplot = 1:4
    sesnow = examses(iplot);
    idnow = examid(iplot);
    load([inpath1,'tsneS',num2str(sesnow),'.mat'],'tsne','lbl')
    ax = subplot('position',pos{iplot});
    hold(ax,'on')
    set(ax,'fontsize',txtfonts(3));
    tsneplottool(ax,tsne,lbl,idnow,cmap1,cmap2,othertc,highmsize)
    if iplot == 1
        colormap(ax,cmap1)
        cbpos =  pos{iplot};
        cbpos(1) = cbpos(1) + 0.005;
        cbpos(3) = cbpos(3) - 0.01;
        cbpos(2) = cbpos(2) + cbpos(4) + 0.006;
        cbpos(4) = 0.01;
        cb = colorbar(ax,'Position',cbpos,'Location','northoutside','fontsize',txtfonts(4));
        cb.Label.String = 'Actual';
        cb.Ticks = [0 1];
        cb.TickLabels = {'Start','End'};
        cb.Label.Position(2) = cb.Label.Position(2) - 1;
    end
    if iplot == 2
        colormap(ax,cmap2)
        cbpos =  pos{iplot};
        cbpos(1) = cbpos(1) + 0.005;
        cbpos(3) = cbpos(3) - 0.01;
        cbpos(2) = cbpos(2) + cbpos(4) + 0.006;
        cbpos(4) = 0.01;
        cb = colorbar(ax,'Position',cbpos,'Location','northoutside','fontsize',txtfonts(4));
        cb.Label.String = 'Prediction';
        cb.Ticks = [0 1];
        cb.TickLabels = {'Start','End'};
        cb.Label.Position(2) = cb.Label.Position(2) - 1;
    end
    if iplot == 3
        hd = plot(ax,nan,nan,'o','MarkerFaceColor',[0.7 0.7 0.7],...
            'markeredgecolor',[0.7 0.7 0.7],'markersize',4,'LineStyle','none');
        cbpos =  pos{iplot};
        cbpos(1) = cbpos(1) + 0.005;
        cbpos(3) = cbpos(3) - 0.01;
        cbpos(2) = cbpos(2) + cbpos(4) + 0.006;
        cbpos(4) = 0.01;
        [lgd,icon] = legend(hd,'Other REMS');
        icon(1).Position(1) = icon(1).Position(1) - 0.2;
        lgd.Position = cbpos;
        legend boxoff
    end
    
end

examses = [1,2,3,3];
examid = [7,3,4,16];

for iplot = 1:4
    sesnow = examses(iplot);
    idnow = examid(iplot);
    load([inpath2,'tsneS',num2str(sesnow),'.mat'],'tsne','lbl')
    ax = subplot('position',pos{iplot+4});
    hold(ax,'on')
    tsneplottool(ax,tsne,lbl,idnow,cmap1,cmap2,othertc,highmsize)
    if iplot == 1
        ylb = ylabel(ax,'TSNE2','fontsize',txtfonts(3));
    end
end


examses = [1,1,2,2];
examid = [18,1,2,12];
for iplot = 1:4
    sesnow = examses(iplot);
    idnow = examid(iplot);
    load([inpath3,'tsneS',num2str(sesnow),'.mat'],'tsne','lbl')
    ax = subplot('position',pos{iplot+8});
    hold(ax,'on')
    tsneplottool(ax,tsne,lbl,idnow,cmap1,cmap2,othertc,highmsize)
    if iplot == 2
        xlb = xlabel(ax,'TSNE1','fontsize',txtfonts(3));
        xlb.Position(1) = xlb.Position(1) + 115;
        
    end
end

end

function PlotREMPredictRunTSNE(pos,allpath)
global txtfonts sescolors sesstr cmap
cmap1 = getNCLColmap('MPL_PuRd.rgb',100);
cmap1 = flipud(cmap1);
cmap1 = cmap1(1:60,:);

cmap2 = getNCLColmap('MPL_GnBu.rgb',100);
cmap2 = flipud(cmap2);
cmap2 = cmap2(1:60,:);
% othertc = [0.94 0.91 0.90];
% otherpc = [0.90 0.91 0.94];
othertc = [0.8 0.8 0.8];
otherpc = [0.8 0.8 0.8];

examses = [1,1,1,1];
examid = [12,12,1,15];
examdir = [1,2,1,1];
ratidx = [1,1,2,2];
highmsize = 2.5;
laplim = 50;
for iplot = 1:4
    sesnow = examses(iplot);
    idnow = examid(iplot);
    ridnow = ratidx(iplot);
    dirnow = examdir(iplot);
    pathnow = allpath{ridnow};
    load([pathnow,'tsneS',num2str(sesnow),'D',num2str(dirnow),'.mat'],'tsne','lbl')
    ax = subplot('position',pos{iplot});
    hold(ax,'on')
    set(ax,'fontsize',txtfonts(3));
    tsneplottool(ax,tsne,lbl,idnow,cmap1,cmap2,othertc,highmsize,laplim)
    if iplot == 1
        ylb = ylabel(ax,{'','Run1'},'fontsize',txtfonts(3));
%         ylb.Position(2) = ylb.Position(2) - 30;
        colormap(ax,cmap1)
        cbpos =  pos{iplot};
        cbpos(1) = cbpos(1) + 0.005;
        cbpos(3) = cbpos(3) - 0.01;
        cbpos(2) = cbpos(2) + cbpos(4) + 0.006;
        cbpos(4) = 0.01;
        cb = colorbar(ax,'Position',cbpos,'Location','northoutside','fontsize',txtfonts(4));
        cb.Label.String = 'Actual';
        cb.Ticks = [0 1];
        cb.TickLabels = {'Start','End'};
        cb.Label.Position(2) = cb.Label.Position(2) - 1;
    end
    if iplot == 2
        colormap(ax,cmap2)
        cbpos =  pos{iplot};
        cbpos(1) = cbpos(1) + 0.005;
        cbpos(3) = cbpos(3) - 0.01;
        cbpos(2) = cbpos(2) + cbpos(4) + 0.006;
        cbpos(4) = 0.01;
        cb = colorbar(ax,'Position',cbpos,'Location','northoutside','fontsize',txtfonts(4));
        cb.Label.String = 'Prediction';
        cb.Ticks = [0 1];
        cb.TickLabels = {'Start','End'};
        cb.Label.Position(2) = cb.Label.Position(2) - 1;
    end
    if iplot == 3
        hd = plot(ax,nan,nan,'o','MarkerFaceColor',[0.7 0.7 0.7],...
            'markeredgecolor',[0.7 0.7 0.7],'markersize',4,'LineStyle','none');
        cbpos =  pos{iplot};
        cbpos(1) = cbpos(1) + 0.005;
        cbpos(3) = cbpos(3) - 0.01;
        cbpos(2) = cbpos(2) + cbpos(4) + 0.006;
        cbpos(4) = 0.01;
        [lgd,icon] = legend(hd,'Other Laps');
        icon(1).Position(1) = icon(1).Position(1) - 0.2;
        lgd.Position = cbpos;
        legend boxoff
    end
    
end

examses = [2,2,2,2];
examid = [4,14,12,3];
examdir = [2,2,1,1];
ratidx = [1,1,2,2];
laplim = 20;
for iplot = 1:4
    sesnow = examses(iplot);
    idnow = examid(iplot);
    ridnow = ratidx(iplot);
    dirnow = examdir(iplot);
    pathnow = allpath{ridnow};
    load([pathnow,'tsneS',num2str(sesnow),'D',num2str(dirnow),'.mat'],'tsne','lbl')
    ax = subplot('position',pos{iplot+4});
    hold(ax,'on')
    set(ax,'fontsize',txtfonts(3));
    tsneplottool(ax,tsne,lbl,idnow,cmap1,cmap2,othertc,highmsize,laplim)
    if iplot == 1
        ylb = ylabel(ax,{'TSNE2','Run2Bi'},'fontsize',txtfonts(3));        
    end
    
end


examses = [3,3,3,3];
examid = [18,12,18,25];
examdir = [1,2,2,2];
ratidx = [1,1,2,2];
laplim = 50;
for iplot = 1:4
    sesnow = examses(iplot);
    idnow = examid(iplot);
    ridnow = ratidx(iplot);
    dirnow = examdir(iplot);
    pathnow = allpath{ridnow};
    load([pathnow,'tsneS',num2str(sesnow),'D',num2str(dirnow),'.mat'],'tsne','lbl')
    ax = subplot('position',pos{iplot+8});
    hold(ax,'on')
    set(ax,'fontsize',txtfonts(3));
    tsneplottool(ax,tsne,lbl,idnow,cmap1,cmap2,othertc,highmsize,laplim)
    if iplot == 1
        ylb = ylabel(ax,{'','Run3'},'fontsize',txtfonts(3));        
    end
    
    if iplot == 2
        xlb = xlabel(ax,'TSNE1','fontsize',txtfonts(3));
        xlb.Position(1) = xlb.Position(1) + 45;
    end
end


end

function PlotREMPredictErr(pos,inpath1)
global txtfonts sescolors sesstr cmap

cmap2 = getNCLColmap('MPL_GnBu.rgb',100);
cmap2 = flipud(cmap2);
cmap2 = cmap2(1:60,:);
load([inpath1,'remsteperr.mat'],'remsteperr')

nsess = max(remsteperr(:,1));
ndata = max(remsteperr(:,end));
plotd = nan(nsess*ndata,2);
for is = 1:nsess
    for id = 1:ndata
        sidx = remsteperr(:,1) == is;
        ratidx = remsteperr(:,end) == id;
        allidx = sidx & ratidx;
        plotd(is+(id-1)*nsess,1) = nanmedian(remsteperr(allidx,2));
        plotd(is+(id-1)*nsess,2) = nanmedian(remsteperr(allidx,3));

    end
end
ax = subplot('position',pos);
set(ax,'fontsize',txtfonts(3))
MyPairSEMBarPlot(plotd,ax,1:2,1,'sigplot',0,'color',{cmap2(20,:),[0.6 0.6 0.6]},'CapSize',8)
[g,p] = AllGroupSignRank(plotd);
myoverlapsigstar_v3(ax,g,p,'starty',4.7,'starsize',7)
xticks(ax,1:2)
xlim(ax,[0.1 2.9])
xticklabels(ax,{'Data','Shuffle'})
ylabel(ax,'LDS One Step Error')

end

function PlotREMPredictRunErr(pos,inpath1)
global txtfonts sescolors sesstr cmap
kstd = 2;
runvcri = 10;
tlencri = 1;
discri = 50;
tbin = 0.2;
runpara = ['Tbin',num2str(tbin),'Kstd',num2str(kstd),'Vcri',num2str(runvcri),'Tlen',num2str(tlencri),'Discri',num2str(discri)];

load([inpath1,'laperr',runpara,'.mat'],'laperr')

nsess = max(laperr(:,1));
ndata = max(laperr(:,end));
plotd = nan(2*ndata,nsess);
for is = 1:nsess
    sesidx = KetSesMap(is);
    for id = 1:ndata
        for idir = 1:2
            sidx = laperr(:,1) == sesidx;
            ratidx = laperr(:,end) == id;
            idiridx = laperr(:,2) == idir;
            allidx = sidx & ratidx & idiridx;
            plotd(idir+(id-1)*2,is) = nanmedian(laperr(allidx,3));
        end
    end
end
ax = subplot('position',pos);
set(ax,'fontsize',txtfonts(3))
MyPairSEMBarPlot(plotd,ax,1:nsess,1,'colors',sescolors,'sigplot',0,'CapSize',6,'linewidth',1.5)
[g,p] = AllGroupSignRank(plotd);
myoverlapsigstar_v3(ax,g([1,4,5]),p([1,4,5]),'starty',3.8,'starsize',7)
xticks(ax,1:nsess)
xticklabels(ax,sesstr)
xlim(ax,[0.3 4.7])
ylabel(ax,'LDS One Step Error')

end

function stars = sigstrfromp(p)
if p<=1E-3
    stars='***';
elseif p<=1E-2
    stars='**';
elseif p<=0.05
    stars='*';
elseif isnan(p)
    stars='n.s.';
else
    stars='n.s.';
end
end

function tsneplottool(ax,tsne,lbl,ilap,cmap1,cmap2,othertc,highmsize,laplim)

if nargin < 9
    laplim = 100;    
end

lapidx = lbl(:,1) <= laplim;
plot(ax,tsne(lapidx,1),tsne(lapidx,2),'o','MarkerFaceColor',othertc,...
    'markeredgecolor',othertc,'markersize',0.5,'LineStyle','none')
% plot(ax,tsne(~trueidx,1),tsne(~trueidx,2),'o','MarkerFaceColor',otherpc,...
%     'markeredgecolor',otherpc,'markersize',0.5,'LineStyle','none')
actidx = lbl(:,1) == ilap-1 & lbl(:,2) == 0;
predidx = lbl(:,1) == ilap-1 & lbl(:,2) == 1;

tsnet = tsne(actidx,:);
tsenp = tsne(predidx,:);
for istep = 1:size(tsnet,1)
    cidx = round(istep/size(tsnet,1)*60);
    cidx = min([cidx,60]);
    cidx = max([cidx,1]);
    plot(ax,tsnet(istep,1),tsnet(istep,2),'o','MarkerFaceColor',cmap1(cidx,:),...
        'markeredgecolor',cmap1(cidx,:),'markersize',highmsize,'LineStyle','none')
    plot(ax,tsenp(istep,1),tsenp(istep,2),'o','MarkerFaceColor',cmap2(cidx,:),...
        'markeredgecolor',cmap2(cidx,:),'markersize',highmsize,'LineStyle','none')
end
box on
xticklabels(ax,[])
yticklabels(ax,[])

xrange = [min(tsne(:,1)),max(tsne(:,1))];
yrange = [min(tsne(:,2)),max(tsne(:,2))];

xlim(ax,[xrange(1)-diff(xrange)*0.1, xrange(2)+diff(xrange)*0.1])
ylim(ax,[yrange(1)-diff(yrange)*0.1, yrange(2)+diff(yrange)*0.1])
end
