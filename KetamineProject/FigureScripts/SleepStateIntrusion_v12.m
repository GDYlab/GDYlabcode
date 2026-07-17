clear;clc;close all
dbstop if error

%% define subpanel positions
% Original positions scaled: y_new = y_old*0.70 + 0.28, h_new = h_old*0.70
testpos = 0;

% low/high vel theta asymmetry
[~,pos11] = tight_subplot(2,2,[0.03 0.023],[0.459 0.396],[0.06 0.795],0);
% REM theta asymmetry
pos12 = [0.235 0.494 0.065 0.067];
% theta asymmetry stats
[~,pos13] = tight_subplot(1,2,[0.035 0.045],[0.310 0.575],[0.06 0.685],0);

% amplitude phase coupling across run and REM
[~,pos31] = tight_subplot(1,3,[0.035 0.015],[0.504 0.396],[0.41 0.35],0);
% Gamma amplitude vs theta phase and phase lag
[~,pos32] = tight_subplot(1,2,[0.035 0.045],[0.315 0.545],[0.405 0.35],0);

% spike phase example
[~,pos21] = tight_subplot(2,3,[0.02 0.015],[0.465 0.390],[0.72 0.01],0);
% kl div between run and rem spike phase distribution
[~,pos22] = tight_subplot(1,2,[0.035 0.04],[0.320 0.563],[0.735 0.01],0);

% rate correlation in run and preREM
[~,pos41] = tight_subplot(2,2,[0.025 0.045],[0.816 0.028],[0.065 0.72],0);
pos42 = [0.09 0.670 0.18 0.10];

% cofiring similarity in run and preREM
[~,pos51] = tight_subplot(2,3,[0.03 0.002],[0.805 0.028],[0.38 0.34],0);
pos52 = [0.44 0.670 0.18 0.10];

% theta cycle spike ppl similarity with preREM
[~,pos61] = tight_subplot(2,2,[0.03 0.025],[0.805 0.028],[0.74 0.05],0);
pos62 = [0.77 0.670 0.18 0.10];

% new bottom row: 4 subplots
pos_pca3d   = [0.06 0.05 0.15 0.185];
pos_projidx = [0.29 0.05 0.17 0.185];
pos_featraj = [0.54 0.05 0.17 0.185];
pos_sliding = [0.79 0.05 0.17 0.175];

if testpos
    f = figure_letter(0.68);
    testsubpos(pos11); testsubpos(pos12); testsubpos(pos13);
    testsubpos(pos21); testsubpos(pos22);
    testsubpos(pos31); testsubpos(pos32);
    testsubpos(pos41); testsubpos(pos42);
    testsubpos(pos51); testsubpos(pos52);
    testsubpos(pos61); testsubpos(pos62);
    testsubpos(pos_pca3d); testsubpos(pos_projidx);
    testsubpos(pos_featraj); testsubpos(pos_sliding);
end
close all
f = figure_letter(0.64);

%% define parameters
subpanelfont = 12;
labelfont    = 7;
legendfont   = 5;
smallfont    = 4.5;

global txtfonts sescolors sesstr pdfcolors cmap
txtfonts  = [subpanelfont, labelfont, legendfont, smallfont];
sescolors = Color_Ketamine_v2;
pdfcolors = {[0.68 0.4 0.1],[0.1 0.4 0.68]};
sesstr    = {'Run1','Run2Bi','Run2woBi','Run3'};
cmap      = getNCLColmap('WhiteYellowOrangeRed.rgb',100);

%% add subpanel annotations
addsubp = 1;
if addsubp
    panellbpos = {[0 1],[0.355 1],[0.69 1],[0 0.65],[0.355 0.65],[0.69 0.65],...
                  [0 0.265],[0.235 0.265],[0.485 0.265],[0.738 0.265]};
    panellb    = {'a','b','c','d','e','f','g','h','i','j'};
    for isub = 1:length(panellbpos)
        posstr = panellbpos{isub};
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
setsname = 'KTMUActNew';

Inpath = cell(10,6);

Inpath{1,1} = '/media/yuchen/data14/CAResults/Solver/KTM_onMazeAll/URat1/kd1_2016-01-30_12-23-12/Decoding/MazeAllLaps_TAlignLposC_Lfit_2cm_byVel_KetImpactRun2/Twin300Step120/WinLen0.04Vel3   10   30  100/';
Inpath{1,2} = '/media/yuchen/data14/CAResults/Solver/KTM_LFPSpec/URat1/kd1_2016-01-30_12-23-12/Run/RunBic_BestTT_KetImpactRun2/Twin300Step120/Vbins0    3   10   30  100/Step3/';
Inpath{1,3} = [getDataRoot(),'CAResults/Preprocess/KTM_SleepLFP/URat1/kd1_2016-01-30_12-23-12/REMThetaLFP/BrainState_EMG_Vcri2/REMLen10/'];
Inpath{1,4} = [getDataRoot(),'CAResults/Solver/KTM_LFPSpec/URat1/kd1_2016-01-30_12-23-12/Sleep/REMBic/BrainState_EMG_Vcri2/REMLen10/'];
Inpath{1,5} = [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/',setsname,'/XRunSleep/RunREMThetaAsym_KetImpactRun2/Twin300Step120/Vbins0    3   10   30  100BrainState_Vcri2_EMG_REMLen10/'];

Inpath{2,1} = [getDataRoot(),'CAResults/Solver/KTM_XRunSleep/URat1/kd1_2016-01-30_12-23-12/SpikePhase/Spike_Local_GPhase_KetImpactRun2/Twin300Step120/BrainState_EMG_Vcri2_Slp10/'];
Inpath{2,2} = [getDataRoot(),'CAResults/Sets_Ave/KTM_XRunSleep/',setsname,'/SpikePhase/Spike_Local_GPhase_KetImpactRun2/Twin300Step120/BrainState_EMG_Vcri2_Slp10/'];

Inpath{3,1} = '/media/yuchen/data14/CAResults/Solver/KTM_LFPSpec/URat2/kd1_2016-03-16_11-03-24/Run/PAC_KetImpactRun2/Twin300Step120/';
Inpath{3,2} = '/media/yuchen/data14/CAResults/Solver/KTM_LFPSpec/URat2/kd1_2016-03-16_11-03-24/Sleep/PAC/BrainState_EMG_Vcri2_REMLen10/';
Inpath{3,3} = [getDataRoot(),'CAResults/Solver/KTM_LFPSpec/URat4/kd1_2024-05-23_16-33-37/Run/GAmp_TPhase_byDir_KetImpactRun2/Twin300Step120/'];
Inpath{3,4} = [getDataRoot(),'CAResults/Solver/KTM_LFPSpec/URat4/kd1_2024-05-23_16-33-37/Sleep/GAmp_TPhase/BrainState_EMG_Vcri2_REMLen10/'];
Inpath{3,5} = [getDataRoot(),'CAResults/Sets_Ave/KTM_LFPSpec/',setsname,'/XRunSleep/GAmp_TPhase_byDir_KetImpactRun2/Twin300Step120/BrainState_EMG_Vcri2_REMLen10/'];

Inpath{4,1} = '/media/yuchen/data14/CAResults/Sets_Ave/KTM_XRunSleep/KTMUActNew/FiringRate/RunLaps_PreREMSlp_Rate_KetImpactRun2/Twin300Step120/EMG_SlpLen10/';

Inpath{5,1} = [getDataRoot(),'CAResults/Solver/KTM_XRunSleep/XRat1/kd1_2024-05-23_14-28-50/Cofiring/PreSlpREM_PyrInt_CellVDiff_KetImpactRun2/Twin300Step120/BrainState_EMG_Vcri2_SlpLen10/'];
Inpath{5,2} = [getDataRoot(),'CAResults/Sets_Ave/KTM_XRunSleep/',setsname,'/Cofiring/PreSlpREM_PyrInt_CellVDiff_KetImpactRun2/Twin300Step120/BrainState_EMG_Vcri2_SlpLen10/'];

Inpath{6,1} = [getDataRoot(),'CAResults/Solver/KTM_XRunSleep/URat2/kd1_2016-03-16_11-03-24/PplVSimi/Ppl_RunPreREMTbinSpk_KetImpactRun2/Twin300Step120/BrainState_EMG_REMLen10/'];
Inpath{6,2} = [getDataRoot(),'CAResults/Sets_Ave/KTM_XRunSleep/',setsname,'/PplVSimi/Ppl_RunPreREMTbinSpk_KetImpactRun2/Twin300Step120/BrainState_EMG_REMLen10/'];

% new subplots data paths
Inpath_intrus  = [getDataRoot(),'CAResults/Sets_Ave/KTM_XRunSleep/',setsname,...
    '/UMAP/REMIntrusionIndex_Dir_v1_KetImpactRun2/Twin300Step120/SlpLen10/'];
Inpath_sliding = [getDataRoot(),'CAResults/Sets_Ave/KTM_XRunSleep/',setsname,...
    '/UMAP/Run2SlidingFeatures_v2_AllREM/Twin600Step150/'];

Outpath     = '/media/yuchen/data14/CAResults/Figures/Ketamine/';
figurename  = 'SleepStateIntrusion_v12';
mkdir(Outpath)

%% make subplots

PlotCofExample(pos51, Inpath{5,1})
PlotThetaPplExample(pos61, Inpath{6,1})

PlotPCA3D(pos_pca3d, Inpath_intrus)
PlotProjectionIndex(pos_projidx, Inpath_intrus)
PlotFeatureMeanTrajectory(pos_featraj, Inpath_intrus)
PlotRun2SlidingProjection(pos_sliding, Inpath_sliding)


PlotSpikeRateCorr(pos41, pos42, Inpath{4,1})
PlotCofPreREMSimilarity(pos52, Inpath{5,2})
PlotThetaPplPreREMSimilarity(pos62, Inpath{6,2})
PlotPCPhaseLag(pos32, Inpath{3,5})
PlotSpikePhaseExample(pos21, Inpath{2,1})
PlotSpikePhaseKLDiv(pos22, Inpath{2,2})
PlotLFPThetaAsymRunREM(pos13, Inpath{1,5})

PlotPCPhaseLagExample(pos32, Inpath{3,3}, Inpath{3,4})
PlotPACRunExample(pos31, Inpath{3,1})
PlotPACREMExample(pos31, Inpath{3,2})
PlotLFPThetaAsymREM(pos12, Inpath{1,4})
PlotLFPThetaAsymXRun(pos11, Inpath{1,2})
PlotLFPThetaExampleXRun(pos11, Inpath{1,1})
PlotLFPThetaExampleREM(pos12, Inpath{1,3})


%% save
File_Path = strcat(Outpath, figurename, '.fig');
saveas(gcf, File_Path);
File_Path = strcat(Outpath, figurename, '.png');
print(gcf, File_Path, '-dpng', '-r300');
print('-painters', gcf, '-dpdf', [Outpath, figurename, '_AI.pdf'])
close all

%% ===================== new plot functions ================================

function PlotPCA3D(pos, Inpath)
global txtfonts
green_light = [0.65, 0.88, 0.70];
green_med   = [0.25, 0.65, 0.42];
green_dark  = [0.05, 0.42, 0.22];
col_run1    = [125,191,194]/255;
col_run2bi  = [140,12,5]/255;
col_run2wo  = [195,185,185]/255;
col_run3    = [120,122,183]/255;
plot_colors = {green_light, col_run1, green_med, col_run2bi, col_run2wo, green_dark, col_run3};
ses_labels  = {'REMS1','Run1','REMS2','Run2Bi','Run2woBi','REMS3','Run3'};
nses_total  = 7;
dotsize = 6.5;
mksize = 6;
falpha = 0.25;

rmpath(genpath([getDataRoot(),'Rcode',filesep,'20221130_CellAss_v2_0',filesep,'00_Libs/TimeseriesAna']))
rmpath(genpath([getDataRoot(),'Rcode',filesep,'20221130_CellAss_v2_0',filesep,'00_Libs/LMT-master']))

load([Inpath, 'REM_intrusion_index_Dir.mat'], 'all_z', 'sample_info');
nfeat = size(all_z, 2);

flat_z_imp = all_z;
for ifeat = 1:nfeat
    col_vals = flat_z_imp(:, ifeat);
    nanmask  = isnan(col_vals);
    if any(~nanmask)
        col_vals(nanmask) = nanmean(col_vals);
    end
    flat_z_imp(:, ifeat) = col_vals;
end

try
    [~, scores, ~, ~, explained] = pca(flat_z_imp);
catch
    disp('PlotPCA3D: PCA failed, skipping.'); return;
end

centroids_pc = nan(nses_total, 3);
for ises = 1:nses_total
    idx = sample_info(:, 2) == ises;
    if any(idx)
        centroids_pc(ises, :) = mean(scores(idx, 1:3), 1);
    end
end

ax = subplot('position', pos);
hold(ax, 'on');
set(ax, 'fontsize', txtfonts(3));

for ises = 1:nses_total
    col  = plot_colors{ises};
    dcol = 1 - falpha * (1 - col);
    idx  = sample_info(:, 2) == ises;
    scatter3(ax, scores(idx,1), scores(idx,2), scores(idx,3), dotsize, ...
        'markerfacecolor', dcol, 'markeredgecolor', dcol);
end

for ises = 1:nses_total
    col = plot_colors{ises};
    if ~any(isnan(centroids_pc(ises,:)))
        plot3(ax, centroids_pc(ises,1), centroids_pc(ises,2), centroids_pc(ises,3), 's', ...
            'markersize', mksize, 'markerfacecolor', col*1, ...
            'markeredgecolor', [0.2,0.2,0.2], 'linewidth', 0.1);
    end
end

xlabel(ax, sprintf('PC1 (%.1f%%)', explained(1)), 'fontsize', txtfonts(4));
ylabel(ax, sprintf('PC2 (%.1f%%)', explained(2)), 'fontsize', txtfonts(4));
zlabel(ax, sprintf('PC3 (%.1f%%)', explained(3)), 'fontsize', txtfonts(4));
xlim(ax,[-4 4])
ylim(ax,[-2 2.5])
zlim(ax,[-1.5 2])
view(ax, [-12.5, 54]);

%% manual legend — 2D overlay axis at the top of the subplot
rem_cols = {green_light, green_med, green_dark,col_run2bi};
rem_lbls = {'REMS1','REMS2','REMS3','Run2Bi'};
run_cols = {col_run1, col_run2wo, col_run3};
run_lbls = {'Run1','Run2woBi','Run3'};

leg_h   = pos(4) * 0.45;   % legend overlay takes top 30% of subplot height
yoff    = 0.03;
leg_pos = [pos(1), pos(2) + pos(4) - leg_h + yoff, pos(3), leg_h];
lag_ax  = axes('Position', leg_pos, 'Color', 'none');
set(lag_ax, 'Visible', 'off');
hold(lag_ax, 'on');
xlim(lag_ax, [0 1]);  ylim(lag_ax, [0 1.5]);

rw = 0.06;  rh = 0.18;  x_gap = 0.015;  % rectangle dimensions (normalized)

% row 1 (top): REM sessions — 3 items
y1 = 0.99;
x1 = linspace(0.05, 0.2, 3);
for ii = 1:3
    xs = x1(ii);
    fill(lag_ax, [xs,xs+rw,xs+rw,xs], [y1-rh/2,y1-rh/2,y1+rh/2,y1+rh/2], ...
        rem_cols{ii}, 'EdgeColor',[0.2,0.2,0.2]);
end

xs = 0.62;
fill(lag_ax, [xs,xs+rw,xs+rw,xs], [y1-rh/2,y1-rh/2,y1+rh/2,y1+rh/2], ...
    rem_cols{end}, 'EdgeColor',[0.2,0.2,0.2]);
text(lag_ax, 0.3, y1, 'REMS1,2,3', ...
    'fontsize',txtfonts(4), 'verticalalignment','middle');
text(lag_ax, 0.7, y1, rem_lbls{end}, ...
    'fontsize',txtfonts(4), 'verticalalignment','middle');

% row 2 (bottom): Run sessions — 4 items
y2 = 0.65;
x2 = [0.05 0.28 0.62];
for ii = 1:3
    xs = x2(ii);
    fill(lag_ax, [xs,xs+rw,xs+rw,xs], [y2-rh/2,y2-rh/2,y2+rh/2,y2+rh/2], ...
        run_cols{ii}, 'EdgeColor',[0.2,0.2,0.2]);
    text(lag_ax, xs+rw+x_gap, y2, run_lbls{ii}, ...
        'fontsize',txtfonts(4), 'verticalalignment','middle');
end
% title(ax, 'Feature space (PCA)', 'fontsize', txtfonts(2), 'fontweight', 'normal');
end

function PlotProjectionIndex(pos, Inpath)
global txtfonts
green_light = [0.65, 0.88, 0.70];
green_med   = [0.25, 0.65, 0.42];
green_dark  = [0.05, 0.42, 0.22];
col_run1    = [125,191,194]/255;
col_run2bi  = [140,12,5]/255;
col_run2wo  = [195,185,185]/255;
col_run3    = [120,122,183]/255;
plot_colors = {green_light, col_run1, green_med, col_run2bi, col_run2wo, green_dark, col_run3};
ses_labels  = {'REMS1','Run1','REMS2','Run2Bi','Run2woBi','REMS3','Run3'};
nses_total  = 7;
ndir        = 2;
key_pairs   = {[4, 2], [4, 7]};
comparisons = { [4, 2], 'Run2Bi vs Run1';
                [4, 7], 'Run2Bi vs Run3';
                [4, 5], 'Run2Bi vs Run2woBi' };

load([Inpath, 'REM_intrusion_index_Dir.mat'], 'all_proj', 'sample_info', 'Datasets');
nanimals = length(Datasets);

proj_by_ses = cell(1, nses_total);
for ises = 1:nses_total
    rows = sample_info(:,2) == ises;
    proj_by_ses{ises} = all_proj(rows);
end

pvals_proj  = nan(1, size(comparisons, 1));
groups_plot = cell(1, size(comparisons, 1));
for ic = 1:size(comparisons, 1)
    p1 = comparisons{ic,1}(1);  p2 = comparisons{ic,1}(2);
    d1 = zeros(0,1);             d2 = zeros(0,1);
    for id = 1:nanimals
        for idir = 1:ndir
            r1 = (sample_info(:,1)==id) & (sample_info(:,2)==p1) & (sample_info(:,3)==idir);
            r2 = (sample_info(:,1)==id) & (sample_info(:,2)==p2) & (sample_info(:,3)==idir);
            if any(r1) && any(r2)
                v1 = all_proj(r1);  v2 = all_proj(r2);
                if ~isnan(v1) && ~isnan(v2)
                    d1(end+1) = v1; %#ok<AGROW>
                    d2(end+1) = v2; %#ok<AGROW>
                end
            end
        end
    end
    if length(d1) >= 2
        try; pvals_proj(ic) = signrank(d1 - d2); catch; pvals_proj(ic) = NaN; end
    else
        pvals_proj(ic) = NaN;
    end
    groups_plot{ic} = [p1, p2];
end

ax = subplot('position', pos);
hold(ax, 'on');
set(ax, 'fontsize', txtfonts(3));

for iel = 1:numel(proj_by_ses)
    proj_by_ses{iel} = (proj_by_ses{iel} - 0.5)*2;
end

MyBarSampleSEMPlot_v2(proj_by_ses, ax, ...
    'plotx', 1:nses_total, 'colors', plot_colors, ...
    'markersize', 2.5, 'linewidth', 1.5, 'barwidth', 0.4, 'baralpha', 0.55);

sig_groups = {};  sig_pvals = [];
for ic = 1:length(groups_plot)
    gp     = groups_plot{ic};
    is_key = any(cellfun(@(k) isequal(k,gp)||isequal(k,fliplr(gp)), key_pairs));
    if is_key || (~isnan(pvals_proj(ic)) && pvals_proj(ic) <= 0.05)
        sig_groups{end+1} = gp;              %#ok<AGROW>
        sig_pvals(end+1)  = pvals_proj(ic);  %#ok<AGROW>
    end
end
ylim(ax, [-1.8, 2.2]);
if ~isempty(sig_groups)
    myoverlapsigstar_v3(ax, sig_groups, sig_pvals, 'starty', 1.85, 'ygap', 0.06, 'starsize', 6.5);
end

p4_sgn = signrank(proj_by_ses{4},0);
myoverlapsigstar_v3(ax, {4}, p4_sgn, 'starty', 1.55, 'ygap', 0.06, 'starsize', 6.5);

set(ax, 'XTick', 1:nses_total, 'XTickLabel', ses_labels);
xtickangle(ax, 40);
ylabel(ax, 'Projection distance', 'fontsize', txtfonts(3));
xlim(ax, [0.3, nses_total + 0.7]);
yticks(ax, [-1 0 1]);  yticklabels(ax, {'REMS1','','Run1'});
set(ax, 'YGrid', 'on', 'XGrid', 'off');
% title(ax, 'REM intrusion index', 'fontsize', txtfonts(2), 'fontweight', 'normal');
end

function PlotFeatureMeanTrajectory(pos, Inpath)
global txtfonts
green_light = [0.65, 0.88, 0.70];
green_med   = [0.25, 0.65, 0.42];
green_dark  = [0.05, 0.42, 0.22];
col_run1    = [125,191,194]/255;
col_run2bi  = [140,12,5]/255;
col_run2wo  = [195,185,185]/255;
col_run3    = [120,122,183]/255;
plot_colors = {green_light, col_run1, green_med, col_run2bi, col_run2wo, green_dark, col_run3};
ses_labels  = {'REMS1','Run1','REMS2','Run2Bi','Run2woBi','REMS3','Run3'};
nses_total  = 7;
feat_names  = {'Cofiring','P→PMI','P→IMI','SpkDim','Theta Asy','T-G Lag'};
feat_markers = {'o','s','d','^','v','p'};
line_gray   = [0.65 0.65 0.65];

load([Inpath, 'REM_intrusion_index_Dir.mat'], 'all_z', 'sample_info');
nfeat = size(all_z, 2);

feat_means = nan(nses_total, nfeat);
for ises = 1:nses_total
    rows = sample_info(:,2) == ises;
    if any(rows)
        feat_means(ises,:) = nanmean(all_z(rows,:), 1);
    end
end

feat_norm = feat_means;
for ifeat = 1:nfeat
    v_rem1 = feat_means(1, ifeat);
    v_run1 = feat_means(2, ifeat);
    denom  = v_run1 - v_rem1;
    if ~isnan(denom) && denom ~= 0
        feat_norm(:, ifeat) = (feat_means(:, ifeat) - v_rem1) / denom;
    end
end

ax = subplot('position', pos);
hold(ax, 'on');
set(ax, 'fontsize', txtfonts(3));

feat_norm = (feat_norm  - 0.5)*2;

for ifeat = 1:nfeat
    y_vals = feat_norm(:, ifeat);
%     y_vals = (y_vals - 0.5)*2;
    for ises = 1:nses_total-1
        if ~isnan(y_vals(ises)) && ~isnan(y_vals(ises+1))
            plot(ax, [ises, ises+1], [y_vals(ises), y_vals(ises+1)], '-', ...
                'color', line_gray, 'linewidth', 0.8);
        end
    end
    for ises = 1:nses_total
        if ~isnan(y_vals(ises))
            plot(ax, ises, y_vals(ises), feat_markers{ifeat}, ...
                'markersize', 5, ...
                'markerfacecolor', plot_colors{ises}, ...
                'markeredgecolor', [0.2 0.2 0.2], 'linewidth', 0.8);
        end
    end
end

g4 = {4};
sig_groups{1} = [2,4];
sig_groups{2} = [4,5];
sig_groups{3} = [4,7];
p4 = signrank(feat_norm(4,:),0);
sig_pvals(1) = signrank(feat_norm(2,:),feat_norm(4,:));
sig_pvals(2) = signrank(feat_norm(4,:),feat_norm(5,:));
sig_pvals(3) = signrank(feat_norm(4,:),feat_norm(7,:));

myoverlapsigstar_v3(ax, sig_groups, sig_pvals, 'starty',2.5, 'ygap', 0.15, 'starsize', 6.5);
myoverlapsigstar_v3(ax, g4, p4, 'starty', 1.8, 'ygap', 0.06, 'starsize', 6.5);


xlim(ax, [0.3, nses_total + 0.7]);
% plot(ax, xlim(ax), [0, 0], '--', 'color', [0.4 0.4 0.4], 'linewidth', 0.8);
% plot(ax, xlim(ax), [1, 1], '--', 'color', [0.4 0.4 0.4], 'linewidth', 0.8);
set(ax, 'XTick', 1:nses_total, 'XTickLabel', ses_labels);
xtickangle(ax, 40);
ylabel(ax, 'Normalized feature value', 'fontsize', txtfonts(3));
yticks(ax, [-1 0 1]);  yticklabels(ax, {'REMS1','','Run1'});
set(ax, 'YGrid', 'on', 'XGrid', 'off');

ylim(ax,[-2.8 2.8])
cur_ylim = ylim(ax);
leg_y1   = cur_ylim(2) + 0.22 * diff(cur_ylim);   % upper row: features 1-3
leg_y2   = cur_ylim(2) + 0.11 * diff(cur_ylim);   % lower row: features 4-6
ylim(ax, [cur_ylim(1), cur_ylim(2) + 0.32 * diff(cur_ylim)]);


n_row1 = 3;
leg_x1 = linspace(1, nses_total - 1.2, n_row1);
leg_x2 = linspace(1, nses_total - 1.2, nfeat - n_row1);

leg_y1 = leg_y1 + 0.5;
leg_y2 = leg_y2 + 0.5;

for ifeat = 1:n_row1
    plot(ax, leg_x1(ifeat), leg_y1, feat_markers{ifeat}, ...
        'markersize', 5, 'markerfacecolor', [0.6 0.6 0.6], ...
        'markeredgecolor', [0.2 0.2 0.2], 'linewidth', 0.8);
    text(ax, leg_x1(ifeat) + 0.28, leg_y1, feat_names{ifeat}, ...
        'fontsize', txtfonts(4), 'verticalalignment', 'middle', 'clipping', 'off');
end
for ifeat = 1:(nfeat - n_row1)
    plot(ax, leg_x2(ifeat), leg_y2, feat_markers{n_row1 + ifeat}, ...
        'markersize', 5, 'markerfacecolor', [0.6 0.6 0.6], ...
        'markeredgecolor', [0.2 0.2 0.2], 'linewidth', 0.8);
    text(ax, leg_x2(ifeat) + 0.28, leg_y2, feat_names{n_row1 + ifeat}, ...
        'fontsize', txtfonts(4), 'verticalalignment', 'middle', 'clipping', 'off');
end
% title(ax, 'Feature dynamics across sessions', 'fontsize', txtfonts(2), 'fontweight', 'normal');
end

function PlotRun2SlidingProjection(pos, Inpath)
global txtfonts  sesstr
col_bi   = [140,12,5]/255;
col_wobi = [0.55 0.55 0.55];
lw_trace = 0.5;
ms       = 4.5;

load([Inpath, 'run2sliding_group.mat'], 'all_trel', 'all_proj', 'all_isbi', 'Datasets');
nanimals = length(Datasets);

nonempty = ~cellfun(@isempty, all_trel);
if ~any(nonempty), disp('PlotRun2SlidingProjection: no data'); return; end
xmax = max(cellfun(@(x) x(end), all_trel(nonempty)));

t_all = [];  y_all = [];
for id = 1:nanimals
    if isempty(all_trel{id}), continue; end
    t = all_trel{id};
    for idir = 1:2
        yv  = all_proj{id}(:, idir);
        ok  = ~isnan(yv);
        tmp_t = t(ok);  tmp_y = yv(ok);
        t_all = [t_all; tmp_t(:)]; %#ok<AGROW>
        y_all = [y_all; tmp_y(:)]; %#ok<AGROW>
    end
end
if length(t_all) >= 3
    [Rmat, Pmat] = corrcoef(t_all, y_all);
    corr_str = sprintf('R=%.2f, P=%.3f', Rmat(1,2), Pmat(1,2));
else
    corr_str = 'r=NA';
end

ax = subplot('position', pos);
hold(ax, 'on');
set(ax, 'fontsize', txtfonts(3));

for id = 1:nanimals
    if isempty(all_trel{id}), continue; end
    t     = all_trel{id};
    is_bi = all_isbi{id};
    for idir = 1:2
        yv    = all_proj{id}(:, idir);
        yv    = (yv - 0.5)*2;
        plot(ax, t, yv, '-', 'color', [0.85 0.85 0.85], 'linewidth', lw_trace);
        bi_ok = is_bi & ~isnan(yv);
        wo_ok = ~is_bi & ~isnan(yv);
        if any(bi_ok)
            plot(ax, t(bi_ok), yv(bi_ok), 'o', 'markersize', ms, ...
                'markerfacecolor', col_bi, 'markeredgecolor', 'none');
        end
        if any(wo_ok)
            plot(ax, t(wo_ok), yv(wo_ok), 'o', 'markersize', ms, ...
                'markerfacecolor', col_wobi, 'markeredgecolor', 'none');
        end
    end
end

xlim(ax, [0, xmax * 1.02]);
xl = xlim(ax);
ylim(ax,[-1.4 1.4])
% plot(ax, xl, [0 0], '--', 'color', [0.4 0.4 0.4], 'linewidth', 0.8);
% plot(ax, xl, [1 1], '--', 'color', [0.4 0.4 0.4], 'linewidth', 0.8);
yticks(ax, [-1 0 1]);  yticklabels(ax, {'REMS1','','Run1'});

[lgd,icon] = legendnan(ax,{col_bi,col_wobi},{sesstr{2}, sesstr{3}},'marker','o','markersize',3);
icon(1).Position(1) = icon(1).Position(1) - 0.19;
icon(2).Position(1) = icon(2).Position(1) - 0.19;
lgd.Position(1) = lgd.Position(1) - 0.065;
lgd.Position(2) = lgd.Position(2) + 0.025;
legend boxoff
ylabel(ax, 'Projection distance (a.u.)', 'fontsize', txtfonts(3));
xlabel(ax, 'Time in Run2 (min)', 'fontsize', txtfonts(3),'HorizontalAlignment','center');
set(ax, 'YGrid', 'on', 'XGrid', 'off');
title(ax,corr_str, 'fontsize', txtfonts(2), 'fontweight', 'normal');
end

%% =================== original plot functions (unchanged) =================

function PlotThetaPplExample(pos,Inpath)
actcri = 4;
parastr = ['ActCellCri',num2str(actcri)];
global txtfonts sescolors sesstr cmap
load([Inpath,'preslpcorr',parastr,'.mat'],'preslpcorr','preslpcorrbinm','preslpcorrdir','oridforplot')
nsess = 4;
for is = 1:nsess
    sid = KetSesMap(is);
    tmpcorr = oridforplot{sid,2};
    ax = subplot('position',pos{is});
    set(ax,'fontsize',txtfonts(4))
    hold(ax,'on')
    imagesc(ax,tmpcorr)
    colormap(ax,cmap)
    title(ax,sesstr(is),'FontWeight','normal')
    axis tight
    if is == 3
        xlb = xlabel(ax,'pre-Run REMS theta cycles');
        xlb.Position(1) = xlb.Position(1) + 650;
        xlb.Position(2) = xlb.Position(2) + 40;
    end
    if is == 1
        ylb = ylabel(ax,'Run theta cycles');
        ylb.Position(2) = ylb.Position(2) - 45;
        ylb.Position(1) = ylb.Position(1) + 150;
    end
    yticks(ax,[]); xticks(ax,[])
    caxis(ax,[0 0.8])
    if is == nsess
        cb = addcolorbar;
        cb.Label.String = 'Cosine similarity';
    end
end
end

function PlotThetaPplPreREMSimilarity(pos,Inpath)
global txtfonts sescolors sesstr
actcri = 4;
parastr = ['ActCellCri',num2str(actcri)];
load([Inpath,'preslpcorr',parastr,'.mat'],'preslpcorr','preslpcorrbinm','meanprecorrdir')
nsess = 4;
ax = subplot('position',pos);
hold(ax,'on')
set(ax,'fontsize',txtfonts(4))
plotd = cell(1,nsess);
im = 2;
for is = 1:nsess
    sid = KetSesMap(is);
    plotd{is} = preslpcorrbinm{sid,im};
end
MyBoxPlot_v3(plotd,ax,'sigplot',0,'color',sescolors)
xticks(ax,1:nsess); xticklabels(ax,sesstr)
ylabel(ax,'Similarity vs. pre-Run REMS')
ylim(ax,[0 0.5])
g{1}=[1,2]; g{2}=[2,3]; g{3}=[2,4];
p(1)=ranksum(plotd{1},plotd{2}); p(2)=ranksum(plotd{2},plotd{3}); p(3)=ranksum(plotd{2},plotd{4});
myoverlapsigstar_v3(ax,g,p,'starty',0.45,'starsize',6.5)
set(ax,'YGrid','on','XGrid','off')
end

function PlotCofExample(pos,Inpath)
cmap = getNCLColmap('MPL_rainbow.rgb',100);
pyronly = 0;
if pyronly; pyrstr = '_pyr'; else; pyrstr = []; end
global txtfonts sescolors sesstr
load([Inpath,'rslpccorr',pyrstr,'.mat'],'rslpccorr','runcof','slpcof')
nsessslp = length(slpcof);
cri = {[-0.16,0.16,0.026],[-0.12,0.12,0.02],[-0.27,0.27,0.026]};
for is = 1:nsessslp
    ax = subplot('position',pos{is});
    set(ax,'fontsize',txtfonts(4)); hold(ax,'on')
    plotcof(ax,squeeze(slpcof(is).corr(:,:)),cmap,cri{is})
    title(ax,['pre-Run',num2str(is),' REMS'],'FontWeight','normal')
end
cri = {[-0.12,0.12,0.02],[-0.18,0.18,0.04],[-0.12,0.12,0.02]};
runstr = sesstr([1,2,4]); idir = 1;
for is = 1:3
    ax = subplot('position',pos{is+3});
    set(ax,'fontsize',txtfonts(4)); hold(ax,'on')
    plotcof(ax,squeeze(runcof(idir,is).corr(:,:)),cmap,cri{is})
    title(ax,runstr{is},'FontWeight','normal')
end
end

function PlotCofPreREMSimilarity(pos,Inpath)
pyrstr = [];
global txtfonts sescolors sesstr
nsess = 4;
load([Inpath,'allcorr',pyrstr,'.mat'],'allcorr')
ax = subplot('position',pos);
hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
plotd = cell(1,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    d1 = squeeze(allcorr(1,sid,:)); d2 = squeeze(allcorr(2,sid,:));
    dtmp = cat(1,d1(~isnan(d1)),d2(~isnan(d2)));
    plotd{is} = dtmp;
end
MyBoxPlot_v3(plotd,ax,'sigplot',0,'color',sescolors)
xticks(ax,1:nsess); xticklabels(ax,sesstr)
ylabel(ax,'Similarity vs. pre-Run REMS'); ylim(ax,[-0.5 1.4])
g{1}=[1,2]; g{2}=[2,3]; g{3}=[2,4];
p(1)=ranksum(plotd{1},plotd{2}); p(2)=ranksum(plotd{2},plotd{3}); p(3)=ranksum(plotd{2},plotd{4});
myoverlapsigstar_v3(ax,g,p,'starty',1.1,'starsize',6.5)
set(ax,'YGrid','on','XGrid','off')
end

function PlotPCPhaseLagExample(pos,Inpath,Inpath2)
thetaedge = -180:10:180;
thetabc = (thetaedge(1:end-1)+thetaedge(2:end))/2;
global txtfonts sescolors sesstr
Gband = [50,100]; vcri = 10;
load([Inpath,'cpsdGband',num2str(Gband),'Vcri',num2str(vcri),'.mat'],'cpsdses','dofses',...
    'gampmsem','lfpses','sesgampthetabin','seslfpthetabin','xsesmod','pdiffmaxcoh','maxcoh','cpsdsescomp')
dir = 2;
ax = subplot('position',pos{1}); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
nsess = 4;
for is = 1:nsess
    sid = KetSesMap(is);
    mnow = squeeze(gampmsem(sid,1,:,dir));
    semnow = squeeze(gampmsem(sid,2,:,dir));
    PlotMeanSemShade(ax,thetabc,mnow,semnow,sescolors{is},0.25,'linewidth',0.5)
    hd(is) = plot(ax,nan,nan,'-','linewidth',1,'color',sescolors{is});
    lgdstr{is} = sesstr{is};
end
load([Inpath2,'cpsdGband',num2str(Gband),'.mat'],'cpsdses','dofses',...
    'gampmsem','lfpses','sesgampthetabin','seslfpthetabin','xsesmod','pdiffmaxcoh','maxcoh','cpsdsescomp')
mnow = squeeze(gampmsem(1,1,:)); semnow = squeeze(gampmsem(1,2,:));
PlotMeanSemShade(ax,thetabc,mnow,semnow,sescolors{nsess+1},0.25,'linewidth',0.5)
hd(nsess+1) = plot(ax,nan,nan,'-','linewidth',1,'color',sescolors{nsess+1});
lgdstr{nsess+1} = 'REMS';
[lgd,icon] = legend(hd,lgdstr);
for ii = 1:5; icon((ii-1)*2+6).XData(1) = icon((ii-1)*2+6).XData(1)+0.45; end
legend boxoff
set(lgd,'Position',[0.435 0.410 0.089 0.06],'fontsize',txtfonts(4)-2);
xlim(ax,[-180 180]); xlabel(ax,'Theta phase (deg)'); ylabel(ax,'Gamma amplitude (mV)')
set(ax,'Yscale','log'); ylim(ax,[1.6e-5 5e-5])
end

function PlotPACRunExample(pos,Inpath)
global txtfonts sescolors sesstr cmap
load([Inpath,'wlpac.mat'],'wlpac','ff','frecoh')
sesplot = [1,2]; vbin = 3;
for is = 1:length(sesplot)
    ax = subplot('position',pos{is}); hold(ax,'on'); set(ax,'fontsize',txtfonts(4)-1)
    dnow = wlpac{sesplot(is),vbin};
    pcolor(ax,frecoh,ff,dnow); xlim(ax,[4 60]); ylim(ax,[4 150])
    shading flat; shading interp
    if is==2; xlabel(ax,'Frequency phase (Hz)'); end
    if is==1; ylabel(ax,'Frequency amplitude (Hz)'); else; yticks(ax,[]); end
    title(ax,sesstr{sesplot(is)},'FontWeight','normal')
    caxis(ax,[0 0.5]); colormap(ax,cmap)
end
end

function PlotPACREMExample(pos,Inpath)
global txtfonts sescolors sesstr cmap
load([Inpath,'wlpac.mat'],'wlpac','ff','frecoh')
remses = 3;
ax = subplot('position',pos{3}); hold(ax,'on'); set(ax,'fontsize',txtfonts(4)-1)
dnow = wlpac{remses,1};
pcolor(ax,frecoh,ff,dnow); xlim(ax,[4 60]); ylim(ax,[4 150])
shading flat; shading interp; yticks(ax,[])
c = addcolorbar; c.Label.String = 'Coherence';
title(ax,'REMS','FontWeight','normal')
caxis(ax,[0 0.5]); colormap(ax,cmap)
end

function PlotPCPhaseLag(pos,Inpath)
global txtfonts sescolors sesstr
gphase = [50,100]; vcri = 15;
load([Inpath,'runcpsdGband',num2str(gphase),'Vcri',num2str(vcri),'.mat'],...
    'runsesgampthetabin','runseslfpthetabin','runpdiff','runcpsdsescomp','rundof')
load([Inpath,'slpcpsdGband',num2str(gphase),'.mat'],...
    'sleepsesgampthetabin','sleepseslfpthetabin','sleeppdiff','sleepdof','sleepcpsdsescomp')
nsess = 4;
s1s3slp = sleeppdiff(:);
pairs1s3 = sleeppdiff([1,3],:)';
tmpplotd = nan(10,nsess);
for is = 1:nsess
    sid = KetSesMap(is);
    tmpd = runpdiff(sid,:); tmpd = tmpd(:);
    tmpplotd(1:length(tmpd),is) = tmpd;
end
ax = subplot('position',pos{2}); set(ax,'fontsize',txtfonts(4))
MyPairSEMPlot_v2(tmpplotd,ax,1:nsess,1,'sigplot',0,'CapSize',5,'linewidth',1.5,'markersize',2,'colors',sescolors);
plotdnow = s1s3slp; postketidx = [3:3:15]; plotdnow(postketidx) = [];
plotx = nsess+1+(rand(size(plotdnow))-0.5)/10;
plot(ax,plotx,plotdnow,'o','markersize',2,...
    'markerfacecolor',1-0.6*(1-sescolors{5}),'MarkerEdgeColor',1-0.6*(1-sescolors{5}))
errorbar(ax,nsess+1,nanmean(plotdnow),sem(plotdnow),'-','linewidth',1.5,'CapSize',5,'color',0.5*sescolors{5});
ylim(ax,[-50 180])
g{1}=[1,5]; g{2}=[2,5]; g{3}=[3,5]; g{4}=[4,5]; g{5}=[1,2];
p(1)=signrank(runpdiff(KetSesMap(1),:)-berow(pairs1s3));
p(2)=signrank(runpdiff(KetSesMap(2),:)-berow(pairs1s3));
p(3)=signrank(runpdiff(KetSesMap(3),:)-berow(pairs1s3));
p(4)=signrank(runpdiff(KetSesMap(4),:)-berow(pairs1s3));
p(5)=signrank(runpdiff(KetSesMap(1),:)-runpdiff(KetSesMap(2),:));
myoverlapsigstar_v3(ax,g,p,'starty',160,'ygap',5,'starsize',6.5,'nsgap',0.85)
xticks(ax,1:(nsess+1)); xlbstr = sesstr; xlbstr{nsess+1} = 'REMS';
xticklabels(ax,xlbstr); ylabel(ax,'Phase lag'); xtickangle(ax,45)
ylim(ax,[-50 220]); set(ax,'YGrid','on','XGrid','off')
end

function PlotSpikePhaseExample(pos,Inpath)
global txtfonts sescolors sesstr
spklim=60; pdfbinlen=4; gband=[50,100]; stdcri=1; vcri=10; bplot=20;
parastr=['Gband',num2str(gband)];
parastrnow=[parastr,'STDCri',num2str(stdcri),'Vel',num2str(vcri)];
parastrsave=[parastrnow,'Spklim',num2str(spklim),'binlen',num2str(pdfbinlen)];
load([Inpath,'pyrphase',parastrsave,'.mat'],'pyrphase','prefpyrphase','allcluphase','intphase');
clexmple=[32,22]; dirs=[1,2]; sesplot=[1,2,5];
ttstr=sesstr(1:2); ttstr{3}='REMS';
for is = 1:length(sesplot)
    for icl = 1:length(clexmple)
        plotidx = is+(icl-1)*3;
        nowspk = allcluphase{sesplot(is),clexmple(icl),dirs(icl)};
        nowspk = wrapTo2Pi(nowspk);
        sumph  = angle(sum(exp(1i.*nowspk))); sumph = wrapTo2Pi(sumph);
        subplot('position',pos{plotidx});
        polarhistogram(nowspk,360/bplot,'FaceColor',sescolors{sesplot(is)},'FaceAlpha',.5,'EdgeColor',sescolors{sesplot(is)});
        ax = gca; hold(ax,'on'); set(ax,'fontsize',txtfonts(4)-2)
        degphase=rad2deg(nowspk); degbin=0:bplot:360; degbinc=(degbin(1:end-1)+degbin(2:end))/2;
        n=histcounts(degphase,degbin); [maxcount,degidx]=max(n); degmod=degbinc(degidx);
        polarplot(ax,[degmod,degmod]./180*pi,[0 maxcount],'-','linewidth',1.25,'color',[0.05 0.05 0.05])
        if icl==1; title(ax,ttstr{is},'fontweight','normal','fontsize',txtfonts(4)); end
        hold(ax,'off')
    end
end
end

function PlotSpikePhaseKLDiv(pos,Inpath)
global txtfonts sescolors sesstr
spklim=60; pdfbinlen=4; stdcri=1; gband=[50 100]; vcri=10;
parastr=['Gband',num2str(gband)];
parastrnow=[parastr,'STDCri',num2str(stdcri),'Vel',num2str(vcri)];
parastrsave=[parastrnow,'Spklim',num2str(spklim),'binlen',num2str(pdfbinlen)];
thetabin=30; thetaedge=0:thetabin:360; thetabc=(thetaedge(1:end-1)+thetaedge(2:end))/2;
load([Inpath,'intrunrempk',parastrsave,'.mat'],'intrunrempk','intphase','spkpdf','mintdiv');
nsess=4;
ax=subplot('position',pos{2}); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
orises=1:nsess; mapses=KetSesMap(orises); intrunrempk=log(intrunrempk);
MyPairSEMPlot_v2(intrunrempk(:,mapses),ax,1:nsess,1,'colors',sescolors,'CapSize',6,'markersize',2)
[g,p]=AllGroupSignRank(intrunrempk(:,mapses));
plotidx=[1,4,5];
myoverlapsigstar_v3(ax,g(plotidx),p(plotidx),'starty',-1.6,'ygap',0.14,'starsize',6.5,'nsgap',0.8)
ylim(ax,[-6.2 -0.7]); ylabel(ax,'Log KLDiv vs. REMS')
xticks(ax,1:nsess); xticklabels(ax,sesstr); xlim(ax,[0.2 4.8]); xtickangle(ax,45)
set(ax,'YGrid','on','XGrid','off')
ax=subplot('position',pos{1}); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
for is = 1:nsess
    sidx=KetSesMap(is); spkpnow=intphase{sidx}; spkpnow=rad2deg(spkpnow);
    [ncount]=histcounts(spkpnow,thetaedge,'Normalization','pdf');
    plot(ax,thetabc,ncount,'linewidth',1.5,'color',sescolors{is})
    hd(is)=plot(nan,nan,'-','linewidth',1,'color',sescolors{is}); lgd{is}=sesstr{is};
end
spkpnow=intphase{nsess+1}; spkpnow=rad2deg(spkpnow);
[ncount]=histcounts(spkpnow,thetaedge,'Normalization','pdf');
plot(ax,thetabc,ncount,'linewidth',1.5,'color',sescolors{nsess+1})
hd(nsess+1)=plot(nan,nan,'-','linewidth',1,'color',sescolors{nsess+1}); lgd{nsess+1}='REMS';
[lgdh,icon]=legend(hd,lgd,'fontsize',txtfonts(4)-1.5);
for it=1:5; icon((it-1)*2+6).XData(1)=icon((it-1)*2+6).XData(1)+0.45; end
set(lgdh,'Position',[0.70 0.392 0.08 0.045]); legend boxoff
xlabel(ax,'Gamma phase (deg)','HorizontalAlignment','center'); ylabel(ax,'Int spike PDF')
xlim(ax,[0 360]); xticks(ax,[0 180 360]); ylim(ax,[2e-3 4e-3])
end

function PlotLFPThetaAsymXRun(pos,Inpath)
global txtfonts sescolors sesstr
cmap1=getNCLColmap('WhiteBlueGreenYellowRed.rgb',256);
cmap2=getNCLColmap('temp_19lev.rgb',100);
load([Inpath,'bic.mat'],'bses','a1ses','a2ses','dofses','ff')
F_range=[0 150]; vels=[0 3 10 30 100]; sesplot=[1,2]; vbin=[2,4];
for is=1:length(sesplot)
    for iv=1:length(vbin)
        plotidx=is+(iv-1)*2;
        ax=subplot('position',pos{plotidx}); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
        allb=bses{sesplot(is),vbin(iv)}; alla1=a1ses{sesplot(is),vbin(iv)};
        alla2=a2ses{sesplot(is),vbin(iv)}; dofall=dofses{sesplot(is),vbin(iv)};
        bic=allb./sqrt((alla1).*(alla2));
        bicplottool(ax,bic,ff,F_range,cmap1,cmap2,dofall,3,0);
        xlim(ax,[0 15]); ylim(ax,[0 15])
        if iv==2; xlabel(ax,'Frequency (Hz)','HorizontalAlignment','center'); else; xticklabels(ax,[]); end
        if is==1; ylabel(ax,{['Vel:',num2str(vels(vbin(iv))),'-',num2str(vels(vbin(iv)+1))],'Frequency (Hz)'}); else; yticklabels(ax,[]); end
        hold(ax,'off')
    end
end
end

function PlotLFPThetaExampleXRun(pos,Inpath)
global txtfonts sescolors sesstr
cmap=getNCLColmap('MPL_BuPu.rgb',100);
load([Inpath,'thetapbDS.mat'],'thetapbDS','thetalfpDS','meshsize','actextend','cttmesh')
sesplot=[1,2]; vbin=[1,3];
for is=1:length(sesplot)
    for iv=1:length(vbin)
        plotidx=is+(iv-1)*2;
        postmp=pos{plotidx}; postmp(2)=postmp(2)+postmp(4)+0.005; postmp(4)=0.015;
        ax=subplot('position',postmp); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
        lfptmp=thetalfpDS{sesplot(is),1,1,vbin(iv)}; lfpnow=[];
        for iseg=1:length(lfptmp)
            if isempty(lfpnow); lfpnow=lfptmp{iseg}; else; lfpnow=cat(2,lfpnow,lfptmp{iseg}); end
        end
        meanlfp=nanmean(lfpnow,2);
        plot(ax,linspace(-actextend,actextend,length(meanlfp)),meanlfp*(5./max(meanlfp)),'-k','linewidth',1)
        axis tight; axis off
        if iv==1; title(ax,sesstr{sesplot(is)},'FontWeight','normal'); end
        hold(ax,'off')
    end
end
end

function PlotLFPThetaAsymREM(pos,Inpath)
global txtfonts sescolors sesstr
cmap1=getNCLColmap('WhiteBlueGreenYellowRed.rgb',256);
cmap2=getNCLColmap('temp_19lev.rgb',100);
load([Inpath,'allslpbic.mat'],'bic','dofall','ff')
F_range=[0 150];
ax=subplot('position',pos); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
bicplottool(ax,bic,ff,F_range,cmap1,cmap2,dofall,3,1);
xlim(ax,[0 15]); ylim(ax,[0 15]); xlabel(ax,'Frequency (Hz)','HorizontalAlignment','center')
hold(ax,'off')
end

function PlotLFPThetaExampleREM(pos,Inpath)
global txtfonts sescolors sesstr
thetaext=0.16;
load([Inpath,'alltlfp.mat'],'alltlfp','sestlfp');
postmp=pos; postmp(2)=postmp(2)+postmp(4)+0.005; postmp(4)=0.015;
ax=subplot('position',postmp); hold(ax,'on'); set(ax,'fontsize',txtfonts(4))
meanlfp=nanmean(alltlfp,1);
plot(ax,linspace(-thetaext,thetaext,length(meanlfp)),meanlfp*(5./max(meanlfp)),'-k','linewidth',1)
axis tight; axis off; title(ax,'REMS','FontWeight','normal')
hold(ax,'off')
end

function PlotLFPThetaAsymRunREM(pos,InPath)
global txtfonts sescolors sesstr pdfcolors
load([InPath,'asysesvel.mat'],'asysesvel','allslpasy','eachslpasy','actvel')
vbin=[0,3,10,30,100];
ax=subplot('position',pos{1}); hold(ax,'on'); set(ax,'FontSize',txtfonts(4))
vbinmean=[2,6,20,60]; vbinmean=log(vbinmean);
vbinstr=cell(size(vbinmean));
nrat=size(asysesvel,3); nsess=size(asysesvel,1); nvbin=length(vbin)-1;
sesnum=becolumn(1:nsess); vbinnum=berow(1:nvbin);
sesnum=repmat(sesnum,[1,nvbin,nrat]); vbinnum=repmat(vbinnum,[nsess,1,nrat]);
for iv=1:length(vbinmean); vbinstr{iv}=num2str(exp(vbinmean(iv))); end
for irat=1:nrat
    for is=1:nsess
        sidx=KetSesMap(is);
        plot(ax,squeeze(log(actvel(sidx,:,irat))),squeeze(asysesvel(sidx,:,irat)),...
            ':o','color',1-0.6*(1-sescolors{is}),'markerfacecolor',1-0.75*(1-sescolors{is}),'markersize',2)
    end
end
for iv=1:nvbin
    plotx=nan(1,nsess); plotd=nan(nsess,2); plotc=cell(1,nsess);
    for is=1:nsess
        sidx=KetSesMap(is);
        plotx(is)=log(nanmean(squeeze(actvel(sidx,iv,:))));
        plotd(is,1)=nanmean(squeeze(asysesvel(sidx,iv,:)));
        plotd(is,2)=nansem(squeeze(asysesvel(sidx,iv,:)));
        plotc{is}=0.9*sescolors{is};
    end
    MySEMPlot_v2(plotd,ax,plotc,plotx,0.4,'linewidth',1.2,'markersize',0)
end
sestr=KetSesStr(nsess);
for is=1:nsess
    sidx=KetSesMap(is);
    hd(is)=plot(ax,nan,nan,'s','markerfacecolor',sescolors{is},'markeredgecolor',sescolors{is},'markersize',2.5);
    dtmp=squeeze(asysesvel(sidx,:,:)); velbin=squeeze(vbinnum(sidx,:,:));
    dxx=dtmp(:); dyy=velbin(:); [R,p]=nancorr(dxx,dyy);
    stars=sigstrfromp(p); lgd{is}=[sestr{sidx},' ',stars];
end
[lgd,icon]=legend(hd,lgd,'fontsize',txtfonts(4)-1.5,'location','northwest');
for ii=1:4; icon(ii).Position(1)=icon(ii).Position(1)-0.275; end
legend boxoff
set(lgd,'Position',[0.0244 0.308 0.144 0.06],'FontSize',3.5);
ylim(ax,[-0.4 0.1]); xlabel(ax,'Velocity (cm/s)','HorizontalAlignment','center'); ylabel(ax,'Theta asymmetry')
xlim(ax,[log(0.4),log(120)]); xticks(ax,vbinmean); xticklabels(ax,vbinstr)
set(ax,'FontSize',txtfonts(4))
plotd=cell(1,nsess+1); tmpplotd=nan(20,nsess);
for is=1:nsess
    sidx=KetSesMap(is); tmpd=asysesvel(sidx,:,:); tmpd=tmpd(:);
    plotd{is}=tmpd(:); tmpplotd(1:length(tmpd),is)=tmpd;
end
plotd{nsess+1}=eachslpasy;
ax=subplot('position',pos{2}); hold(ax,'on'); set(ax,'FontSize',txtfonts(4))
MyPairSEMPlot_v2(tmpplotd,ax,1:nsess,1,'sigplot',0,'CapSize',5,'linewidth',1.5,'markersize',2,'colors',sescolors);
plotdnow=eachslpasy; postketidx=[3:3:15]; plotdnow(postketidx)=[];
plotx=nsess+1+(rand(size(plotdnow))-0.5)/10;
plot(ax,plotx,plotdnow,'o','markersize',2,...
    'markerfacecolor',1-0.6*(1-sescolors{5}),'MarkerEdgeColor',1-0.6*(1-sescolors{5}))
errorbar(ax,nsess+1,nanmean(plotdnow),sem(plotdnow),'-','linewidth',1.5,'CapSize',5,'color',0.5*sescolors{5});
ylim(ax,[-0.3 0.26])
g{1}=[1,5]; g{2}=[2,5]; g{3}=[3,5]; g{4}=[4,5]; g{5}=[1,2];
p(1)=ranksum(plotd{1},plotd{5}); p(2)=ranksum(plotd{2},plotd{5});
p(3)=ranksum(plotd{3},plotd{5}); p(4)=ranksum(plotd{4},plotd{5});
p(5)=ranksum(plotd{1},plotd{2}); 
myoverlapsigstar_v3(ax,g,p,'starty',0.08,'ygap',0.012,'starsize',6.5,'nsgap',0.8)
xticks(ax,1:(nsess+1)); xtickangle(ax,45)
xlbstr=sesstr; xlbstr{nsess+1}='REMS'; xticklabels(ax,xlbstr)
set(ax,'YGrid','on','XGrid','off')
end

function PlotSpikeRateCorr(pos,pos2,Inpath)
global txtfonts sescolors sesstr
ratecri = 0; 
parastr = ['Rate',num2str(ratecri)];
load([Inpath,'rate',parastr,'.mat'],'rateppl','ratec','rateraw')
nsess=4;
for is=1:nsess
    ax=subplot('position',pos{is}); sid=KetSesMap(is); dnow=rateraw{sid,1};
    plot2vcorr(ax,dnow,'linewidth',1,'marker','.','rpval',0,'color',sescolors{is},'markersize',1)
    px=dnow(:,1); py=dnow(:,2); goodind=~isnan(px)&~isnan(py);
    px=px(goodind); py=py(goodind); [r,p]=corr(becolumn(px),becolumn(py));
    text(ax,11,3.75,{sesstr{is},['R=',num2str(r,'%.2f')],['P=',num2str(p,'%.3f')]},'fontsize',txtfonts(4))
    xticks(ax,0:4:12); yticks(ax,0:2:10); xlim(ax,[-0.1 14]); ylim(ax,[-0.1 4])
    set(ax,'fontsize',txtfonts(4))
    if is<=2; xticklabels(ax,[]); end
    if is==3; xlb=xlabel(ax,'Run rate (Hz)','HorizontalAlignment','center'); xlb.Position(1)=xlb.Position(1)+10; end
    if rem(is,2)~=1; yticklabels(ax,[]); end
    if is==1; ylb=ylabel(ax,'pre-Run REMS rate (Hz)'); ylb.Position(2)=ylb.Position(2)-3; end
end
ax=subplot('position',pos2); set(ax,'fontsize',txtfonts(4))
plotd=cell(nsess,1);
for is=1:nsess
    sid=KetSesMap(is);
    for idir=1:2; plotd{is}=cat(1,plotd{is},becolumn(ratec{sid,idir})); end
end
MyBoxPlot_v3(plotd,ax,'sigplot',0,'color',sescolors)
xticks(ax,1:nsess); xticklabels(ax,sesstr); ylabel(ax,'Similarity vs. pre-Run REMS'); ylim(ax,[0 1.18])
g{1}=[1,2]; g{2}=[2,3]; g{3}=[2,4];
p(1)=ranksum(plotd{1},plotd{2}); p(2)=ranksum(plotd{2},plotd{3}); p(3)=ranksum(plotd{2},plotd{4});
myoverlapsigstar_v3(ax,g,p,'starty',1.03,'starsize',6.5)
set(ax,'YGrid','on','XGrid','off')
end

%% =================== shared helper functions =============================

function stars = sigstrfromp(p)
if p<=1E-3; stars='***'; elseif p<=1E-2; stars='**'; elseif p<=0.05; stars='*';
elseif isnan(p); stars='n.s.'; else; stars='n.s.'; end
end

function bicplottool(ax,bic,ff,F_range,cmap1,cmap2,dof,mod,addcb)
trg=[5 10];
ss=bic; ind=find(ff<=F_range(2)&ff>=F_range(1)); ss=ss(ind,ind); ff=ff(ind);
ss=triu(ss); ss=fliplr(ss); ss=triu(ss); ss=fliplr(ss);
switch mod; case 1; splot=abs(ss); case 2; splot=real(ss); case 3; splot=imag(ss); end
imagesc(ax,ff([1 end]),ff([1 end]),splot); set(ax,'Ydir','normal')
switch mod; case 1; caxis(ax,[0 0.4]); shading interp;
    case 2; caxis(ax,[-0.8 0.8]); shading interp;
    case 3; caxis(ax,[-0.8 0.8]); shading interp; end
plot(ax,[F_range(1),mean(F_range),F_range(2)],[F_range(1),mean(F_range),F_range(1)],'-k','linewidth',0.75)
xlim(ax,[F_range(1),F_range(2)]); ylim([F_range(1),mean(F_range)]);
switch mod; case 1; colormap(ax,cmap1); case 2; colormap(ax,cmap2); case 3; colormap(ax,cmap2); end
if addcb
    c=addcolorbar;
    switch mod; case 1; c.Label.String='Bicoherence'; case 2; c.Label.String='Skewness'; case 3; c.Label.String='Asymmetry'; end
end
plot(ax,[trg(1),trg(2),trg(2),trg(1)],[trg(1),trg(1),trg(2),trg(1)],'color',[210,160,13]/255,'linewidth',1.5)
end

function plotcof(ax,cofm,cmap,cri)
nscale=1.02; axis(ax,'square'); axis(ax,'off'); hold(ax,'on')
nc=size(cofm,1); cradius=1; diffrad=2*pi/nc;
leftrads=0:diffrad:2*pi-diffrad;
leftx=cradius.*cos(leftrads); lefty=cradius.*sin(leftrads);
for ic=1:nc
    plot(ax,nscale*leftx(ic),nscale*lefty(ic),'o',...
        'markerfacecolor',[0.3 0.3 0.3],'markeredgecolor',[0.3 0.3 0.3],'markersize',1)
end
allx=leftx; ally=lefty; allplot=[];
for ic=1:nc-1
    for jc=ic+1:nc
        cofv=cofm(ic,jc);
        if isnan(cofv); continue; end
        if abs(cofv)<cri(3); continue; end
        plotcolor=round(interp1([cri(1),cri(2)],[1,100],cofv));
        plotcolor=min([plotcolor,100]); plotcolor=max([plotcolor,1]);
        idxnow=[ic,jc,plotcolor,1]; allplot=cat(1,allplot,idxnow);
    end
end
[~,newidx]=sort(allplot(:,3),'ascend'); nowplot=allplot(newidx,:);
for iplot=1:size(nowplot,1)
    ic=nowplot(iplot,1); jc=nowplot(iplot,2); plotcolor=nowplot(iplot,3);
    plot(ax,[allx(ic),allx(jc)],[ally(ic),ally(jc)],'-','color',cmap(plotcolor,:),'linewidth',nowplot(iplot,4))
end
axis(ax,'tight')
end
