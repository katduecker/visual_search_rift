%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition
%% GLM Spectrum analysis - Interactions
% Plot results of cluster-based test for RT main effect, RT*guided, RT*set32 interactions
% (c), Katharina Duecker
% last edited, May-02-2026

clear all; close all; clc

%% paths
addpath('RainCloudPlots/tutorial_matlab/')
rmpath(genpath('/rds/projects/2018/jenseno-entrainment/fieldtrip'))

suf = '_interactions_piv';

addpath('/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';
load(fullfile(pth,'matlab_scripts/',"preproc_meg/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha','pow');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');
plotpth = fullfile(pth,'results','meg','9 GLM', 'fig','GLM_spec_interactions');
mkdir(plotpth)

% Set the end_time used in the main analysis (e.g., 0.4 -> 400)
end_time = 0.4;

% Load main effect / interaction stats and jackknife (now both at end_time*1000)
load(fullfile(outpth,['stat_GLM_spec_occi_sens',suf,'_',num2str(end_time*1000),'.mat']))
% This loads: stat_all, stat_occi_RT, stat_occi_guiRT, stat_occi_set32RT, occi_grad,
% interval_jack, frequency_jack

addpath(fullfile(pth,'matlab_scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');

%% Load template subject for SPEC structure
load(fullfile(alphapth, subj{1}, 'data_fourier_winl_10.mat'), 'TFR_avg')

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_avg);

cfg = [];
cfg.latency = [-1.5 0.5];
SPEC = ft_selectdata(cfg,TFR_alpha);

clear TFR_avg TFR_alpha
SPEC = rmfield(SPEC,'cfg');

%% Find significant negative cluster for the RT main effect
% (used to define the time-frequency window for the topoplots)
if ~isempty(stat_occi_RT.negclusters)
    neg_cluster_pvals = [stat_occi_RT.negclusters(:).prob];
    neg_cluster_id = find(neg_cluster_pvals < 0.05);
    neg_pos = ismember(stat_occi_RT.negclusterslabelmat, neg_cluster_id);
else
    neg_pos = false(size(stat_occi_RT.mask));
end

%% Jackknife plot: time and frequency bounds
[t1_count, t1] = groupcounts(interval_jack(:,1));
[t2_count, t2] = groupcounts(interval_jack(:,2));
[f1_count, f1] = groupcounts(frequency_jack(:,1));
[f2_count, f2] = groupcounts(frequency_jack(:,2));

[~, p1] = max(t1_count);
[~, p2] = max(t2_count);
glm_time_sig = [t1(p1), t2(p2)];

[~, p1] = max(f1_count);
[~, p2] = max(f2_count);
glm_freq_sig = [f1(p1), f2(p2)];

col_time = [27, 158, 119; 217, 95, 2]./255;
col_freq = [117, 112, 179; 231, 41, 138]./255;

fig = figure('Position',[0 0 1940/2.5 1080/2]);
subplot(121)
scatter(ones(1,length(t1)), t1, t1_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_time(1,:))
hold on
scatter(ones(1,length(t2)).*1.5, t2, t2_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_time(2,:))
xlim([0.75, 1.75])
xticks([1, 1.5])
% ylim([-1.1 0.05])
% yticks([-1:0.25:0])
xticklabels({'onset', 'offset'})
ylabel('time (s)')

subplot(122)
scatter(ones(1,length(f1)), f1, f1_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_freq(1,:))
hold on
scatter(ones(1,length(f2)).*1.5, f2, f2_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_freq(2,:))
xlim([0.75, 1.75])
xticks([1, 1.5])
% yticks([8:2:22])
xticklabels({'lower bound', 'upper bound'})
ylabel('frequency (Hz)')

print(fig,fullfile(plotpth,['jackknife_time_freq',suf]),'-dpng')
print(fig,fullfile(plotpth,['jackknife_time_freq',suf]),'-dsvg')
close all

%% Identify representative channels (from RT main effect cluster)
f_idx = logical(sum(sum(neg_pos,3)));
t_idx = find(logical(sum(sum(neg_pos,2))));

frep_idx = logical(sum(squeeze(sum(neg_pos(:,:,t_idx),3))));
f_rep = stat_occi_RT.freq(frep_idx);
chan_rep = stat_occi_RT.label(logical(sum(sum(neg_pos(:,frep_idx,t_idx),3),2)));

%% Compute Cohen's f² grand average from individual subject files
effsize_cohF = cell(1,length(subj));
for s = 1:length(subj)
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'CohensF_rt')
    SPEC_tmp = SPEC;
    SPEC_tmp.powspctrm(:,:,:) = CohensF_rt;
    effsize_cohF{s} = SPEC_tmp;
end

cfg = [];
cfg.parameter = 'powspctrm';
effsize = ft_freqgrandaverage(cfg, effsize_cohF{:});

% cfg = [];
% cfg.latency = [-1 0];
% effsize = ft_selectdata(cfg, effsize);

%% Helper: build outline arrays from a stat mask
% Returns mask_line_vert (horizontal-edge segments) and mask_line_horz (vertical-edge segments)
build_outline = @(mask) deal( ...
    [diff(mask); zeros(1,size(mask,2))], ...
    [zeros(size(mask,1),1), diff(mask, [], 2)]);

% RT main
mask_RT = squeeze(logical(mean(stat_occi_RT.mask(ismember(stat_occi_RT.label,chan_rep),:,:))));
[mask_line_vert_RT, mask_line_horz_RT] = build_outline(mask_RT);
stat_T_RT = squeeze(mean(stat_occi_RT.stat(ismember(stat_occi_RT.label,chan_rep),:,:)));

% guided × RT
mask_gui = squeeze(logical(mean(stat_occi_guiRT.mask(ismember(stat_occi_guiRT.label,chan_rep),:,:))));
[mask_line_vert_gui, mask_line_horz_gui] = build_outline(mask_gui);
stat_T_gui = squeeze(mean(stat_occi_guiRT.stat(ismember(stat_occi_guiRT.label,chan_rep),:,:)));

% set size × RT
mask_set = squeeze(logical(mean(stat_occi_set32RT.mask(ismember(stat_occi_set32RT.label,chan_rep),:,:))));
[mask_line_vert_set, mask_line_horz_set] = build_outline(mask_set);
stat_T_set = squeeze(mean(stat_occi_set32RT.stat(ismember(stat_occi_set32RT.label,chan_rep),:,:)));

%% Helper: draw outline on the current axes
draw_outline = @(time_ax, freq_ax, m_vert, m_horz) draw_outline_fn(time_ax, freq_ax, m_vert, m_horz);

%% Main figure: 3 effects (RT main, gui×RT, set32×RT) + Cohen's f²
% Each effect takes 1 row: [topoplot] [TFR with outline (spans 2 cols)]
% Cohen's f² is the last row.
% Grid: 4 rows × 3 cols; topo in col 1, TFR in cols 2-3.

fig = figure('Position',[0 0 1940/1.6 1080]);
set(gcf, 'DefaultAxesFontSize', 11);

% ===== Row 1: RT main effect =====
subplot(4,3,1)
cfg = [];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [glm_time_sig(1), glm_time_sig(end)];
cfg.ylim = [glm_freq_sig(1) glm_freq_sig(end)];
cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg, stat_occi_RT)
colormap(cm)
title('RT main')

subplot(4,3,2:3)
imagesc(stat_occi_RT.time, stat_occi_RT.freq, stat_T_RT)
hold on
draw_outline_fn(stat_occi_RT.time, stat_occi_RT.freq, mask_line_vert_RT, mask_line_horz_RT)
axis xy
xlabel('time (s)'); xticks(-1:0.5:0)
ylabel('frequency (Hz)'); yticks(10:10:30)
title('RT main (T-values)')
cb = colorbar;
caxis([-max(abs(stat_T_RT(:))), max(abs(stat_T_RT(:)))])
cb.Label.String = 'T-value';
box off
colormap(cm)

% ===== Row 2: guided × RT =====
subplot(4,3,4)
cfg = [];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [glm_time_sig(1), glm_time_sig(end)];
cfg.ylim = [glm_freq_sig(1) glm_freq_sig(end)];
cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg, stat_occi_guiRT)
colormap(cm)
title('guided × RT')

subplot(4,3,5:6)
imagesc(stat_occi_guiRT.time, stat_occi_guiRT.freq, stat_T_gui)
hold on
draw_outline_fn(stat_occi_guiRT.time, stat_occi_guiRT.freq, mask_line_vert_gui, mask_line_horz_gui)
axis xy
xlabel('time (s)'); xticks(-1:0.5:0)
ylabel('frequency (Hz)'); yticks(10:10:30)
title('guided × RT (T-values)')
cb = colorbar;
if max(abs(stat_T_gui(:))) > 0
    caxis([-max(abs(stat_T_gui(:))), max(abs(stat_T_gui(:)))])
end
cb.Label.String = 'T-value';
box off
colormap(cm)

% ===== Row 3: set size × RT =====
subplot(4,3,7)
cfg = [];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [glm_time_sig(1), glm_time_sig(end)];
cfg.ylim = [glm_freq_sig(1) glm_freq_sig(end)];
cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg, stat_occi_set32RT)
colormap(cm)
title('set size × RT')

subplot(4,3,8:9)
imagesc(stat_occi_set32RT.time, stat_occi_set32RT.freq, stat_T_set)
hold on
draw_outline_fn(stat_occi_set32RT.time, stat_occi_set32RT.freq, mask_line_vert_set, mask_line_horz_set)
axis xy
xlabel('time (s)'); xticks(-1:0.5:0)
ylabel('frequency (Hz)'); yticks(10:10:30)
title('set size × RT (T-values)')
cb = colorbar;
if max(abs(stat_T_set(:))) > 0
    caxis([-max(abs(stat_T_set(:))), max(abs(stat_T_set(:)))])
end
cb.Label.String = 'T-value';
box off
colormap(cm)

% ===== Row 4: Cohen's f² =====
ax4 = subplot(4,3,10);q
cfg = [];
cfg.parameter = 'powspctrm';
%cfg.zlim = [0, 2e-3];
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [glm_time_sig(1) glm_time_sig(end)];
cfg.ylim = [glm_freq_sig(1) glm_freq_sig(end)];
cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg, effsize)
colormap(ax4, cm(length(cm)/2:end,:))
% cb = colorbar;
% title('Cohen''s f² (RT contribution)')

ax4 = subplot(4,3,11:12);
imagesc(effsize.time, effsize.freq, ...
        squeeze(mean(effsize.powspctrm(ismember(effsize.label, chan_rep),:,:),1)))
axis xy
xlabel('time (s)'); xticks(-1:0.5:0)
ylabel('frequency (Hz)'); yticks(10:10:30)
title('Cohen''s f² (TFR)')
cb = colorbar;
%caxis([0, 2e-3])
cb.Label.String = 'Cohen''s f²';
box off
colormap(ax4, cm(length(cm)/2:end,:))


print(fig, fullfile(plotpth, ['main_fig_interactions', suf]), '-dsvg')
print(fig, fullfile(plotpth, ['main_fig_interactions', suf]), '-dpng')

%% Local function for outline drawing
function draw_outline_fn(time_ax, freq_ax, m_vert, m_horz)
% Horizontal segments (transitions in frequency direction)
for i = 1:size(m_vert,1)
    for m = 2:size(m_vert,2)
        if m_vert(i,m)
            plot(time_ax(m-1:m)+0.025, ones(1,2)*freq_ax(i)+0.45, ...
                 'Color','k','Linewidth',1.5)
        end
    end
end
% Vertical segments (transitions in time direction)
for i = 2:size(m_horz,1)
    for m = 2:size(m_horz,2)
        if m_horz(i,m)
            plot(time_ax(m-1)*ones(1,2)+0.025, freq_ax(i-1:i)+0.45, ...
                 'Color','k','Linewidth',1.5)
        end
    end
end
end