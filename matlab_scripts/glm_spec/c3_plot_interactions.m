%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition
%% GLM Spectrum analysis - Interactions
% Plot results of cluster-based test for RT main effect, RT*guided, RT*set32 interactions
% (c), Katharina Duecker
% last edited, May-03-2026

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

end_time = 0.4;

% Load main effect / interaction stats and jackknife
load(fullfile(outpth,['stat_GLM_spec_occi_sens',suf,'_',num2str(end_time*1000),'.mat']))
% Loaded: stat_all, stat_occi_RT, stat_occi_guiRT, stat_occi_set32RT, occi_grad,
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
if ~isempty(stat_occi_RT.negclusters)
    neg_cluster_pvals = [stat_occi_RT.negclusters(:).prob];
    neg_cluster_id = find(neg_cluster_pvals < 0.05);
    neg_pos = ismember(stat_occi_RT.negclusterslabelmat, neg_cluster_id);
else
    neg_pos = false(size(stat_occi_RT.mask));
end

%% Compute jackknife time/frequency bounds (used by all later panels)
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

%% Identify representative channels (group-level, from RT main effect cluster)
f_idx = logical(sum(sum(neg_pos,3)));
t_idx = find(logical(sum(sum(neg_pos,2))));

frep_idx = logical(sum(squeeze(sum(neg_pos(:,:,t_idx),3))));
f_rep = stat_occi_RT.freq(frep_idx);
chan_rep = stat_occi_RT.label(logical(sum(sum(neg_pos(:,frep_idx,t_idx),3),2)));

%% Identify subject-level SOI (from each subject's z-map within the cluster window)
load(fullfile(pth, 'matlab_scripts', 'preproc_meg','occi_sens.mat'))
occi_grad = zeros(size(SPEC.label));
for cc = 1:length(occi_soi)
    occi_grad = occi_grad + cell2mat(cellfun(@(x) ~isempty(x), regexp(SPEC.label,occi_soi{cc}),'UniformOutput',false));
end
labels = SPEC.label(logical(occi_grad));

subj_soi = cell(1,length(subj));
for s = 1:length(subj)
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'z_score_T')
    SPEC_T = SPEC;
    SPEC_T.powspctrm = squeeze(z_score_T(5,:,:,:));   % RT main effect

    cfg = [];
    cfg.channel = labels;
    cfg.frequency = [glm_freq_sig(1), glm_freq_sig(end)];
    cfg.latency = [glm_time_sig(1), glm_time_sig(end)];
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    SPEC_T = ft_selectdata(cfg, SPEC_T);
    subj_soi{s} = labels(SPEC_T.powspctrm < 0);
end
subj_soi(cell2mat(cellfun(@isempty, subj_soi, 'UniformOutput', false))) = {chan_rep};

%% Compute Cohen's f² grand average
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

cfg = [];
cfg.latency = [-1 0];
effsize_lim = ft_selectdata(cfg, effsize);

% Cohen's f² TFR averaged over each subject's SOI
cohF_subj = zeros(length(subj), length(SPEC.freq), length(SPEC.time));
for s = 1:length(subj)
    chan_idx = ismember(effsize_cohF{s}.label, subj_soi{s});
    if sum(chan_idx) > 1
        cohF_subj(s,:,:) = squeeze(mean(effsize_cohF{s}.powspctrm(chan_idx,:,:),1));
    else
        cohF_subj(s,:,:) = squeeze(effsize_cohF{s}.powspctrm(chan_idx,:,:));
    end
end
cohF_avg_TFR = squeeze(mean(cohF_subj, 1));

%% ============================================================
%% FIGURE 1: Jackknife scatter + Cohen's f²
%% ============================================================

col_time = [27, 158, 119; 217, 95, 2]./255;
col_freq = [117, 112, 179; 231, 41, 138]./255;

fig = figure('Position',[0 0 1940/1.5 1080/2]);
set(gcf, 'DefaultAxesFontSize', 11);

% --- Jackknife: time bounds ---
subplot(2,4,1)
scatter(ones(1,length(t1)), t1, t1_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_time(1,:))
hold on
scatter(ones(1,length(t2)).*1.5, t2, t2_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_time(2,:))
xlim([0.75, 1.75])
xticks([1, 1.5])
xticklabels({'onset', 'offset'})
ylabel('time (s)')
title('jackknife: time')

% --- Jackknife: frequency bounds ---
subplot(2,4,2)
scatter(ones(1,length(f1)), f1, f1_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_freq(1,:))
hold on
scatter(ones(1,length(f2)).*1.5, f2, f2_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_freq(2,:))
xlim([0.75, 1.75])
xticks([1, 1.5])
xticklabels({'lower bound', 'upper bound'})
ylabel('frequency (Hz)')
title('jackknife: frequency')

% --- Cohen's f² topo ---
ax_cohF_topo = subplot(2,4,3);
cfg = [];
cfg.parameter = 'powspctrm';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [glm_time_sig(1) glm_time_sig(end)];
cfg.ylim = [glm_freq_sig(1) glm_freq_sig(end)];
cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg, effsize_lim)
colormap(ax_cohF_topo, cm(length(cm)/2:end,:))
title('Cohen''s f² (RT)')

% --- Cohen's f² TFR (subject SOI averaged) ---
ax_cohF_tfr = subplot(2,4,4:8);
imagesc(SPEC.time, SPEC.freq, cohF_avg_TFR)
axis xy
xlim([-1 0])
xlabel('time (s)'); xticks(-1:0.5:0)
ylabel('frequency (Hz)'); yticks(10:10:30)
title('Cohen''s f² TFR (subject SOI avg)')
cb = colorbar;
cb.Label.String = 'Cohen''s f²';
box off
colormap(ax_cohF_tfr, cm(length(cm)/2:end,:))

print(fig, fullfile(plotpth, ['fig1_jackknife_cohensF', suf]), '-dsvg')
print(fig, fullfile(plotpth, ['fig1_jackknife_cohensF', suf]), '-dpng')

%% Compute projected spectra averaged over each subject's own SOI
[~, t_proj_start] = min(abs(SPEC.time - glm_time_sig(1)));
[~, t_proj_end]   = min(abs(SPEC.time - glm_time_sig(end)));

proj_min_subj = zeros(length(subj), length(SPEC.freq));
proj_max_subj = zeros(length(subj), length(SPEC.freq));

for s = 1:length(subj)
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'proj_spec_min_rt','proj_spec_max_rt')
    chan_idx = ismember(SPEC.label, subj_soi{s});

    if sum(chan_idx) > 1
        proj_min_subj(s,:) = squeeze(mean(mean(proj_spec_min_rt(chan_idx,:,t_proj_start:t_proj_end),3),1));
        proj_max_subj(s,:) = squeeze(mean(mean(proj_spec_max_rt(chan_idx,:,t_proj_start:t_proj_end),3),1));
    else
        proj_min_subj(s,:) = squeeze(mean(proj_spec_min_rt(chan_idx,:,t_proj_start:t_proj_end),3));
        proj_max_subj(s,:) = squeeze(mean(proj_spec_max_rt(chan_idx,:,t_proj_start:t_proj_end),3));
    end
end

proj_min_avg = mean(proj_min_subj, 1);
proj_max_avg = mean(proj_max_subj, 1);
proj_diff_avg = mean(proj_max_subj - proj_min_subj, 1);

%% Build outline arrays for the RT main TFR (averaged over chan_rep)
build_outline = @(mask) deal( ...
    [diff(mask); zeros(1,size(mask,2))], ...
    [zeros(size(mask,1),1), diff(mask, [], 2)]);

mask_RT = squeeze(logical(mean(stat_occi_RT.mask(ismember(stat_occi_RT.label,chan_rep),:,:))));
[mask_line_vert_RT, mask_line_horz_RT] = build_outline(mask_RT);
stat_T_RT = squeeze(mean(stat_occi_RT.stat(ismember(stat_occi_RT.label,chan_rep),:,:)));

%% Compute time-on-task TFR (regressor 6) averaged over each subject's SOI
% Also build a SPEC structure for the topoplot via grand average
ToT_subj = zeros(length(subj), length(SPEC.freq), length(SPEC.time));
ToT_freq_struct = cell(1, length(subj));
for s = 1:length(subj)
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'model_T')
    SPEC_tmp = SPEC;
    SPEC_tmp.powspctrm(:,:,:) = squeeze(model_T(6,:,:,:));   % time-on-task
    ToT_freq_struct{s} = SPEC_tmp;

    chan_idx = ismember(SPEC.label, subj_soi{s});
    if sum(chan_idx) > 1
        ToT_subj(s,:,:) = squeeze(mean(SPEC_tmp.powspctrm(chan_idx,:,:),1));
    else
        ToT_subj(s,:,:) = squeeze(SPEC_tmp.powspctrm(chan_idx,:,:));
    end
end
ToT_avg_TFR = squeeze(mean(ToT_subj, 1));

cfg = [];
cfg.parameter = 'powspctrm';
ToT_grand = ft_freqgrandaverage(cfg, ToT_freq_struct{:});

%% ============================================================
%% FIGURE 2: RT main effect | Projected spectra | Time-on-task
%% ============================================================

fig = figure('Position',[0 0 1940/1.5 1080]);
set(gcf, 'DefaultAxesFontSize', 11);

% ===== Row 1: RT main effect =====
subplot(3,4,1)
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

subplot(3,4,2:4)
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

% ===== Row 2: Projected spectra (line plots) =====
subplot(3,4,5:6)
plot(SPEC.freq, proj_min_avg, 'LineWidth', 1.5, 'Color', [8, 115, 59]./255)
hold on
plot(SPEC.freq, proj_max_avg, 'LineWidth', 1.5, 'Color', 'k')
xlim([4 30])
ylabel('magnitude')
xlabel('frequency (Hz)')
legend('fast RT', 'slow RT', 'Location', 'best')
box off
title('projected spectra (subject SOI avg)')

subplot(3,4,7:8)
plot(SPEC.freq, proj_diff_avg, 'LineWidth', 1.5, 'Color', 'k')
ylabel('slow - fast')
xlabel('frequency (Hz)')
xlim([4 30])
box off
title('difference (slow - fast)')

% ===== Row 3: Time-on-task =====
subplot(3,4,9)
cfg = [];
cfg.parameter = 'powspctrm';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [glm_time_sig(1) glm_time_sig(end)];
cfg.ylim = [glm_freq_sig(1) glm_freq_sig(end)];
cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg, ToT_grand)
colormap(cm)
title('time-on-task')

subplot(3,4,10:12)
imagesc(SPEC.time, SPEC.freq, ToT_avg_TFR)
axis xy
xlim([-1 0])
xlabel('time (s)'); xticks(-1:0.5:0)
ylabel('frequency (Hz)'); yticks(10:10:30)
title('time-on-task (T-values, subject SOI avg)')
cb = colorbar;
if max(abs(ToT_avg_TFR(:))) > 0
    caxis([-max(abs(ToT_avg_TFR(:))), max(abs(ToT_avg_TFR(:)))])
end
cb.Label.String = 'T-value';
box off
colormap(cm)

print(fig, fullfile(plotpth, ['fig2_main_rt_proj_tot', suf]), '-dsvg')
print(fig, fullfile(plotpth, ['fig2_main_rt_proj_tot', suf]), '-dpng')

%% ============================================================
%% SUPPLEMENT: Other regressors (constant, conditions, slouch, interactions)
%% ============================================================
close all;

regr_names = {'constant','guided','set size','target present','rt','tot','slouch','guided × RT','set size × RT'};
regr_idx_to_plot = [1, 2, 3, 4, 7, 8, 9];   % skip RT (5) and ToT (6) - both shown elsewhere

IAFall = cell(length(regr_names), length(subj));
for s = 1:length(subj)
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'model_T')
    for i = 1:length(regr_names)
        SPEC_tmp = SPEC;
        SPEC_tmp.powspctrm(:,:,:) = squeeze(model_T(i,:,:,:));
        IAFall{i,s} = SPEC_tmp;
    end
end

n_regr = length(regr_idx_to_plot);
n_per_row = 4;
n_rows = ceil(n_regr / n_per_row);
n_cols = n_per_row * 3;

fig = figure('Position',[0 0 1940 1080/3*n_rows]);
set(gcf, 'DefaultAxesFontSize', 10);

for k = 1:n_regr
    i = regr_idx_to_plot(k);

    cfg = [];
    cfg.parameter = 'powspctrm';
    grand_freq = ft_freqgrandaverage(cfg, IAFall{i,:});

    tfr_subj = zeros(length(subj), length(SPEC.freq), length(SPEC.time));
    for s = 1:length(subj)
        chan_idx = ismember(IAFall{i,s}.label, subj_soi{s});
        if sum(chan_idx) > 1
            tfr_subj(s,:,:) = squeeze(mean(IAFall{i,s}.powspctrm(chan_idx,:,:),1));
        else
            tfr_subj(s,:,:) = squeeze(IAFall{i,s}.powspctrm(chan_idx,:,:));
        end
    end
    tfr_avg = squeeze(mean(tfr_subj, 1));

    row = ceil(k/n_per_row);
    col_in_row = mod(k-1, n_per_row);
    base_col = col_in_row * 3 + 1;

    subplot(n_rows, n_cols, (row-1)*n_cols + base_col)
    cfg = [];
    cfg.zlim = 'maxabs';
    cfg.marker = 'off';
    cfg.layout = 'neuromag306cmb_helmet.mat';
    cfg.xlim = [glm_time_sig(1) glm_time_sig(end)];
    cfg.ylim = [f_rep(1) f_rep(end)];
    cfg.figure = 'gca';
    cfg.comment = 'no';
    ft_topoplotTFR(cfg, grand_freq)
    colormap(cm)
    title(regr_names{i})

    subplot(n_rows, n_cols, (row-1)*n_cols + (base_col+1):(row-1)*n_cols + (base_col+2))
    imagesc(SPEC.time, SPEC.freq, tfr_avg)
    axis xy
    xlabel('time (s)'); xticks(-1:0.5:0)
    ylabel('frequency (Hz)'); yticks(10:10:30)
    title('')
    cb = colorbar;
    if max(abs(tfr_avg(:))) > 0
        caxis([-max(abs(tfr_avg(:))), max(abs(tfr_avg(:)))])
    end
    box off
    colormap(cm)
end

print(fig, fullfile(plotpth, ['supp_other_reg_interactions', suf]), '-dpng')
print(fig, fullfile(plotpth, ['supp_other_reg_interactions', suf]), '-dsvg')

%% Local function for outline drawing
function draw_outline_fn(time_ax, freq_ax, m_vert, m_horz)
for i = 1:size(m_vert,1)
    for m = 2:size(m_vert,2)
        if m_vert(i,m)
            plot(time_ax(m-1:m)+0.025, ones(1,2)*freq_ax(i)+0.45, ...
                 'Color','k','Linewidth',1.5)
        end
    end
end
for i = 2:size(m_horz,1)
    for m = 2:size(m_horz,2)
        if m_horz(i,m)
            plot(time_ax(m-1)*ones(1,2)+0.025, freq_ax(i-1:i)+0.45, ...
                 'Color','k','Linewidth',1.5)
        end
    end
end
end