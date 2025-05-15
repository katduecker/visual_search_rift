%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Blanket inhibition

%% GLM Spectrum analysis
% Plot results of cluster-based test

% (c), Katharina Duecker
% last edited, Nov-29-2024


clear all; close all; clc

%% paths

addpath('RainCloudPlots/tutorial_matlab/')

suf = '_fourier_noz';
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

alphapth = fullfile(pth,'results','meg','6 Alpha','pow');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');
plotpth = fullfile(pth,'results','meg','9 GLM', 'fig','GLM_spec_results');
mkdir(plotpth)
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

load(fullfile(outpth,['stat_GLM_spec_occi_sens',suf,'.mat']))

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

%% load template subject

load(fullfile(alphapth, subj{1}, 'data_fourier_winl_10.mat'), 'TFR_avg')

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_avg);

cfg = [];
cfg.latency = [-1.5 0.5];
SPEC = ft_selectdata(cfg,TFR_alpha);

%
clear TFR_avg TFR_alpha


SPEC = rmfield(SPEC,'cfg');

%% Find negative clusters

if ~isempty(stat_occi.negclusters)
    neg_cluster_pvals = [stat_occi.negclusters(:).prob];
    neg_cluster_id = find(neg_cluster_pvals < 0.05);
    neg_pos     = ismember(stat_occi.negclusterslabelmat, neg_cluster_id);
    neg_cluster = sum(neg_pos(:));
else
    neg_cluster = 0;
end

fig = figure('Position',[0 0 1940/1.5 1080/2]);
subplot(221)
h1 = histogram(interval_jack(:,1), 10);
t1_bins = h1.BinEdges;
t1_counts = h1.BinCounts;
ylabel('count')
xlabel('onset')
xlim([-1 0])
subplot(222)
h2 = histogram(interval_jack(:,2), 10);
xlabel('offset')
xlim([-1 0])
t2_bins = h2.BinEdges;
t2_counts = h2.BinCounts;
subplot(223)
h3 = histogram(frequency_jack(:,1), 10);
ylabel('count')
xlabel('min frequency')
f1_bins = h3.BinEdges;
f1_counts = h3.BinCounts;
subplot(224)
h4 = histogram(frequency_jack(:,2), 10);
xlabel('max frequency')
f2_bins = h4.BinEdges;
f2_counts = h4.BinCounts;

[t1_count, t1] = groupcounts(interval_jack(:,1));
[t2_count, t2] = groupcounts(interval_jack(:,2));

[f1_count, f1] = groupcounts(frequency_jack(:,1));
[f2_count, f2] = groupcounts(frequency_jack(:,2));

col_time = [27, 158, 119; 217, 95, 2]./255;
col_freq = [117, 112, 179; 231, 41, 138]./255;
fig = figure('Position',[0 0 1940/2.5 1080/2]);
subplot(121)
scatter(ones(1,length(t1)), t1, t1_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_time(1,:))
hold on
scatter(ones(1,length(t2)).*1.5, t2, t2_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_time(2,:))
xlim([0.75, 1.75])
xticks([1, 1.5])
ylim([-1 0.05])
yticks([-1:0.25:0])
xticklabels({'onset', 'offset'})
ylabel ('time (s)')

subplot(122)
scatter(ones(1,length(f1)), f1, f1_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_freq(1,:))
hold on
scatter(ones(1,length(f2)).*1.5, f2, f2_count.*30, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerFaceColor', col_freq(2,:))
xlim([0.75, 1.75])
xticks([1, 1.5])
yticks([8:2:22])
xticklabels({'lower bound', 'upper bound'})
ylabel('frequency (Hz)')

print(fig,fullfile(plotpth,['jackknife_time_freq',suf]),'-dpng')
print(fig,fullfile(plotpth,['jackknife_time_freq',suf]),'-dsvg')

close all

%% Plot 1: Significant channels over time points

fig = figure('Position',[0 0 1940 1080]);
f_idx = logical(sum(sum(neg_pos,3)));
freq = stat_occi.freq(f_idx);

t_idx = find(logical(sum(sum(neg_pos,2))));
time_oi = stat_occi.time(t_idx);
time_oi_all = [-0.9 0];
% size of subplot grid
c = 5;
r = ceil(length(t_idx)/c);

h = 0;
for t = t_idx'
    h = h + 1;
    subplot(r,c,h)
    cfg = [];
    cfg.parameter = 'stat';
    cfg.zlim = [-2, 2];
    cfg.highlight = 'on';
    cfg.marker = 'off';
    cfg.layout = 'neuromag306cmb_helmet.mat';
    cfg.xlim = [stat_occi.time(t) stat_occi.time(t)];
    cfg.ylim = [freq(1) freq(end)];
    cfg.highlightchannel = stat_occi.label(logical(sum(neg_pos(:,f_idx,t),2)));
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    cfg.zlim = 'maxabs';
    cfg.figure = 'gca';
    ft_topoplotTFR(cfg,stat_occi)
    colormap(cm)
end

print(fig,fullfile(plotpth,['topo_timepoints',suf]),'-dpng')
print(fig,fullfile(plotpth,['topo_timepoints',suf]),'-dsvg')



frep_idx = logical(sum(squeeze(sum(neg_pos(:,:,t_idx),3))));
f_rep = stat_occi.freq(frep_idx);
chan_rep = stat_occi.label(logical(sum(sum(neg_pos(:,frep_idx,t_idx),3),2)));


% find the combined planars belonging to the occipital sensors
% load occipital sensors
load(fullfile(pth, 'matlab scripts', 'preprocessing MEG','occi_sens.mat'))
occi_grad = zeros(size(SPEC.label));

for c = 1:length(occi_soi)
    occi_grad = occi_grad + cell2mat(cellfun(@(x) ~isempty(x), regexp(SPEC.label,occi_soi{c}),'UniformOutput',false));
end

%% Get fitted GLM for each participant
close all
% get the occipital sensors that show the effect
labels = SPEC.label(logical(occi_grad));
subj_soi = cell(1,length(subj));

fig = figure('Position',[0 0 1940/2 1080]);

for s = 1:length(subj)
    
   
    % identify sensors of interest
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'model_T', 'T_perm')

    SPEC_T = SPEC;
    SPEC_T.powspctrm = squeeze(model_T(5,:,:,:));
    SPEC_perm = SPEC;
    SPEC_perm.powspctrm = T_perm;
    SPEC_perm.dimord = 'subj_chan_freq_time';
    SPEC_perm_std = SPEC_perm;

    % get subject SOI
    cfg = [];
    cfg.channel = labels;
    cfg.frequency = [f_rep(1), f_rep(end)];
    cfg.latency = [time_oi(1), time_oi(end)];
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    SPEC_T = ft_selectdata(cfg, SPEC_T);
    SPEC_perm = ft_selectdata(cfg, SPEC_perm); 
    SPEC_perm_std = ft_selectdata(cfg, SPEC_perm_std); 
  
    subj_soi{s} = labels((SPEC_T.powspctrm - squeeze(mean(SPEC_perm.powspctrm))')./squeeze(std(SPEC_perm_std.powspctrm))' < -1.96);

    subplot(8,4,s)
    cfg = [];
    cfg.zlim = 'maxabs';
    cfg.highlight = 'on';
    cfg.marker = 'off';
    cfg.layout = 'neuromag306cmb_helmet.mat';

    cfg.highlightchannel = subj_soi{s};
    cfg.zlim = 'maxabs';
    cfg.figure = 'gca';
    cfg.comment = 'no';
    ft_topoplotTFR(cfg,SPEC_T)
    colormap(cm)
    title('')
end

print(fig,fullfile(plotpth,'rt_regr_subj_zsmaller-196'),'-dsvg')

save(fullfile(outpth,['glm_rt_chan',suf,'.mat']),'chan_rep', 'f_rep', 'subj_soi', 'time_oi', 'time_oi_all')

%% Main figure: average over time points + show TFR of RT regressor

stat_T = mean(stat_occi.stat(ismember(stat_occi.label,chan_rep),:,:));
mask = squeeze(logical(mean(stat_occi.mask(ismember(stat_occi.label,chan_rep),:,:))));

mask_line_vert = diff(mask);
mask_line_vert = [mask_line_vert;zeros(1,size(mask_line_vert,2))];

mask_line_horz = diff(mask, [], 2);
mask_line_horz = [zeros(size(mask_line_horz,1),1),mask_line_horz];

fig = figure('Position',[0 0 1940/2 1080/3]);
subplot(141)
cfg = [];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [time_oi(1), time_oi(end)];
cfg.ylim = [f_rep(1) f_rep(end)];
cfg.highlightchannel = chan_rep;
cfg.zlim = 'maxabs';
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg,stat_occi)
colormap(cm)
title('')


subplot(1,4,2:4)
imagesc(stat_occi.time,stat_occi.freq,squeeze(stat_T))
hold on

for i = 1:size(mask_line_vert,1)
    for m = 1:size(mask_line_vert,2)
        if mask_line_vert(i,m)
            plot(stat_occi.time(m-1:m)+0.025, ones(1,2)*stat_occi.freq(i)+0.45, 'Color', 'k', 'Linewidth',1.5)
        end
    end
end

for i = 1:size(mask_line_horz,1)
    for m = 1:size(mask_line_horz,2)
        if mask_line_horz(i,m)
            plot(stat_occi.time(m-1)*ones(1,2)+0.025, stat_occi.freq(i-1:i)+0.45, 'Color', 'k', 'Linewidth',1.5)
        end
    end
end

axis xy
xlabel('time (s)')
xticks(-1:0.5:0)
yticks(10:10:30)
ylabel('frequency (Hz)')
title('')
cb = colorbar;
cb.Limits = [-max(abs(stat_T(:))) max(abs(stat_T(:)))];
caxis(round(cb.Limits,1))
cb.Ticks = cb.Limits(1):cb.Limits(end):cb.Limits(end);
cb.Label.String = 'T-value';
box off

print(fig,fullfile(plotpth,['main_fig_rt_regressor_outline',suf]),'-dsvg')


%% Plot effect size
% load regressors

effsize_cohF = cell(1,length(subj));

fig = figure('Position',[0 0 1940 1080]);
cfg = [];
cfg.parameter = 'powspctrm';
cfg.zlim = 'zeromax';
%cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [time_oi(1) time_oi(end)];
cfg.ylim = [f_rep(1) f_rep(end)];
%cfg.highlightchannel = chan_rep;
cfg.figure = 'gca';
cfg.comment = 'no';
for s = 1:length(subj)

    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'CohensF_rt')

    SPEC.powspctrm(:,:,:) = CohensF_rt;

    effsize_cohF{s} = SPEC;
    subplot(4,8,s)
    
    %cfg.zlim = 'zeromax';
    
    ft_topoplotTFR(cfg,effsize_cohF{s})
    colormap(cm)
    colorbar
    title('')

end

effsize = ft_freqgrandaverage([],effsize_cohF{:});

cfg = [];
cfg.latency = [-1 0];
effsize = ft_selectdata(cfg, effsize);

% plot average effect size

fig = figure('Position',[0 0 1940/2 1080/3]);
subplot(141)
cfg = [];
cfg.parameter = 'powspctrm';
cfg.zlim = [0, 2e-3];
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [time_oi(1) time_oi(end)];
cfg.ylim = [f_rep(1) f_rep(end)];
cfg.highlightchannel = chan_rep;
%cfg.zlim = 'zeromax';
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg,effsize)
colormap(cm)
cb = colorbar;
title('')


subplot(1,4,2:4)
imagesc(effsize.time,effsize.freq,squeeze(mean(effsize.powspctrm(ismember(effsize.label, chan_rep),:,:),1)))

hold on

for i = 1:size(mask_line_vert,1)
    for m = 1:size(mask_line_vert,2)
        if mask_line_vert(i,m)
            plot(stat_occi.time(m-1:m)+0.025, ones(1,2)*stat_occi.freq(i)+0.45, 'Color', 'k', 'Linewidth',1.5)
        end
    end
end

for i = 1:size(mask_line_horz,1)
    for m = 1:size(mask_line_horz,2)
        if mask_line_horz(i,m)
            plot(stat_occi.time(m-1)*ones(1,2)+0.025, stat_occi.freq(i-1:i)+0.45, 'Color', 'k', 'Linewidth',1.5)
        end
    end
end

axis xy
xlabel('time (s)')
xticks(-1:0.5:0)
yticks(10:10:30)
ylabel('frequency (Hz)')
title('')

cb = colorbar;

caxis([0, 2e-3])
cb.Ticks = [0,2e-3];
cb.Label.String = 'Cohen''s F';
box off

print(fig,fullfile(plotpth,['cohensF_rt_regr_outline',suf]),'-dsvg')

%% plot projected spectra
close all;
% load projected spectra
GA_proj_spc_min = zeros(length(subj),size(neg_pos,1),size(neg_pos,2),length(SPEC.time));
GA_proj_spc_max = zeros(length(subj),size(neg_pos,1),size(neg_pos,2),length(SPEC.time));

fig = figure('Position',[0 0 1940/2 1080]);

[~, t1] = min(abs(SPEC.time+0.3));
[~, t2] = min(abs(SPEC.time));


f_peak = zeros(length(subj),1);
for s = 1:length(subj)

    chan_idx = ismember(SPEC.label, subj_soi{s});
    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'proj_spec_min_rt','proj_spec_max_rt')
    
    bsl = mean(proj_spec_max_rt(:,:,t1:t2),3);

    GA_proj_spc_min(s,:,:,:) = proj_spec_min_rt./max(bsl,[],2);
    GA_proj_spc_max(s,:,:,:) = proj_spec_max_rt./max(bsl,[],2);
    
    if sum(chan_idx) > 1

        proj_slow = mean(mean(proj_spec_max_rt(chan_idx,:,t1:t2),3));
        proj_fast = mean(mean(proj_spec_min_rt(chan_idx,:,t1:t2),3));
    else
        proj_slow = mean(proj_spec_max_rt(chan_idx,:,t1:t2),3);
        proj_fast = mean(proj_spec_min_rt(chan_idx,:,t1:t2),3);
    end
    
    [peaks, locs] = findpeaks(proj_fast);
    [~, p] = max(peaks);
    
    
    f_peak(s) = SPEC.freq(locs(p));
    subplot(8,4,s)
    plot(SPEC.freq,proj_slow.*10e11, 'Color', 'k')
    hold on
    plot(SPEC.freq,proj_fast.*10e11,'Color',[8, 115,59]./255)
    xlim([4,30])
    %plot(SPEC.freq(locs(p)),peaks(p),'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k', 'MarkerSize', 5)
end
legend('slow', 'fast')
print(fig,fullfile(plotpth,['rt_glm_proj_spec_all_subj',suf]),'-dsvg')

save(fullfile(outpth,['glm_rt_chan',suf,'.mat']),'f_peak', '-append')

close all
diff_proj_spc = GA_proj_spc_max./GA_proj_spc_min-1;

% average
PROJ_min = SPEC;
PROJ_min.powspctrm = squeeze(mean(GA_proj_spc_min));
PROJ_max = SPEC;
PROJ_max.powspctrm = squeeze(mean(GA_proj_spc_max));
PROJ_diff = SPEC;
PROJ_diff.powspctrm = squeeze(mean(diff_proj_spc));
PROJ_diff_std = SPEC;
PROJ_diff_std.powpsctrm = squeeze(std(diff_proj_spc));

close all;
fig = figure('Position',[0 0 1940 1080/3.5]);
set(gcf, 'DefaultAxesFontSize', 12); % Set default font size for axes in the current figure
subplot(171)
cfg = [];
cfg.parameter = 'powspctrm';
cfg.zlim = 'maxabs';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [time_oi(1) time_oi(end)];
cfg.ylim = [f_rep(1) f_rep(end)];
cfg.zlim = 'maxmin';
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg,PROJ_diff)
colormap(cm)
colorbar
title('relative difference fast vs slow')


subplot(1,7,2:3)
cfg = [];
cfg.parameter = 'powspctrm';
cfg.zlim = 'maxabs';
cfg.channel = chan_rep;
cfg.figure = 'gca';
cfg.layout = 'neuromag306cmb_helmet.mat';
ft_singleplotTFR(cfg,PROJ_diff)

xlabel('time (s)')
xticks(-1:0.5:0)
yticks(10:10:30)
ylabel('frequency (Hz)')
title('')
cb = colorbar;
cb.Ticks = cb.Limits(1):cb.Limits(end):cb.Limits(end);
cb.Label.String = 'difference log power slow - fast';
box off


cfg = [];
cfg.latency = [time_oi(1) time_oi(end)];
cfg.channel = chan_rep;
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
SPEC_min = ft_selectdata(cfg,PROJ_min);
SPEC_max = ft_selectdata(cfg,PROJ_max);
SPEC_diff = ft_selectdata(cfg,PROJ_diff);
SPEC_diff_std = ft_selectdata(cfg,PROJ_diff_std);

subplot(1,7,4:5)
plot(SPEC_min.freq,SPEC_min.powspctrm,'LineWidth',2, 'Color', [8, 115,59]./255)
hold on
plot(SPEC_max.freq,SPEC_max.powspctrm,'LineWidth',2, 'Color', 'k')
xlim([4 30])

ylabel('magnitude')
xlabel('frequency (Hz)')
legend('fast RT', 'slow RT')
box off

subplot(1,7,6:7)
plot(SPEC_diff.freq,SPEC_diff.powspctrm,'LineWidth',2,'Color','k')
ylabel('difference log power slow - fast')
xlabel('frequency (Hz)')
xlim([4 30])
% ylim([-1*10^(-14) 1*10^(-14)])
% yticks(-1*10^(-14):1*10^(-14):1*10^(-14))
box off

print(fig,fullfile(plotpth,['rt_glm_proj_spec',suf]),'-dsvg')


fig = figure('Position',[0 0 1940/2 1080/3.5]);
% projected spectra
subplot(1,2,1)
plot(SPEC_min.freq,SPEC_min.powspctrm,'LineWidth',1.5, 'Color', [8, 115,59]./255)
hold on
plot(SPEC_max.freq,SPEC_max.powspctrm,'LineWidth',1.5, 'Color', 'k')
xlim([4 30])

ylabel('magnitude')
xlabel('frequency (Hz)')
legend('fast RT', 'slow RT')
box off

subplot(1,2,2)
plot(SPEC_diff.freq,SPEC_diff.powspctrm,'LineWidth',1.5,'Color','k')
ylabel('difference log power slow - fast')
xlabel('frequency (Hz)')
xlim([4 30])
% ylim([-1*10^(-14) 1*10^(-14)])
% yticks(-1*10^(-14):1*10^(-14):1*10^(-14))
box off
print(fig,fullfile(plotpth,['rt_glm_proj_spec_small',suf]),'-dsvg')


%% Supplement plot: Other regressors in the model
close all;
h = 1;

regr_names = {'constant','guided','set size','target present','rt','tot','slouch'};

IAFall = cell(length(regr_names),length(subj));

% load regressors
for s = 1:length(subj)

    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'model_T')

    for i = 1:length(regr_names)

        % get T-value
        SPEC.powspctrm(:,:,:) = squeeze(model_T(i,:,:,:));

        IAFall{i,s} = SPEC;

    end

end

fig = figure('Position',[0 0 1940 1080/3]);

for i = 1:length(regr_names)

    if strcmp(regr_names{i},'rt')
        continue
    end
    cfg = [];
    cfg.parameter = 'powspctrm';

    grand_freq = ft_freqgrandaverage(cfg,IAFall{i,:});


    subplot(2,9,h)
    cfg = [];
    cfg.zlim = 'maxabs';
    cfg.marker = 'off';
    cfg.layout = 'neuromag306cmb_helmet.mat';
    cfg.xlim = [time_oi(1) time_oi(end)];
    cfg.ylim = [f_rep(1) f_rep(end)];
    cfg.zlim = 'maxabs';
    cfg.figure = 'gca';
    cfg.comment = 'no';
    ft_topoplotTFR(cfg,grand_freq)
    colormap(cm)
    title(regr_names{i})

    h = h+1;
    subplot(2,9,h:h+1)
    cfg = [];
    cfg.zlim = 'maxabs';
    cfg.channel = chan_rep;
    cfg.figure = 'gca';
    cfg.layout = 'neuromag306cmb_helmet.mat';
    ft_singleplotTFR(cfg,grand_freq)
    xlabel('time (s)')
    xticks(-1:0.5:0)
    yticks(10:10:30)
    ylabel('frequency (Hz)')
    title('')
    cb = colorbar;
    caxis(round(cb.Limits,1))
    %cb.Limits = round(cb.Limits,11);
    cb.Ticks = cb.Limits(1):cb.Limits(end):cb.Limits(end);
    %cb.Label.String = 'T value';
    box off

    h = h+2;
end

print(fig,fullfile(plotpth,['supp_other_reg_fourier',suf]),'-dpng')
print(fig,fullfile(plotpth,['supp_other_reg_fourier',suf]),'-dsvg')

