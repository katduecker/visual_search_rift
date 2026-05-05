%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Blanket inhibition

%% Balanced median split RIFT ~ alpha
% Plot results of cluster-based test

% (c), Katharina Duecker
% last edited, Nov-29-2024

clear all; close all; clc

set(0,'defaultAxesFontSize',16,'defaultAxesFontName','Arial')
col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;

%rmpath(genpath('/rds/projects/2018/jenseno-entrainment/fieldtrip'))
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';

addpath(fullfile(pth,'matlab_scripts','alpha'))
addpath(fullfile(pth,'raacampbell-shadedErrorBar-19cf3fe/'))
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','balanced_split');
plotpth = fullfile(pth,'results','meg','5 COH hilb','fig');
mkdir(plotpth)
load(fullfile(cohpth,'RIFT_balanced_split_glm_H1.mat'))


lwidth = 3;
timevec = linspace(-0.5,0.5,1001);
fig = figure('Position',[0 -100 1920/2.5 1080]);

subplot(321)
shadedErrorBar(timevec,squeeze(avg_cohT_high(:,2,:)),{@mean, @std},'lineProps',{'Color',col_palette(1,:),'markerfacecolor',col_palette(1,:)});
hold on
shadedErrorBar(timevec,squeeze(avg_cohT_low(:,2,:)),{@mean, @std},'lineProps',{'Color',col_palette(2,:),'markerfacecolor',col_palette(2,:)});
shadedErrorBar(timevec,squeeze(avg_cohD_high(:,2,:)),{@mean, @std},'lineProps',{'Color',col_palette(3,:),'markerfacecolor',col_palette(3,:)});
shadedErrorBar(timevec,squeeze(avg_cohD_low(:,2,:)),{@mean, @std},'lineProps',{'Color',col_palette(4,:),'markerfacecolor',col_palette(4,:)});

xlim([-0.1 0.5])
xticks([-0.5:0.5:0.5])
ylim([0.00 0.055])
yticks(0:0.025:0.05)
xlabel('time (s)')
box off
subplot(322)
ung_high = squeeze((avg_cohT_high(:,1,:) + avg_cohD_high(:,1,:))./2);
ung_low = squeeze((avg_cohT_low(:,1,:) + avg_cohD_low(:,1,:))./2);

shadedErrorBar(timevec,ung_high,{@mean, @std},'lineProps',{'Color',col_palette(5,:),'markerfacecolor',col_palette(5,:)});
hold on
shadedErrorBar(timevec,ung_low,{@mean, @std},'lineProps',{'Color',col_palette(6,:),'markerfacecolor',col_palette(6,:)});

xlim([-0.1 0.5])
xticks([-0.5:0.5:0.5])
ylim([0.00 0.055])
yticks(0:0.025:0.05)
xlabel('time (s)')
box off

subplot(323)
shadedErrorBar(timevec,squeeze(avg_cohT_high(:,4,:)),{@mean, @std},'lineProps',{'Color',col_palette(1,:),'markerfacecolor',col_palette(1,:)});
hold on
shadedErrorBar(timevec,squeeze(avg_cohT_low(:,4,:)),{@mean, @std},'lineProps',{'Color',col_palette(2,:),'markerfacecolor',col_palette(2,:)});
shadedErrorBar(timevec,squeeze(avg_cohD_high(:,4,:)),{@mean, @std},'lineProps',{'Color',col_palette(3,:),'markerfacecolor',col_palette(3,:)});
shadedErrorBar(timevec,squeeze(avg_cohD_low(:,4,:)),{@mean, @std},'lineProps',{'Color',col_palette(4,:),'markerfacecolor',col_palette(4,:)});

plot(stat.time(stat.mask),stat.mask(stat.mask) .* 0.025,'Color','k','LineWidth',6)
xlim([-0.1 0.5])
xticks([-0.5:0.5:0.5])
ylim([0.00 0.055])
yticks(0:0.025:0.05)
xlabel('time (s)')
ylabel('Rapid Frequency Tagging (coherence')
box off

subplot(324)

ung_high = squeeze((avg_cohT_high(:,3,:) + avg_cohD_high(:,3,:))./2);
ung_low = squeeze((avg_cohT_low(:,3,:) + avg_cohD_low(:,3,:))./2);

shadedErrorBar(timevec,ung_high,{@mean, @std},'lineProps',{'Color',col_palette(5,:),'markerfacecolor',col_palette(5,:)});
hold on
shadedErrorBar(timevec,ung_low,{@mean, @std},'lineProps',{'Color',col_palette(6,:),'markerfacecolor',col_palette(6,:)});

xlim([-0.1 0.5])
xticks([-0.5:0.5:0.5])
ylim([0.00 0.055])
yticks(0:0.025:0.05)
xlabel('time (s)')
box off

subplot(325)
plot(cohensd.time, cohensd.cohensd, 'Color', 'k', 'LineWidth', lwidth)
xlim([-0.1 0.5])
xticks([-0.5:0.5:0.5])
ylim([-0.5, 0.5])
yticks([-0.5,0, 0.5])
xlabel('time (s)')
ylabel('Cohen''s d')
box off

subplot(326)
plot(cohensd.time, cohensdung.cohensd, 'Color', 'k', 'LineWidth', lwidth)
xlim([-0.1 0.5])
xticks([-0.5:0.5:0.5])
ylim([-0.5, 0.5])
yticks([-0.5,0, 0.5])
xlabel('time (s)')
box off
print(fig,fullfile(plotpth,'GLM_balanced_split_H1_longtoi'),'-dsvg')
print(fig,fullfile(plotpth,'GLM_balanced_split_H1_longtoi'),'-dpng')


load(fullfile(pth,'matlab_scripts/',"preproc_meg/",'idx_subjoi.mat'));

subj_hits = zeros(length(subj),2);
for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},'balanced_split_-1000_to_400.mat'), 'hits_high', 'hits_low')
    subj_hits(s,1) = sum(hits_high)/length(hits_high);
    subj_hits(s,2) = sum(hits_low)/length(hits_low);
end

addpath(fullfile(pth,'Violinplot-Matlab-master'))

fig = figure('Position',[0 0 1920/6 1080/3]);
violinplot(subj_hits,{'high', 'low'},'ViolinColor',[[229, 44, 135]./255; [0,0,0]],'ShowMean',true);
ylabel('hit rate')

ytickformat('%.1f')
print(fig,fullfile(plotpth,'hit_rate_high_vs_low'),'-dsvg')
print(fig,fullfile(plotpth,'hit_rate_high_vs_low'),'-dpng')
[h, p, ci, stats] = ttest(subj_hits(:,1), subj_hits(:,2))