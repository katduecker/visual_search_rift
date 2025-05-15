%% Blanket inhibition test


clear all; close all; clc;

which_set = 'gui';

rmpath(genpath('/rds/projects/2018/jenseno-entrainment/fieldtrip'))
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT')
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
plotpth = fullfile(pth,'results','meg','9 GLM', 'fig','GLM_rift_results');
mkdir(plotpth)
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_rift');

cohpth = fullfile(pth,'results','meg','8 COH single trl');

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

% template
load(fullfile(cohpth,subj{1},'coh_single_trial_pearson.mat'),'corrT')
corrT.avg = repmat(corrT.avg,1,1,2);
corrT.time = 1:2;

cfg = [];
cfg.method = 'sum';
corrT = ft_combineplanar(cfg,corrT);

cfg = [];
%cfg.latency = [1 1];
cfg.avgoverrpt = 'yes';
corrT = ft_selectdata(cfg,corrT);

corrT.avg = corrT.trial;
corrT = rmfield(corrT,'trial');

glm_subj_alpha = cell(length(subj),1);
glm_subj_alpha_tot = cell(length(subj),1);

soi_alpha_tot = cell(1,length(subj));
fig = figure('Position',[0 0 1980 1080]);
cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.parameter = 'avg';
cfg.zlim = 'maxabs';
cfg.marker = 'off';
cfg.xlim = [1 1];
cfg.highlight = 'on';
cfg.comment = 'no';
cfg.highlightsize = 15;
for s = 1:length(subj)
    load(fullfile(outpth,subj{s},append('glm_coh_distractor',which_set,'.mat')))
    
    glm_alpha = corrT; 
    glm_alpha.avg(:,1) = T_alpha_z;
    glm_alpha_tot = corrT;
    glm_alpha_tot.avg(:,1) = T_alpha_tot_z;
    glm_subj_alpha{s} = glm_alpha;
    glm_subj_alpha_tot{s} = glm_alpha_tot;
    
    soi_alpha_tot{s} = glm_alpha_tot.label(T_alpha_tot_z < 0);
    cfg.highlightchannel = soi_alpha_tot{s};
    
    cfg.figure = 'gca';
    cfg.zlim = [-2 2];
    subplot(4,8,s)
    ft_topoplotER(cfg,glm_alpha_tot);
    colormap(cm)
    glm_subj_alpha{s} = glm_alpha;
    glm_subj_alpha_tot{s} = glm_alpha_tot;
    
    clear glm_alpha glm_alpha_tot
end

print(fig,fullfile(plotpth,['distractor_glm_single_subj_',which_set]),'-dsvg')
print(fig,fullfile(plotpth,['distractor_glm_single_subj_', which_set]),'-dpng')

cfg = [];
cfg.keepindividual = 'yes';
grand_alpha = ft_timelockgrandaverage(cfg,glm_subj_alpha{:});
grand_alpha_tot = ft_timelockgrandaverage(cfg,glm_subj_alpha_tot{:});

null_hyp = grand_alpha;
null_hyp.individual = zeros(size(null_hyp.individual));

cfg = [];
cfg.feedback = 'no';
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';
neighbours = ft_prepare_neighbours(cfg);

% load occipital sensors
load(fullfile(pth, 'matlab scripts', 'preprocessing MEG','occi_sens.mat'))

cfg = [];
cfg.neighbours       = neighbours;                 % fieldtrip template                    
cfg.method           = 'montecarlo';
cfg.minnbchan        = 1;
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.tail             = -1;
cfg.clustertail      = -1;
cfg.numrandomization = 5000;
cfg.alpha            = 0.05;
cfg.clusteralpha = 0.05;
cfg.latency = [1 1];
design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

stat_alpha = ft_timelockstatistics(cfg,grand_alpha,null_hyp);
stat_alpha_tot = ft_timelockstatistics(cfg,grand_alpha_tot,null_hyp);

soi_alpha_tot = stat_alpha.label(stat_alpha.mask);
fig = figure('Position',[0 0 600 200]);
cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.parameter = 'avg';
cfg.zlim = 'maxabs';
cfg.marker = 'off';
cfg.xlim = [1 1];
cfg.highlight = 'on';
cfg.parameter = 'stat';
cfg.comment = 'no';
cfg.highlightchannel = stat_alpha.label(stat_alpha.mask);
cfg.highlightsize = 15;


cfg.figure = 'gca';
cfg.zlim = [-3 3];
subplot(121)
ft_topoplotER(cfg,stat_alpha);
colormap(cm)
cb = colorbar;
cb.Ticks = -3:3:3;
title('RIFTavg - alpha');
subplot(122)
cfg.highlightchannel = stat_alpha_tot.label(stat_alpha_tot.mask);

ft_topoplotER(cfg,stat_alpha_tot);
colormap(cm)
cb = colorbar;
cb.Ticks = -3:3:3;
cb.Label.String = 'zscore';
title('RIFTavg - alpha + tot');

print(fig,fullfile(plotpth,['distractor_glm_coh_MEGGRAD_', which_set]),'-dsvg')
print(fig,fullfile(plotpth,['distractor_glm_coh_MEGGRAD_', which_set]),'-dpng')
