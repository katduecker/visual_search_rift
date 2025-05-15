%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 1: Guided Visual Search is associated with a feature-based priority
% map in early visual cortex

%% GLM Spectrum analysis
% Target boosting and distractor suppression using a single-trial measure
% of coherence based on correlations

% (c), Katharina Duecker
% last edited, Nov-29-2024

% Cluster based permutation test on GLM factors

clear all; close all; clc;

z_corr = 0;
%%
% path settings
rmpath(genpath('/rds/projects/2018/jenseno-entrainment/'))
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT')
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

outpth = fullfile(pth,'results','meg','9 GLM', 'glm spec quinn et al','T boost D supp');

% colormap
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);

set32 = 0;
which_corr = 'pearson';

%% Load individual GLMs
glm_subj = cell(length(subj),1);

for s = 1:length(subj)
    
    if set32
        load(fullfile(outpth,subj{s},['glm_TvsD_set32_',which_corr,'_firws.mat']),'coh_RIFT')
    else
        load(fullfile(outpth,subj{s},['glm_TvsD_',which_corr,'_firws.mat']),'coh_RIFT')
    end

    glm_subj{s} = coh_RIFT;
end

cfg = [];
cfg.keepindividual = 'yes';
grand_coh = ft_timelockgrandaverage(cfg,glm_subj{:});

% null distribution
null_hyp = grand_coh;
null_hyp.individual = zeros(size(null_hyp.individual));  


clear glm_subj*


%% Cluster-based test
cfg = [];
cfg.feedback = 'no';
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';
neighbours = ft_prepare_neighbours(cfg,grand_coh);

cfg = [];
cfg.neighbours       = neighbours;                 % fieldtrip template                    
cfg.method           = 'montecarlo';
cfg.minnbchan        = 1;
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.numrandomization = 2000;
cfg.alpha            = 0.05;
cfg.clusteralpha = 0.05;
cfg.latency = [1 1];
design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

% Target boosting
statT = ft_timelockstatistics(cfg,grand_coh,null_hyp);


% Disttractor suppression
cfg.latency = [2 2];
cfg.tail             = -1;
cfg.clustertail      = -1;
statD = ft_timelockstatistics(cfg,grand_coh,null_hyp);

if z_corr
    save(fullfile(outpth,['tboost_dsuppr_cluster_',which_corr,'_zcorr.mat']),'statT','statD')
else
    save(fullfile(outpth,['tboost_dsuppr_cluster_',which_corr,'.mat']),'statT','statD')
end

%% Plot
fig = figure;
cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.parameter = 'avg';
cfg.zlim = 'maxabs';
cfg.marker = 'off';
cfg.highlight = 'on';
cfg.parameter = 'stat';
cfg.highlightchannel = statT.label(statT.mask(:,1));
cfg.highlightsize = 15;


cfg.comment = 'Target';
cfg.commentpos = 'title';
cfg.figure = 'gca';
cfg.zlim = [-2 2];
subplot(121)
ft_topoplotER(cfg,statT);
colormap(cm)
colorbar

subplot(122)
cfg.highlightchannel = statD.label(statD.mask(:,1));
cfg.highlightsize = 15;
cfg.comment = 'Distractor';
cfg.commentpos = 'title';
ft_topoplotER(cfg,statD);
colormap(cm)
colorbar

if set32
   print(fig, fullfile(pth,'results','meg','9 GLM', 'glm spec quinn et al',['tboost_dsuppr_zcorr_set32_', which_corr]),'-dsvg')
else
    print(fig, fullfile(pth,'results','meg','9 GLM', 'glm spec quinn et al',['tboost_dsuppr_zcorr_', which_corr]),'-dsvg')
end



