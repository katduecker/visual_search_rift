%% Visual Search + RFT
clear all; close all; clc;

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
glmpth = fullfile(pth,'results','meg','9 GLM');

cohpth = fullfile(pth,'results','meg','8 COH single trl');


%% load template subject

s = 1;
d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];


load(fullfile(inpth,subj{s},file{1}))

% TFR
winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEGGRAD';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'yes';
cfg.trials = 1;

TFR_alpha = ft_freqanalysis(cfg,data);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);

cfg = [];
cfg.latency = [-1.5 0.5];
cfg.avgoverrpt = 'yes';
IAF = ft_selectdata(cfg,TFR_alpha);

% 
clear data_trig

% clear TFR_alpha

IAF = rmfield(IAF,'cfg');



%% Coefficients

IAF_GLM = cell(6,length(subj));

% check alpha effect
for c = 1:6
    
    for s = 1:length(subj)
        
        load(fullfile(glmpth,subj{s},'GLM_RT_RIFT.mat'),'model_beta')
        
        IAF_GLM{c,s} = IAF;
        
        IAF_GLM{c,s}.powspctrm = squeeze(model_beta(c,:,:,:));
        
        
    end
end

GLM = cell(1,6);
GLM_avg = cell(1,6);
cfg = [];
for c = 1:6
    cfg.keepindividual = 'yes';
    GLM{c} = ft_freqgrandaverage(cfg,IAF_GLM{c,:});
    cfg.keepindividual = 'no';
    GLM_avg{c} = ft_freqgrandaverage(cfg,IAF_GLM{c,:});
end


regr_names = {'constant','RT','tot','slouch','target','distractor'};

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);
figure;
cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [12 16];
cfg.ylim = [-1 0];
cfg.zlim = 'maxabs';
cfg.commentpos = 'title';
for c = 1:6
    subplot(2,3,c)
    cfg.comment = regr_names{c};
    
    ft_topoplotTFR(cfg,GLM{c});
    colorbar
    colormap(cm)
end

null_hyp = GLM{1};
null_hyp.powspctrm = zeros(size(GLM{1}.powspctrm));

cfg = [];
cfg.neighbours       = 'neuromag306cmb_neighb.mat';                 % fieldtrip template                    
cfg.method           = 'montecarlo';
cfg.minnbchan        = 2;
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.numrandomization = 2000;
cfg.alpha            = 0.025;
cfg.clusteralpha = 0.05;
cfg.latency          = [-1.25 0.5];
cfg.parameter = 'powspctrm';
design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

for c = 1:7
    stat{c} = ft_freqstatistics(cfg,GLM{c},null_hyp);
end
