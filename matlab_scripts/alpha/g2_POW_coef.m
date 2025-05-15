%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d: prep: alpha vs time-on-task!

% Input
% -s : subject index


% Output
% csv file per participant, containing RT and alpha power for each trial 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 5 Aug 2023

clear all; close all; clc; beep off;

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

alphapth = fullfile(pth,'results','meg','6 Alpha');
outpth = fullfile(alphapth,'GLM of TFR int');
soipth = fullfile(alphapth,'iaf_soi');
plotpth = fullfile(pth,'results','exploratory fig');


addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

%% load template subject

s = 1;
d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];

load(fullfile(soipth,subj{s},'iaf_soi.mat'))

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

TFR_alpha = ft_freqanalysis(cfg,data_trig);

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


%% all Channel


regr_names = {'rt','tot','slouch','rt*tot','tot*slouch','guided','set32','gui*set32'};

IAFall = cell(length(regr_names),length(subj));

% check alpha effect
for s = 1:length(subj)
    
    load(fullfile(outpth,subj{s},'Power_GLM_full.mat'),'model_T')
    
    for i = 1:length(regr_names)
        IAF.powspctrm(:,:,:) = squeeze(model_T(i+1,:,:,:));
        
        IAFall{i,s} = IAF;
        
    end
    
end


%% Cluster test
% RT

% get grand average

for i = 1:length(regr_names)
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'powspctrm';
    
    grand_freq = ft_freqgrandaverage(cfg,IAFall{i,:});
    
    % copy in grad struct from one subject
    grand_freq.grad = IAFall{i,15}.grad;
    
    null_hyp = grand_freq;
    null_hyp.powspctrm = zeros(size(null_hyp.powspctrm));
    
    cfg = [];
    cfg.feedback = 'no';
    cfg.method = 'template';
    cfg.template = 'neuromag306cmb_neighb.mat';
    cfg.neighbours = ft_prepare_neighbours(cfg,grand_freq);
    
    cfg.method           = 'montecarlo';
    %cfg.minnbchan        = 2;
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.correcttail     = 'alpha';
    cfg.numrandomization = 2000;
    cfg.alpha            = 0.01;
    cfg.clusteralpha = 0.01;
    cfg.latency          = [-1.25 0.5];
    cfg.parameter = 'powspctrm';
    design = ones(2,2*length(subj));
    design(1,:) = [1:length(subj), 1:length(subj)];
    design(2,length(subj)+1:2*length(subj)) = 2;
    
    

    
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    
    stat{i} = ft_freqstatistics(cfg,grand_freq,null_hyp);
end


save(fullfile(plotpth,'stat_GLM_full_model_template_wcm.mat'),'stat','regr_names')

% 
% %% all Channel
% 
% 
% regr_names = {'rt','tot','slouch','tot*slouch'};
% 
% IAFall = cell(length(regr_names),length(subj));
% 
% % check alpha effect
% for s = 1:length(subj)
%     
%     load(fullfile(outpth,subj{s},'Power_GLM_rt.mat'),'model_T')
%     
%     for i = 1:size(model_T,1)-1
%         IAF.powspctrm(:,:,:) = squeeze(model_T(i+1,:,:,:));
%         
%         IAFall{i,s} = IAF;
%         
%     end
%     
% end
% 
% 
% 
% %% Cluster test
% % RT
% 
% % get grand average
% 
% for i = 1:length(regr_names)
%     cfg = [];
%     cfg.keepindividual = 'yes';
%     cfg.parameter = 'powspctrm';
%     
%     grand_freq = ft_freqgrandaverage(cfg,IAFall{i,:});
%     
%     null_hyp = grand_freq;
%     null_hyp.powspctrm = zeros(size(null_hyp.powspctrm));
%     
%     cfg = [];
%     cfg.neighbours       = 'neuromag306cmb_neighb.mat';                 % fieldtrip template
%     cfg.method           = 'montecarlo';
%     %cfg.minnbchan        = 2;
%     cfg.statistic        = 'ft_statfun_depsamplesT';
%     cfg.correctm         = 'cluster';
%     cfg.clusterstatistic = 'maxsum';
%     cfg.tail             = 0;
%     cfg.clustertail      = 0;
%     cfg.numrandomization = 2000;
%     cfg.alpha            = 0.025;
%     cfg.clusteralpha = 0.05;
%     cfg.latency          = [-1.25 0.5];
%     cfg.parameter = 'powspctrm';
%     design = ones(2,2*length(subj));
%     design(1,:) = [1:length(subj), 1:length(subj)];
%     design(2,length(subj)+1:2*length(subj)) = 2;
%     
%     cfg.design   = design;
%     cfg.uvar     = 1;
%     cfg.ivar     = 2;
%     
%     stat{i} = ft_freqstatistics(cfg,grand_freq,null_hyp);
% end
% 
% 
% save(fullfile(plotpth,'stat_GLM_rt.mat'),'stat','regr_names')
