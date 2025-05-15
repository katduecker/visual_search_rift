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


addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('qual','Dark2',8);

%% Check multicollinearity first

reg_short = {'gui','set','tp','rt','tot','slouch'};
regr_names = {'guided','set size','target present','rt','tot','slouch'};
% fig = figure('Position',[0 0 1940 1080/2.5]);
% multi_coll_all = zeros(length(reg_short)+1,length(reg_short)+1,31);
% idx_regr = 2:length(reg_short)+1;
% for s = 1:length(subj)
%     
%     load(fullfile(outpth,subj{s},'Power_GLM_allset_condi.mat'),'multi_coll')
%     
%     multi_coll_all(:,:,s) = multi_coll;
%     
%     for i = 2:size(multi_coll,1)
%         subplot(2,4,i-1)
%         sctr = scatter(1:length(reg_short)-1,multi_coll(i,setdiff(idx_regr,i)),'filled','MarkerFaceAlpha',0.5,'MarkerFaceColor',cm(i-1,:));
%         hold on
%         xticks(1:length(reg_short))
%         xlim([0 length(reg_short)])
%         xticklabels(reg_short(setdiff(idx_regr,i)-1))
%         ylim([-1 1])
%         title(['corr ',regr_names{i-1},' with'])
%         xtickangle(45)
% 
%     end
%     
% end
% 
% print(fig,fullfile(plotpth,'multi_collin_condi_model'),'-dpng')

%% load template subject
close all;
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


TFRall = cell(length(regr_names),length(subj));

for s = 1:length(subj)
    
    load(fullfile(outpth,subj{s},'Power_GLM_allset_condi.mat'),'model_T')
    
    for i = 1:length(regr_names)
        IAF.powspctrm(:,:,:) = squeeze(model_T(i+1,:,:,:));
        
        TFRall{i,s} = IAF;
        
    end
end

% average T's
for i = 1:length(regr_names)
    Tregr{i} = ft_freqgrandaverage([],TFRall{i,:});
end


addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

% plot alpha band
fig = figure;
cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [-1 0];
cfg.ylim = [8 12];
cfg.commentpos = 'title';
for i = 1:length(regr_names)
    subplot(2,4,i)
    cfg.comment = regr_names{i};
    ft_topoplotTFR(cfg,Tregr{i});
    colormap(cm)
end


print(fig,fullfile(plotpth,'GLM_condi_model'),'-dpng')


%% Cluster-based test
close all;

for i = 1:length(regr_names)
    cfg = [];
    cfg.keepindividual = 'yes';
    cfg.parameter = 'powspctrm';
    
    grand_freq = ft_freqgrandaverage(cfg,TFRall{i,:});
    
    % copy in grad struct from one subject
    grand_freq.grad = TFRall{i,15}.grad;
    
    null_hyp = grand_freq;
    null_hyp.powspctrm = zeros(size(null_hyp.powspctrm));
    
    cfg = [];
    cfg.feedback = 'no';
    cfg.method = 'template';
    cfg.template = 'neuromag306cmb_neighb.mat';
    cfg.neighbours = ft_prepare_neighbours(cfg,grand_freq);
    
    cfg.method           = 'montecarlo';
    cfg.minnbchan        = 1;
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.numrandomization = 5000;
    cfg.alpha            = 0.005;
    cfg.clusteralpha     = 0.01;
    cfg.latency          = [-1.2 0.4];
    cfg.parameter = 'powspctrm';
    design = ones(2,2*length(subj));
    design(1,:) = [1:length(subj), 1:length(subj)];
    design(2,length(subj)+1:2*length(subj)) = 2;
    
    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;
    
    stat{i} = ft_freqstatistics(cfg,grand_freq,null_hyp);
end



save(fullfile(plotpth,'GLM_stat_allset_condi.mat'),'stat','Tregr','regr_names')