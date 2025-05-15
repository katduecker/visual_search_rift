%% VS + RFT
% PhD project 2

% e. alpha power for fast vs slow trials

% Inputs
% - s: subject index
% - winl: window length in seconds
% - c_idx: condition index

% Output
% - TFR contrasts between conditions

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 28/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions 

clear all; close all; clc; beep off;
condi = {{'ti','16t'}, {'ni','16t'},{'ti','32t'},{'ni','32t'}};

%% settings
pth = 'Z:\Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha','RT');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi/');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphafigpth = fullfile(pth,'results','meg','6 Alpha','fig');

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

addpath('Z:\fieldtrip')
ft_defaults;

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
% load template TFR
load(fullfile(outpth,subj{1},[strjoin(condi{1},'_'),'_5_RT.mat']),'TFR_cond')
TFR_cond.label = {'MEG2012+2013'};

pow_fast = cell(1,length(subj));
pow_slow = cell(1,length(subj));
pow_contr = cell(1,length(subj));
fig = figure('Position',[0 0 1920/1.5 1080]);

for c = 1:length(condi)

    for s = 1:length(subj)
        load(fullfile(outpth,subj{s},[strjoin(condi{c},'_'),'_RTsplit.mat']))

        pow_fast{s} = TFR_cond;
        pow_slow{s} = TFR_cond;

        pow_fast{s}.powspctrm = iaf_fast;
        pow_slow{s}.powspctrm = iaf_slow;

        cfg = [];
        cfg.latency = [-1.25 0.5];
        pow_fast{s} = ft_selectdata(cfg,pow_fast{s});
        pow_slow{s} = ft_selectdata(cfg,pow_slow{s});

        pow_contr{s} = pow_fast{s};

        pow_contr{s}.powspctrm = pow_fast{s}.powspctrm./pow_slow{s}.powspctrm-1;
    end

    cfg = [];
    cfg.channel          = 'MEG2012+2013';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = 1;
    cfg.clustertail      = 1;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 500;
    cfg.neighbours = [];
    % specifies with which sensors other sensors can form clusters

    design = ones(2,2*length(subj));
    design(1,:) = [1:length(subj), 1:length(subj)];
    design(2,length(subj)+1:2*length(subj)) = 2;

    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;

    [stat] = ft_freqstatistics(cfg, pow_fast{:}, pow_slow{:});

    GA_contr = ft_freqgrandaverage([],pow_contr{:});

    GA_contr.mask = stat.mask;

    cfg = [];
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
   %cfg.parameter = 'stat';
    cfg.colormap = cm;
    cfg.zlim = [-0.15 0.15];
    cfg.title = strjoin(condi{c},' ');
    subplot(2,2,c)
    ft_singleplotTFR(cfg,GA_contr);
    xlabel('time (s)')
    xticks([-1:0.5:0.5])
    cb = colorbar;
    cb.FontName = 'Arial';
    cb.Ticks = -0.15:0.15:0.15;
    
end

print(fig,fullfile(alphafigpth,'alpha_TFR_fast_vs_slow'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'alpha_TFR_fast_vs_slow'),'-dpng','-r0')

