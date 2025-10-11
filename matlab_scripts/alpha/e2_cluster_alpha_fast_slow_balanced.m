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
condi = {{'ti','16ta'}, {'ni','16ta'},{'ti','32ta'},{'ni','32ta'},{'ti','16tp'}, {'ni','16tp'},{'ti','32tp'},{'ni','32tp'}};


select_soi = 'occi';
%% settings
pth = 'Z:\Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha','RT balanced split');
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

load(fullfile(pth,'results','meg','6 Alpha','pow',subj{1},'data_winl_5.mat'),'TFR_alpha_avg')
% 
cfg = [];
cfg.method = 'sum';
TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

cfg = [];
cfg.latency = [-1.25 0.5];
cfg.frequency = [4 TFR_alpha_avg.freq(12)];
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);

pow_fast = cell(1,length(subj));
pow_slow = cell(1,length(subj));
pow_contr = cell(1,length(subj));

fig = figure('Position',[0 0 1920/1.5 1080]);


for s = 1:length(subj)


    TFR_fast = cell(1,length(condi));
    TFR_slow = cell(1,length(condi));

    for c = 1:length(condi)

        load(fullfile(outpth,subj{s},[strjoin(condi{c},'_'),select_soi,'_RTsplit.mat']))

       TFR_fast{c} = IAFfast;
       TFR_slow{c} = IAFslow;

        clear IAFfast IAFslow
    end

    TFR_fast = ft_appendfreq([],TFR_fast{:});
    TFR_slow = ft_appendfreq([],TFR_slow{:});

    cfg = [];
    cfg.avgoverrpt = 'yes';
    TFR_fast = ft_selectdata(cfg,TFR_fast);
    TFR_slow = ft_selectdata(cfg,TFR_slow);

    cfg = [];
    cfg.method = 'sum';
    TFR_fast = ft_combineplanar(cfg,TFR_fast);
    TFR_slow = ft_combineplanar(cfg,TFR_slow);

    % average over baseline and extract IAF power
    cfg = [];
    cfg.latency = [-1.25 0.5];
    TFR_fast = ft_selectdata(cfg,TFR_fast);
    TFR_slow = ft_selectdata(cfg,TFR_slow);

    cfg = [];
    cfg.avgoverchan = 'yes';

    IAF_fast = ft_selectdata(cfg,TFR_fast);
    IAF_slow = ft_selectdata(cfg,TFR_slow);

    pow_fast{s} = TFR_alpha_avg;
    pow_slow{s} = TFR_alpha_avg;

    pow_fast{s}.powspctrm(1,:,:) = IAF_fast.powspctrm;
    pow_slow{s}.powspctrm(1,:,:) = IAF_slow.powspctrm;

    pow_contr{s} = pow_fast{s};

    pow_contr{s}.powspctrm = pow_fast{s}.powspctrm./pow_slow{s}.powspctrm-1;

end

cfg = [];
cfg.method           = 'montecarlo';
cfg.channel          = pow_fast{s}.label{1};
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 1;
cfg.clustertail      = 1;
cfg.alpha            = 0.05;
cfg.frequency = [4 8];
cfg.numrandomization = 10000;
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
cfg = [];
cfg.frequency = [4 8];
GA_contr = ft_selectdata(cfg,GA_contr)
GA_contr.freq = IAF_fast.freq(1:3);
GA_contr.mask = stat.mask;

cfg = [];
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';
%cfg.parameter = 'stat';
cfg.colormap = cm;
%cfg.zlim = [-0.15 0.15];
cfg.channel          = pow_fast{s}.label{1};
ft_singleplotTFR(cfg,GA_contr);
xlabel('time (s)')
xticks([-1:0.5:0.5])
cb = colorbar;
cb.FontName = 'Arial';
% cb.Ticks = -0.15:0.15:0.15;



print(fig,fullfile(alphafigpth,'alpha_balanced_split_TFR_fast_vs_slow'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'alpha_balanced_split_TFR_fast_vs_slow'),'-dpng','-r0')

