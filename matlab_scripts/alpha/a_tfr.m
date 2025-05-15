%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a. TFR of alpha over all trials

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs
% - s: subject index
% - c_idx: condition index (it's just 1 because we are looking at all data)
% - toi: start and end time point the sliding window is centered over
% - winl: window length
% Output
% - TFR of power with window length 1 sec, Hanning window, moved in 50 sec
% steps

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 23/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

function a_tfr(s)

toi = -1.75:0.05:0.5;
winl = 1;
%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
outpth = fullfile(pth,'results','meg','6 Alpha', 'pow');

rmpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

outpth = fullfile(outpth,subj{s});
mkdir(outpth)

% list condition files that contain 'condi'
d = dir(fullfile(inpth,subj{s}));

d = {d.name};

files = d(end);

load(fullfile(inpth,subj{s},files{1}));

%% TFR
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEGGRAD';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = toi;
cfg.keeptrials = 'yes';
TFR = ft_freqanalysis(cfg,data);

cfg = [];
cfg.avgoverrpt = 'yes';
TFR_avg = ft_selectdata(cfg,TFR);

TFR = rmfield(TFR,'cfg');
TFR_avg = rmfield(TFR_avg,'cfg');

save(fullfile(outpth,['data_winl_',num2str(winl*10),'.mat']),'TFR','TFR_avg','-v7.3')