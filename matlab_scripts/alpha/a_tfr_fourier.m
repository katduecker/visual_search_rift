%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a. TFR of alpha over all trials

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 4 May 2026


function a_tfr_fourier(s)

% init spectrogram 
toi = -1.75:0.05:0.5;
winl = 0.5;
%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
outpth = fullfile(pth,'results','meg','6 Alpha', 'pow');

rmpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')
addpath('/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/fieldtrip')

load(fullfile(pth,'matlab_scripts/',"preproc_meg/",'idx_subjoi.mat'));

outpth = fullfile(outpth,subj{s});
mkdir(outpth)

%% list all files
d = dir(fullfile(inpth,subj{s}));

d = {d.name};

files = d(end);

load(fullfile(inpth,subj{s},files{1}));

%% TFR over all trials
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEGGRAD';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = toi;
cfg.keeptrials = 'yes';
cfg.output = 'fourier';
cfg.pad = 'nextpow2';
TFR = ft_freqanalysis(cfg,data);

TFR.powspctrm = abs(TFR.fourierspctrm);

TFR = rmfield(TFR,'fourierspctrm');

cfg = [];
cfg.avgoverrpt = 'yes';
TFR_avg = ft_selectdata(cfg,TFR);

TFR = rmfield(TFR,'cfg');
TFR_avg = rmfield(TFR_avg,'cfg');

save(fullfile(outpth,['data_fourier_winl_',num2str(winl*10),'.mat']),'TFR','TFR_avg','-v7.3')