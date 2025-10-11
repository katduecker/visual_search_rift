%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c2. align frequencies of TFR to IAF 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow


function c2_align_TFR(s,winl)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','6 Alpha','pow');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi/');
outpth = fullfile(pth,'results','meg','6 Alpha','pow align iaf');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

d = dir(inpth);
d = {d.name};
folds = d(strncmp(d,'202',3));


load('alpha_align_vec.mat')

load(fullfile(inpth,folds{s},['data_winl_',num2str(winl*10),'.mat']))
load(fullfile(iafpth,folds{s},'iaf_soi.mat'),'iaf_grad')

TFR_alpha.freq = round(TFR_alpha.freq - iaf_grad);
TFR_alpha_avg.freq = round(TFR_alpha_avg.freq - iaf_grad);


cfg = [];
cfg.frequency = [max_minf min_maxf];
TFR_alpha = ft_selectdata(cfg,TFR_alpha);
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);

mkdir(fullfile(outpth,folds{s}))
save(fullfile(outpth,folds{s},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha','TFR_alpha_avg','-v7.3')

