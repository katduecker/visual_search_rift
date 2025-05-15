%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c3. calculate grandaverage alpha TFR (Fig 4 a)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

function c3_alpha_GA()

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
outpth = fullfile(pth,'results','meg','6 Alpha','pow align iaf');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi/');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));


% load in first subject to get correct time and frequency vector

s = 1;

load(fullfile(iafpth,subj{s},'iaf_soi_not_align.mat'))
load(fullfile(outpth,subj{s},'data_winl_5.mat'),'TFR_alpha_avg')


% combine planar
cfg = [];
cfg.method = 'sum';
TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

% log power
TFR_alpha_avg.powspctrm = log(TFR_alpha_avg.powspctrm);

% average
cfg = [];
cfg.channel = soi_grad_cmb;
cfg.avgoverchan = 'yes';
cfg.latency = [-1.5 0.75];
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);


timevec = TFR_alpha_avg.time;
freq = TFR_alpha_avg.freq;
tfr_ga = zeros(length(subj),length(freq),length(timevec));

tfr_ga(s,:,:) = TFR_alpha_avg.powspctrm;

clear TFR_alpha_avg

% fill in all subjects

for s = 2:length(subj)
    
    load(fullfile(iafpth,subj{s},'iaf_soi_not_align.mat'))
    load(fullfile(outpth,subj{s},'data_winl_5.mat'),'TFR_alpha_avg')
    
    cfg = [];
    cfg.method = 'sum';
    TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);
    
    % log power
    TFR_alpha_avg.powspctrm = log(TFR_alpha_avg.powspctrm);
    
    cfg = [];
    cfg.channel = soi_grad_cmb;
    cfg.latency = [-1.5 0.75];
    cfg.avgoverchan = 'yes';
    TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);
    
    tfr_ga(s,:,:) = TFR_alpha_avg.powspctrm;
    
    clear TFR_alpha_avg

end

save(fullfile(outpth,'tfr_ga.mat'),'timevec','freq','tfr_ga')
% plot
% 
% fig = figure;
% imagesc(timevec,freq,squeeze(mean(tfr_ga,1)))
% xlabel('time (s)')
% ylabel('IAF + (Hz)')
% colormap(cm)
% cb = colorbar;
% cb.Label.String = 'power (rel.chan.)';
