%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c2. Beamformer (DICS) for alpha during search (estimate a new spatial
% filter)

% Input
% - s: subject index

% Output:
% - source localized iaf power high vs low (during search; here actually
% +-2HZ)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Source analysis
% a: align digitized headshape to T1 -> realign Tq
% b: Forward model/Lead field
% c: DICS beamformer
function c2_alpha_during(s)

%% settings & paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path

toi = [-2.5 2];

toi_coh = [0.15 0.5];

dtpth = fullfile(pth, 'data'); % data path
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft','fieldtrip'))
ft_defaults;

% location where fieldtrip is installed
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template', 'sourcemodel');
% load 4 mmm sourcemodel
template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm'));
templatedir = fullfile(ftdir, 'external','spm8','templates');
templmri = ft_read_mri(fullfile(templatedir,'T1.nii'));

% paths
dtpth = fullfile(pth,'data');                                       % raw data
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');

mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
exppth = fullfile(pth,'experiment');
outpth = fullfile(pth,'results','meg','7 Beamformer');
alphapth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
addpath(genpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/'));

% list subj
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));


% load in trial structure (information on samples defining each trial)
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))

fs =1000;
% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)+fs*toi(1),meginfo.alltrl_bl{p}(:,3)+toi(2)*fs,zeros(length(meginfo.alltrl_bl{p}),1)+toi(1)*fs];
end

% list clean data files
d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];

load(fullfile(inpth,subj{s},files{1}));
data_load = data_trig;
perf_coh = perf_cur;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,subj{s},files{f}));
    data_load = ft_appenddata([],data_load,data_trig);
    perf_coh = [perf_coh;perf_cur];
end

data = data_load;

clear data_load

% find kappa (rank of maxfiltered data)
cfg = [];
cfg.channel = 'MEG';
meg = ft_selectdata(cfg,data);


cfg = [];
cfg.covariance = 'yes';
cov = ft_timelockanalysis(cfg,meg);

[~,sv,~] = svd(cov.cov);
d = -diff(log10(diag(sv)));
d = d./std(d);
kappa= find(d>4,1,'first');

clear cov d sv meg



%% Split into alpha high low during search

load(fullfile(alphapth,subj{s},'iaf_soi.mat'))

cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.channel = 'MEGGRAD';
cfg.foi = 4:2:30;
cfg.toi = -2.25:0.05:1.75;
cfg.t_ftimwin = ones(length(cfg.foi)).*0.5;
cfg.keeptrials = 'yes';
TFR_alpha = ft_freqanalysis(cfg,data);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);

cfg = [];
cfg.latency = [0.25 0.5];
cfg.frequency = iaf_grad;
cfg.channel = soi_grad_cmb;
cfg.avgovertime = 'yes';
cfg.avgoverfrequency = 'yes';
cfg.avgoverchan = 'yes';
TFR_alpha = ft_selectdata(cfg,TFR_alpha);


iaf_high = TFR_alpha.powspctrm>median(TFR_alpha.powspctrm);
iaf_low = TFR_alpha.powspctrm<median(TFR_alpha.powspctrm);

clear TFR_alpha

%% ALPHA

% spatial filter

cfg = [];
cfg.toilim = [0.25 0.5];
cfg.minlength = 'maxperlen';
data_search = ft_redefinetrial(cfg,data);
cfg.toilim = [-1 0];
data_bsl = ft_redefinetrial(cfg,data);

cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = [iaf_grad-2 iaf_grad+2];
cfg.keeptrials = 'yes';
cfg.channel = {'MEGGRAD'};
freq_csd_iaf = ft_freqanalysis(cfg,data_search);

load(fullfile(outpth,subj{s},'head_leadf_avg.mat'))
cfg = [];
cfg.method            = 'dics';
cfg.keeptrials        = 'no';
cfg.dics.keepfilter   = 'yes';
cfg.dics.projectnoise = 'yes';
cfg.dics.kappa        = kappa;
cfg.senstype          = 'meg';
cfg.grad              = mGrad;
cfg.sourcemodel       = leadf;
cfg.headmodel         = headmodel;
cfg.projectmom   = 'no';
iaf_dics       = ft_sourceanalysis(cfg, freq_csd_iaf);

clear freq_csd_iaf

iaf_dics.dim = template.sourcemodel.dim;
iaf_dics.pos = template.sourcemodel.pos;

save(fullfile(outpth,subj{s},'dics_filt_iaf_search.mat'),'iaf_dics')


% high vs low
cfg = [];
cfg.trials = iaf_high;
data_high = ft_selectdata(cfg,data_search);
cfg.trials = iaf_low;
data_low = ft_selectdata(cfg,data_search);


cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = iaf_grad;
cfg.keeptrials = 'yes';
cfg.channel = {'MEGGRAD'};
freq_csd_high = ft_freqanalysis(cfg,data_high);
freq_csd_low = ft_freqanalysis(cfg,data_low);


cfg = [];
cfg.method            = 'dics';
cfg.keeptrials        = 'no';
cfg.dics.keepfilter   = 'yes';
cfg.dics.projectnoise = 'yes';
cfg.dics.kappa        = kappa;
cfg.senstype          = 'meg';
cfg.grad              = mGrad;
cfg.sourcemodel       = leadf;
cfg.headmodel         = headmodel;
cfg.projectmom   = 'yes';
cfg.sourcemodel.filter = iaf_dics.avg.filter;
cfg.sourcemodel.label = iaf_dics.avg.label;

iaf_high_dics = ft_sourceanalysis(cfg,freq_csd_high);
iaf_low_dics = ft_sourceanalysis(cfg,freq_csd_low);


iaf_high_dics = iaf_high_dics.avg.pow;
iaf_low_dics = iaf_low_dics.avg.pow;

save(fullfile(outpth,subj{s},'dics_iaf_high_low_search.mat'),'iaf_high_dics','iaf_low_dics')


clear freq_* iaf_high_dics iaf_low_dics data_high data_low
