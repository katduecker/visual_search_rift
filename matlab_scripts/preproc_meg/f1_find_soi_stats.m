%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% f1. Identify sensors w/ significant RIFT response (see Maris & Oostenveld,
% 2007; Pan et al. 2021)
% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Output
% - soi_stat: sensors with significant tagging response

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions

function f1_find_soi_stats(s)


% settings to identify sensors as descirbed in manuscript
tbsl = [-.5 -0.2];
start_bsl = -2.5;
end_trl = 2;
tstim = [0.2 .5];
rft_freq = [60,67];
phshft = [0, 0.65];                     % frequency tagging phase shift for 60 vs 67 Hz as used in experiment
% rsmpfrq = 1000;                          % resampling frequency
noise = 0;



% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
%pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
rmpath(genpath('/rds/projects/2018/jenseno_entrainment/fieldtrip'))
% addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA', '1 all subj');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb', 'soi','sinusoid');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

mkdir(fullfile(cohpth,subjfolds{s}))
% 
% % don't run if soi have been identified for this subject
% if exist(fullfile(cohpth,subjfolds{s},'soi_stat.mat'))
%     error(['soi already found for ',subjfolds{s}])
% end

% load clean trial structure

load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'meginfo')

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)+start_bsl*fs,meginfo.alltrl_bl{p}(:,3)+end_trl*fs,repmat(start_bsl*fs,size(meginfo.alltrl_bl{p},1),1)];
end

% list fif files
d = dir(fullfile(dtpth,subjfolds{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load in maxfiltered data using trial structure

for p = 1:length(f)
    cfg = [];
    cfg.dataset = fullfile(dtpth,subjfolds{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = {'MEG','MISC004', 'MISC005'};
    
    % load data and diode separately (to be able to reject ICA)
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);
    
end

data = ft_appenddata([],dtprt{:});
data.hdr = dtprt{1}.hdr;
clear dtprt
meginfo.alltrl_list = meginfo.alltrl_list(logical(meginfo.keeptrl_all),:);

% separate MEG sensors and diodes
cfg = [];
cfg.channel = 'MEG';
cfg.trials =find(meginfo.keeptrl_all);                   % select clean trials
datameg = ft_selectdata(cfg,data);
cfg.channel = {'MISC004','MISC005'};
diodes = ft_selectdata(cfg,data);
clear data

% fix sampleinfo
trl_end_samp = cell2mat(cellfun(@length,datameg.trial,'UniformOutput',false)).*[1:length(datameg.trial)];
datameg.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];
diodes.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];



%% Reject ICA comps

load(fullfile(icapth, [subjfolds{s},'_ica.mat']))

% reject comps
cfg = [];
cfg.component = badcomps; % to be removed component(s)
dataclean = ft_rejectcomponent(cfg, dataICA, datameg);
clear dataICA dtprt diodes_trl datameg

% average grad positions
d = dir(fullfile(dtpth,subjfolds{s}));
files = {d.name};
files = files(cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'.fif'),'UniformOutput',false)));
% 
grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(dtpth,subjfolds{s},files{fl}))];
end
% 
% mean channel position
mGrad = grad(1);
% average grad structure
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
mGrad.coilpos = mGrad.coilpos + grad(g).coilpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
mGrad.coilpos = mGrad.coilpos./length(grad);
clear grad

% concatenate cleaned data and diodes
data = ft_appenddata([],dataclean,diodes);
clear dataclean diodes
data.grad = mGrad;

% relabel diode: 60 Hz vs 67 Hz
load(fullfile(pth,'experiment','trigdef.mat'))


% trials in which yellow (MISC004) was tagged at 60 Hz

% find trigger
misc4_60_trig = cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'126067'),'UniformOutput',false))...
    + cell2mat(cellfun(@(x) ~isempty(x),strfind(trigdef(:,2),'216760'),'UniformOutput',false));

% find trial indices with these triggers
misc4_60_trl = ismember(meginfo.alltrl_list(:,1),[trigdef{logical(misc4_60_trig),1}]);
% trials in which cyan (MISC005) was tagged at 60 Hz
misc5_60_trl = ~misc4_60_trl;

% check if numbers make sense
sum(misc4_60_trl)
sum(misc5_60_trl)

% select data
cfg = [];
cfg.trials = misc4_60_trl;
dat_misc4_60 = ft_selectdata(cfg,data);
cfg.trials = misc5_60_trl;
dat_misc5_60 = ft_selectdata(cfg,data);

clear data

% select diodes and MEG
cfg = [];
cfg.channel = 'MEG';
dat_misc4_60_meg = ft_selectdata(cfg,dat_misc4_60);
dat_misc5_60_meg = ft_selectdata(cfg,dat_misc5_60);
cfg.channel = 'MISC004';
dat_misc4_60_m4 = ft_selectdata(cfg,dat_misc4_60);
dat_misc5_60_m4 = ft_selectdata(cfg,dat_misc5_60);
cfg.channel = 'MISC005';
dat_misc4_60_m5 = ft_selectdata(cfg,dat_misc4_60);
dat_misc5_60_m5 = ft_selectdata(cfg,dat_misc5_60);

clear dat_misc4_60 dat_misc5_60

% relabel based on tagging frequency
dat_misc4_60_m4.label = {'diode 60'};
dat_misc4_60_m5.label = {'diode 67'};
dat_misc5_60_m5.label = {'diode 60'};
dat_misc5_60_m4.label = {'diode 67'};

% append
dat_misc4_60_newlbl = ft_appenddata([],dat_misc4_60_meg,dat_misc4_60_m4,dat_misc4_60_m5);
dat_misc5_60_newlbl = ft_appenddata([],dat_misc5_60_meg,dat_misc5_60_m5,dat_misc5_60_m4);

clear dat_misc4_60_m4 dat_misc4_60_m5 dat_misc4_60_meg dat_misc5_60_m4 dat_misc5_60_m5 dat_misc5_60_meg

% append
data = ft_appenddata([],dat_misc4_60_newlbl,dat_misc5_60_newlbl);
clear data_misc*
diode_idx = find(ismember(data.label,{'diode 60','diode 67'}));

% replace diode with perfect sinusoid
data = kd_replace_diode_sinu(data,diode_idx,rft_freq,phshft,start_bsl,fs,noise);

% % resample
% cfg = [];
% cfg.resamplefs = rsmpfrq;
% data = ft_resampledata(cfg,data);


% separate baseline and flicker
cfg = [];
cfg.latency = tbsl;
data_bsl = ft_selectdata(cfg,data);
cfg.latency = tstim;
data_rft = ft_selectdata(cfg,data);

clear data
% fourier spectra
cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 60;
%cfg.pad          = 'nextpow2';
fourier          = ft_freqanalysis(cfg,data_rft);
fourier_bl       = ft_freqanalysis(cfg,data_bsl);
clear data*

% statistical test for each occipital sensor
load(fullfile(pth,'matlab scripts','preprocessing MEG','occi_sens.mat'))
nchancom = length(occi_soi);
stat_mask_60 = zeros(nchancom,1);

for i = 1:nchancom
    cfg            = [];
    cfg.channel    = {occi_soi{i},'diode 60'};
    fourier_tmp    = ft_selectdata(cfg, fourier);
    fourier_bl_tmp = ft_selectdata(cfg, fourier_bl);
    
    cfg                  = [];
    cfg.parameter        = 'fourierspctrm';
    cfg.statistic        = 'ft_statfun_indepsamplesZcoh';  %%%% take fourierspctrm as input, so no time domain information
    cfg.method           = 'montecarlo';
    cfg.frequency        = 60;
    cfg.tail             = 1; %% right sided, grp1 is bigger than grp2
    cfg.alpha            = 0.01;
    cfg.numrandomization = 10000;
    ntrl_1 = size(fourier.fourierspctrm,1);
    ntrl_2 = size(fourier_bl.fourierspctrm,1);
    design = zeros(1, ntrl_1 + ntrl_2);
    design(1,1:ntrl_1) = 1;
    design(1,(ntrl_1 +1):(ntrl_1 + ntrl_2))= 2;
    cfg.design = design;
    cfg.ivar   = 1;
    cfg.design = design;
    stat = ft_freqstatistics(cfg, fourier_tmp, fourier_bl_tmp);
    stat_mask_60(i) = stat.mask;
    clear fourier_tmp fourier_bl_tmp
end

% save sois for all p

%soi_stat = cell(60,1);

%load(fullfile(pth,'matlab scripts', 'preprocessing','soi_stat_allp.mat'))
soi_stat = occi_soi(logical(stat_mask_60));

save(fullfile(cohpth,subjfolds{s},'soi_stat_consv.mat'),'soi_stat','noise','tbsl','tstim')



