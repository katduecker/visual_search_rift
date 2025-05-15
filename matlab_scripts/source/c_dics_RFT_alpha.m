%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c. Beamformer (DICS)
% For each tagging frequency -> beamform tagging signal for each P
% For alpha before the search -> beamform for each P

% Input
% - s: subject index

% Output:
% - source localized rft vs baseline
% - source localized iaf power high vs low

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Source analysis
% a: align digitized headshape to T1 -> realign Tq
% b: Forward model/Lead field
% c: DICS beamformer



function c_dics_RFT_alpha(s)

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


%% rename diodes to Diode 60 and Diode 67

load(fullfile(pth, 'experiment','trigdef.mat'))

% trials T=60; D=67

x = cell2mat(cellfun(@(x) ~isempty(x), regexp(trigdef(:,2),'6067'), 'UniformOutput',false));
trig6067 = [trigdef{x,1}];
trl6067 = ismember([perf_coh{:,1}],trig6067);

cfg = [];
cfg.trials = find(trl6067);
data6067 = ft_selectdata(cfg,data);

data6067.label{end-1} = 'diode 60';
data6067.label{end} = 'diode 67';

% trials T=67; D=60

x = cell2mat(cellfun(@(x) ~isempty(x), regexp(trigdef(:,2),'6760'), 'UniformOutput',false));
trig6760 = [trigdef{x,1}];
trl6760 = ismember([perf_coh{:,1}],trig6760);

cfg = [];
cfg.trials = find(trl6760);
data6760 = ft_selectdata(cfg,data);

data6760.label{end-1} = 'diode 67';
data6760.label{end} = 'diode 60';

cfg = [];
cfg.channel = 'MEG';
data6067meg = ft_selectdata(cfg,data6067);
data6760meg = ft_selectdata(cfg,data6760);

cfg.channel = 'diode 60';
data6067d60 = ft_selectdata(cfg,data6067);
data6760d60 = ft_selectdata(cfg,data6760);

cfg.channel = 'diode 67';
data6067d67 = ft_selectdata(cfg,data6067);
data6760d67 = ft_selectdata(cfg,data6760);

data6067 = ft_appenddata([],data6067meg,data6067d60,data6067d67);
data6760 = ft_appenddata([],data6760meg,data6760d60,data6760d67);

data_rename = ft_appenddata([],data6067,data6760);

%% Frequency analysis

cfg = [];
cfg.toilim = [0 .5];
data_rft = ft_redefinetrial(cfg,data_rename);

cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = 60;
cfg.pad = 1;
cfg.keeptrials = 'yes';
cfg.channel = {'MEGGRAD' 'diode 60'};
cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode 60'};
freq_csd_60 = ft_freqanalysis(cfg,data_rft);
cfg.foi = 67;
cfg.channel = {'MEGGRAD' 'diode 67'};
cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode 67'};
freq_csd_67 = ft_freqanalysis(cfg,data_rft);

clear data6* data_rft 


%% Estimate spatial filter
% 
% estimate dics filter 
% save filter and apply to data
% downproject after beamforming

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
cfg.refchan           = 'diode 60';
rft_dics{1}       = ft_sourceanalysis(cfg, freq_csd_60);
cfg.refchan           = 'diode 67';
rft_dics{2}       = ft_sourceanalysis(cfg, freq_csd_67);

clear freq_csd*


% change frequency by hand to be able to average
rft_dics{1}.freq = 0;
rft_dics{2}.freq = 0;

% adjust dimension and position vectors (use template model
rft_dics{1}.dim = template.sourcemodel.dim;
rft_dics{1}.pos = template.sourcemodel.pos;
rft_dics{2}.dim = template.sourcemodel.dim;
rft_dics{2}.pos = template.sourcemodel.pos;

% close all

save(fullfile(outpth,subj{s},'dics_filt_rft.mat'),'rft_dics')


clear data_rename src*
%% Beamform data for each condition separately

% use only set size 32 because better signal
condi = {{'ni','32t','6067'},{'ti','32t','6067'},{'ni','32t','6760'},{'ti','32t','6760'}};

rft_freq = [60,67];
for c = 1:length(condi)
    % select current condition
    c_idx = (cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{1}),'UniformOutput',false)) + ...
        cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{2}),'UniformOutput',false)) + ...
        cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{3}),'UniformOutput',false))) == 3;
    
    trig_condi = [trigdef{c_idx,1}];
    
    trl_condi = ismember([perf_coh{:,1}],trig_condi);
    
    cfg = [];
    cfg.trials = trl_condi;
    data_cond = ft_selectdata(cfg,data);
    
    cfg = [];
    cfg.toilim = [0.0 .5];
    cfg.minlength = 'maxperlen';
    data_rft = ft_redefinetrial(cfg,data_cond);
    
    cfg.toilim = [-.5 0];
    cfg.minlength = 'maxperlen';
    bsl_rft = ft_redefinetrial(cfg,data_cond);

    cfg = [];
    cfg.output = 'powandcsd';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';
    cfg.pad = 1;

    % coherence to Target diode
    cfg.foi = [str2num(condi{c}{3}(1:2))];
    cfg.channel = {'MEGGRAD' 'diode T'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode T'};
    freq_csdT = ft_freqanalysis(cfg,data_rft);
    freq_csdT_bsl = ft_freqanalysis(cfg,bsl_rft);
    
    % coherence to distractor diode
    cfg.channel = {'MEGGRAD' 'diode D'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode D'};
    cfg.foi = [str2num(condi{c}{3}(3:4))];
    freq_csdD = ft_freqanalysis(cfg,data_rft);
    freq_csdD_bsl = ft_freqanalysis(cfg,bsl_rft);
    
    f = find(str2num(condi{c}{3}(1:2)) == rft_freq);
    cfg = [];
    cfg.method            = 'dics';
    cfg.grad              = mGrad;
    cfg.sourcemodel       = leadf;
    cfg.dics.projectnoise = 'yes';
    cfg.dics.kappa        = kappa;
    cfg.senstype          = 'meg';
    cfg.projectmom        = 'yes';
    cfg.sourcemodel.filter = rft_dics{f}.avg.filter;
    cfg.sourcemodel.label = rft_dics{f}.avg.label;
    cfg.headmodel         = headmodel;
    
    % Target coherence
    cfg.refchan = 'diode T';
    
    % extract source localized vectors to be able to save
    x = ft_sourceanalysis(cfg,freq_csdT);
    rft_src{c,1}.coh = x.avg.coh;
    rft_src{c,1}.pow = x.avg.pow;
    
    x = ft_sourceanalysis(cfg,freq_csdT_bsl);
    bsl_src{c,1}.coh = x.avg.coh;
    bsl_src{c,1}.pow = x.avg.pow;
    
    
    f = find(str2num(condi{c}{3}(3:4)) == rft_freq);
    cfg.sourcemodel.filter = rft_dics{f}.avg.filter;
    cfg.sourcemodel.label = rft_dics{f}.avg.label;
    
    % Distractor coherence
    cfg.refchan = 'diode D';

    x = ft_sourceanalysis(cfg,freq_csdD);
    rft_src{c,2}.coh = x.avg.coh;
    rft_src{c,2}.pow = x.avg.pow;
    
    x = ft_sourceanalysis(cfg,freq_csdD_bsl);
    bsl_src{c,2}.coh = x.avg.coh;
    bsl_src{c,2}.pow = x.avg.pow;
    
    clear data_cond freq_c* trig_condi trl_condi data_rft bsl_rft
end


save(fullfile(outpth,subj{s},'dics_rft_condi.mat'),'bsl_src','rft_src','condi','-v7.3')


%% Beamformer for ALPHA

% Frequency analysis
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
cfg.latency = [-1 0];
cfg.frequency = iaf_grad;
cfg.channel = soi_grad_cmb;
cfg.avgovertime = 'yes';
cfg.avgoverfrequency = 'yes';
cfg.avgoverchan = 'yes';
TFR_alpha = ft_selectdata(cfg,TFR_alpha);


iaf_high = TFR_alpha.powspctrm>median(TFR_alpha.powspctrm);
iaf_low = TFR_alpha.powspctrm<median(TFR_alpha.powspctrm);


condi = {{'ti','32t','6067'},{'ti','32t','6760'}};

rft_freq = [60,67];

for c = 1:length(condi)
    
    %% high alpha
    % select current condition
    c_idx = (cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{1}),'UniformOutput',false)) + ...
        cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{2}),'UniformOutput',false)) + ...
        cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{3}),'UniformOutput',false))) == 3;
    
    trig_condi = [trigdef{c_idx,1}];
    
    trl_condi_high = find((ismember([perf_coh{:,1}],trig_condi) + iaf_high') == 2);
    trl_condi_low = find((ismember([perf_coh{:,1}],trig_condi) + iaf_low') == 2);
    
    if length(trl_condi_high) > length(trl_condi_low)
        x = randperm(length(trl_condi_high));
        trl_condi_high = sort(trl_condi_high(x(1:length(trl_condi_low))));
        
    elseif length(trl_condi_high) < length(trl_condi_low)
        x = randperm(length(trl_condi_low));
        trl_condi_low = sort(trl_condi_low(x(1:length(trl_condi_high))));
    end

    cfg = [];
    cfg.trials = trl_condi_high;
    data_cond = ft_selectdata(cfg,data);
    
    cfg = [];
    cfg.toilim = [0.0 .5];
    cfg.minlength = 'maxperlen';
    data_rft = ft_redefinetrial(cfg,data_cond);
    
    cfg.toilim = [-.5 0];
    cfg.minlength = 'maxperlen';
    bsl_rft = ft_redefinetrial(cfg,data_cond);

    cfg = [];
    cfg.output = 'powandcsd';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';

    % coherence to Target diode
    cfg.foi = [str2num(condi{c}{3}(1:2))];
    cfg.channel = {'MEGGRAD' 'diode T'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode T'};
    freq_csdT = ft_freqanalysis(cfg,data_rft);
    freq_csdT_bsl = ft_freqanalysis(cfg,bsl_rft);
    
    % coherence to distractor diode
    cfg.channel = {'MEGGRAD' 'diode D'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode D'};
    cfg.foi = [str2num(condi{c}{3}(3:4))];
    freq_csdD = ft_freqanalysis(cfg,data_rft);
    freq_csdD_bsl = ft_freqanalysis(cfg,bsl_rft);
    
    f = find(str2num(condi{c}{3}(1:2)) == rft_freq);
    cfg = [];
    cfg.method            = 'dics';
    cfg.grad              = mGrad;
    cfg.sourcemodel       = leadf;
    cfg.dics.projectnoise = 'yes';
    cfg.dics.kappa        = kappa;
    cfg.senstype          = 'meg';
    cfg.projectmom        = 'yes';
    cfg.sourcemodel.filter = rft_dics{f}.avg.filter;
    cfg.sourcemodel.label = rft_dics{f}.avg.label;
    cfg.headmodel         = headmodel;
    
    % Target coherence
    cfg.refchan = 'diode T';
    
    % extract source localized vectors to be able to save
    x = ft_sourceanalysis(cfg,freq_csdT);
    rft_src_high{c,1}.coh = x.avg.coh;
    rft_src_high{c,1}.pow = x.avg.pow;
    
    x = ft_sourceanalysis(cfg,freq_csdT_bsl);
    bsl_src_high{c,1}.coh = x.avg.coh;
    bsl_src_high{c,1}.pow = x.avg.pow;
    
    
    f = find(str2num(condi{c}{3}(3:4)) == rft_freq);
    cfg.sourcemodel.filter = rft_dics{f}.avg.filter;
    cfg.sourcemodel.label = rft_dics{f}.avg.label;
    
    % Distractor coherence
    cfg.refchan = 'diode D';

    x = ft_sourceanalysis(cfg,freq_csdD);
    rft_src_high{c,2}.coh = x.avg.coh;
    rft_src_high{c,2}.pow = x.avg.pow;
    
    x = ft_sourceanalysis(cfg,freq_csdD_bsl);
    bsl_src_high{c,2}.coh = x.avg.coh;
    bsl_src_high{c,2}.pow = x.avg.pow;
    
    %% Low alpha
    
    % select current condition
       
    cfg = [];
    cfg.trials = trl_condi_low;
    data_cond = ft_selectdata(cfg,data);
    
    cfg = [];
    cfg.toilim = [0.0 .5];
    cfg.minlength = 'maxperlen';
    data_rft = ft_redefinetrial(cfg,data_cond);
    
    cfg.toilim = [-.5 0];
    cfg.minlength = 'maxperlen';
    bsl_rft = ft_redefinetrial(cfg,data_cond);

    cfg = [];
    cfg.output = 'powandcsd';
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.keeptrials = 'yes';

    % coherence to Target diode
    cfg.foi = [str2num(condi{c}{3}(1:2))];
    cfg.channel = {'MEGGRAD' 'diode T'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode T'};
    freq_csdT = ft_freqanalysis(cfg,data_rft);
    freq_csdT_bsl = ft_freqanalysis(cfg,bsl_rft);
    
    % coherence to distractor diode
    cfg.channel = {'MEGGRAD' 'diode D'};
    cfg.channelcmb = {'MEGGRAD' 'MEGGRAD';'MEGGRAD' 'diode D'};
    cfg.foi = [str2num(condi{c}{3}(3:4))];
    freq_csdD = ft_freqanalysis(cfg,data_rft);
    freq_csdD_bsl = ft_freqanalysis(cfg,bsl_rft);
    
    f = find(str2num(condi{c}{3}(1:2)) == rft_freq);
    cfg = [];
    cfg.method            = 'dics';
    cfg.grad              = mGrad;
    cfg.sourcemodel       = leadf;
    cfg.dics.projectnoise = 'yes';
    cfg.dics.kappa        = kappa;
    cfg.senstype          = 'meg';
    cfg.projectmom        = 'yes';
    cfg.sourcemodel.filter = rft_dics{f}.avg.filter;
    cfg.sourcemodel.label = rft_dics{f}.avg.label;
    cfg.headmodel         = headmodel;
    
    % Target coherence
    cfg.refchan = 'diode T';
    
    % extract source localized vectors to be able to save
    x = ft_sourceanalysis(cfg,freq_csdT);
    rft_src_low{c,1}.coh = x.avg.coh;
    rft_src_low{c,1}.pow = x.avg.pow;
    
    x = ft_sourceanalysis(cfg,freq_csdT_bsl);
    bsl_src_low{c,1}.coh = x.avg.coh;
    bsl_src_low{c,1}.pow = x.avg.pow;
    
    
    f = find(str2num(condi{c}{3}(3:4)) == rft_freq);
    cfg.sourcemodel.filter = rft_dics{f}.avg.filter;
    cfg.sourcemodel.label = rft_dics{f}.avg.label;
    
    % Distractor coherence
    cfg.refchan = 'diode D';

    x = ft_sourceanalysis(cfg,freq_csdD);
    rft_src_low{c,2}.coh = x.avg.coh;
    rft_src_low{c,2}.pow = x.avg.pow;
    
    x = ft_sourceanalysis(cfg,freq_csdD_bsl);
    bsl_src_low{c,2}.coh = x.avg.coh;
    bsl_src_low{c,2}.pow = x.avg.pow;
    clear data_cond freq_c* trig_condi trl_condi data_rft bsl_rft
end

save(fullfile(outpth,subj{s},'dics_rft_condi_alpha_high_low.mat'),'rft_src_high','rft_src_low','condi','iaf_high','-v7.3')


clear rft_* bsl* trig* trl* data_rename data_trig


%% ALPHA

% spatial filter

cfg = [];
cfg.toilim = [-1 0];
cfg.minlength = 'maxperlen';
data_bsl = ft_redefinetrial(cfg,data);


%% Estimate CSD over all trials
cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = [iaf_grad-2 iaf_grad+2];
cfg.keeptrials = 'yes';
cfg.channel = {'MEGGRAD'};
freq_csd_iaf = ft_freqanalysis(cfg,data_bsl);

%% Estimate spatial filter

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

save(fullfile(outpth,subj{s},'dics_filt_iaf.mat'),'iaf_dics')


% high vs low
cfg = [];
cfg.trials = iaf_high;
data_high = ft_selectdata(cfg,data_bsl);
cfg.trials = iaf_low;
data_low = ft_selectdata(cfg,data_bsl);

%% Estimate CSD for high vs low
cfg = [];
cfg.output = 'powandcsd';
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = iaf_grad;
cfg.keeptrials = 'yes';
cfg.channel = {'MEGGRAD'};
freq_csd_high = ft_freqanalysis(cfg,data_high);
freq_csd_low = ft_freqanalysis(cfg,data_low);


%% Beamform using spatial filter over all trials
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

save(fullfile(outpth,subj{s},'dics_iaf_high_low.mat'),'iaf_high_dics','iaf_low_dics')


clear freq_* iaf_high_dics iaf_low_dics data_high data_low





