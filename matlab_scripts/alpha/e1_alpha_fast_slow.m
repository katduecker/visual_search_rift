%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e1. calculate TFRs for fast vs slow trials (Fig. 4 d)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs
% - s: subject index
% - c_idx: condition index
% - data_trim (bool): discard trials w/ RT +- 3*std?
% - split_ta_tp (bool): split into target absent/present?

% Output
% - soi_grad: gradiometers with high alpha power
% iaf_grad: identified IAF in gradiometers

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 23/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow


function e1_alpha_fast_slow(s,c_idx,data_trim,split_ta_tp)



condi_all = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

condi = condi_all{c_idx};


if split_ta_tp

    split_suf = '_ta_tp';
else
    split_suf = '';
end

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha','RT');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphasoipth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_not_align.mat'));

%% load data

d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];


% split into Target absent/Target present?

if split_ta_tp
    % if yes, separate trials for ta/tp
    ta_tp = {'ta','tp'};
else
    % if not, just concatenate a t to the set size (in filenames) which
    % looks for ta/tp
    ta_tp = {'t'};
end
for ti = 1:length(ta_tp)

% find relevant data files
condi_files = zeros(length(files),1);
for c = 1:length(condi)
    
    condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c}),'UniformOutput',false))';
    
end
condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,[condi{2}(end-2:end-1),ta_tp{ti}]),'UniformOutput',false))';

condi_files = condi_files == length(condi)+1;   

c_files = files(condi_files);

load(fullfile(inpth,subj{s},c_files{1}));
trl_idx = find(trlcur);

% performance in first trials
trl_perf = perf_cur;

% load data
data_load = data_trig;
% data_load = rmfield(data_load,'sampleinfo');
% load & append data
for f = 2:length(c_files)
    clear data_trig trlcur perf_cur
    load(fullfile(inpth,subj{s},c_files{f}));
    % append data
    data_load = ft_appenddata([],data_load,data_trig);
    %append performance
    trl_perf = [trl_perf;perf_cur];
end

data = data_load;

clear data_load


%% Median split fast slow


rt_all = [trl_perf{:,3}];


% trim data
if data_trim
    m = mean(rt_all);
    std_rt = std(rt_all);
    
    trl_idx = logical(([trl_perf{:,3}]< m-3*std_rt) + ([trl_perf{:,3}] > m+3*std_rt));

    cfg = [];
    cfg.trials = ~trl_idx;
    data = ft_selectdata(cfg,data);

    trl_perf(trl_idx,:) = [];
 
end

rt_all = [trl_perf{:,3}];

mRT = median(rt_all);

trl_fast = [trl_perf{:,3}] < mRT;

trl_slow = [trl_perf{:,3}] > mRT;


% select trials
cfg = [];
cfg.trials = find(trl_fast);
data_fast = ft_selectdata(cfg,data);

cfg.trials = find(trl_slow);
data_slow = ft_selectdata(cfg,data);

clear data


%% Alpha power in these trials
load(fullfile(alphasoipth,subj{s},'iaf_soi_not_align.mat'))

load('alpha_align_vec_not_aligned.mat')


winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEG';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'no';
TFRfast = ft_freqanalysis(cfg,data_fast);
TFRslow = ft_freqanalysis(cfg,data_slow);

cfg = [];
cfg.method = 'sum';
TFRfast = ft_combineplanar(cfg,TFRfast);
TFRslow = ft_combineplanar(cfg,TFRslow);

% average over baseline and extract IAF power
cfg = [];
cfg.channel = soi_grad_cmb;
cfg.avgoverchan = 'yes';
IAFfast = ft_selectdata(cfg,TFRfast);
IAFslow = ft_selectdata(cfg,TFRslow);


IAFfast.freq = IAFfast.freq-iaf_grad;
IAFslow.freq = IAFslow.freq-iaf_grad;

cfg = [];
cfg.frequency = [max_minf min_maxf];
IAFfast = ft_selectdata(cfg,IAFfast);
IAFslow = ft_selectdata(cfg,IAFslow);

iaf_fast = squeeze(IAFfast.powspctrm);
iaf_slow = squeeze(IAFslow.powspctrm);

freqvec = IAFfast.freq;

mkdir(fullfile(outpth,subj{s}))
condname = strjoin(condi,'_');
save(fullfile(outpth,subj{s},[condname(1:end-1),ta_tp{ti},'_RTsplit.mat']),'iaf_fast','iaf_slow','freqvec')

end

end
