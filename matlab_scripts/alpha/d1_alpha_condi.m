%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d1. Alpha power per condition/subject (Supplementary Fig. 5)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs
% - s: subject index
% - c_idx: condition index (it's just 1 because we are looking at all data)
% - rt_trim (bool): trim according to reaction time? (3 std above mean);
% can be 0 for this step

% Output
% - iaf_search: alpha power during search
% - iaf_bsl: alpha power pre-search

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

function d1_alpha_condi(s,c_idx,rt_trim)

condi_all = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},...
    {'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};

condi = condi_all{c_idx};

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
alphapth = fullfile(pth,'results','meg','6 Alpha');
outpth = fullfile(alphapth,'pow align iaf');
soipth = fullfile(alphapth,'iaf_soi');


addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

load(fullfile(soipth,subj{s},'iaf_soi.mat'))

%% load data

d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];



% find relevant data files
condi_files = zeros(length(files),1);
for c = 1:length(condi)

    condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c}),'UniformOutput',false))';

end

condi_files = condi_files == length(condi);

% load files
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

if rt_trim
    m_rt = mean([trl_perf{:,3}]);
    std_rt = std([trl_perf{:,3}]);

    trl_rej = [[trl_perf{:,3}] > m_rt + 3*std_rt]';
    trl_keep = (strcmp(trl_perf(:,2),'h') + ~trl_rej) == 2;
else
    trl_keep = strcmp(trl_perf(:,2),'h');
end

hit_rate = sum(strcmp(trl_perf(:,2),'h'))/size(trl_perf,1);

%% alpha in condition

winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEG';
cfg.taper = 'hanning';
cfg.foi = iaf_grad;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'no';
cfg.trials = trl_keep;
TFR_alpha = ft_freqanalysis(cfg,data_load);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);



cfg.channel = soi_grad_cmb;

cfg.avgoverchan = 'yes';

IAF = ft_selectdata(cfg,TFR_alpha);


% iaf power over time
iaf_time = IAF.powspctrm;

% average over time
cfg = [];
cfg.avgovertime = 'yes';

% baseline power
cfg.latency = [-1 0];
IAFbsl = ft_selectdata(cfg,IAF);
iaf_bsl = IAFbsl.powspctrm;

% power during search
cfg.latency = [0.25 0.5];
IAFs = ft_selectdata(cfg,IAF);
iaf_search = IAFs.powspctrm;


rt = mean([trl_perf{trl_keep,3}]);

condname = strjoin(condi,'_');
mkdir(fullfile(outpth,subj{s}))
save(fullfile(outpth,subj{s},[condname,'.mat']),'iaf_time','iaf_bsl','iaf_search','rt','hit_rate')

