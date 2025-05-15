%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b1. coherence per condition (Fig. 3 c & d)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% b1. coherence per condition (guided, unguided; set size)


% Input:
% - s: subject index
% - c_idx: condition index
% - foi: frequencies of interest (for which coherence will be computed,
% could be [60,67] or [55:75]
% - fwdth: bandwith of bpfilter applied before Hilbert transform
% (recommended: 5 Hz)
% - filttype: {'but','twopass'}

% Output
% struct containing coherence, PSD, CSD for gradiometers

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)


function b1_rift_condition(s,c_idx,foi,fwdth,filttype)


% all conditions
condi_all = {{'ni','16t','6067'},{'ti','16t','6067'}, {'ni','32t','6067'},{'ti','32t','6067'},...
    {'ni','16t','6760'},{'ti','16t','6760'}, {'ni','32t','6760'},{'ti','32t','6760'}};
    
% select current condition
condi = condi_all{c_idx};

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
outpth = fullfile(pth,'results','meg','5 COH hilb','coh','conditions', 'no noise diode');
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid', 'no noise diode');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

% list condition files that contain 'condi'
d = dir(fullfile(inpth,subj{s}));
d = {d.name};

for c = 1:length(condi)
    cond_idx(c,:) = cell2mat(cellfun(@(x) ~isempty(x),regexp(d,condi(c)),'UniformOutput',false));
end
cond_idx = sum(cond_idx,1) == length(condi);
files = d(cond_idx);

load(fullfile(inpth,subj{s},files{1}));
data_load = data_trig;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,subj{s},files{f}));
    data_load = ft_appenddata([],data_load,data_trig);
end

data = data_load;

clear data_load


%% Coherence

% calculate coherence


% "{'MEGGRAD'} can be replaced with cell containing sensors of interest
% gradiometers
[coh.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data,'diode T', {'MEGGRAD'},foi, fwdth,filttype);
% coherence distractor
[coh.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data,'diode D', {'MEGGRAD'}, foi, fwdth,filttype);

% number of trials
coh.ntrial = length(data.trial);


% % magnetometers
% [coh.cohTmag, ~, ~, ~] = kd_coh_hilb_fun(data,'diode T', 'MEGMAG',foi, fwdth,filttype);
% % coherence distractor
% [coh.cohDmag, ~, ~, ~] = kd_coh_hilb_fun(data,'diode D', 'MEGMAG', foi, fwdth,filttype);

mkdir(fullfile(outpth,subj{s}))
save(fullfile(outpth,subj{s},strcat('coh_',strjoin(condi,'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat')),'coh')


