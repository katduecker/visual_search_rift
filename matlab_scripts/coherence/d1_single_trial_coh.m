%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 1: Guided Visual Search is associated with a feature-based priority
% map in early visual cortex

%% Single trial coherence
% Estimate single-trial coherence using Welch's method (mscohere

% (c), Katharina Duecker
% last edited, Nov-29-2024

function d1_single_trial_coh(s)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
outpth = fullfile(pth,'results','meg','5 COH hilb','coh','conditions');
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
cohsoipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid/');
cohpth = fullfile(pth,'results','meg','8 COH single trl');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

% list files
d = dir(fullfile(inpth,subj{s}));
d = {d.name};

files = d(3:end);

load(fullfile(inpth,subj{s},files{1}));
behav_array = perf_cur;
data_load = data_trig;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,subj{s},files{f}));
    data_load = ft_appenddata([],data_load,data_trig);
    behav_array = [behav_array;perf_cur];
end

data = data_load;

clear data_load

tot = [behav_array{:,4}];
rt = [behav_array{:,3}];

load(fullfile(pth,'experiment','trigdef.mat'))

condi_specs = {'ti','32t','tp','6067'};

condi = zeros(length(behav_array),3);


% find trigger condition
% store: 1: guided; 1: set 32; 1: target present; Target60 Distractor 67

for c = 1:length(condi_specs)
    trig = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x), strfind(trigdef(:,2),condi_specs{c}),'UniformOutput',false)),1));
       
    idx = ismember([behav_array{:,1}],trig);
    condi(:,c) = idx;
end


cohT_bsl = zeros(length(condi), 204);
cohD_bsl = cohT_bsl;
cohT = cohT_bsl;
cohD = cohD_bsl;

% cfg = [];
% cfg.hilbert     = 'complex';
% cfg.keeptrials  = 'yes';
% % bp filter
% cfg.bpfilter = 'yes';
% cfg.bpfreq = [60-3.5 60+3.5];
% cfg.bpfilttype = 'firws';
% cfg.bpfiltdir = 'twopass';
% cfg.binstabilityfix = 'reduce';
% 
% cfg.channel = {'MEGGRAD'};
% data60 = ft_preprocessing(cfg,data);
% cfg.bpfreq = [67-3.5 67+3.5];
% data67 = ft_preprocessing(cfg,data);
% 
% cfg.channel = 'diode T';
% cfg.bpfreq = [60-3.5 60+3.5];
% 
% diodeT60 = ft_preprocessing(cfg,data);
% 
% cfg.bpfreq = [67-3.5 67+3.5];
% 
% diodeT67 = ft_preprocessing(cfg,data);
% 
% cfg.channel = 'diode D';
% cfg.bpfreq = [60-3.5 60+3.5];
% diodeD60 = ft_preprocessing(cfg,data);
% 
% cfg.bpfreq = [67-3.5 67+3.5];
% diodeD67 = ft_preprocessing(cfg,data);
% 
% 
% 
% cfg = [];
% cfg.latency = [0.1, 0.5];
% data60 = ft_selectdata(cfg, data60);
% data67 = ft_selectdata(cfg, data67);
% 
% diodeT60 = ft_selectdata(cfg, diodeT60);
% diodeT67 = ft_selectdata(cfg, diodeT67);
% 
% diodeD60 = ft_selectdata(cfg, diodeD60);
% diodeD67 = ft_selectdata(cfg, diodeD67);



cfg = [];
%cfg.hilbert     = 'complex';
cfg.keeptrials  = 'yes';
% bp filter
cfg.bpfilter = 'yes';
cfg.bpfreq = [30 80];
% cfg.bpfilttype = 'firws';
cfg.bpfiltdir = 'twopass';
%cfg.binstabilityfix = 'reduce';

cfg.channel = {'MEGGRAD', 'diode T', 'diode D'};
data = ft_preprocessing(cfg,data);
cfg.channel = 'diode T';
diodeT = ft_preprocessing(cfg, data);
cfg.channel = 'diode D';
diodeD = ft_preprocessing(cfg, data);
% clear data


cfg = [];
cfg.latency = [0.2, 0.5];
cfg.channel = 'MEGGRAD';
data_filt = ft_selectdata(cfg, data);
cfg.channel = 'diode T';
diodeT = ft_selectdata(cfg, diodeT);
cfg.channel = 'diode D';
diodeD = ft_selectdata(cfg, diodeD);
cfg.latency = [-0.5, -0.2];
cfg.channel = 'MEGGRAD';
data_bsl = ft_selectdata(cfg, data);

%% Coherence
Fs = 1000;
window_length = 100;
noverlap = 75;
nfft = 512;

hann_window = hann(window_length);

trial_idx6067 = find(condi(:,4));
clear cohT
for idx = 1:length(trial_idx6067)
    
    % Target
    t = trial_idx6067(idx);
    meg = data_filt.trial{t};
    meg_bsl = data_bsl.trial{t};
    diode = repmat(diodeT.trial{t},length(data_filt.label),1);
    
    [coh, f] = mscohere(meg', diode', hann_window, noverlap, nfft, Fs);
    [coh_bsl, f] = mscohere(meg_bsl', diode', hann_window, noverlap, nfft, Fs);
    [~, fi] = min(abs(f-60));
    cohT(t,:) = coh(fi,:);
    cohT_bsl(t,:) = coh(fi,:)-coh_bsl(fi,:);
    
    % Distractor
    diode = repmat(diodeD.trial{t},length(data_filt.label),1);
    
    [coh, f] = mscohere(meg', diode', hann_window, noverlap, nfft, Fs);
    [coh_bsl, f] = mscohere(meg_bsl', diode', hann_window, noverlap, nfft, Fs);

    [~, fi] = min(abs(f-67));
    cohD(t,:) = coh(fi,:);
    cohD_bsl(t,:) = coh(fi,:)-coh_bsl(fi,:);
        
end

trial_idx6760 = find(~condi(:,4));

for idx = 1:length(trial_idx6760)
    
    % Target
    t = trial_idx6760(idx);
    meg = data_filt.trial{t};
    meg_bsl = data_bsl.trial{t};
    diode = repmat(diodeT.trial{t},length(data_filt.label),1);
    
    [coh, f] = mscohere(meg', diode', hann_window, noverlap, nfft, Fs);
    [coh_bsl, f] = mscohere(meg_bsl', diode', hann_window, noverlap, nfft, Fs);
    [~, fi] = min(abs(f-67));
    cohT(t,:) = coh(fi,:);
    cohT_bsl(t,:) = coh(fi,:)-coh_bsl(fi,:);
    
    % Distractor
    diode = repmat(diodeD.trial{t},length(data_filt.label),1);
    
    [coh, f] = mscohere(meg', diode', hann_window, noverlap, nfft, Fs);
    [coh_bsl, f] = mscohere(meg_bsl', diode', hann_window, noverlap, nfft, Fs);
    [~, fi] = min(abs(f-60));
    cohD(t,:) = coh(fi,:);
    cohD_bsl(t,:) = coh(fi,:)-coh_bsl(fi,:);
end
% % % 
% load(fullfile(cohsoipth, subj{s}, 'soi_stat.mat'))
% idx = ismember(data_filt.label, soigrad);
% close all
% plot(1:204, squeeze(mean(cohT_bsl)), 'o')
% hold on
% plot(find(idx), squeeze(mean(cohT_bsl(:,idx))), 'o')


save(fullfile(cohpth,subj{s},'coh_single_trial_mcohere.mat'),'cohT','cohD','cohT_bsl', 'cohD_bsl','condi','tot', 'rt')
