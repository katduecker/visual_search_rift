%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b2. create one matrix containing data from all participants (Fig. 3 c & d)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023
 
% Input:
% - foi: frequencies of interest (for which coherence will be computed,
% could be [60,67] or [55:75]
% - fwdth: bandwith of bpfilter applied before Hilbert transform
% (recommended: 5 Hz)

% Output
% arrays with coherence for each condition, separately for target and
% distractor

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)

%% settings & paths


clear all; close all; clc;

fwdth = 5;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

set(0,'defaultAxesFontSize',16,'defaultAxesFontName','Arial')
col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;
addpath(fullfile(pth,'raacampbell-shadedErrorBar-19cf3fe/'))


cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','conditions');
cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
mkdir(cohfigpth)
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath('/rds/projects/j/jenseno-visual-search-rft/shadederror')
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/RT')

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

toi = [-2.5,2];                         % trial time window

fs = 1000;                              % sampling rate

condi = {'32t'};                 % conditions

% list subjects
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

% list maxfiltered data
d = dir(fullfile(maxfpth,subj{1}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

%% load example data to get labels
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))
% read header of subject 1 to get labels
trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)+fs*toi(1),meginfo.alltrl_bl{1}(:,3)+toi(2)*fs,zeros(length(meginfo.alltrl_bl{1}),1)+toi(1)*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(maxfpth,subj{1},f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = {'MEGGRAD','MISC004','MISC005'};
% load in data for this part
data = ft_preprocessing(cfg);

labels = data.label;
clear data trlstruct *info

load(fullfile(cohpth,['cohspct_subj_fwdth_',num2str(fwdth),'.mat']), 'coh_subj_ni32_T','coh_subj_ni32_D','coh_subj_ti32_T','coh_subj_ti32_D')

min_rt = size(coh_subj_ni32_T,2)/fs-0.5;                        % minimum RT
timevec = linspace(-0.5,min_rt,size(coh_subj_ni32_T,2));        % time vector

T_diff = zeros(size(coh_subj_ti32_T));
D_diff = T_diff;

for s = 1:length(subj)

    % guided
    T = coh_subj_ti32_T(s,:);
    D = coh_subj_ti32_D(s,:);

    % unguided
    U = (coh_subj_ni32_T(s,:) + coh_subj_ni32_D(s,:))./2;

    T_diff(s,:) = T-U;
    D_diff(s,:) = D-U;

end

[half_maxT, posT] = max(mean(T_diff));
[half_minD, posD] = min(mean(D_diff));

fig = figure;

shadedErrorBar(timevec,T_diff,{@mean, @std},'lineProps',{'Color',col_palette(1,:),'markerfacecolor',col_palette(5,:)});
hold on
shadedErrorBar(timevec,D_diff,{@mean, @std},'lineProps',{'Color',col_palette(3,:),'markerfacecolor',col_palette(6,:)});
xlabel('time (s)')
ylabel('coherence guided - unguided')
print(fig,fullfile(cohfigpth, 'latency_guided_unguided'),'-dpng','-r300')
print(fig,fullfile(cohfigpth, 'latency_guided_unguided'),'-dsvg','-r300')

