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

function b2_store_rift_condition(fwdth,foi)


%% settings & paths

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','conditions', 'no noise diode');
cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
mkdir(cohfigpth)
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath('/rds/projects/j/jenseno-visual-search-rft/shadederror')
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/preprocessing MEG')

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

toi = [-2.5,2];                         % trial time window
filttype = {'firws','twopass'};

fs = 1000;                              % sampling rate

condi = {'16t','32t'};                 % conditions

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


%% prepare extraction of the data

% minimum reaction time (how much of the trial should be cut out?
min_rt = kd_find_minrt(mergepth,subj);
min_rt_min = round(min(min_rt),3);



% find frequencies of interest in foi
f1 = foi==60;
f2 = foi==67;

% fill in coherence over time for each condition
coh_subj_ni16_T = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ni16_D = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ti16_T = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ti16_D = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ni32_T = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ni32_D = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ti32_T = zeros(length(subj),(0.5+min_rt_min)*fs);
coh_subj_ti32_D = zeros(length(subj),(0.5+min_rt_min)*fs);

for s = 1:length(subj)
    
    % subject SOI
    load(fullfile(soipth,subj{s},'soi_stat.mat'))
    
    % separate magnetometers and gradiometers
    soimag = soi_stat(logical(cell2mat(cellfun(@(x) strcmp(x(end),'1'),soi_stat,'UniformOutput',false))));
    soigrad = soi_stat(~ismember(soi_stat,soimag));
    
    soi_idx = ismember(labels,soigrad);
    
    %% unguided, set size 16
    
    
    % target 60 Hz
    curcond = ['coh_',strjoin({'ni','16t','6067'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));

    % cut out from 2 seconds (0.5 ms before display onset) to minimum
    % reaction time
    
    
    % baseline correct
    bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
    coh.cohTgrad = coh.cohTgrad-bsl;
    % frequency 1=60 Hz
    
    if size(coh.cohTgrad,2) > length(foi)
        error(['old coh version subj ',num2str(s),' ', subj{s}])
    end
    avg_soi_T1 = squeeze(mean(coh.cohTgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 2=67
    
    avg_soi_D1 = squeeze(mean(coh.cohDgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));

    % target 67 Hz
    curcond = ['coh_',strjoin({'ni','16t','6760'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));
    % baseline correct
    bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
    coh.cohTgrad = coh.cohTgrad-bsl;
    % Target frequency 2=67 Hz
    avg_soi_T2 = squeeze(mean(coh.cohTgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 60 Hz
    avg_soi_D2 = squeeze(mean(coh.cohDgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % store
    
    % Target
    coh_subj_ni16_T(s,:) = (avg_soi_T1+avg_soi_T2)./2;

    
    % Distractor
    coh_subj_ni16_D(s,:) = (avg_soi_D1+avg_soi_D2)./2;
    
    clear coh avg_soi*
    
    %% guided, set size 16
    
    curcond = ['coh_',strjoin({'ti','16t','6067'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));

    % cut out from 2 seconds (0.5 ms before display onset) to minimum
    % reaction time
    % baseline correct
    bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
    coh.cohTgrad = coh.cohTgrad-bsl;
    % frequency 1=60 Hz
    avg_soi_T1 = squeeze(mean(coh.cohTgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 2=67
    avg_soi_D1 = squeeze(mean(coh.cohDgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
     % target 67 Hz
    curcond = ['coh_',strjoin({'ti','16t','6760'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));
    
%     % baseline correct
%     bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
%     coh.cohTgrad = coh.cohTgrad-bsl;
    % Target frequency 2=67 Hz
    avg_soi_T2 = squeeze(mean(coh.cohTgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 60 Hz
    avg_soi_D2 = squeeze(mean(coh.cohDgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % Target
    coh_subj_ti16_T(s,:) = (avg_soi_T1+avg_soi_T2)./2;
    
    
    % Distractor
    coh_subj_ti16_D(s,:) = (avg_soi_D1+avg_soi_D2)./2;
    clear coh avg_soi*
    
    %% unguided, set size 32
    % target 60 Hz
    curcond = ['coh_',strjoin({'ni','32t','6067'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));

%     cut out from 2 seconds (0.5 ms before display onset) to minimum
%     reaction time
    
    % baseline correct
    bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
    coh.cohTgrad = coh.cohTgrad-bsl;   
    % frequency 1=60 Hz
    avg_soi_T1 = squeeze(mean(coh.cohTgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 2=67
    avg_soi_D1 = squeeze(mean(coh.cohDgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));

    % target 67 Hz
    curcond = ['coh_',strjoin({'ni','32t','6760'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));
    % baseline correct
    bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
    coh.cohTgrad = coh.cohTgrad-bsl;
    % Target frequency 2=67 Hz
    avg_soi_T2 = squeeze(mean(coh.cohTgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 60 Hz
    avg_soi_D2 = squeeze(mean(coh.cohDgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % Target
    coh_subj_ni32_T(s,:) = (avg_soi_T1+avg_soi_T2)./2;

    
    % Distractor
    coh_subj_ni32_D(s,:) = (avg_soi_D1+avg_soi_D2)./2;
    
    clear coh avg_soi*
    
    %% guided, set size 32
    
    curcond = ['coh_',strjoin({'ti','32t','6067'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));

    % cut out from 2 seconds (0.5 ms before display onset) to minimum
    % reaction time
    % baseline correct
    bsl = mean(coh.cohTgrad(:,:,1*fs:2.5*fs),3);
    coh.cohTgrad = coh.cohTgrad-bsl;
    % frequency 1=60 Hz
    avg_soi_T1 = squeeze(mean(coh.cohTgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 2=67
    avg_soi_D1 = squeeze(mean(coh.cohDgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
     % target 67 Hz
    curcond = ['coh_',strjoin({'ti','32t','6760'},'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'.mat'];
    load(fullfile(cohpth,subj{s},curcond));

    avg_soi_T2 = squeeze(mean(coh.cohTgrad(soi_idx,f2,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % baseline correct
    bsl = mean(coh.cohDgrad(:,:,1*fs:2.5*fs),3);
    coh.cohDgrad = coh.cohDgrad-bsl;
    % Distractor frequency 60 Hz
    avg_soi_D2 = squeeze(mean(coh.cohDgrad(soi_idx,f1,2*fs:round(2.5*fs+fs*min_rt_min)-1),1));
    
    
    % Target
    coh_subj_ti32_T(s,:) = (avg_soi_T1+avg_soi_T2)./2;
    
    
    % Distractor
    coh_subj_ti32_D(s,:) = (avg_soi_D1+avg_soi_D2)./2;
    clear coh avg_soi*
end

save(fullfile(cohpth,['cohspct_subj_fwdth_',num2str(fwdth),strjoin(filttype,'_'),'_bslcor.mat']),'coh_subj_ni16_T','coh_subj_ni16_D','coh_subj_ti16_T','coh_subj_ti16_D',...
    'coh_subj_ni32_T','coh_subj_ni32_D','coh_subj_ti32_T','coh_subj_ti32_D')
