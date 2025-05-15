%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a1. coherence for each subject, averaged over trials (-> Supplementary
% Fig. 1)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023


% Input
% - s: subject index
% - foi: frequency range for coherence ([50:75])
% - fwdth: width of bandpass filter (5 Hz works well)
% - filttype: {'but','twopass'} (two-pass butterworth filter)
% - senstype: 'MEGGRAD' or 'MEGMAG'

% Output
% - coh60; coh67: coherence over all trials at 60 and 67 Hz

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)


function a1_rift_SNR(s,foi,fwdth,filttype,senstype)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','5 COH hilb','SNR');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

folds = subj;
clear subj

% list condition files that contain 'condi'
d = dir(fullfile(inpth,folds{s}));
files = {d.name};

condi = {'6067','6760'};

for c = 1:length(condi)
    c_files = files(cell2mat(cellfun(@(x) ~isempty(x), regexp(files,condi{c}),'UniformOutput',false)));
    
    load(fullfile(inpth,folds{s},c_files{1}));
    data_load = data_trig;
    % load & append data
    for f = 2:length(c_files)
        load(fullfile(inpth,folds{s},c_files{f}));
        data_load = ft_appenddata([],data_load,data_trig);
    end
    
    data = data_load;
    
    clear data_load
    
    %% Coherence
    
    % calculate coherence

    % gradiometers
    [cohT(c,:,:,:),~, ~, ~] = kd_coh_hilb_fun(data,'diode T', {senstype},foi, fwdth,filttype);
    % coherence distractor
    [cohD(c,:,:,:), ~, ~, ~] = kd_coh_hilb_fun(data,'diode D', {senstype}, foi, fwdth,filttype);
    
    clear data
end


mkdir(fullfile(outpth,folds{s}))

% average coherence over Target and Distractors
coh60 = squeeze(((cohT(1,:,:,:)+cohD(2,:,:,:))./2));            % coherence with 60 Hz diode
coh67 = squeeze(((cohT(2,:,:,:)+cohD(1,:,:,:))./2));            % coherence with 67 Hz diode



save(fullfile(outpth,folds{s},['SNR_',senstype,'_fwdth_',num2str(fwdth),'_filttype_', filttype{1},'.mat']),'coh60','coh67')