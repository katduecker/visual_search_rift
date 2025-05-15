%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c2: store coherence per condtion, median split into fast and slow trials


% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Input
% - split_ta_tp: (0 or 1): split into target absent/present?
% - bslcor: baseline corrected coherence?

% Output:
% coh_condiT: cell of size number of conditions x 3, storing condition
% name, RIFT coherence to Target for fast trials, RIFT coherence to Target
% for slow trials
% coh_condiD: same for Distractor

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)

function c2_store_rift_RTsplit(split_ta_tp,bslcor)


% conditions based on set size
if split_ta_tp

    condi = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},{'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};
    split_suf = '_ta_tp';
else
    condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
    split_suf = '';
end


pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','RT');
alphapth = fullfile(pth,'results','meg','6 Alpha','not align','pow align iaf');

cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath('/rds/projects/j/jenseno-visual-search-rft/shadederror')
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/RT')

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

% parameters
toi = [-2.5,2];                 % length of trial in seconds (-2.5 to 2 s)
fs = 1000;                      % sampling rate MEG
filttype = {'but','twopass'};
% list subjects
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

% time interval of interest
timevec = linspace(toi(1),toi(2),diff(toi)*fs);      % timevector based on g_split_trl.m in preprocessing

min_rt = .5;  % minimum reaction time for each participant
min_rt_min = round(min(min_rt),3);      % minimum reaction time

% for average over time
[~,startsamp] = min(abs(timevec+0.5)); % start sample
[~,endsamp] = min(abs(timevec-min_rt_min)); % end sample

% read first file to get channel labels
cfg = [];
cfg.dataset = fullfile(maxfpth,subj{1},'part1_sss.fif');
cfg.channel = {'MEGGRAD','MISC004','MISC005'};

% load in data for this part
data = ft_preprocessing(cfg);
label_grad = data.label;
clear data



% target coherence
coh_condiT = cell(length(condi),3);
coh_condiD = cell(length(condi),3);

for c = 1:length(condi)

    % Target
    coh_condiT{c,1} = strjoin(condi{c},'_');
    % store coherence fast trials
    coh_condiT{c,2} = zeros(length(subj),endsamp-startsamp+1);
    % store coherence slow trials
    coh_condiT{c,3} = zeros(length(subj),endsamp-startsamp+1);

    % Distractor
    coh_condiD{c,1} = strjoin(condi{c},'_');
    % store coherence fast trials
    coh_condiD{c,2} = zeros(length(subj),endsamp-startsamp+1);
    % store coherence slow trials
    coh_condiD{c,3} = zeros(length(subj),endsamp-startsamp+1);
    for s = 1:length(subj)
        load(fullfile(soipth,subj{s},'soi_stat.mat'))
        label_idx = ismember(label_grad,soigrad);
        
        load(fullfile(cohpth,subj{s},strjoin({strjoin(condi{c},'_'),'RTsplit.mat'},'_')),'coh6067','coh6760')
        
        % baseline corrected?
        if bslcor
            % Target
            x1 = squeeze(mean(coh6067.fast.bslcor.cohTgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.fast.bslcor.cohTgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiT{c,2}(s,:) = (x1+x2)./2;
            
            x1 = squeeze(mean(coh6067.slow.bslcor.cohTgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.slow.bslcor.cohTgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiT{c,3}(s,:) = (x1+x2)./2;
            
            
            x1 = squeeze(mean(coh6067.fast.bslcor.cohDgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.fast.bslcor.cohDgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiD{c,2}(s,:) = (x1+x2)./2;
            
            x1 = squeeze(mean(coh6067.slow.bslcor.cohDgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.slow.bslcor.cohDgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiD{c,3}(s,:) = (x1+x2)./2;
        else
            
            % Target
            x1 = squeeze(mean(coh6067.fast.cohTgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.fast.cohTgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiT{c,2}(s,:) = (x1+x2)./2;
            
            x1 = squeeze(mean(coh6067.slow.cohTgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.slow.cohTgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiT{c,3}(s,:) = (x1+x2)./2;
            
            
            x1 = squeeze(mean(coh6067.fast.cohDgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.fast.cohDgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiD{c,2}(s,:) = (x1+x2)./2;
            
            x1 = squeeze(mean(coh6067.slow.cohDgrad(label_idx,:,startsamp:endsamp),1));
            x2 = squeeze(mean(coh6760.slow.cohDgrad(label_idx,:,startsamp:endsamp),1));
            
            coh_condiD{c,3}(s,:) = (x1+x2)./2;
        end
        
        clear coh6067 coh6760
    end
    
end

if bslcor
    save(fullfile(cohpth,['coh_RT_bslcor',split_suf,'.mat']),'coh_condiT','coh_condiD','-v7.3')

else
save(fullfile(cohpth,['coh_RT',split_suf,'.mat']),'coh_condiT','coh_condiD','-v7.3')
end
