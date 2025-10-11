%% VS RFT
% PhD project 2

% e3: matrix alpha high vs low

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 2/06/2022
% katharina.duecker@gmail.com



function e1_matrix_alpha_high_low()

%% set up paths

per_suf = '';               % peri_search suffix? (set to '_peri_search') - this should be an input but I don't think anyone will use it (no effects)
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/RT')
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions','alpha RFT');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphapth = fullfile(pth,'results','meg','6 Alpha');
alphapowpth = fullfile(alphapth,'pow');

% list subjects
d = dir(cohpth);
d = {d.name};
folds = d(strncmp(d,'202',3));

% time intervals to be averaged
load(fullfile(alphapowpth,folds{1},['ni_16t',per_suf,'.mat']),'timevec') % load timevec first subject  
[~,p1] = min(abs(timevec+1.25));
[~,p2] = min(abs(timevec+0.25));
pre_stim = p1:p2;                % samples pre-stimulus -750 to -250 ms (trials: pre-trial -1 sec, baseline -1.5)

% during trial
min_rt = kd_find_minrt(mergepth,folds); % minimum RT
min_rt = min(min_rt);
% select timepoint +250 ms to minimum RT
[~,p1] = min(abs(timevec-0.25));
[~,p2] = min(abs(timevec-min_rt));
peri_stim = p1:p2;

% for each condition, store average pre- and per-stim alpha
% size: subject x pre/peri 
ni16_high = zeros(length(folds),2);
ni16_low = zeros(length(folds),2);
ni32_high = zeros(length(folds),2);
ni32_low = zeros(length(folds),2);

ti16_high = zeros(length(folds),2);
ti16_low = zeros(length(folds),2);
ti32_high = zeros(length(folds),2);
ti32_low = zeros(length(folds),2);

for s = 1:length(folds)

    % unguided, 16
    load(fullfile(alphapowpth,folds{s},['ni_16t',per_suf,'.mat']))

    ni16_high(s,1) = mean(IAFpow_high(pre_stim));
    ni16_high(s,2) = mean(IAFpow_high(peri_stim));

    ni16_low(s,1) = mean(IAFpow_low(pre_stim));
    ni16_low(s,2) = mean(IAFpow_low(peri_stim));

    % guided, 16
    load(fullfile(alphapowpth,folds{s},['ti_16t',per_suf,'.mat']))

    ti16_high(s,1) = mean(IAFpow_high(pre_stim));
    ti16_high(s,2) = mean(IAFpow_high(peri_stim));

    ti16_low(s,1) = mean(IAFpow_low(pre_stim));
    ti16_low(s,2) = mean(IAFpow_low(peri_stim));

    % unguided, 32
    load(fullfile(alphapowpth,folds{s},['ni_32t',per_suf,'.mat']))

    ni32_high(s,1) = mean(IAFpow_high(pre_stim));
    ni32_high(s,2) = mean(IAFpow_high(peri_stim));

    ni32_low(s,1) = mean(IAFpow_low(pre_stim));
    ni32_low(s,2) = mean(IAFpow_low(peri_stim));

    % guided, 32
    load(fullfile(alphapowpth,folds{s},['ti_32t',per_suf,'.mat']))

    ti32_high(s,1) = mean(IAFpow_high(pre_stim));
    ti32_high(s,2) = mean(IAFpow_high(peri_stim));

    ti32_low(s,1) = mean(IAFpow_low(pre_stim));
    ti32_low(s,2) = mean(IAFpow_low(peri_stim));
end

save(fullfile(alphapowpth,['iaf_pow_high_low',per_suf,'.mat']),'ni16_high','ni16_low','ti16_low','ti16_high','ni32_low','ni32_high','ti32_low','ti32_high')