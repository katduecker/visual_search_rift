%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e2. Manually identify bad components for each subject

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023


% Output
% - adds variable "badcomp" to [subj_id]_ica.mat

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions

clear all; close all; clc; beep off
% define paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter/');
% ICA projections
icapth = fullfile(pth,'results','meg', '3 ICA', '1 all subj');
addpath('/rds/projects/2018/jenseno-entrainment/fieldtrip')            % fieldtrip
ft_defaults;

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

for s = 1:length(subjfolds)
    subj_id = subjfolds{s};
    load(fullfile(icapth,[subj_id,'_ica.mat']))
    
    cfg = [];
    cfg.layout = 'neuromag306all.lay'; % specify the layout file that should be used for plotting
    cfg.viewmode = 'component';
    ft_databrowser(cfg, dataICA)
    shg
    pause
    q = input('continue? y/n','s');
    close all
    while ~strcmp(q,'y')
        q = input('continue? y/n','s');
        pause
    end
%     
    badcomps = input('bad components? []');
    save(fullfile(icapth,[subj_id,'_ica.mat']), 'badcomps','-append')

    clear dataICA
    close all
end