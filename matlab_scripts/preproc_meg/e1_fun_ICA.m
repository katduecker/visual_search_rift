%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e. Run ICA

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Input:
% - s: subject index

% Output
% - result of ft_componentanalysis for each subject (saved as matfile)

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions

function e1_fun_ICA(s)

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path
dtpth = fullfile(pth, 'data'); % data path
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))
ft_defaults;
fs = 1000;
dtpth = fullfile(pth,'data');                                       % raw data
maxfpth = fullfile(pth,'results','meg', '1 maxfilter');             % max filter
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');   % scripts
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
icapth = fullfile(pth,'results','meg', '3 ICA', '1 all subj');
mkdir(icapth)

% list subj
d = dir(dtpth);
folds = {d.name};
subj = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% load(fullfile(mergepth,'docu_merge.mat'))
% subj2 = subj(strcmp(mergesubj(:,3),'adjusted'));
% d = dir(icapth);
% folds = {d.name};
% subj3 = folds(strncmp(folds,'202',3));
% subj3 = cellfun(@(x) x(1:13),subj3,'UniformOutput',false);
% 
% mergesubj(find(ismember(subj,setdiff(subj2,subj3))),:)
% load in trial structure
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*2.5,meginfo.alltrl_bl{p}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];
    trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;
end

% list maxfiltered data
d = dir(fullfile(maxfpth,subj{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load data
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(maxfpth,subj{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = 'MEG';
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);   
end

data = ft_appenddata([],dtprt{:});

trl_idx = 1:length(data.trial);

% discard artefactual trials
cfg = [];
cfg.trials = trl_idx(meginfo.keeptrl_all);
data = ft_selectdata(cfg,data);

% sampleinfo
trl_end_samp = cell2mat(cellfun(@length,data.trial,'UniformOutput',false)).*[1:length(data.trial)];
data.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];

% resample to 250 Hz (1/4 of sampling rate)
disp(['resample subj ', subj{s}])
cfg = [];
cfg.resamplefs = 250;
dataresamp = ft_resampledata(cfg,data);

%% ICA per subject
disp(['ICA ', subj{s}])
cfg = [];
cfg.method = 'runica';
cfg.numcomponent = 68;                      % identified based on first 20 subjects
dataICA = ft_componentanalysis(cfg,dataresamp);
save(fullfile(icapth, [subj{s},'_ica.mat']),'dataICA')