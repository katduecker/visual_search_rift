%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a12. find occipital sensors to be used in RIFT soi algorithm (only once!)
% a1 needs to be done

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement -> just sanity check, kept trials w/ saccades
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RIFT response
% g. Split trials into conditions

clear all; close all; clc; beep off

% add fieldtrip path
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

%% folder structure
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/';

% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter');

ft_defaults;

% list subject filders
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
clear d folds

% sampling rate
fs = 1000;


%% Load data first subject
s = 1;
% load trial structure
load(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'alltrl_bl')

% trial structure to load in trl
for p = 1:length(alltrl_bl)
    trlstruct{p} = [alltrl_bl{p}(:,2)-fs,alltrl_bl{p}(:,2)+fs*2,repmat(-2.5*fs,size(alltrl_bl{p},1),1)];
end

% list files in subject folder
d = dir(fullfile(dtpth,subjfolds{s}));
f = {d.name};

% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% load in data from one file
cfg = [];
cfg.dataset = fullfile(dtpth,subjfolds{s},f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1};
cfg.channel = {'MEG','MISC004','MISC005'};
dtprt= ft_preprocessing(cfg);

% average
ERF = ft_timelockanalysis(cfg,dtprt);

cfg = [];
cfg.layout = 'neuromag306all_helmet.mat';
helm_lay = ft_prepare_layout(cfg);

% plot
cfg = [];
cfg.layout = helm_lay;
ft_multiplotER(cfg,ERF)

% occi soi -> these were selected by hand and copied and pasted
occi_soi = {'MEG2231', 'MEG2021', 'MEG2011', 'MEG1841', 'MEG1811', 'MEG2211', 'MEG2241',...
    'MEG1831', 'MEG1821', 'MEG0731', 'MEG0741', 'MEG2311', 'MEG2031', 'MEG2041', 'MEG1911',...
    'MEG1631', 'MEG2341', 'MEG2111', 'MEG1921', 'MEG1941', 'MEG2331', 'MEG2121', 'MEG1931', ...
    'MEG1731', 'MEG2131', 'MEG2141', 'MEG1741', 'MEG2232', 'MEG2022', 'MEG2012', 'MEG1842', ...
    'MEG2212', 'MEG2242', 'MEG1832', 'MEG1822', 'MEG0732', 'MEG0742', 'MEG2312', 'MEG2032', ...
    'MEG2042','MEG1912', 'MEG2322', 'MEG2342', 'MEG2112', 'MEG1922', 'MEG1942', 'MEG2512', ...
    'MEG2332','MEG2122', 'MEG1932', 'MEG1732', 'MEG2132', 'MEG2142', 'MEG1742', 'MEG2233', 'MEG2023',...
    'MEG2013', 'MEG1843', 'MEG2213', 'MEG2243', 'MEG1833', 'MEG1823', 'MEG0733', 'MEG0743',...
    'MEG2313', 'MEG2033', 'MEG2043', 'MEG1913', 'MEG2323', 'MEG2343', 'MEG2113', 'MEG1923',...
    'MEG1943', 'MEG2513', 'MEG2333', 'MEG2123', 'MEG1933', 'MEG1733', 'MEG2133', 'MEG2143',...
    'MEG1743'};

% save
save('occi_sens.mat','occi_soi')