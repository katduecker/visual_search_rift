%% VS + RFT
% PhD project 2

% b. Semi-automatic artefact rejection

% [c] Katharina Duecker

% Inputs
% - s: subject id
% - pth: project path
% - fs: sampling rate

% Outputs:
% - sets meginfo.keeptrl_all = 0 for artefactual trials

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement -> just sanity check, kept trials w/ saccades
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RIFT response
% g. Split trials into conditions

clear all; close all; clc; beep off;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                 % main path
fs = 1000;

for s=1:48
    semi_auto_artef_rej(s,pth)
end

function meginfo = semi_auto_artef_rej(s,pth)
close all
dtpth = fullfile(pth,'data');                                       % raw data
maxfpth = fullfile(pth,'results','meg', '1 maxfilter');             % max filter
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');   % scripts
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
% list subj
d = dir(dtpth);
folds = {d.name};
subj = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

% load in trial structure
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))
    

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*2.5,meginfo.alltrl_bl{p}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];
end

% find trials with wrong number of samples
trlstruct_list = vertcat(trlstruct{:});

% find trials starting before recording & reject
idx_normal_startsamp = (trlstruct_list(:,1)>0);
meginfo.keeptrl_all(~idx_normal_startsamp) = 0;
trlstruct_list = trlstruct_list(idx_normal_startsamp,:);

% reject abnormal stuff
if find(~idx_normal_startsamp)
    % trial structure to load in trl
    for p = 1:length(meginfo.alltrl_bl)
        trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*2.5,meginfo.alltrl_bl{p}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];

        % getting rid of trials starting before recording (as set to 0
        % above!)
        trlstruct{p}(find(trlstruct{p}(:,1) < 0),:) = [];               % trial start < 0 samples
        meginfo.alltrl_bl{p}(find(trlstruct{p}(:,1) < 0),:) = [];
    end
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

% sampleinfo
len_trl = cell2mat(cellfun(@length,data.trial,'UniformOutput',false));

trl_end_samp = len_trl.*[1:length(data.trial)];
data.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];

%% Artefact identification

% split into grad and mag
cfg = [];
cfg.channel = 'MEGGRAD';
graddata = ft_selectdata(cfg,data);
cfg.channel = 'MEGMAG';
magdata = ft_selectdata(cfg,data);

% reject visual
cfg = [];
cfg.method = 'summary';
cfg.layout = 'neuromag306planar.lay';
cfg.latency = 'all';
grad_rej = ft_rejectvisual(cfg,graddata);

% get idx of to be kept trials
[~, keep_trl_grad] = intersect(graddata.sampleinfo(:,1),grad_rej.sampleinfo(:,1));

cfg.layout = 'neuromag306mag.lay';
mag_rej = ft_rejectvisual(cfg,magdata);
[~, keep_trl_mag] = intersect(magdata.sampleinfo(:,1),mag_rej.sampleinfo(:,1));

keep_trl_idx = intersect(keep_trl_grad,keep_trl_mag);

trl_idx = 1:length(data.trial);
trl_idx_rej = setdiff(trl_idx,keep_trl_idx);

% change keeptrl_all
meginfo.keeptrl_all(trl_idx_rej) = 0;

save(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'),'meginfo','-append')
end