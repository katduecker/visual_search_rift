%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Blanket inhibition

%% Balanced median split RIFT ~ alpha
% Cluster-based test

% (c), Katharina Duecker
% last edited, Nov-29-2024


function e2_rift_balanced_split_cluster()

% paths
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

addpath(fullfile(pth,'matlab scripts','alpha'))


inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');

alphapth = fullfile(pth,'results','meg','6 Alpha');
soipth = fullfile(alphapth,'iaf_soi');

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','balanced_split');
cohsoipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

avg_cohT_high = zeros(length(subj),4,1001);
avg_cohT_low = zeros(length(subj),4,1001);
avg_cohD_high = zeros(length(subj),4,1001);
avg_cohD_low = zeros(length(subj),4,1001);

timevec = linspace(-0.5,0.5,1001);

for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},'balanced_split_glm_-300_0_4_blocksfirws_twopass.mat'))
    
    avg_cohT_high(s,:,:) = cohT_high(:,2000:3000);
    avg_cohD_high(s,:,:) = cohD_high(:,2000:3000);
    avg_cohT_low(s,:,:) = cohT_low(:,2000:3000);
    avg_cohD_low(s,:,:) = cohD_low(:,2000:3000);
    
    clear cohT_* cohD_*
end


%% Paste into fieldtrip structure
% load in example subject to get ft structure
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
fs = 1000;
% list subjects
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

d = dir(fullfile(maxfpth,subj{1}));
fast = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),fast,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
fast = fast(idxx);
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(maxfpth,subj{1},fast{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = 'MEGGRAD';
% load in data for this part
data = ft_preprocessing(cfg);

% timelock analysis: get ERF strutcure

ERF = ft_timelockanalysis([],data);
ERF.time = timevec;
ERF.avg = zeros(204,1001);
ERF.var = ERF.var(1:204,1:length(timevec));
ERF.dof = ERF.dof(1:204,1:length(timevec));


% guided search, 32

COH_high = cell(1,length(subj));
COH_low = cell(1,length(subj));

for s = 1:length(subj)
    COH_high{s} = ERF;
    COH_high{s}.avg(1,:) = squeeze((avg_cohT_high(s,4,:,:)+avg_cohD_high(s,4,:,:)./2));
    
    COH_low{s} = ERF;
    COH_low{s}.avg(1,:) = squeeze((avg_cohT_low(s,4,:,:)+avg_cohD_low(s,4,:,:)./2));
end


%% cluster-based test
% Main Effect
cfg                  = [];
cfg.method           = 'montecarlo';
cfg.channel          = ERF.label{1};
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
cfg.minnbchan        = 1;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
cfg.latency = [0.1 0.5];

design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

% compare coherence high low

cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = -1;

stat = ft_timelockstatistics(cfg,COH_high{:},COH_low{:});

% unguided search, 32

COH_high = cell(1,length(subj));
COH_low = cell(1,length(subj));

for s = 1:length(subj)
    COH_high{s} = ERF;
    COH_high{s}.avg(1,:) = squeeze((avg_cohT_high(s,3,:,:)+avg_cohD_high(s,3,:,:))./2);
    
    COH_low{s} = ERF;
    COH_low{s}.avg(1,:) = squeeze((avg_cohT_low(s,3,:,:)+avg_cohD_low(s,3,:,:))./2);
end


% cluster-based test

cfg                  = [];
cfg.method           = 'montecarlo';
cfg.channel          = ERF.label{1};
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; 
% sample-specific t-value
cfg.minnbchan        = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
cfg.latency = [0.1 0.5];

design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

% compare coherence high low

cfg.tail             = 0;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = 0;

%statung = ft_timelockstatistics(cfg,COH_high{:},COH_low{:});

%% Compare individually for Target and Distractor
% not included in paper

% Target

COHT_high = cell(1,length(subj));
COHT_low = cell(1,length(subj));
cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = -1;
cfg.alpha = 0.05;
for s = 1:length(subj)
    COHT_high{s} = ERF;
    COHT_high{s}.avg(1,:) = avg_cohT_high(s,4,:,:);
    
    COHT_low{s} = ERF;
    COHT_low{s}.avg(1,:) = avg_cohT_low(s,4,:,:);
end

statT = ft_timelockstatistics(cfg,COHT_high{:},COHT_low{:});

cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = -1;
% Distractor
COHD_high = cell(1,length(subj));
COHD_low = cell(1,length(subj));

for s = 1:length(subj)
    COHD_high{s} = ERF;
    COHD_high{s}.avg(1,:) = avg_cohD_high(s,4,:,:);
    
    COHD_low{s} = ERF;
    COHD_low{s}.avg(1,:) = avg_cohD_low(s,4,:,:);
end

statD = ft_timelockstatistics(cfg,COHD_high{:},COHD_low{:});

save(fullfile(cohpth,'RIFT_balanced_split_glm_H1.mat'),'stat','avg_cohT_high','avg_cohT_low','avg_cohD_high','avg_cohD_low','statT','statD')
%% High minus low

% guided search, 32

COH_high = cell(1,length(subj));
COH_low = cell(1,length(subj));

for s = 1:length(subj)
    COH_high{s} = ERF;
    COH_high{s}.avg(1,:) = squeeze((avg_cohT_high(s,4,:,:)-avg_cohD_high(s,4,:,:)));
    
    COH_low{s} = ERF;
    COH_low{s}.avg(1,:) = squeeze((avg_cohT_low(s,4,:,:)-avg_cohD_low(s,4,:,:)));
end

cfg                  = [];
cfg.method           = 'montecarlo';
cfg.channel          = ERF.label{1};
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; 
% sample-specific t-value
cfg.minnbchan        = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 5000;
cfg.latency = [0.1 0.5];

design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

% compare coherence high low

cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = 1;


stat = ft_timelockstatistics(cfg,COH_high{:},COH_low{:});

save(fullfile(cohpth,'RIFT_balanced_split_glm_H2.mat'),'stat','avg_cohT_high','avg_cohT_low','avg_cohD_high','avg_cohD_low')
