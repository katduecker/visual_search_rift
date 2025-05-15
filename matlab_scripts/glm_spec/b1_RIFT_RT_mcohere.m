%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 1: Guided Visual Search is associated with a feature-based priority
% map in early visual cortex

%% GLM analysis
% Target boosting and distractor suppression using a single-trial measure
% of coherence based on correlations

% (c), Katharina Duecker
% last edited, Nov-29-2024

% Code based on Quinn et al., The GLM-spectrum: [...], 2024, Imaging
% Neuroscience

% Fit the GLM for each individual participant.

function d1_RIFT_RT_mcohere(s)

% Inputs:
% - s: subject index
% - set32: (bool) only use set size 32 trials

% paths
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha');
soipth = fullfile(alphapth,'iaf_soi');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm spec quinn et al','T boost D supp');
mkdir(fullfile(outpth,subj{s}))

cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');

cohpth = fullfile(pth,'results','meg','8 COH single trl');

%% Fill in template structure

datpth = fullfile(pth,'data');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat'); 
d = dir(fullfile(datpth,subj{1},'meg'));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);
fs = 1000;
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(datpth,subj{1},'meg',f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = {'MEGGRAD'};
% load in data for this part
data = ft_preprocessing(cfg);

ERP = ft_timelockanalysis([], data);

%% GLM Y
% load single trial correlations
load(fullfile(cohpth,subj{s},['coh_single_trial_mcohere.mat']))

condi = logical(condi);

coh_T = ERP;
coh_T.avg = cohT;
coh_T.var = cohT;
coh_T.dof = ones(size(coh_T.avg));

coh_D = ERP;
coh_D.avg = cohD;
coh_D.var = cohD;
coh_D.dof = ones(size(coh_D.avg));

% change layout to please fieldtrip
coh_T.avg = repmat(coh_T.avg,1,1,2);
coh_T.time = [1, 2];
coh_D.avg = repmat(coh_D.avg,1,1,2);
coh_D.time = [1, 2];



cfg = [];
cfg.method = 'sum';
coh_T = ft_combineplanar(cfg,coh_T);
coh_D = ft_combineplanar(cfg,coh_D);


% select guided only
cfg = [];
cfg.trials = condi(:,1);
cfg.latency = [1, 1];
coh_T = ft_selectdata(cfg, coh_T);
coh_D = ft_selectdata(cfg, coh_D);

Y_T = coh_T.trial;
Y_D = coh_D.trial;

rt = rt(condi(:,1)) - mean(rt(condi(:,1)));
tot = tot(condi(:,1))-mean(tot(condi(:,1)));

%% Set up Design matrix

X = [ones(length(Y_T),1), tot', rt'];
% pseudoinverse
Xp = pinv(X);

% cope
Tcontr = [0, 0, 1];        % target vs unguided contrast
Dcontr = [0, 0, 1];

% fit model with pseudoinverse matrix
betaT = Xp*Y_T;
betaD = Xp*Y_D;
copeT = Tcontr*betaT;
copeD = Dcontr*betaD;


%% Target T
% residuals
res = Y_T - (betaT'*X')';
% variance of residuals
res_dot = diag(res'*res);

% degrees of freedom
dof_error = length(Y_T) - rank(X);

varres = res_dot/dof_error;

% varcope
residue_matrix = pinv(X'*X);
varcopeT = diag(Tcontr*residue_matrix*Tcontr')*varres;


%% Distractor T
% residuals
res = Y_D - (betaD'*X')';
% variance of residuals
res_dot = diag(res'*res);

% degrees of freedom
dof_error = length(Y_D) - rank(X);

varres = res_dot/dof_error;

% varcope
residue_matrix = pinv(X'*X);
varcopeD = diag(Dcontr*residue_matrix*Dcontr')*varres;

% T stats
T = [copeT'./sqrt(varcopeT),copeD'./sqrt(varcopeD)];



coh_beta = ERP;
coh_beta.label = coh_T.label;
coh_beta.time = 1:size(X,2);
coh_beta.avg = betaT;
coh_beta.dimord = 'chan_time';

coh_RIFT = coh_beta;

% Organize T-values in fieldtrip format
cfg = [];
cfg.latency = [1 2];
coh_RIFT = ft_selectdata(cfg,coh_RIFT);
coh_RIFT.avg = T; 


save(fullfile(outpth,subj{s},append('glm_coh_RT.mat')),'coh_beta', 'coh_RIFT')
