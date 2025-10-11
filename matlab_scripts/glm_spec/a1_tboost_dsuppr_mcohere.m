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

function a1_tboost_dsuppr_mcohere(s, which_set)

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

% get set size 32 trials
if strcmp(which_set, 'set32')
    tot = tot(condi(:,2));
    cohT = cohT(condi(:,2),:);
    cohD = cohD(condi(:,2),:);
    condi = condi(condi(:,2),:);
elseif strcmp(which_set, 'set16')
    tot = tot(~condi(:,2));
    cohT = cohT(~condi(:,2),:);
    cohD = cohD(~condi(:,2),:);
    condi = condi(~condi(:,2),:);
end


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

Y = [coh_T.trial(:,:,1);coh_D.trial(:,:,1)];
%Y = [coh_T.avg(:,:,1);coh_D.avg(:,:,1)]

tot = tot-mean(tot);
tot = repmat(tot',2,1);
%% Set up Design matrix

% Guided Target
gui = repmat(condi(:,1),2,1);
TvsD = [ones(length(condi),1);zeros(length(condi),1)];

guiT = gui.*TvsD;
%unguiT = ~gui.*TvsD;

% guided distractor
DvsT = [zeros(length(condi),1);ones(length(condi),1)];
guiD = gui.*DvsT;
%unguiD = ~gui.*DvsT;

X = [guiT,guiD,~gui, tot];

% pseudoinverse
Xp = pinv(X);

% cope
Tcontr = [1 0 -1, 0];        % target vs unguided contrast
Dcontr = [0 1 -1, 0];

% fit model with pseudoinverse matrix
beta = Xp*Y;

copeT = Tcontr*beta;
copeD = Dcontr*beta;

% residuals
res = Y - (beta'*X')';
% variance of residuals
res_dot = diag(res'*res);

% degrees of freedom
dof_error = length(Y) - rank(X);

varres = res_dot/dof_error;

% varcope
residue_matrix = pinv(X'*X);
varcopeT = diag(Tcontr*residue_matrix*Tcontr')*varres;
varcopeD = diag(Dcontr*residue_matrix*Dcontr')*varres;

% T stats
T = [copeT'./sqrt(varcopeT),copeD'./sqrt(varcopeD)];



coh_beta = ERP;
coh_beta.label = coh_T.label;
coh_beta.time = 1:size(X,2);
coh_beta.avg = beta;
coh_beta.dimord = 'chan_time';

coh_RIFT = coh_beta;

% Organize T-values in fieldtrip format
cfg = [];
cfg.latency = [1 2];
coh_RIFT = ft_selectdata(cfg,coh_RIFT);
coh_RIFT.avg = T; 


save(fullfile(outpth,subj{s},append('glm_coh_TvsD_',which_set,'.mat')),'coh_beta', 'coh_RIFT')
