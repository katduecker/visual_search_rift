%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition

%% GLM Spectrum analysis
% Fit a GLM to the fourier Spectrum

% (c), Katharina Duecker
% last edited, 4 May 2026

% Code based on Quinn et al., The GLM-spectrum: [...], 2024, Imaging
% Neuroscience

% GLM fitted for each participant

function c1_GLM_spectrum_interactions(s)

% Input:
% -s: subject index

winl=0.5;        % window length

%% paths
rmpath(genpath('/rds/projects/2018/jenseno-entrainment'))
addpath('/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
load(fullfile(pth,'matlab_scripts/',"preproc_meg",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha','pow');

cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');
mkdir(outpth,subj{s})

load(fullfile(pth,'matlab_scripts','coherence','occi_grad.mat'))        % load occipital sensors

%% Load trial data 


d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];

% load RT, tot and slouch (all trials)
load(fullfile(inpth,subj{s},file{1}),'behav_array','dist_z_array')

% extract specs condition
load(fullfile(pth,'experiment','trigdef.mat'))

condi_specs = {'ti','32t','tp'};

condi = zeros(length(behav_array),3);

% find trigger condition
% store: 1: guided; 1: set 32; 1: target present

for c = 1:length(condi_specs)
    trig = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x), strfind(trigdef(:,2),condi_specs{c}),'UniformOutput',false)),1));
       
    idx = ismember([behav_array{:,1}],trig);
    condi(:,c) = idx;
end

%% Code conditions and prepare factors for design matrix

% combination conditions
condi_comb = [[0,0,0];[0,0,1];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[1,1,0];[1,1,1]];

rt = [behav_array{:,3}]';


% demean tot and distance along z axis
tot = [behav_array{:,4}];
tot = tot'-mean(tot);

% change sign -> large positive should be large distance from initial
% position
dist_z_array = dist_z_array .* sign(dist_z_array);
dist_z = dist_z_array - mean(dist_z_array);

% zscore RT within condition to address multicollinearity

for c = 1:length(condi_comb)
    cur_cond = zeros(size(condi,1),1);

    for ti = 1:size(condi,1)
        cur_cond(ti) = isequal(condi(ti,:),condi_comb(c,:));
    end
    
    rt(logical(cur_cond)) = zscore(rt(logical(cur_cond)));
    
end

% distinguish between blocks for randomization
condi_comb = [[0,0];[0,1];[1,0];[1,1]];

condi_idx = zeros(length(condi),1);

for c = 1:length(condi_comb)
    cur_cond = zeros(size(condi,1),1);

    for ti = 1:size(condi,1)
        cur_cond(ti) = isequal(condi(ti,1:2),condi_comb(c,:));
        
    end
    
    condi_idx(logical(cur_cond)) = c;
    
end


%% Load TFR of fourier

load(fullfile(alphapth,subj{s},['data_fourier_winl_',num2str(winl*10),'.mat']),'TFR')

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);
cfg = [];
cfg.latency = [-1.5 .5];
IAF = ft_selectdata(cfg,TFR);

%% GLM

% prepare design matrix
X = ones(length(rt),9);

% guided
X(:,2) = condi(:,1)*2-1;
% set size 32
X(:,3) = condi(:,2)*2-1;
% target
X(:,4) = condi(:,3)*2-1;
% rt
X(:,5) = rt;

% tot
X(:,6) = tot;

% distance to parietal sensors
X(:,7) = dist_z;

% interactions

% guided*RT
X(:,8) = X(:,2).*X(:,5);

% guided * set size
X(:,9) = X(:,3).*X(:,5);

% VIF
X_vif = X(:,2:end);
X_vif = X_vif./std(X_vif, 0, 1);
R = corrcoef(X_vif);
VIF = diag(inv(R));

% sanity
p = size(X_vif,2);
VIF2 = zeros(p,1);
for j = 1:p
    others =setdiff(1:p,j);
    b=X_vif(:,others)\X_vif(:,j);
    resid=X_vif(:,j)-X_vif(:,others)*b;
    R2 = 1-sum(resid.^2)/sum(X_vif(:,j).^2);
    VIF2(j)=1/(1-R2);
end

if all(~(abs(VIF-VIF2)<1e-4)')
    error('VIF methods in disagreement!')
end


% for effect size, estimate model without any regressors related to RT;
X_no_rt = X(:,[1:4,6:7]);

% prepare empty matrices
model_beta = zeros(size(X,2),length(IAF.label),length(IAF.freq),length(IAF.time)); % regressor
model_T = model_beta; % t-values
proj_spec_max_rt = zeros(length(IAF.label),length(IAF.freq),length(IAF.time));    % GLM-projected spectrum
proj_spec_min_rt = zeros(length(IAF.label),length(IAF.freq),length(IAF.time));
CohensF_rt = proj_spec_min_rt; % effect size

varcope = model_T;    % variance of contrast

cov_forming_matrix = pinv(X'*X); % covariance forming matrix

Y = IAF.powspctrm;      % fourier spectrum

Xp = pinv(X);
Xp_no_rt = pinv(X_no_rt);
for c = 1:length(IAF.freq)*length(IAF.time)*length(IAF.label)

    Yc = Y(:,c);
    % fit model
    model_beta(:,c) = Xp*Yc;

    % sum of squared errors
    var_forming_matrix = diag(cov_forming_matrix);
    RSS = sum((Yc - (X*model_beta(:,c))).^2);
    % degrees of freedom
    dof_error = length(Y) - rank(X);

    varres = RSS/dof_error;

    % varcope
    varcope(:,c) = var_forming_matrix*varres;

    % max and min spectra for plotting
    proj_spec_min_rt(c) = min(rt)*model_beta(5,c)+model_beta(1,c);
    proj_spec_max_rt(c) = max(rt)*model_beta(5,c)+model_beta(1,c);

    %% effect size:
    
    % model without RT
    beta_nort = Xp_no_rt*Yc;
    RSS2 = sum((Yc - (X_no_rt*beta_nort)).^2);

    % Variance explained by each model
    TSS = sum((Yc - mean(Yc)).^2);        % total sum of squares
    R2_full = 1-(RSS/TSS);
    R2_nort = 1-(RSS2/TSS);

    CohensF_rt(c) = (R2_full-R2_nort)/(1-R2_full);

end

% T stats
model_T = model_beta./sqrt(varcope);


%% Null distribution

% shuffle blocks

block_boundaries = [1];

block_size = 0;

for c = 1:length(condi_idx)
    if c > 1 && (condi_idx(c) ~= condi_idx(c-1) || block_size ==40)
       block_boundaries = [block_boundaries;c];
       block_size = 0;
    end
    block_size = block_size + 1;
end

block_boundaries = [block_boundaries;length(condi)+1];


num_blocks = length(block_boundaries)-1;

% permute over blocks
num_perm = 500;
model_beta_perm = model_beta;

T_perm = zeros(num_perm,size(X,2),length(IAF.label),length(IAF.freq),length(IAF.time));

tic
for n =1:num_perm
    disp(['perm: ', num2str(n)])
    
    shuffled_order = randperm(num_blocks);
    
    shuf_Y = zeros(size(Y));
    cnt = 1;
    for i = 1:num_blocks
        start_index = block_boundaries(shuffled_order(i));
        end_index = block_boundaries(shuffled_order(i)+1)-1;
        block = Y(start_index:end_index,:,:,:);
        shuf_Y(cnt:cnt+(end_index-start_index),:,:,:) = block;
        cnt = cnt+(end_index-start_index)+1;
    end
    
    if size(Y,1) ~= (cnt-1)
        error('number of trials not accounted for')
    end
    
    for c = 1:length(IAF.freq)*length(IAF.time)*length(IAF.label)
        
        Yc = shuf_Y(:,c);
        % fit model
        model_beta_perm(:,c) = Xp*Yc;
        
        % sum of squared errors
        var_forming_matrix = diag(cov_forming_matrix);
        RSS = sum((Yc - (X*model_beta_perm(:,c))).^2);
        % degrees of freedom
        dof_error = length(Yc) - rank(X);
        
        varres = RSS/dof_error;
        
        % varcope
        varcope(:,c) = var_forming_matrix*varres;
        
        %clear RSS varres var_forming_matrix Yc
        
    end
    
    % T stats
    model_T_perm = model_beta_perm./sqrt(varcope);
    
    T_perm(n,:,:,:,:) = model_T_perm;
    
    % make sure there weren't any indexing errors
     if n>1
        assert(~any(T_perm(n,:) == 0), sprintf('Permutation %d has all-zero T-values', n))
    end
    
end
toc
z_score_T = (model_T - squeeze(mean(T_perm,1)))./squeeze(std(T_perm,1));




save(fullfile(outpth,subj{s},'glm_spec_rt_interactions_piv.mat'),'model_beta','model_T','proj_spec_max_rt','proj_spec_min_rt','CohensF_rt','z_score_T', 'T_perm', 'VIF', '-v7.3')
