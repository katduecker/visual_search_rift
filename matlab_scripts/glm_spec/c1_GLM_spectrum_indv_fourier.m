%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Blanket inhibition

%% GLM Spectrum analysis
% Fit a GLM to the fourier Spectrum

% (c), Katharina Duecker
% last edited, Nov-29-2024

% Code based on Quinn et al., The GLM-spectrum: [...], 2024, Imaging
% Neuroscience

% GLM fitted for each participant

function b1_GLM_spectrum_indv_fourier(s)

% Input:
% -s: subject index

%% paths
rmpath(genpath('/rds/projects/2018/jenseno-entrainment'))
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha','pow');

cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');
mkdir(outpth,subj{s})

load(fullfile(pth,'matlab scripts','coherence','occi_grad.mat'))        % load occipital sensors

%% Load trial data 
d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];

% load RT, tot and slouch
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
condi_comb = [[0,0,0];[0,0,1];[0,1,0];[0,1,1];[1,0,0];[1,0,1];[1,1,0];[1,1,1]];

rt = [behav_array{:,3}]';


% demean tot and distance along z axis
tot = [behav_array{:,4}];
tot = tot'-mean(tot);

% change sign -> large positive should be large distance from initial
% position
dist_z_array = dist_z_array .* sign(dist_z_array);
dist_z = dist_z_array - mean(dist_z_array);

%% zscore RT within condition to address multicollinearity

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

load(fullfile(alphapth,subj{s},'data_fourier_winl_10.mat'),'TFR')

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);

% select time of interest
cfg = [];
cfg.latency = [-1.5 .5];
IAF = ft_selectdata(cfg,TFR);

%% GLM

% prepare design matrix
X = ones(length(rt),7);

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


% pseudoinverse (works well in case of multicollinearity - doesn't change
% results here
Xp = pinv(X);

% for effect size, estimate model without RT;
X_no_rt = X(:,[1:4,6:7]);
Xp_no_rt = pinv(X_no_rt);

% prepare empty matrices
model_beta = zeros(size(X,2),length(IAF.label),length(IAF.freq),length(IAF.time)); % regressor
model_T = model_beta; % t-values
proj_spec_max_rt = zeros(length(IAF.label),length(IAF.freq),length(IAF.time));    % GLM-projected spectrum
proj_spec_min_rt = zeros(length(IAF.label),length(IAF.freq),length(IAF.time));
CohensF_rt = proj_spec_min_rt; % effect size

varcope = model_T;    % variance of contrast

residue_matrix = pinv(X'*X); % residual matrix

Y = IAF.powspctrm;      % fourier

for c = 1:length(IAF.freq)*length(IAF.time)*length(IAF.label)

    Yc = zscore(Y(:,c));
    % fit model
    model_beta(:,c) = X\Yc;

    % sum of squared errors
    var_forming_matrix = diag(residue_matrix);
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
    beta_nort = X_no_rt\Yc;
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

T_perm = zeros(num_perm,length(IAF.label),length(IAF.freq),length(IAF.time));

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
        
        Yc = zscore(shuf_Y(:,c));
        % fit model
        model_beta_perm(:,c) = X\Yc;
        
        % sum of squared errors
        var_forming_matrix = diag(residue_matrix);
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
    
    T_perm(i,:,:,:) = model_T_perm(5,:,:,:);
    
end
toc
z_score_T = (model_T(5,:,:,:) - mean(T_perm))./std(T_perm);


save(fullfile(outpth,subj{s},'glm_spec_rt_fourier.mat'),'model_beta','model_T','proj_spec_max_rt','proj_spec_min_rt','CohensF_rt','z_score_T', 'T_perm', '-v7.3')
