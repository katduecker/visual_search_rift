%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition

%% GLM RIFT ~ alpha power
% Fit a GLM to the single-trial RIFT data

% (c), Katharina Duecker
% last edited, 5 May 2026

% Inputs
% -s: subject ID
% -which_set: '': all trials, 'set32': only set size 32, 'gui': only guided
% -time_int (s): time interval for alpha power average, can be [nan, nan] if
% using time points identified by RT GLM analysis, replace nan for fixed
% time windows

function d1_RIFT_alpha_coh_distractor(s, which_set, time_init)

winl = 0.5;
% Inputs
% - s: subject index
% - which_set (str), which trials to include


pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/';
addpath(fullfile(pth, 'fieldtrip'))
ft_defaults;

load(fullfile(pth,'matlab_scripts/',"preproc_meg",'idx_subjoi.mat'));


alphapth = fullfile(pth,'results','meg','6 Alpha','pow');
soipth = fullfile(alphapth,'iaf_soi');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_rift');
mkdir(fullfile(outpth,subj{s}))
glmrtpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');

load(fullfile(glmrtpth, 'glm_RT_soi_iaf_subj'))
soi_glm = subj_soi{s};
soi_iaf = iaf_glm(s); 

cohpth = fullfile(pth,'results','meg','8 COH single trl');

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');

% load GLM results
load(fullfile(glmrtpth, 'glm_RT_soi_iaf_subj'), 'subj_soi', 'glm_time_sig', 'iaf_glm')
soi = subj_soi{s};

% if pre-defined time_oi, take that; otherwise use GLM time points
time_oi = zeros(1,2);
time_oi(isnan(time_init)) = glm_time_sig(isnan(time_init));
time_oi(~isnan(time_init)) = time_init(~isnan(time_init));
iaf = iaf_glm(s);

%% Init template structure

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

% load single RIFT
load(fullfile(cohpth,subj{s},['coh_single_trial_mcohere.mat']))

% select conditions based on "which_set"
condi = logical(condi);



% get set size 32 trials & guided (always guided because this is distractor
% suppression

if strcmp(which_set, 'set32')
    trl_oi = (condi(:,1) + condi(:,2)) == 2;
    tot = tot(trl_oi);
    cohD = cohD(trl_oi);
    condi = condi(trl_oi);
    
elseif strcmp(which_set, 'set16')
    trl_oi = (condi(:,1) + ~condi(:,2)) == 2;
    tot = tot(trl_oi);
    cohD = cohD(trl_oi);
    condi = condi(trl_oi);
else
    % select only guided
    trl_oi = condi(:,1);
    tot = tot(condi(:,1));
    cohT = cohT(condi(:,1),:);
    cohD = cohD(condi(:,1),:);
    condi = condi(condi(:,1),:);
    
    
end

tot = (tot - mean(tot))';           % demean tot


coh_D = ERP;
coh_D.avg = cohD;
coh_D.var = cohD;
coh_D.dof = ones(size(coh_D.avg));

% change layout to please fieldtrip (otherwise combineplanar won't work)
coh_D.avg = repmat(coh_D.avg,1,1,2);
coh_D.time = [1, 2];


% combine planar
cfg = [];
cfg.method = 'sum';
coh_D = ft_combineplanar(cfg,coh_D);

% select just one time point
cfg = [];
cfg.latency = [1 1];
coh_D = ft_selectdata(cfg, coh_D);

%% alpha power

% load data in same order as coherence/d1_single_trial_coh

d = dir(fullfile(inpth,subj{s}));
d = {d.name};

files = d(3:end);

load(fullfile(inpth,subj{s},files{1}));
behav_array = perf_cur;
data_load = data_trig;
% load & append data
for f = 2:length(files)
    load(fullfile(inpth,subj{s},files{f}));
    data_load = ft_appenddata([],data_load,data_trig);
    behav_array = [behav_array;perf_cur];
end

data = data_load;

clear data_load

load(fullfile(pth,'experiment','trigdef.mat'))

condi_specs = {'ti','32t','tp','6067'};

condi_check = zeros(length(behav_array),3);


%% Split data into conditions 
% create array that stores the condition keys
% store: 1: guided; 1: set 32; 1: target present; Target60 Distractor 67

for c = 1:length(condi_specs)
    trig = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x), strfind(trigdef(:,2),condi_specs{c}),'UniformOutput',false)),1));
       
    idx = ismember([behav_array{:,1}],trig);
    condi_check(:,c) = idx;
end

condi_check = logical(condi_check(trl_oi,:));

if ~isequal(condi, condi_check)
    error('order of trials not the same!')
end


%% TFR of alpha power (ensures coherence and alpha power are in the same order)
% uncombine soi to speed things up

soi_uncmb = {};
for i = 1:length(soi)
    soi_uncmb = [soi_uncmb; {soi{i}(1:7)}; {['MEG',soi{i}(9:end)]}];
end

cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = soi_uncmb;
cfg.taper = 'hanning';
cfg.foi = iaf;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:0.5;
cfg.keeptrials = 'yes';
cfg.output = 'fourier';
cfg.pad = 'nextpow2';
cfg.trials = trl_oi;

TFR = ft_freqanalysis(cfg,data);
TFR.powspctrm = abs(TFR.fourierspctrm);
TFR = rmfield(TFR,'fourierspctrm');
% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);

% select time of interest
cfg = [];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.latency = time_oi;
IAF = ft_selectdata(cfg,TFR);

alpha_pow = IAF.powspctrm;

%% GLM: Time-On-Task

const = ones(length(condi),1);  

Xtot = [const,tot];% cope
tot_contr = [0 1];


% H1: blanket inhibition
Y = coh_D.trial;

% T stats
T_tot = totGLM(Xtot,tot_contr, Y);


%% alpha model

% H1: blanket inhibition
Y = coh_D.trial;

% T stats
T_alpha = alphaGLM(Y,alpha_pow);


%% Tot and alpha
% H1: blanket
T_alpha_tot = totalphaGLM(alpha_pow,tot,Y);


%% Permutations
block_boundaries = [1];

block_size = 0;


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

for c = 1:length(condi_idx)
    if c > 1 && (condi_idx(c) ~= condi_idx(c-1) || block_size ==40)
       block_boundaries = [block_boundaries;c];
       block_size = 0;
    end
    block_size = block_size + 1;
end

block_boundaries = [block_boundaries;length(condi)+1];

num_blocks = length(block_boundaries)-1;

num_perm = 1000;
T_alpha_perm = zeros(num_perm,length(coh_D.label));

T_alpha_tot_perm = T_alpha_perm;

for n =1:num_perm
    disp(['perm: ', num2str(n)])
    
    shuffled_order = randperm(num_blocks);
    
    shuf_D = zeros(size(coh_D.trial,1),102);
    
    cnt = 1;
    for i = 1:num_blocks
        start_index = block_boundaries(shuffled_order(i));
        end_index = block_boundaries(shuffled_order(i)+1)-1;
       
        
        % cut out D
        block = coh_D.trial(start_index:end_index,:);
        shuf_D(cnt:cnt+(end_index-start_index),:) = block;
        cnt = cnt+(end_index-start_index)+1;
    end
    
    if size(coh_D.trial,1) ~= (cnt-1)
        error('number of trials not accounted for')
    end
    
    Y  = shuf_D;
    T_alpha_perm(n,:) = alphaGLM(Y,alpha_pow);
    
    T_alpha_tot_perm(n,:) = totalphaGLM(alpha_pow,tot,Y);
    
 
    
    if n>1
        assert(~all(T_alpha_perm(n,:) == 0), sprintf('Permutation %d has all-zero T-values', n))
        assert(~all(T_alpha_perm(n,:) == 0), sprintf('Permutation %d has all-zero T-values', n))
        

    end
end


T_alpha_z = (T_alpha' - mean(T_alpha_perm))./std(T_alpha_perm);

T_alpha_tot_z =(T_alpha_tot' - mean(T_alpha_tot_perm))./std(T_alpha_tot_perm);

save(fullfile(outpth,subj{s},append('glm_coh_distractor',which_set,'_start_',num2str(time_oi(1)*1000),'_end_',num2str(time_oi(2)*1000),'.mat')),'T_tot', 'T_alpha_z','T_alpha_tot_z')  

end

function T_tot = totGLM(Xtot,tot_contr, Y)

    % pseudoinverse
    Xptot = pinv(Xtot);

    
    % fit model with pseudoinverse matrix
    beta_tot = Xptot*Y;

    cope_tot = tot_contr*beta_tot;

    % residuals
    res = Y - (beta_tot'*Xtot')';
    % variance of residuals
    res_dot = diag(res'*res);

    % degrees of freedom
    dof_error = size(Y,1) - rank(Xtot);

    varres = res_dot/dof_error;

    % varcope
    residue_matrix = pinv(Xtot'*Xtot);
    varcope_tot = diag(tot_contr*residue_matrix*tot_contr')*varres;
    
    T_tot = cope_tot'./sqrt(varcope_tot);
    
end

function T_alpha = alphaGLM(Y,alpha_pow)
    const = ones(size(Y,1),1);
    T_alpha = zeros(size(Y,2),1);
    alpha_z = zscore(alpha_pow(:));


    X_alpha = [const, alpha_z];
    alpha_contr = [0 1];



    % pseudoinverse
    Xpalpha = pinv(X_alpha);
    for c = 1:size(Y,2)
       
        Yc = Y(:,c);

        % fit model with pseudoinverse matrix
        beta_alpha = Xpalpha*Yc;

        cope_alpha = alpha_contr*beta_alpha;

        % residuals
        res = Yc - (beta_alpha'*X_alpha')';
        % variance of residuals
        res_dot = diag(res'*res);

        % degrees of freedom
        dof_error = size(Yc,1) - rank(X_alpha);

        varres = res_dot/dof_error;

        % varcope
        residue_matrix = pinv(X_alpha'*X_alpha);
        varcope_alpha = diag(alpha_contr*residue_matrix*alpha_contr')*varres;

        % T stats
        T_alpha(c) = cope_alpha./sqrt(varcope_alpha);

    end
end


function T_tot_alpha = totalphaGLM(alpha_pow,tot,Y)

    const = ones(size(Y,1),1);
    T_tot_alpha = zeros(size(Y,2),1);
    alpha_z = zscore(alpha_pow(:));
    

    X_tot_alpha = [const, tot,alpha_z];
    % cope
    tot_alpha_contr = [0 0 1];        % target vs unguided contrast

    
    % pseudoinverse
    Xptotalpha = pinv(X_tot_alpha);
    
    for c = 1:size(Y,2)
        Yc = Y(:,c);

        % fit model with pseudoinverse matrix
        beta_tot_alpha = Xptotalpha*Yc;

        cope_tot_alpha = tot_alpha_contr*beta_tot_alpha;

        % residuals
        res = Yc - (beta_tot_alpha'*X_tot_alpha')';
        % variance of residuals
        res_dot = diag(res'*res);

        % degrees of freedom
        dof_error = size(Yc,1) - rank(X_tot_alpha);

        varres = res_dot/dof_error;

        % varcope
        residue_matrix = pinv(X_tot_alpha'*X_tot_alpha);
        varcope_tot_alpha = diag(tot_alpha_contr*residue_matrix*tot_alpha_contr')*varres;

        % T stats
        T_tot_alpha(c) = cope_tot_alpha'./sqrt(varcope_tot_alpha);
    end
end