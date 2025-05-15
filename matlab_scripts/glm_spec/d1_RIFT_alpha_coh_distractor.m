%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Blanket inhibition

%% GLM RIFT ~ alpha power
% Fit a GLM to the single-trial coherence analysis

% (c), Katharina Duecker
% last edited, Nov-29-2024

% GLM fitted for each participant


function c1_RIFT_alpha_coh_distractor(s, which_set)

% Inputs
% - s: subject index
% - set32 (bool) only set size 32 trials

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha','pow');
soipth = fullfile(alphapth,'iaf_soi');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_rift');
mkdir(fullfile(outpth,subj{s}))
glmrtpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');


cohpth = fullfile(pth,'results','meg','8 COH single trl');

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');

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
    trl_oi = condi(:,2);
    tot = tot(condi(:,2));
    cohT = cohT(condi(:,2),:);
    cohD = cohD(condi(:,2),:);
    condi = condi(condi(:,2),:);
    
elseif strcmp(which_set, 'set16')
    trl_oi = ~condi(:,2);
    tot = tot(~condi(:,2));
    cohT = cohT(~condi(:,2),:);
    cohD = cohD(~condi(:,2),:);
    condi = condi(~condi(:,2),:);

elseif strcmp(which_set, 'gui')
    trl_oi = condi(:,1);
    tot = tot(condi(:,1));
    cohT = cohT(condi(:,1),:);
    cohD = cohD(condi(:,1),:);
    condi = condi(condi(:,1),:);
    
elseif strcmp(which_set, 'ungui')
    trl_oi = ~condi(:,1);
    tot = tot(~condi(:,1));
    cohT = cohT(~condi(:,1),:);
    cohD = cohD(~condi(:,1),:);
    condi = condi(~condi(:,1),:);
    
else
    trl_oi = logical(ones(length(condi),1));
    
end

tot = (tot - mean(tot))';

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

cfg = [];
cfg.latency = [1 1];
coh_T = ft_selectdata(cfg, coh_T);
coh_D = ft_selectdata(cfg, coh_D);

%% alpha power

% load data in same order as f3_single_trial_coh

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


% find trigger condition
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

load(fullfile(glmrtpth,'glm_rt_chan_fourier_noz.mat'))

% estimate TFR again because order of trials in coh_T is not consecutive!

% uncombine soi to speed things up

soi_uncmb = {};
if ~isempty(subj_soi{s})
for i = 1:length(subj_soi{s})
    soi_uncmb = [soi_uncmb; {subj_soi{s}{i}(1:7)}; {['MEG',subj_soi{s}{i}(9:end)]}];
end
else
    for i = 1:length(can_rep)
        soi_uncmb = [soi_uncmb; {chan_rep{i}(1:7)}; {['MEG',chan_rep{i}(9:end)]}];
    end
end

toi = -1.75:0.05:0.5;
winl = 1;

cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = soi_uncmb;
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = toi;
cfg.keeptrials = 'yes';
% cfg.output = 'fourier';
cfg.trials = trl_oi;

TFR = ft_freqanalysis(cfg,data);

% TFR.powspctrm = abs(TFR.fourierspctrm);
% 
% TFR = rmfield(TFR,'fourierspctrm');


% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);

% select time of interest
cfg = [];
cfg.avgoverchan = 'yes';
cfg.latency = [time_oi(1), time_oi(end)];
cfg.avgovertime = 'yes';
cfg.frequency = [f_rep(1), f_rep(end)];             
cfg.avgoverfreq = 'yes';
IAF = ft_selectdata(cfg,TFR);

alpha_pow = IAF.powspctrm;

%% time on task model

const = ones(length(condi),1);
% Guided Target
gui = condi(:,1);

Xtot = [const,gui,tot];

% distractor inhibition
Y = coh_D.trial;

% T stats
T_tot = totGLM(Xtot,Y);


%% alpha model

% H1: blanket inhibition
Y = coh_D.trial;

% T stats
T_alpha = alphaGLM(Y,alpha_pow,gui);


%% Tot and alpha
% H1: blanket
T_alpha_tot = totalphaGLM(gui,alpha_pow,tot,Y);


%% Permutations
block_boundaries = [1];

block_size = 0;

for c = 1:length(condi)
    if c > 1 && (condi(c) ~= condi(c-1) || block_size ==40)
       block_boundaries = [block_boundaries;c];
       block_size = 0;
    end
    block_size = block_size + 1;
end

block_boundaries = [block_boundaries;length(condi)];

num_blocks = length(block_boundaries)-1;

num_perm = 1000;
T_alpha_perm = zeros(num_perm,length(coh_T.label));

T_alpha_tot_perm = T_alpha_perm;

for n =1:num_perm
    disp(['perm: ', num2str(n)])
    
    shuffled_order = randperm(num_blocks);
    
    shuf_T = zeros(size(coh_T.trial,1),102);
    shuf_D = shuf_T;
    
    cnt = 1;
    for i = 1:num_blocks
        start_index = block_boundaries(shuffled_order(i));
        end_index = block_boundaries(shuffled_order(i)+1)-1;
        
        % cut out T 
        block = coh_T.trial(start_index:end_index,:);
        shuf_T(cnt:cnt+(end_index-start_index),:) = block;
        
        % cut out D
        block = coh_D.trial(start_index:end_index,:);
        shuf_D(cnt:cnt+(end_index-start_index),:) = block;
        cnt = cnt+(end_index-start_index)+1;
    end
    
    if size(coh_T.trial,1) ~= (cnt)
        error('number of trials not accounted for')
    end
    
    % H1: blanket
    Y  = (shuf_T+shuf_D)./2;
    T_alpha_perm(n,:) = alphaGLM(Y,alpha_pow,gui);
    
    T_alpha_tot_perm(n,:) = totalphaGLM(gui,alpha_pow,tot,Y);
    
    
end


T_alpha_z = (T_alpha' - mean(T_alpha_perm))./std(T_alpha_perm);

T_alpha_tot_z =(T_alpha_tot' - mean(T_alpha_tot_perm))./std(T_alpha_tot_perm);

save(fullfile(outpth,subj{s},append('glm_coh_distractor',which_set,'.mat')),'T_tot', 'T_alpha_z','T_alpha_tot_z')  

end

function T_tot = totGLM(Xtot,Y)

    % pseudoinverse
    Xptot = pinv(Xtot);

    % cope
    tot_contr = [0 0 1];        % target vs unguided contrast
    % fit model with pseudoinverse matrix
    beta_tot = Xptot*Y;

    cope_tot = tot_contr*beta_tot;

    % residuals
    res = Y - (beta_tot'*Xtot')';
    % variance of residuals
    res_dot = diag(res'*res);

    % degrees of freedom
    dof_error = length(Y) - rank(Xtot);

    varres = res_dot/dof_error;

    % varcope
    residue_matrix = pinv(Xtot'*Xtot);
    varcope_tot = diag(tot_contr*residue_matrix*tot_contr')*varres;
    
    T_tot = cope_tot'./sqrt(varcope_tot);
    
end

function T_alpha = alphaGLM(Y,alpha_pow,gui)
   const = ones(size(Y,1),1);
    T_alpha = zeros(size(Y,2),1);
    for c = 1:size(Y,2)

        alpha_z = zscore(alpha_pow(:));
        Yc = Y(:,c);
        X_alpha = [const, gui,alpha_z];

        % pseudoinverse
        Xpalpha = pinv(X_alpha);

        % cope
        alpha_contr = [0 0 1];        % target vs unguided contrast

        % fit model with pseudoinverse matrix
        beta_alpha = Xpalpha*Yc;

        cope_alpha = alpha_contr*beta_alpha;

        % residuals
        res = Yc - (beta_alpha'*X_alpha')';
        % variance of residuals
        res_dot = diag(res'*res);

        % degrees of freedom
        dof_error = length(Yc) - rank(X_alpha);

        varres = res_dot/dof_error;

        % varcope
        residue_matrix = pinv(X_alpha'*X_alpha);
        varcope_alpha = diag(alpha_contr*residue_matrix*alpha_contr')*varres;

        % T stats
        T_alpha(c) = cope_alpha./sqrt(varcope_alpha);

    end
end


function T_tot_alpha = totalphaGLM(gui,alpha_pow,tot,Y)

    const = ones(size(Y,1),1);
    T_tot_alpha = zeros(size(Y,2),1);

    for c = 1:size(Y,2)
        alpha_z = zscore(alpha_pow(:));
        Yc = Y(:,c);
        X_tot_alpha = [const, gui,tot,alpha_z];

        % pseudoinverse
        Xptotalpha = pinv(X_tot_alpha);

        % cope
        tot_alpha_contr = [0 0 0 1];        % target vs unguided contrast

        % fit model with pseudoinverse matrix
        beta_tot_alpha = Xptotalpha*Yc;

        cope_tot_alpha = tot_alpha_contr*beta_tot_alpha;

        % residuals
        res = Yc - (beta_tot_alpha'*X_tot_alpha')';
        % variance of residuals
        res_dot = diag(res'*res);

        % degrees of freedom
        dof_error = length(Y) - rank(X_tot_alpha);

        varres = res_dot/dof_error;

        % varcope
        residue_matrix = pinv(X_tot_alpha'*X_tot_alpha);
        varcope_tot_alpha = diag(tot_alpha_contr*residue_matrix*tot_alpha_contr')*varres;

        % T stats
        T_tot_alpha(c) = cope_tot_alpha'./sqrt(varcope_tot_alpha);
    end
end