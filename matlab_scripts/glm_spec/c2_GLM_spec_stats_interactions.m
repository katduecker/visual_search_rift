%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition

%% GLM Spectrum analysis
% Test the GLM RT coefficients in occipital sensors against zero

% (c), Katharina Duecker
% last edited, Oct-05-2025

% Cluster-based permutation test

function c2_GLM_spec_stats_interactions(end_time)

suf = '_interactions_piv';
%% paths
addpath('/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';
load(fullfile(pth,'matlab_scripts/',"preproc_meg",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');
plotpth = fullfile(outpth,'fig');
mkdir(plotpth)
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
load(fullfile(pth,'matlab_scripts/',"preproc_meg/",'idx_subjoi.mat'));

% load occipital sensors
load(fullfile(pth, 'matlab_scripts', 'preproc_meg','occi_sens.mat'))

%% load template subject

s = 1;
d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];

load(fullfile(inpth,subj{s},file{1}))

% TFR
winl=1;
cfg = [];
cfg.method = 'mtmconvol';
cfg.channel = 'MEGGRAD';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'yes';
cfg.trials = 1;

TFR_alpha = ft_freqanalysis(cfg,data_trig);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);

cfg = [];
cfg.latency = [-1.5 0.5];
cfg.avgoverrpt = 'yes';
SPEC = ft_selectdata(cfg,TFR_alpha);

clear data_trig

clear TFR_alpha

SPEC = rmfield(SPEC,'cfg');

% find the combined planars belonging to the occipital sensors
occi_grad = zeros(size(SPEC.label));

for c = 1:length(occi_soi)
    occi_grad = occi_grad + cell2mat(cellfun(@(x) ~isempty(x), regexp(SPEC.label,occi_soi{c}),'UniformOutput',false));
end

%% Get fitted GLM for each participant

% get the occipital sensors that show the effect
labels = SPEC.label(logical(occi_grad));
subj_soi = cell(1,length(subj));


Z_RT_main = cell(1,length(subj));
Z_guiRT = cell(1,length(subj));
Z_set32RT = cell(1,length(subj));

for s = 1:length(subj)

    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'z_score_T')

    % get T-value
    SPEC.powspctrm(:,:,:) = squeeze(z_score_T(5,:,:,:));
    
    Z_RT_main{s} = SPEC;
    
    SPEC.powspctrm(:,:,:) = squeeze(z_score_T(8,:,:,:));
    
    Z_guiRT{s} = SPEC;
    
    SPEC.powspctrm(:,:,:) = squeeze(z_score_T(9,:,:,:));
    
    Z_set32RT{s} = SPEC;
end

%% Cluster based test

% data format
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = 'powspctrm';

grand_freq_RTmain = ft_freqgrandaverage(cfg,Z_RT_main{:});
grand_freq_guiRT = ft_freqgrandaverage(cfg,Z_guiRT{:});
grand_freq_set32RT = ft_freqgrandaverage(cfg,Z_set32RT{:});

% copy in grad struct from one subject
grand_freq_RTmain.grad = Z_RT_main{15}.grad;

% H0: z < 0
null_hyp = grand_freq_RTmain;
null_hyp.powspctrm =  zeros(size(null_hyp.powspctrm));


% Neighbours
cfg = [];
cfg.feedback = 'no';
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';
neighbours = ft_prepare_neighbours(cfg,grand_freq_RTmain);


cfg.neighbours = neighbours;
% cluster-based one-sided test on the baseline interval (RT regressor only)
cfg.method           = 'montecarlo';
cfg.minnbchan        = 2;
cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.tail             = -1;
cfg.clustertail      = -1;
cfg.numrandomization = 5000;
cfg.alpha            = 0.05;
cfg.clusteralpha = 0.05;
cfg.latency          = [-1, end_time];
cfg.parameter = 'powspctrm';
design = ones(2,2*length(subj));
design(1,:) = [1:length(subj), 1:length(subj)];
design(2,length(subj)+1:2*length(subj)) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

stat_all = ft_freqstatistics(cfg,grand_freq_RTmain,null_hyp);

cfg.channel = grand_freq_RTmain.label(logical(occi_grad));

stat_occi_test_RT_main = ft_freqstatistics(cfg,grand_freq_RTmain,null_hyp);
stat_occi_test_guiRT = ft_freqstatistics(cfg,grand_freq_guiRT,null_hyp);
stat_occi_test_set32RT = ft_freqstatistics(cfg,grand_freq_set32RT,null_hyp);

% fill in stat_occi -> mask with only tested channels
stat_occi_RT = stat_all;
stat_occi_RT.prob = ones(size(stat_occi_RT.prob));
stat_occi_RT.negclusters = stat_occi_test_RT_main.negclusters;
stat_occi_RT.negclusterslabelmat = zeros(size(stat_occi_RT.negclusterslabelmat));
stat_occi_RT.mask = zeros(size(stat_occi_RT.mask));
stat_occi_RT.stat = zeros(size(stat_occi_RT.stat));
stat_occi_RT.ref = zeros(size(stat_occi_RT.stat));

stat_occi_RT.prob(logical(occi_grad),:,:) = stat_occi_test_RT_main.prob;
stat_occi_RT.negclusterslabelmat(logical(occi_grad),:,:) = stat_occi_test_RT_main.negclusterslabelmat;
stat_occi_RT.mask(logical(occi_grad),:,:) = stat_occi_test_RT_main.mask;
stat_occi_RT.stat(logical(occi_grad),:,:) = stat_occi_test_RT_main.stat;
stat_occi_RT.ref(logical(occi_grad),:,:) = stat_occi_test_RT_main.ref;


% fill in stat_occi -> mask with only tested channels
stat_occi_guiRT = stat_all;
stat_occi_guiRT.prob = ones(size(stat_occi_guiRT.prob));
stat_occi_guiRT.negclusters = stat_occi_test_guiRT.negclusters;
stat_occi_guiRT.negclusterslabelmat = zeros(size(stat_occi_guiRT.negclusterslabelmat));
stat_occi_guiRT.mask = zeros(size(stat_occi_guiRT.mask));
stat_occi_guiRT.stat = zeros(size(stat_occi_guiRT.stat));
stat_occi_guiRT.ref = zeros(size(stat_occi_guiRT.stat));

stat_occi_guiRT.prob(logical(occi_grad),:,:) = stat_occi_test_guiRT.prob;
stat_occi_guiRT.negclusterslabelmat(logical(occi_grad),:,:) = stat_occi_test_guiRT.negclusterslabelmat;
stat_occi_guiRT.mask(logical(occi_grad),:,:) = stat_occi_test_guiRT.mask;
stat_occi_guiRT.stat(logical(occi_grad),:,:) = stat_occi_test_guiRT.stat;
stat_occi_guiRT.ref(logical(occi_grad),:,:) = stat_occi_test_guiRT.ref;

% fill in stat_occi -> mask with only tested channels
stat_occi_set32RT = stat_all;
stat_occi_set32RT.prob = ones(size(stat_occi_set32RT.prob));
stat_occi_set32RT.negclusters = stat_occi_test_set32RT.negclusters;
stat_occi_set32RT.negclusterslabelmat = zeros(size(stat_occi_set32RT.negclusterslabelmat));
stat_occi_set32RT.mask = zeros(size(stat_occi_set32RT.mask));
stat_occi_set32RT.stat = zeros(size(stat_occi_set32RT.stat));
stat_occi_set32RT.ref = zeros(size(stat_occi_set32RT.stat));

stat_occi_set32RT.prob(logical(occi_grad),:,:) = stat_occi_test_set32RT.prob;
stat_occi_set32RT.negclusterslabelmat(logical(occi_grad),:,:) = stat_occi_test_set32RT.negclusterslabelmat;
stat_occi_set32RT.mask(logical(occi_grad),:,:) = stat_occi_test_set32RT.mask;
stat_occi_set32RT.stat(logical(occi_grad),:,:) = stat_occi_test_set32RT.stat;
stat_occi_set32RT.ref(logical(occi_grad),:,:) = stat_occi_test_set32RT.ref;

save(fullfile(outpth,['stat_GLM_spec_occi_sens',suf,'_', num2str(end_time*100), '.mat']),'stat_all','stat_occi_RT','stat_occi_set32RT','stat_occi_guiRT','occi_grad')




%% Jackknife to find onset of time and frequency

null_hyp.powspctrm = null_hyp.powspctrm(2:end,:,:,:);

interval_jack = zeros(length(subj),2);
frequency_jack = zeros(length(subj),2);

for s = 1:length(subj)
    
    idx = setdiff(1:length(subj),s);
    
    jack_spec = grand_freq_RTmain;
    jack_spec.powspctrm = grand_freq_RTmain.powspctrm(idx,:,:,:);
    
   
    
    cfg = [];
    cfg.neighbours = neighbours;
    % cluster-based one-sided test on the baseline interval (RT regressor only)
    cfg.method           = 'montecarlo';
    cfg.minnbchan        = 2;
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusterstatistic = 'maxsum';
    cfg.tail             = -1;
    cfg.clustertail      = -1;
    cfg.numrandomization = 5000;
    cfg.alpha            = 0.05;
    cfg.clusteralpha = 0.05;
    cfg.latency          = [-1.0 0.0];
    cfg.parameter = 'powspctrm';
    design = ones(2,2*length(subj)-2);
    design(1,:) = [1:length(subj)-1, 1:length(subj)-1];
    design(2,length(subj):2*length(subj)-2) = 2;

    cfg.design   = design;
    cfg.uvar     = 1;
    cfg.ivar     = 2;


    cfg.channel = grand_freq_RTmain.label(logical(occi_grad));

    stat_occi_test_RT_main = ft_freqstatistics(cfg,jack_spec,null_hyp);
    
    time_idx = find(squeeze(sum(sum(stat_occi_test_RT_main.mask,1),2)),1,'first');
    interval_jack(s,1) = stat_occi_test_RT_main.time(time_idx);
    time_idx = find(squeeze(sum(sum(stat_occi_test_RT_main.mask,1),2)),1,'last');
    interval_jack(s,2) = stat_occi_test_RT_main.time(time_idx);

    freq_idx = find(squeeze(sum(sum(stat_occi_test_RT_main.mask,1),3)),1,'first');
    frequency_jack(s,1) = stat_occi_test_RT_main.freq(freq_idx);
    freq_idx = find(squeeze(sum(sum(stat_occi_test_RT_main.mask,1),3)),1,'last');
    frequency_jack(s,2) = stat_occi_test_RT_main.freq(freq_idx);

end

save(fullfile(outpth,['stat_GLM_spec_occi_sens',suf,'_',num2str(end_time*1000),'_piv.mat']),'interval_jack','frequency_jack','-append')
