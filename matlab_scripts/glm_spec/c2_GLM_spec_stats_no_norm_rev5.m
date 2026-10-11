%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition

%% GLM Spectrum analysis
% Test the GLM RT coefficients in occipital sensors against zero

% (c), Katharina Duecker
% last edited, Oct-05-2025

% Cluster-based permutation test

end_time = 0.4;
suf = '_fourier_piv';
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


SPECall = cell(1,length(subj));

max_T = zeros(size(subj));
max_T_perm = zeros(size(subj));
z_T = zeros(size(subj));
for s = 1:length(subj)

    load(fullfile(outpth,subj{s},['glm_spec_rt',suf,'.mat']),'model_T')
    
    T_rt = squeeze(model_T(5,:,:,:));

    SPEC.powspctrm(:,:,:) =T_rt;
    
    SPECall{s} = SPEC;
     
end

%% Cluster based test

% data format
cfg = [];
cfg.keepindividual = 'yes';
cfg.parameter = 'powspctrm';

grand_freq = ft_freqgrandaverage(cfg,SPECall{:});

% copy in grad struct from one subject
grand_freq.grad = SPECall{15}.grad;

% H0: z < 0
null_hyp = grand_freq;
null_hyp.powspctrm = zeros(size(null_hyp.powspctrm));


% Neighbours
cfg = [];
cfg.feedback = 'no';
cfg.method = 'template';
cfg.template = 'neuromag306cmb_neighb.mat';
neighbours = ft_prepare_neighbours(cfg,grand_freq);


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

stat_all = ft_freqstatistics(cfg,grand_freq,null_hyp);

cfg.channel = grand_freq.label(logical(occi_grad));

stat_occi_test = ft_freqstatistics(cfg,grand_freq,null_hyp);

% fill in stat_occi -> mask with only tested channels
stat_occi = stat_all;
stat_occi.prob = ones(size(stat_occi.prob));
stat_occi.negclusters = stat_occi_test.negclusters;
stat_occi.negclusterslabelmat = zeros(size(stat_occi.negclusterslabelmat));
stat_occi.mask = zeros(size(stat_occi.mask));
stat_occi.stat = zeros(size(stat_occi.stat));
stat_occi.ref = zeros(size(stat_occi.stat));

stat_occi.prob(logical(occi_grad),:,:) = stat_occi_test.prob;
stat_occi.negclusterslabelmat(logical(occi_grad),:,:) = stat_occi_test.negclusterslabelmat;
stat_occi.mask(logical(occi_grad),:,:) = stat_occi_test.mask;
stat_occi.stat(logical(occi_grad),:,:) = stat_occi_test.stat;
stat_occi.ref(logical(occi_grad),:,:) = stat_occi_test.ref;

%% Find negative clusters
addpath(fullfile(pth,'matlab_scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);
if ~isempty(stat_occi.negclusters)
    neg_cluster_pvals = [stat_occi.negclusters(:).prob];
    neg_cluster_id = find(neg_cluster_pvals < 0.05);
    neg_pos     = ismember(stat_occi.negclusterslabelmat, neg_cluster_id);
    neg_cluster = sum(neg_pos(:));
else
    neg_cluster = 0;
end
% Plots
t_idx = find(logical(sum(sum(neg_pos,2))));
time_oi = stat_occi.time(t_idx);
time_oi_all = [-0.9 0];

frep_idx = logical(sum(squeeze(sum(neg_pos(:,:,t_idx),3))));
f_rep = stat_occi.freq(frep_idx);
chan_rep = stat_occi.label(logical(sum(sum(neg_pos(:,frep_idx,t_idx),3),2)));



stat_T = mean(stat_occi.stat(ismember(stat_occi.label,chan_rep),:,:));
mask = squeeze(logical(mean(stat_occi.mask(ismember(stat_occi.label,chan_rep),:,:))));

mask_line_vert = diff(mask);
mask_line_vert = [mask_line_vert;zeros(1,size(mask_line_vert,2))];

mask_line_horz = diff(mask, [], 2);
mask_line_horz = [zeros(size(mask_line_horz,1),1),mask_line_horz];

fig = figure('Position',[0 0 1940/2 1080/3]);
subplot(141)
cfg = [];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [-1, .4];
cfg.ylim = [10 20];
cfg.highlightchannel = chan_rep;
cfg.zlim = 'maxabs';
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg,stat_occi)
colormap(cm)
title('')


subplot(1,4,2:4)
imagesc(stat_occi.time,stat_occi.freq,squeeze(stat_T))
hold on

for i = 1:size(mask_line_vert,1)
    for m = 2:size(mask_line_vert,2)
        if mask_line_vert(i,m)
            plot(stat_occi.time(m-1:m)+0.025, ones(1,2)*stat_occi.freq(i)+0.45, 'Color', 'k', 'Linewidth',1.5)
        end
    end
end

for i = 2:size(mask_line_horz,1)
    for m = 2:size(mask_line_horz,2)
        if mask_line_horz(i,m)
            plot(stat_occi.time(m-1)*ones(1,2)+0.025, stat_occi.freq(i-1:i)+0.45, 'Color', 'k', 'Linewidth',1.5)
        end
    end
end

axis xy
xlabel('time (s)')
xticks(-1:0.5:0)
yticks(10:10:30)
ylabel('frequency (Hz)')
title('')
cb = colorbar;
cb.Limits = [-max(abs(stat_T(:))) max(abs(stat_T(:)))];
caxis(round(cb.Limits,1))
cb.Ticks = cb.Limits(1):cb.Limits(end):cb.Limits(end);
cb.Label.String = 'T-value';
box off

print(fig,fullfile(plotpth,['main_fig_rt_regressor_outline',suf,'_no_norm']),'-dsvg')

