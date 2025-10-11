%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b3. cluster based stats on the results of b1/b2; and plotting

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)

clear all; close all; clc; beep off;

fwdth = 3.5;
filttype = {'firws', 'twopass'};
pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','conditions', 'no noise diode');
alphapth = fullfile(pth,'results','meg','6 Alpha','pow','alpha RT');
cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');
mkdir(cohfigpth)
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter/');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

addpath(fullfile(pth,'Violinplot-Matlab-master'))
ft_defaults;

set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Arial')
col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;

addpath(fullfile(pth,'matlab_scripts','RT'))


% load coherence over time per condition
load(fullfile(cohpth,['cohspct_subj_fwdth_',num2str(fwdth),strjoin(filttype,'_'),'_bslcor.mat']))



fs = 1000;                                                      % sampling rate
min_rt = size(coh_subj_ni32_T,2)/fs-0.5;                        % minimum RT
timevec = linspace(-0.5,min_rt,size(coh_subj_ni32_T,2));        % time vector

% list subjects
load(fullfile(pth,'matlab_scripts/',"preproc_meg/",'idx_subjoi.mat'));


%% load in example subject to get ft structure
d = dir(fullfile(maxfpth,subj{1}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(maxfpth,subj{1},f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = 'MEGGRAD';
% load in data for this part
data = ft_preprocessing(cfg);

% timelock analysis: get ERF strutcure

ERF = ft_timelockanalysis([],data);
ERF.time = timevec;
ERF.avg = repmat(timevec,204,1);
ERF.var = ERF.var(1:204,1:length(timevec));
ERF.dof = ERF.dof(1:204,1:length(timevec));
clear d f idx idxx trlstruct


% save data for data sharing
id = [1:length(subj)]';
variablenames = [{'id'}, arrayfun(@num2str, timevec, 'UniformOutput', false)];

% unguided 16 Target
T = array2table([id, coh_subj_ni16_T], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2d_coherence_set16.xlsx'), 'Sheet', 1)

% guided 16 Target
T = array2table([id, coh_subj_ti16_T], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2d_coherence_set16.xlsx'), 'Sheet', 2)

% unguided 16 Distractor
T = array2table([id, coh_subj_ni16_D], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2d_coherence_set16.xlsx'), 'Sheet', 3)

% guided 16 Distractor
T = array2table([id, coh_subj_ti16_D], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2d_coherence_set16.xlsx'), 'Sheet', 4)


% unguided 32 Target
T = array2table([id, coh_subj_ni32_T], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2e_coherence_set32.xlsx'), 'Sheet', 1)

% guided 32 Target
T = array2table([id, coh_subj_ti32_T], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2e_coherence_set32.xlsx'), 'Sheet', 2)

% unguided 16 Distractor
T = array2table([id, coh_subj_ni32_D], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2e_coherence_set32.xlsx'), 'Sheet', 3)

% guided 16 Distractor
T = array2table([id, coh_subj_ti32_D], 'VariableNames', variablenames);
writetable(T,fullfile(cohpth,'fig2e_coherence_set32.xlsx'), 'Sheet', 4)

%% set size 16

% Prepare cell arrays for Monte Carlo test

T_gui = cell(1,length(subj));
D_gui = cell(1,length(subj));
ung = cell(1,length(subj));



for s = 1:length(subj)


    % guided
    T = coh_subj_ti16_T(s,:);
    D = coh_subj_ti16_D(s,:);

    ERF.avg(1,:) = T;
    T_gui{s} = ERF;

    ERF.avg(1,:) = D;
    D_gui{s} = ERF;

    % unguided
    U = (coh_subj_ni16_T(s,:) + coh_subj_ni16_D(s,:))./2;

    ERF.avg(1,:) = U;
    ung{s}  = ERF;


end

% compare target - unguided & distractor - unguided

cfg         = [];

cfg.method           = 'montecarlo';                                            % cluster-based
cfg.channel          = ERF.label{1};
cfg.statistic        = 'depsamplesT';                                           % within subject
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution   
cfg.minnbchan        = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
cfg.latency = [0.1 min_rt];

Nsubj  = length(subj);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

% compare Target and Distractor fast vs slow

% set size 16
cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = 1;

[stat_16_T] = ft_timelockstatistics(cfg, T_gui{:}, ung{:});

cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = -1;
[stat_16_D] = ft_timelockstatistics(cfg, D_gui{:}, ung{:});




%% repeat for set size 32

T_gui = cell(1,length(subj));
D_gui = cell(1,length(subj));
ung = cell(1,length(subj));



for s = 1:length(subj)

    % guided
    T = coh_subj_ti32_T(s,:);
    D = coh_subj_ti32_D(s,:);

    ERF.avg(1,:) = T;
    T_gui{s} = ERF;

    ERF.avg(1,:) = D;
    D_gui{s} = ERF;

    % unguided
    U = (coh_subj_ni32_T(s,:) + coh_subj_ni32_D(s,:))./2;

    ERF.avg(1,:) = U;
    ung{s}  = ERF;



end

cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = 1;
[stat_32_T] = ft_timelockstatistics(cfg, T_gui{:}, ung{:});
cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
cfg.clustertail = -1;
[stat_32_D] =  ft_timelockstatistics(cfg, D_gui{:}, ung{:});


save(fullfile(cohpth,'coh_cluster_Tboost_Dsupp_new.mat'),'stat_16_T','stat_16_D','stat_32_T','stat_32_D')


%% plot

timevec_stats = cfg.latency(1):1/fs:cfg.latency(2);

GA_ti16_T = mean(coh_subj_ti16_T,1);
GA_ti16_D = mean(coh_subj_ti16_D,1);

GA_ni16_T = mean(coh_subj_ni16_T,1);
GA_ni16_D = mean(coh_subj_ni16_D,1);


GA_ti32_T = mean(coh_subj_ti32_T,1);
GA_ti32_D = mean(coh_subj_ti32_D,1);

GA_ni32_T = mean(coh_subj_ni32_T,1);
GA_ni32_D = mean(coh_subj_ni32_D,1);

GA_ni32_avg =  mean((coh_subj_ni32_T+coh_subj_ni32_D)./2,1);


addpath(fullfile(pth,'raacampbell-shadedErrorBar-19cf3fe/'))
fig = figure('Position',[0 0  1980/1.5 1080/3]);

subplot(121);
plot(timevec,GA_ti16_T,'Color',col_palette(1,:),'LineWidth',2)
hold on
plot(timevec,GA_ti16_D,'Color',col_palette(3,:),'LineWidth',2)
plot(timevec,GA_ni16_T,'Color',col_palette(5,:),'LineWidth',2)
plot(timevec,GA_ni16_D,'Color',col_palette(6,:),'LineWidth',2)
plot(timevec_stats(stat_16_T.mask),stat_16_T.mask(stat_16_T.mask).*0.035,'Color',col_palette(2,:),'LineWidth',2)
plot(timevec_stats(stat_16_D.mask),stat_16_D.mask(stat_16_D.mask).*0.0375,'Color',col_palette(4,:),'LineWidth',2)

xlabel('time (s)')
xlim([-0.1 0.5])
xticks([0:0.5:0.5])
ylim([-0.01 0.04])
yticks([0:0.02:0.04])

yticklabels({})
box off
ytickformat('%.2f')


subplot(122);
plot(timevec,GA_ti32_T,'Color',col_palette(1,:),'LineWidth',2)
hold on
plot(timevec,GA_ti32_D,'Color',col_palette(3,:),'LineWidth',2)
plot(timevec,GA_ni32_T,'Color',col_palette(5,:),'LineWidth',2)
plot(timevec,GA_ni32_D,'Color',col_palette(6,:),'LineWidth',2)
plot(timevec_stats(stat_32_T.mask),stat_32_T.mask(stat_32_T.mask).*0.035,'Color',col_palette(2,:),'LineWidth',2)
plot(timevec_stats(stat_32_D.mask),stat_32_D.mask(stat_32_D.mask).*0.0375,'Color',col_palette(4,:),'LineWidth',2)


ylabel('coherence')
xlabel('time (s)')
ylim([-0.01 0.04])
yticks([0:0.02:0.04])
% ylim([0 0.04])
% yticks([0:0.02:0.04])
xlim([-0.1 0.5])
xticks([0:0.5:0.5])

box off
ytickformat('%.2f')


print(fig,fullfile(cohfigpth,['rift_condition_fwdth_',num2str(fwdth*10),'_bslcor']),'-dpng')
print(fig,fullfile(cohfigpth,['rift_condition_fwdth_',num2str(fwdth*10),'_bslcor']),'-dsvg')



