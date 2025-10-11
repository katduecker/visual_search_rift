%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c3: cluster-based tests 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)


clear all; close all; clc

%% Input!

split_ta_tp = 0;            % split for target present/absent?
bslcor = 0;                 % baseline corrected coherence?
pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';
addpath(fullfile(pth,'Violinplot-Matlab-master'))
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))
ft_defaults;

set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Arial')
col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','RT');
cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig');



if split_ta_tp

    condi = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},{'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};
    split_suf = '_ta_tp';

    % if split: only consider guided search
    condi_ni = {{'ni','16ta'},{'ni','32ta'},{'ni','16tp'}, {'ni','32tp'}};
    condi_ti = {{'ti','16ta'}, {'ti','32ta'},{'ti','16tp'}, {'ti','32tp'}};

    varnames = {'gui 16 absent','gui 32 absent','gui 16 present','gui 32 present'};

else
    condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
    split_suf = '';

    varnames = {'ungu 16', 'gui 16', 'ungui 32', 'gui 32'};
end

if bslcor
    bslcor_suf = '_bslcor';
else
    bslcor_suf = '';
end



load(fullfile(cohpth,['coh_RT',bslcor_suf,split_suf,'.mat']))

% time vector
min_rt = length(coh_condiT{1,2})/1000-0.5;
timevec = linspace(-0.5,min_rt,length(coh_condiT{1,2}))

load(fullfile(pth,'matlab_scripts/',"preproc_meg/",'idx_subjoi.mat'));

id = [1:length(subj)]';
variablenames = [{'id'}, arrayfun(@num2str, timevec, 'UniformOutput', false)];

% gradnaverage
for c = 1:length(condi)
    coh_condiT{c,4} = mean(coh_condiT{c,2},1);
    coh_condiT{c,5} = mean(coh_condiT{c,3},1);

    coh_condiD{c,4} = mean(coh_condiD{c,2},1);
    coh_condiD{c,5} = mean(coh_condiD{c,3},1);
    
    % export data for sharing
    T = array2table([id, coh_condiT{c,2}], 'VariableNames', variablenames);
    writetable(T, fullfile(cohpth, 'Suppfig2_coherence_fast_slow.xlsx'), 'Sheet', [varnames{c},' target fast'])
    T = array2table([id, coh_condiT{c,3}], 'VariableNames', variablenames);
    writetable(T, fullfile(cohpth, 'Suppfig2_coherence_fast_slow.xlsx'), 'Sheet', [varnames{c},' target slow'])
    T = array2table([id, coh_condiD{c,2}], 'VariableNames', variablenames);
    writetable(T, fullfile(cohpth, 'Suppfig2_coherence_fast_slow.xlsx'), 'Sheet', [varnames{c},' distractor fast'])
    T = array2table([id, coh_condiD{c,3}], 'VariableNames', variablenames);
    writetable(T, fullfile(cohpth, 'Suppfig2_coherence_fast_slow.xlsx'), 'Sheet', [varnames{c},' distractor slow'])
end

% if split based on target present/absent only consider guided search from
% here on...
if split_ta_tp
    condi = condi_ti;
    clear condi_ti condi_ni
end

%% Cluster-based permutation t test

% load in example subject to get ft structure
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
fs = 1000;
% list subjects
d = dir(cohpth);
d = {d.name};
subj = d(strncmp(d,'202',3));
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


% store ERF structus in these cells
Tf = cell(1,length(subj));
Df = cell(1,length(subj));
Ts = cell(1,length(subj));
Ds = cell(1,length(subj));

TvDf = cell(1,length(subj));
TvDs = cell(1,length(subj));
fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

for c = 1:length(condi)

    for s = 1:length(subj)


        condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi{c},'_'));

        % fast
        T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));

        ERF.avg(1,:) = T;
        Tf{s} = ERF;

        ERF.avg(1,:) = D;
        Df{s} = ERF;


        % slow
        T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
        D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));

        ERF.avg(1,:) = T;
        Ts{s} = ERF;

        ERF.avg(1,:) = D;
        Ds{s} = ERF;


        %% fast vs slow

        cf = abs(squeeze(coh_condiT{condi_idx_ti,2}(s,:))-squeeze(coh_condiD{condi_idx_ti,2}(s,:)));
        cs = abs(squeeze(coh_condiT{condi_idx_ti,3}(s,:))-squeeze(coh_condiD{condi_idx_ti,3}(s,:)));

        ERF.avg(1,:) = cf;
        TvDf{s} = ERF;

        ERF.avg(1,:) = cs;
        TvDs{s} = ERF;

    end

    % T-test Target vs Distractor

      cfg         = [];

    cfg.method           = 'montecarlo';
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.minnbchan        = 0;
    cfg.alpha            = 0.05;
    cfg.numrandomization = 2000;
    cfg.latency = [0.2 min_rt];
    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare Target and Distractor fast vs slow

    cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail = 1;

    [statT{c}] = ft_timelockstatistics(cfg, Tf{:}, Ts{:});
    timevec_stats = statT{c}.time;

    cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail = -1;
    [statD{c}] = ft_timelockstatistics(cfg, Df{:}, Ds{:});


    cfg.statistic        = 'ft_statfun_depsamplesT';
    cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
    cfg.clustertail = 1;
    %cfg.correcttail  = 'prob';

    design = zeros(2, Nsubj*2);
     design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    [statfvd{c}] = ft_timelockstatistics(cfg, TvDf{:}, TvDs{:});

    % cohens d
    cfg         = [];

    cfg.method           = 'analytic';                                            % cluster-based
    cfg.channel          = ERF.label{1};
    cfg.statistic        = 'cohensd';                                           % within subject
    %cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
    cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution
    cfg.minnbchan        = 0;
    % cfg.neighbours       = neighbours;  % same as defined for the between-trials experiment
    cfg.alpha            = 0.05;
    cfg.numrandomization = 1000;
    cfg.latency = [-.5 min_rt];

    Nsubj  = length(subj);
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;

    % compare coherence for high and low alpha

    % Target
    cfg.clustertail = 1;
    cfg.tail        = 1;                                                       % two-sided test for both clustering and cluster-level statistic

    [statfd{c}] = ft_timelockstatistics(cfg, Tf{:}, Df{:});


    % cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
    % cfg.clustertail = -1;

    [statsd{c}] = ft_timelockstatistics(cfg, Ts{:}, Ds{:});



    %% new plot

    subplot(2,2,c);
    if strcmp(condi{c}{1},'ti')
        plot(timevec,coh_condiT{condi_idx_ti,4},'Color',col_palette(1,:),'LineWidth',3)
        hold on
        plot(timevec,coh_condiT{condi_idx_ti,5},'Color',col_palette(2,:),'LineWidth',3)
        plot(timevec,coh_condiD{condi_idx_ti,4},'Color',col_palette(3,:),'LineWidth',3)
        plot(timevec,coh_condiD{condi_idx_ti,5},'Color',col_palette(4,:),'LineWidth',3)
        
        
        plot(timevec_stats(statT{c}.mask),statT{c}.mask(statT{c}.mask).*0.0275,'Color',col_palette(5,:),'LineWidth',4)
        plot(timevec_stats(statD{c}.mask),statD{c}.mask(statD{c}.mask).*0.03,'Color',col_palette(5,:),'LineWidth',4,'LineStyle','-.')
        
        ylabel('coherence')
        xlabel('time (s)')
        xlim([-0.2 0.5])
        xticks([0:0.5:0.5])
        ylim([.00 0.05])
        yticks([0.0:0.025:0.05])
        title(varnames{c})
        box off
    else
        plot(timevec,(coh_condiT{condi_idx_ti,4}+coh_condiT{condi_idx_ti,5})./2,'Color',col_palette(5,:),'LineWidth',3)
        hold on
        plot(timevec,(coh_condiD{condi_idx_ti,4}+coh_condiD{condi_idx_ti,5})./2,'Color',col_palette(6,:),'LineWidth',3)        
        
        plot(timevec_stats(statT{c}.mask),statT{c}.mask(statT{c}.mask).*0.0275,'Color',col_palette(5,:),'LineWidth',4)
        plot(timevec_stats(statD{c}.mask),statD{c}.mask(statD{c}.mask).*0.03,'Color',col_palette(5,:),'LineWidth',4,'LineStyle','-.')
        
        ylabel('coherence')
        xlabel('time (s)')
        xlim([-0.2 0.5])
        xticks([0:0.5:0.5])
        ylim([0.00 0.05])
        yticks([0.0:0.025:0.05])
        title(varnames{c})
        box off
        
    end

end

print(fig,fullfile(cohfigpth,['rift_RT_cluster',bslcor_suf,split_suf]),'-dpng')
print(fig,fullfile(cohfigpth,['rift_RT_cluster',bslcor_suf,split_suf]),'-dsvg')

fig = figure('Position',[0 0  1980/1.5 1080/1.5]);

for c = 1:length(condi)
 % set size 16
    subplot(2,2,c);
    plot(statfd{c}.time,statfd{c}.cohensd,'Color',col_palette(5,:),'LineWidth',3)
    hold on
    plot(statsd{c}.time,statsd{c}.cohensd,'Color',col_palette(5,:),'LineWidth',3,'LineStyle','-.')
    plot(statfvd{c}.time(statfvd{c}.mask),statfvd{c}.mask(statfvd{c}.mask).*0.8,'Color',col_palette(5,:),'LineWidth',3)

    ylabel('delta RIFT Target vs Distractor Cohen''s d')
    xlabel('time (s)')
    xlim([-0.2 0.5])
    xticks([0 0.5])
    box off
    title(varnames{c})
%     ylim([-0.8 0.8])
%     yticks([-0.8:0.8:0.8])
end

if bslcor
    print(fig,fullfile(cohfigpth,['rift_RT_cohensd',split_suf,'_bslcor']),'-dpng')
    print(fig,fullfile(cohfigpth,['rift_RT_cohensd',split_suf,'_bslcor']),'-dsvg')
else
    print(fig,fullfile(cohfigpth,['rift_RT_cohensd',split_suf]),'-dpng')
    print(fig,fullfile(cohfigpth,['rift_RT_cohensd',split_suf]),'-dsvg')
end
% 
% %% unguided average
% 
% condi = {{'ni','16t'},{'ni','32t'}};
% varnames = {'ung16','ung32'};
% fig = figure('Position',[0 0  1980/1.5 1080/1.5]);
% 
% for c = 1:length(condi)
% 
%     for s = 1:length(subj)
% 
% 
%         condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi{c},'_'));
% 
%         % fast
%         T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
%         D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));
% 
%         ERF.avg(1,:) = (T+D)./2;
%         f_ung{s} = ERF;
% 
%         % slow
%         T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
%         D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));
% 
%         ERF.avg(1,:) = (T+D)./2;
%         s_ung{s} = ERF;
% 
% 
%     end
% 
%     % T-test Target vs Distractor
% 
%     cfg         = [];
% 
%     cfg.method           = 'montecarlo';
%     cfg.channel          = ERF.label{1};
%     cfg.statistic        = 'ft_statfun_depsamplesT';
%     cfg.correctm         = 'cluster';
%     cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
%     cfg.minnbchan        = 0;
%     cfg.alpha            = 0.05;
%     cfg.numrandomization = 2000;
%     cfg.latency = [0.15 min_rt];
%     Nsubj  = length(subj);
%     design = zeros(2, Nsubj*2);
%     design(1,:) = [1:Nsubj 1:Nsubj];
%     design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
% 
%     cfg.design = design;
%     cfg.uvar   = 1;
%     cfg.ivar   = 2;
% 
%     % compare Target and Distractor fast vs slow
% 
%     cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
%     cfg.clustertail =-1;
% 
%     [statf{c}] = ft_timelockstatistics(cfg, f_ung{:}, s_ung{:});
% 
% 
% 
% 
%        %% new plot
% 
% 
%     timevec_stats = timevec(650:end);
% 
%     % set size 16
%     subplot(2,2,c);
%     plot(timevec,(coh_condiT{condi_idx_ti,4}+coh_condiD{condi_idx_ti,4})./2,'Color',col_palette(5,:),'LineWidth',3)
%     hold on
%     plot(timevec,(coh_condiT{condi_idx_ti,5}+coh_condiD{condi_idx_ti,5})./2,'Color',col_palette(6,:),'LineWidth',3)
% 
%     plot(timevec_stats(statf{c}.mask),statf{c}.mask(statf{c}.mask).*0.0275,'Color',col_palette(5,:),'LineWidth',5)
% 
% 
%     ylabel('coherence')
%     xlabel('time (s)')
%     xlim([-0.2 0.5])
%     xticks([0:0.5:0.5])
%     ylim([-0.01 0.03])
%     yticks([0.0:0.03:0.03])
%     title(varnames{c})
%     box off
% 
% end
% 
% print(fig,fullfile(cohfigpth,['rift_rt_cluster_ung',bslcor_suf,split_suf]),'-dpng')
% print(fig,fullfile(cohfigpth,['rift_rt_cluster_ung',bslcor_suf,split_suf]),'-dsvg')
% % fast_TD = cell(1,length(subj));
% % slow_TD = cell(1,length(subj));
% % fig = figure('Position',[0 0  1980/1.5 1080/1.5]);
% % 
% % for c = 1:length(condi_ti)
% % 
% %     for s = 1:length(subj)
% % 
% % 
% %         condi_idx_ti = strcmp(coh_condiT(:,1),strjoin(condi_ti{c},'_'));
% % 
% %         % fast
% %         T = squeeze(coh_condiT{condi_idx_ti,2}(s,:));
% %         D = squeeze(coh_condiD{condi_idx_ti,2}(s,:));
% % 
% %         ERF.avg(1,:) = T-D;
% %         fast_TD{s} = ERF;
% % 
% %         coh_diff_fast(s,:) = T-D;
% % 
% %         % slow
% %         T = squeeze(coh_condiT{condi_idx_ti,3}(s,:));
% %         D = squeeze(coh_condiD{condi_idx_ti,3}(s,:));
% % 
% %         ERF.avg(1,:) = T-D;
% %         slow_TD{s} = ERF;
% %         coh_diff_slow(s,:) = T-D;
% % 
% %     end
% % 
% %     % T-test Target vs Distractor
% % 
% %     cfg         = [];
% % 
% %     cfg.method           = 'montecarlo';                                            % cluster-based
% %     cfg.channel          = ERF.label{1};
% %     cfg.statistic        = 'depsamplesT';                                           % within subject
% %     cfg.correctm         = 'cluster';
% %     cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
% %     cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution
% %     cfg.minnbchan        = 0;
% %     cfg.alpha            = 0.05;
% %     cfg.numrandomization = 1000;
% %     cfg.latency = [0 min_rt];
% %     Nsubj  = length(subj);
% %     design = zeros(2, Nsubj*2);
% %     design(1,:) = [1:Nsubj 1:Nsubj];
% %     design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
% % 
% %     cfg.design = design;
% %     cfg.uvar   = 1;
% %     cfg.ivar   = 2;
% % 
% %     % compare Target and Distractor fast vs slow
% % 
% %     cfg.tail             = 0;                                                       % two-sided test for both clustering and cluter-level statistic
% %     cfg.clustertail = 0;
% % 
% %     [stat_fast_slow{c}] = ft_timelockstatistics(cfg, fast_TD{:}, slow_TD{:});
% % 
% %     subplot(2,2,c)
% %     plot(timevec,mean(coh_diff_fast,1))
% %     hold on 
% %     plot(timevec,mean(coh_diff_slow,1))
% % end
% 
% % close all
% % 
% % scatter(mean(coh_diff_fast(:,500:end),2),mean(coh_diff_slow(:,500:end),2))
% % hold on
% % plot(-0.03:0.001:0.03,-0.03:0.001:0.03)
% 
% % 
% % % set size 32
% % ax2 = subplot(132);
% % plot(t_vec,GA_ti_32_fast(1,:),'Color',col_palette(1,:),'LineWidth',3)
% % hold on
% % plot(t_vec,GA_ti_32_slow(1,:),'Color',col_palette(2,:),'LineWidth',3)
% % plot(t_vec,GA_ti_32_fast(2,:),'Color',col_palette(3,:),'LineWidth',3)
% % plot(t_vec,GA_ti_32_slow(2,:),'Color',col_palette(4,:),'LineWidth',3)
% % xlabel('time (s)')
% % ylabel('coherence')
% % plot(timevec_stats(statT32.mask),statT32.mask(statT32.mask).*0.04,'Color',col_palette(1,:),'LineWidth',2)
% % plot(timevec_stats(statD32.mask),statD32.mask(statD32.mask).*0.0425,'Color',col_palette(4,:),'LineWidth',2)
% % ylabel('coherence')
% % xlabel('time (s)')
% % xlim([-0.2 0.5])
% % xticks([0:0.5:0.5])
% % ylim([0.01 0.045])
% % yticks([0:0.02:0.04])
% % box off
% % legend('guided T fast', 'guided T slow','guided D fast', 'guided D slow',Location='southeast')
% % 
% % 
% % print(fig,fullfile(cohfigpth,'coh_fast_slow'),'-dpng','-r0')
% % print(fig,fullfile(cohfigpth,'coh_fast_slow'),'-dsvg','-r0')
% % 
% % 
% % fig = figure('Position',[0 0  1980 1080/2.5]);
% % 
% % timevec_stats = timevec(500:end);
% % 
% % 
% % %% Unguided
% % % for set size 32
% % u16f = cell(1,length(subj));
% % u16s = cell(1,length(subj));
% % u32f = cell(1,length(subj));
% % u32s = cell(1,length(subj));
% % 
% % for s = 1:length(subj)
% % 
% %     % 16
% % 
% %     % fast
% %     u = squeeze(((niT(1,s,1,:) + niD(1,s,1,:)))./2);
% % 
% %     ERF.avg(1,:) = u;
% %     u16f{s} = ERF;
% % 
% %     % slow
% %     u = squeeze(((niT(1,s,2,:) + niD(1,s,2,:)))./2);
% %     ERF.avg(1,:) = u;
% %     u16s{s} = ERF;
% % 
% % 
% %     % 32
% % 
% %     % fast
% %     u = squeeze(((niT(2,s,1,:) + niD(2,s,1,:)))./2);
% % 
% %     ERF.avg(1,:) = u;
% %     u32f{s} = ERF;
% % 
% %     % slow
% %     u = squeeze(((niT(2,s,2,:) + niD(2,s,2,:)))./2);
% %     ERF.avg(1,:) = D;
% %     u32s{s} = ERF;
% % 
% % 
% % end
% % 
% % cfg         = [];
% % 
% % cfg.method           = 'montecarlo';                                            % cluster-based
% % cfg.channel          = ERF.label{1};
% % cfg.statistic        = 'depsamplesT';                                           % within subject
% % cfg.correctm         = 'cluster';
% % cfg.clusteralpha     = 0.05;                                                    % sample-specific t-value
% % cfg.clusterstatistic = 'maxsum';                                                % maximum t-value will be evaluated under permutation distribution   
% % cfg.minnbchan        = 0;
% % cfg.alpha            = 0.05;
% % cfg.numrandomization = 500;
% % cfg.latency = [0 min_rt];
% % Nsubj  = length(subj);
% % design = zeros(2, Nsubj*2);
% % design(1,:) = [1:Nsubj 1:Nsubj];
% % design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
% % 
% % cfg.design = design;
% % cfg.uvar   = 1;
% % cfg.ivar   = 2;
% % 
% % % compare Target and Distractor fast vs slow
% % 
% % cfg.tail             = 1;                                                       % two-sided test for both clustering and cluter-level statistic
% % cfg.clustertail = 1;
% % 
% % [statu16] = ft_timelockstatistics(cfg, u16f{:}, u16s{:});
% % 
% % cfg.tail             = -1;                                                       % two-sided test for both clustering and cluter-level statistic
% % cfg.clustertail = -1;
% % [statu32] = ft_timelockstatistics(cfg, u32f{:}, u32s{:});
% % 
% % fig = figure('Position',[0 0  1980 1080/2.5]);
% % 
% % % set size 16
% % ax1 = subplot(131);
% % plot(t_vec,mean(GA_ni_16_fast,1),'Color',col_palette(5,:),'LineWidth',3)
% % hold on
% % plot(t_vec,mean(GA_ni_16_slow,1),'Color',col_palette(6,:),'LineWidth',3)
% % 
% % ylabel('coherence')
% % xlabel('time (s)')
% % xlim([-0.1 0.5])
% % xticks([0:0.5:0.5])
% % ylim([0.01 0.05])
% % yticks([0:0.025:0.05])
% % box off
% % 
% % legend('guided D fast', 'guided D slow',Location='southeast')
% % 
% % 
% % % set size 32
% % ax2 = subplot(132);
% % plot(t_vec,mean(GA_ni_32_fast,1),'Color',col_palette(5,:),'LineWidth',3)
% % hold on
% % plot(t_vec,mean(GA_ni_32_slow,1),'Color',col_palette(6,:),'LineWidth',3)
% % 
% 
% ylabel('coherence')
% xlabel('time (s)')
% xlim([-0.1 0.5])
% xticks([0:0.5:0.5])
% ylim([0.01 0.05])
% yticks([0:0.025:0.05])
% box off
% legend('guided D fast', 'guided D slow',Location='southeast')
% 
% 
% print(fig,fullfile(cohfigpth,'coh_fast_slow_ung'),'-dpng','-r0')
% print(fig,fullfile(cohfigpth,'coh_fast_slow_ung'),'-dsvg','-r0')
% 
