%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d3. Alpha topoplots (Supplementary Fig. 3)
% prepare statistical test
% violin plots

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow


clear all; close all; clc


%% settings
pth = 'Z:\Visual Search RFT';
outpth = fullfile(pth,'results','meg','6 Alpha','pow');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
alphafigpth = fullfile(pth,'results','meg','6 Alpha','fig');
col_palette = [20,156,140;0,0,0]./255;
cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');

addpath('Z:\fieldtrip')
addpath(genpath(fullfile(pth,'ScientificColourMaps7')))
ft_defaults;
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
load(fullfile(pth,'matlab scripts/','coherence','occi_grad.mat'))

%% prepare cluster test

% spectrum IAF SOI
spec_iaf_pre = cell(1,length(subj));
% spectrum RIFT soi
spec_rift_pre = cell(1,length(subj));

% spectrum IAF SOI
spec_iaf_post = cell(1,length(subj));
% spectrum RIFT soi
spec_rift_post = cell(1,length(subj));

% pre-stimulus spectrum
spec_high = cell(1,length(subj));
spec_low = cell(1,length(subj));

% post-stimulus spectrum
spec_high_post = cell(1,length(subj));
spec_low_post = cell(1,length(subj));

for s = 1:length(subj)

    load(fullfile(outpth,subj{s},['data_winl_5.mat']),'TFR_alpha')
    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'soi_grad_cmb','iaf_grad')
    % RIFT sensors
    load(fullfile(cohsoipth,subj{s},'soi_stat.mat'))
   

    cfg = [];
    cfg.channel = soi_occi;
    TFR_RIFT = ft_selectdata(cfg,TFR_alpha);
    
    cfg = [];
    cfg.method = 'sum';
    TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha);
    TFR_RIFT_avg = ft_combineplanar(cfg,TFR_RIFT);
    clear TFR_alpha TFR_RIFT

    cfg = [];
    cfg.channel = soi_grad_cmb;
    TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);

    cfg = [];
    % pre-search
    cfg.latency = [-1 0];
    cfg.avgovertime = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.avgoverrpt = 'yes';
    TFR_alpha_avg_pre = ft_selectdata(cfg,TFR_alpha_avg);
    TFR_RIFT_pre = ft_selectdata(cfg,TFR_RIFT_avg);

    % search
    cfg.latency = [0.25 0.5];
    TFR_alpha_avg_post = ft_selectdata(cfg,TFR_alpha_avg);
    TFR_RIFT_post = ft_selectdata(cfg,TFR_RIFT_avg);
    clear TFR_alpha_avg TFR_alpha
    
    spec_iaf_pre{s} = TFR_alpha_avg_pre.powspctrm;
    spec_iaf_post{s} = TFR_alpha_avg_post.powspctrm;

    spec_rift_pre{s} = TFR_RIFT_pre.powspctrm;
    spec_rift_post{s} = TFR_RIFT_post.powspctrm;

end

freqvec = TFR_alpha_avg_post.freq;

%% plot pre-stimulus
fig = figure('Position',[0 0 1900/2 1080]);
for s = 1:length(subj)

    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'iaf_grad')
    % RIFT sensors
    load(fullfile(cohsoipth,subj{s},'soi_stat.mat'),'iaf_RIFT')

    f = find(freqvec==iaf_grad);
    subplot(8,4,s)
    plot(freqvec, spec_iaf_pre{s}.*10^23,'Color',col_palette(2,:),'Marker','*','MarkerIndices',[f f])
    hold on
    [~, f] = min(abs(freqvec - iaf_grad));
    plot(freqvec, spec_rift_pre{s}.*10^23,'Color',col_palette(1,:),'Marker','*','MarkerIndices',[f f])
    ylim([0 max(spec_iaf_pre{s}).*10^23+1])

    yticks([0 max(spec_iaf_pre{s}).*10^23+1])
    ytickformat('%.1f')
    if find(1:4:length(subj) == s)
        ylabel('power (T/m)^2')
    end

    if find(length(subj)-2:length(subj) == s)
        xlabel('frequency (Hz)')
    end
    xlim([4 30])
    xticks([10:10:30])
    box off
    box off
    clear iaf_pow* soi_grad_cmb iaf_grad
end

print(fig,fullfile(alphafigpth,'spectra_alpha_iaf_occi_sens'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'spectra_alpha_iaf_occi_sens'),'-dpng','-r0')


fig = figure('Position',[0 0 1900/2 1080]);
for s = 1:length(subj)

    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'iaf_grad')
    % RIFT sensors
    load(fullfile(cohsoipth,subj{s},'soi_stat.mat'),'iaf_RIFT')

    f = find(freqvec==iaf_grad);
    subplot(8,4,s)
    plot(freqvec, spec_iaf_post{s}.*10^23,'Color',col_palette(2,:),'Marker','*','MarkerIndices',[f f])
    hold on
    [~, f] = min(abs(freqvec - iaf_RIFT));
    plot(freqvec, spec_rift_post{s}.*10^23,'Color',col_palette(1,:),'Marker','*','MarkerIndices',[f f])
    ylim([0 max(spec_iaf_post{s}).*10^23+1])

    yticks([0 max(spec_iaf_post{s}).*10^23+1])
    ytickformat('%.1f')
    if find(1:4:length(subj) == s)
        ylabel('power (T/m)^2')
    end

    if find(length(subj)-2:length(subj) == s)
        xlabel('frequency (Hz)')
    end
    xlim([4 30])
    xticks([10:10:30])
    box off
    box off
    clear iaf_pow* soi_grad_cmb iaf_grad
end


print(fig,fullfile(alphafigpth,'spectra_alpha_iaf_occi_sens_post'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'spectra_alpha_iaf_occi_sens_post'),'-dpng','-r0')