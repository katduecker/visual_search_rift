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
outpth = fullfile(pth,'results','meg','6 Alpha','not align','pow');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
alphafigpth = fullfile(pth,'results','meg','6 Alpha','fig');
col_palette = [20,156,140;0,0,0]./255;

addpath('Z:\fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

%% prepare cluster test
% pre-stimulus spectrum
spec_high = cell(1,length(subj));
spec_low = cell(1,length(subj));

% post-stimulus spectrum
spec_high_post = cell(1,length(subj));
spec_low_post = cell(1,length(subj));

for s = 1:length(subj)

    subplot(8,4,s)
    load(fullfile(outpth,subj{s},['data_winl_5.mat']),'TFR_alpha')
    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'soi_grad_cmb','iaf_grad')

    cfg = [];
    cfg.method = 'sum';
    TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha);

    clear TFR_alpha
    cfg = [];
    % pre-search
    cfg.latency = [-1 0];
    cfg.avgovertime = 'yes';
    TFR_alpha_avg_pre = ft_selectdata(cfg,TFR_alpha_avg);
    % search
    cfg.latency = [0.25 0.5];
    TFR_alpha_avg_post = ft_selectdata(cfg,TFR_alpha_avg);

    clear TFR_alpha_avg TFR_alpha
    cfg = [];
    cfg.sensors = soi_grad_cmb;
    cfg.avgoverchan = 'yes';
    TFR_avg_sens = ft_selectdata(cfg,TFR_alpha_avg_pre);

    % find high alpha trials
    m_iaf_pow = median(squeeze(TFR_avg_sens.powspctrm(:,:,TFR_avg_sens.freq == iaf_grad)));
    h_iaf_trl = squeeze(TFR_avg_sens.powspctrm(:,:,TFR_avg_sens.freq == iaf_grad)) > m_iaf_pow;
    l_iaf_trl = squeeze(TFR_avg_sens.powspctrm(:,:,TFR_avg_sens.freq == iaf_grad)) < m_iaf_pow;

    clear TFR_avg_sens TFR_alpha
    freqvec = TFR_alpha_avg_pre.freq;

    % select the high alpha trials (keep all sensors so that struct is
    % intact for cluster based test!)
    cfg = [];
    cfg.trials = h_iaf_trl;
    cfg.avgoverrpt = 'yes';
    iaf_pow_high = ft_selectdata(cfg,TFR_alpha_avg_pre);
    spec_high{s} = iaf_pow_high;
    spec_high_post{s} = ft_selectdata(cfg,TFR_alpha_avg_post);
    cfg.trials = l_iaf_trl;
    cfg.avgoverrpt = 'yes';
    iaf_pow_low = ft_selectdata(cfg,TFR_alpha_avg_pre);
    spec_low{s} = iaf_pow_low;
    spec_low_post{s} = ft_selectdata(cfg,TFR_alpha_avg_post);

    % average over sensors and replace first sensor with averaged TFR
    cfg = [];
    cfg.channel = soi_grad_cmb;
    cfg.avgoverchan = 'yes';
    iaf_pow_high = ft_selectdata(cfg,iaf_pow_high);
    iaf_pow_low = ft_selectdata(cfg,iaf_pow_low);

    spec_high{s}.powspctrm(1,:) = iaf_pow_high.powspctrm;
    spec_low{s}.powspctrm(1,:) = iaf_pow_low.powspctrm;

    iaf_pow_high = ft_selectdata(cfg,spec_high_post{s});
    iaf_pow_low = ft_selectdata(cfg,spec_low_post{s});
    spec_high_post{s}.powspctrm(1,:) = iaf_pow_high.powspctrm;
    spec_low_post{s}.powspctrm(1,:) = iaf_pow_low.powspctrm;

end


%% plot pre-stimulus
freqvec = spec_high{1}.freq;
fig = figure('Position',[0 0 1900/2 1080]);
for s = 1:length(subj)
    subplot(8,4,s)
    plot(freqvec, spec_high{s}.powspctrm(1,:).*10^23,'Color',col_palette(2,:))
    hold on
    plot(freqvec, spec_low{s}.powspctrm(1,:).*10^23,'Color',col_palette(1,:))
    ylim([min(spec_low{s}.powspctrm(1,:)).*10^23 max(spec_high{s}.powspctrm(1,:)).*10^23])
    yticks([min(spec_low{s}.powspctrm(1,:)).*10^23 max(spec_high{s}.powspctrm(1,:)).*10^23])
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

print(fig,fullfile(alphafigpth,'spectra_alpha_high_vs_low_baseline'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'spectra_alpha_high_vs_low_baseline'),'-dpng','-r0')

close all

fig = figure('Position',[0 0 1900/2 1080]);
for s = 1:length(subj)
    subplot(8,4,s)
%% plot post-stimulus
    plot(freqvec, spec_high_post{s}.powspctrm(1,:).*10^23,'Color',col_palette(2,:))
    hold on
    plot(freqvec, spec_low_post{s}.powspctrm(1,:).*10^23,'Color',col_palette(1,:))
    ylim([min(spec_low_post{s}.powspctrm(1,:).*10^23) max(spec_high_post{s}.powspctrm(1,:).*10^23)])
    yticks([min(spec_low_post{s}.powspctrm(1,:).*10^23) max(spec_high_post{s}.powspctrm(1,:).*10^23)])
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
end

print(fig,fullfile(alphafigpth,'spectra_alpha_high_vs_low_search'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'spectra_alpha_high_vs_low_search'),'-dpng','-r0')

save(fullfile(outpth,'alpha_spec_high_low.mat'),'spec_high_post','spec_high','spec_low','spec_low_post')