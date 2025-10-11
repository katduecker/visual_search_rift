%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c1. prep alignment of TFR to soi: find alpha frequency range in sample
% plot alpha spectra & mark IAF (Supplementary Fig. 4)

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
set(0,'defaultAxesFontName','Arial')

addpath('Z:\fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

winl = 0.5;

iaf = zeros(size(subj,2),1);
minf = zeros(size(subj,2),1);
maxf = zeros(size(subj,2),1);

load(fullfile(outpth,subj{1},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha_avg')

for s = 1:length(subj)
    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'iaf_grad')
    
    freqvec = round(TFR_alpha_avg.freq - iaf_grad);
    
    iaf(s) = iaf_grad;
    minf(s) = freqvec(1);
    maxf(s) = freqvec(end);
    
    clear iaf_grad
end

max_minf = max(minf);
min_maxf = min(maxf);

save('alpha_align_vec.mat','minf','maxf','max_minf','min_maxf')

freqvec = TFR_alpha_avg.freq;

fig = figure('Position',[0 0 1900/2 1080]);
for s = 1:length(subj)

    subplot(8,4,s)
    load(fullfile(outpth,subj{s},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha_avg')
    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'soi_grad_cmb','iaf_grad')

    cfg = [];
    cfg.method = 'sum';
    TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

    cfg = [];
    cfg.channel = soi_grad_cmb;
    cfg.avgoverchan = 'yes';
    cfg.latency = [-1 0];
    cfg.avgovertime = 'yes';
    TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);
    
    plot(freqvec,TFR_alpha_avg.powspctrm.*10^23,'-*','Color','k','MarkerIndices',find(TFR_alpha_avg.freq == iaf_grad),'MarkerSize',10)
    %ylim([min(TFR_alpha_avg.powspctrm.*10^23) max(TFR_alpha_avg.powspctrm.*10^23)])
    yticks([min(TFR_alpha_avg.powspctrm.*10^23) max(TFR_alpha_avg.powspctrm.*10^23)])
    ytickformat('%.1f')
    box off
    clear TFR_alpha_avg soi_grad_cmb iaf_grad
end

% adjust P's

s = 1;

load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'soi_grad_cmb','iaf_grad')
load(fullfile(outpth,subj{s},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha_avg')

cfg = [];
cfg.method = 'sum';
TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

cfg = [];
cfg.channel = soi_grad_cmb;
cfg.avgoverchan = 'yes';
cfg.latency = [-0.5 0];
cfg.avgovertime = 'yes';
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);
plot(freqvec,TFR_alpha_avg.powspctrm,'-+','Color','k','MarkerIndices',find(TFR_alpha_avg.freq == iaf_grad),'MarkerSize',10)
iaf_grad = TFR_alpha_avg.freq(find(TFR_alpha_avg.freq == iaf_grad)-1);

save(fullfile(iafpth,subj{s},'iaf_soi.mat'),'-append','iaf_grad')
close all
s=13;
load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'soi_grad_cmb','iaf_grad')
load(fullfile(outpth,subj{s},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha_avg')

cfg = [];
cfg.method = 'sum';
TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

cfg = [];
cfg.channel = soi_grad_cmb;
cfg.avgoverchan = 'yes';
cfg.latency = [-0.5 0];
cfg.avgovertime = 'yes';
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);
plot(freqvec,TFR_alpha_avg.powspctrm,'-+','Color','k','MarkerIndices',find(TFR_alpha_avg.freq == iaf_grad),'MarkerSize',10)
iaf_grad = TFR_alpha_avg.freq(find(TFR_alpha_avg.freq == iaf_grad)-1);

save(fullfile(iafpth,subj{s},'iaf_soi.mat'),'-append','iaf_grad')

close all

fig = figure('Position',[0 0 1900/2 1080]);
for s = 1:length(subj)

    subplot(8,4,s)
    load(fullfile(outpth,subj{s},['data_winl_',num2str(winl*10),'.mat']),'TFR_alpha_avg')
    load(fullfile(iafpth,subj{s},'iaf_soi.mat'),'soi_grad_cmb','iaf_grad')

    cfg = [];
    cfg.method = 'sum';
    TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

    cfg = [];
    cfg.channel = soi_grad_cmb;
    cfg.avgoverchan = 'yes';
    cfg.latency = [-1 0];
    cfg.avgovertime = 'yes';
    TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);
    
    plot(freqvec,TFR_alpha_avg.powspctrm.*10^23,'-*','Color','k','MarkerIndices',find(TFR_alpha_avg.freq == iaf_grad),'MarkerSize',10)
    %ylim([min(TFR_alpha_avg.powspctrm.*10^23) max(TFR_alpha_avg.powspctrm.*10^23)])
    yticks([min(TFR_alpha_avg.powspctrm.*10^23) max(TFR_alpha_avg.powspctrm.*10^23)])
    ytickformat('%.1f')
    if find(1:4:length(subj) == s)
        ylabel('power (T/m)^2')
    end

    if find(length(subj)-3:length(subj) == s)
        xlabel('frequency (Hz)')
    end
    xlim([4 30])
    xticks([10:10:30])
    box off
    clear TFR_alpha_avg soi_grad_cmb iaf_grad
end


print(fig,fullfile(alphafigpth,'iaf_subj'),'-dsvg','-r0')
print(fig,fullfile(alphafigpth,'iaf_subj'),'-dpng','-r0')