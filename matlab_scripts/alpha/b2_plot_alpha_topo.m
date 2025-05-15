%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b2. Alpha topoplots (Supplementary Fig. 3)
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
% e. control analysis: compare alpha for fast vs slow


clear all; close all; clc


%% settings
pth = 'Z:\Visual Search RFT';
outpth = fullfile(pth,'results','meg','6 Alpha','pow');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
alphafigpth = fullfile(pth,'results','meg','6 Alpha','fig');

addpath('Z:\fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))

cm = cbrewer('seq','Greens',101);

cm(cm>1)=1;
cm(cm<0)=0;

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));



fig = figure('Position',[0 0 920 1080]);
for s = 1:length(subj)
    load(fullfile(outpth,subj{s},'data_winl_5.mat'),'TFR_alpha_avg')
    load(fullfile(iafpth,subj{s},'iaf_soi.mat'))

    cfg = [];
    cfg.method = 'sum';
    TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

    TFR_alpha_avg.powspctrm = TFR_alpha_avg.powspctrm.*10^23;

    subplot(8,4,s)
    cfg = [];
    cfg.marker = 'off';
    cfg.zlim = 'zeromax';
    cfg.comment = 'no';
    cfg.colormap = cm;
    cfg.colorbar = 'yes';
    cfg.layout = 'neuromag306cmb_helmet.mat';

    cfg.ylim = [iaf_grad iaf_grad];
    cfg.xlim = [-1 0];
    cfg.highlight = 'on';
    cfg.highlightsymbol = 'o';
    cfg.highlightcolor = [255,0,255]./255;
    cfg.highlightchannel = find(ismember(TFR_alpha_avg.label,soi_grad_cmb));
    cfg.highlightsize = 4;
    ft_topoplotTFR(cfg,TFR_alpha_avg)
    cb = colorbar;
    cb.Ticks = cb.Limits;
    cb.TickLabels = strsplit(sprintf('%0.1f ',cb.Ticks));

    clear TFR_alpha* soi_*
end

print(fig,fullfile(alphafigpth,'alpha_topo'),'-dsvg','-r600','-painters')

print(fig,fullfile(alphafigpth,'alpha_topo'),'-dpng','-r0')