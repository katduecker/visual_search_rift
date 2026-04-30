%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition

%% GLM Spectrum analysis
% Fit a GLM to the fourier Spectrum

% (c), Katharina Duecker
% last edited, Nov-26-2025
clear all; close all; clc

%% paths
rmpath(genpath('/rds/projects/2018/jenseno-entrainment'))
addpath('/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';

inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
load(fullfile(pth,'matlab_scripts/',"preproc_meg",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha','pow');

cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');
outpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');

load(fullfile(pth,'matlab_scripts','coherence','occi_grad.mat'))        % load occipital sensors
s= 2;
addpath(fullfile(pth,'matlab_scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);
plotpth = fullfile(pth,'results','meg','9 GLM', 'fig','GLM_spec_results');

%% Load TFR of fourier

load(fullfile(alphapth,subj{s},'data_fourier_winl_10.mat'),'TFR_avg')

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR_avg);

fig = figure('Position',[10 0 1800 1080/3]);

subplot(1,2,1)
cfg = [];
cfg.zlim = 'maxabs';
cfg.marker = 'off';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.xlim = [-1 0];
cfg.ylim = [8, 12];
cfg.figure = 'gca';
cfg.comment = 'no';
ft_topoplotTFR(cfg,TFR)
subplot(1,2,2)
cfg = [];
%cfg.channel = soi_occi;
cfg.figure = 'gca';
cfg.layout = 'neuromag306cmb_helmet.mat';
ft_singleplotTFR(cfg,TFR)
xlabel('time (s)')
xlim([-1.5, 0.5])
xticks(-1.5:0.5:0.5)
yticks(10:10:30)
ylabel('frequency (Hz)')
title('')
cb = colorbar;
box off
colormap(cm)
print(fig,fullfile(plotpth,'graph_abstract_alpha'),'-dpng')
print(fig,fullfile(plotpth,'graph_abstract_alpha'),'-dsvg')
print(fig,fullfile(plotpth,'graph_abstract_alpha'),'-dpdf')