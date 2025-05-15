%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d: prep: alpha vs time-on-task!

% Input
% -s : subject index


% Output
% csv file per participant, containing RT and alpha power for each trial 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 5 Aug 2023

clear all; close all; clc; beep off;

addpath('Z:\fieldtrip')
ft_defaults;

pth = 'Z:\Visual Search RFT';
addpath('Z:\Visual Search RFT\Violinplot-Matlab-master')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
plotpth = fullfile(pth,'results','exploratory fig');

col_palette = [228,131,12; 146, 90,20; 12, 146, 237; 20, 87, 132; 0 0 0; 120 120 120]./255;
condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

addpath(genpath('Z:\Visual Search RFT\ScientificColourMaps7\'))
load('batlow10.mat')

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

%% GLM RT ~ alpha

load(fullfile(plotpth,"stat_GLM_pow.mat"))

% find TOI in alpha band
[~,foi] = intersect(round(stat_pow.freq),[8,10,12]);

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_pow.time)

    if ~isempty(find(stat_pow.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_pow.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_pow.label(find(sum(stat_pow.mask(:,foi,toi_idx(t)),2)));
end


cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs'
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(3,4,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    %cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_pow)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end
print(fig,fullfile(plotpth,'GLM_POW_regr_time'),'-dpng')

soi_pow_occi = unique(vertcat(chan_oi{6:8}));






%% TOT - uniform?

% find TOI in alpha band
[~,foi] = intersect(round(stat_tot.freq),[8,10,12]);

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_tot.time)

    if ~isempty(find(stat_tot.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_tot.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_pow.label(find(sum(stat_tot.mask(:,foi,toi_idx(t)),2)));
end

cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs'
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(6,6,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    %cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_pow)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end

soi_tot = unique(vertcat(chan_oi{:}));

%% General plot

fig = figure('Position',[0 0 1920/3 1080/2]);
subplot(231)
cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_pow_occi;
cfg.comment = 'no';
cfg.xlim = [-0.25 -0.1];
ft_topoplotTFR(cfg,stat_pow)


subplot(2,3,2:3)
cfg = [];
cfg.channel = soi_pow_occi;
cfg.parameter = 'stat';
cfg.zlim = [-3 3];
cfg.title = '';
ft_singleplotTFR(cfg,stat_pow)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([5:5:30])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 'power regressor (T statistic)';
cb.Ticks = [-3 0 3];


subplot(234)
cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_tot;
cfg.comment = 'no';
% cfg.xlim = [-0.25 -0.1];
cfg.zlim = [-3 3];
ft_topoplotTFR(cfg,stat_pow)


subplot(2,3,5:6)
cfg = [];
cfg.channel = soi_tot;
cfg.parameter = 'stat';
cfg.zlim = [-3 3];
cfg.title = '';
ft_singleplotTFR(cfg,stat_pow)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([5:5:30])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 'tot regressor (T statistic)';
cb.Ticks = [-3 0 3];
print(fig,fullfile(plotpth,'GLM_pow_tot'),'-dpng')


%% GLM figure Alpha ~ RT

close all
load(fullfile(plotpth,"stat_GLM_rt_slouch_beta.mat"))

[~,foi] = intersect(round(stat_rt.freq),[8,10,12]);

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_rt.time)

    if ~isempty(find(stat_rt.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_rt.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_rt.label(find(sum(stat_rt.mask(:,foi,toi_idx(t)),2)));
end

% save the SOIs
soi_GLM = stat_rt.label(find(sum(logical(sum(stat_rt.mask(:,foi,toi_idx),2)),3) >= 10));


cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs'
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(4,5,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_rt)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end
print(fig,fullfile(plotpth,'GLM_rt_regr_time'),'-dpng')

close all


%% find slouch SOI
toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_rt.time)

    if ~isempty(find(stat_slouch.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_slouch.time(t)];
    end
end

chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_slouch.label(find(sum(stat_slouch.mask(:,foi,toi_idx(t)),2)));
end

cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(6,6,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    %cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_slouch)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end
print(fig,fullfile(plotpth,'GLM_slouch_regr_time'),'-dpng')

soi_GLM_slouch = unique(vertcat(chan_oi{:}));

%% PLOT regressors
% select first 6 time points

toi_bsl = -0.35:0.05:-0.1;
chan_oi_bsl = cell(1,length(toi_bsl));
for t = 1:length(toi_bsl)
    [~,p] = min(abs(stat_rt.time - toi_bsl(t)));
    chan_oi_bsl{t} = stat_rt.label(find(sum(stat_rt.mask(:,foi,p),2)));
end

soi_rt_bsl = unique(vertcat(chan_oi_bsl{:}));

cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
layout_soi = ft_prepare_layout(cfg);

ft_plot_layout(layout_soi)
% RT regressor


% 
fig = figure('Position',[0 0 1920/3 1080/1.7]);
subplot(331)
cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_rt_bsl;
cfg.comment = 'no';
cfg.xlim = [toi_bsl(1) toi_bsl(end)];
ft_topoplotTFR(cfg,stat_rt)


subplot(3,3,2:3)
cfg = [];
cfg.channel = soi_rt_bsl;
cfg.parameter = 'stat';
cfg.zlim = [-2.5 2.5];
cfg.title = '';
ft_singleplotTFR(cfg,stat_rt)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([5:5:30])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 'reaction time regressor (T statistic)';
cb.Ticks = [-2.5 0 2.5];


[~,foi] = intersect(round(stat_tot.freq),[8,10,12]);

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_tot.time)

    if ~isempty(find(stat_tot.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_tot.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_tot.label(find(sum(stat_tot.mask(:,foi,toi_idx(t)),2)));
end
soi_GLM_tot = stat_tot.label(find(sum(logical(sum(stat_tot.mask(:,foi,toi_idx),2)),3) >= 10));


cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_GLM_tot;

subplot(334)
cfg.xlim = [stat_tot.time(1) stat_tot.time(end)];
cfg.title = num2str(toi(t));
cfg.comment = 'no';
ft_topoplotTFR(cfg,stat_tot)


subplot(3,3,5:6)
cfg= [];
cfg.channel = soi_GLM_tot;
cfg.parameter = 'stat';
cfg.title = '';
cfg.zlim = [-5 5];
ft_singleplotTFR(cfg,stat_tot)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([5:5:30])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 't-o-t regressor (T statistic)';
cb.Ticks = [-5 0 5];

cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.highlightchannel = soi_GLM_slouch;
cfg.layout = 'neuromag306cmb_helmet.mat';

subplot(337)
cfg.xlim = [stat_slouch.time(1) stat_slouch.time(end)];
cfg.title = num2str(toi(t));
cfg.comment = 'no';
ft_topoplotTFR(cfg,stat_slouch)


subplot(3,3,8:9)
cfg = [];
cfg.channel = soi_GLM_slouch;
cfg.parameter = 'stat';
cfg.title = '';
cfg.zlim = [-3 3];
ft_singleplotTFR(cfg,stat_slouch)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([5:5:30])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 'slouch regressor (T statistic)';
cb.Ticks = [-3 0 3];

print(fig,fullfile(plotpth,'GLM_main_results_rt_tot_slouch'),'-dpng')






fig = figure;

subplot(121)
cfg = [];
cfg.ylim = [10 10];
cfg.parameter = 'stat';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.zlim = 'maxabs';
ft_topoplotTFR(cfg,stat_int)
cb = colorbar;
cb.Label.String = 'rt*t-o-t regressor (T statistic)';
subplot(122)
ft_topoplotTFR(cfg,stat_const)
cb = colorbar;
cb.Label.String = 'constant (T statistic)';
colormap(cm)
print(fig,fullfile('GLM_inter_constant_tot_z'),'-dpng')



%% interaction - n.s.

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_int.time)

    if ~isempty(find(stat_int.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_rt.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_int.label(find(sum(stat_int.mask(:,foi,toi_idx(t)),2)));
end

% save the SOIs
soi_GLM_int = stat_int.label(find(sum(logical(sum(stat_int.mask(:,foi,toi_idx),2)),3) >= 10));


cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs'
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(4,5,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_int)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end
print(fig,fullfile('GLM_int_regr_time'),'-dpng')




%% Aligned IAF GLM

load('alpha_align_vec.mat')
freqvec = max_minf:2:min_maxf;

load(fullfile(plotpth,"stat_GLM_rt_slouch_beta_align.mat"))

foi = stat_rt.freq==0;

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_rt.time)

    if ~isempty(find(stat_rt.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_rt.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_rt.label(find(sum(stat_rt.mask(:,foi,toi_idx(t)),2)));
end

% save the SOIs
soi_GLM = stat_rt.label(find(sum(logical(sum(stat_rt.mask(:,foi,toi_idx),2)),3) >= 10));


cfg = [];
cfg.ylim = [0 0];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs'
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(2,6,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_rt)
    title(num2str(toi(t)));
    colormap(cm)
%    cb = colorbar;
%     cb .Ticks = [-3 0 3];
end
print(fig,fullfile(plotpth,'GLM_rt_regr_time_align'),'-dpng')

close all


%% find slouch SOI
toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_rt.time)

    if ~isempty(find(stat_slouch.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_slouch.time(t)];
    end
end

chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_slouch.label(find(sum(stat_slouch.mask(:,foi,toi_idx(t)),2)));
end

cfg = [];
cfg.ylim = [-4 18];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(6,6,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    %cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_slouch)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end
print(fig,fullfile(plotpth,'GLM_slouch_regr_time_align'),'-dpng')

soi_GLM_slouch = unique(vertcat(chan_oi{:}));

%% PLOT regressors
% select first 6 time points

toi_bsl = -0.35:0.05:-0.1;
chan_oi_bsl = cell(1,length(toi_bsl));
for t = 1:length(toi_bsl)
    [~,p] = min(abs(stat_rt.time - toi_bsl(t)));
    chan_oi_bsl{t} = stat_rt.label(find(sum(stat_rt.mask(:,foi,p),2)));
end

soi_rt_bsl = unique(vertcat(chan_oi_bsl{:}));

cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
layout_soi = ft_prepare_layout(cfg);

ft_plot_layout(layout_soi)
% RT regressor

close all
% 
fig = figure('Position',[0 0 1920/3 1080/1.7]);
subplot(331)
cfg = [];
cfg.ylim = [0 0];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_rt_bsl;
cfg.comment = 'no';
cfg.xlim = [toi_bsl(1) toi_bsl(end)];
ft_topoplotTFR(cfg,stat_rt)


subplot(3,3,2:3)
cfg = [];
cfg.channel = soi_rt_bsl;
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.title = '';
ft_singleplotTFR(cfg,stat_rt)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([-4:4:18])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 'reaction time regressor (T statistic)';
%cb.Ticks = [-2.5 0 2.5];


[~,foi] = intersect(round(stat_tot.freq),[8,10,12]);

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_tot.time)

    if ~isempty(find(stat_tot.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_tot.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_tot.label(find(sum(stat_tot.mask(:,foi,toi_idx(t)),2)));
end
soi_GLM_tot = stat_tot.label(find(sum(logical(sum(stat_tot.mask(:,foi,toi_idx),2)),3) >= 10));


cfg = [];
cfg.ylim = [0 0];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_GLM_tot;

subplot(334)
cfg.xlim = [stat_tot.time(1) stat_tot.time(end)];
cfg.title = num2str(toi(t));
cfg.comment = 'no';
ft_topoplotTFR(cfg,stat_tot)


subplot(3,3,5:6)
cfg= [];
cfg.channel = soi_GLM_tot;
cfg.parameter = 'stat';
cfg.title = '';
cfg.zlim = 'maxabs';
ft_singleplotTFR(cfg,stat_tot)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([-4:4:18])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 't-o-t regressor (T statistic)';
%cb.Ticks = [-5 0 5];

cfg = [];
cfg.ylim = [0 0];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs';
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.highlightchannel = soi_GLM_slouch;

subplot(337)
cfg.xlim = [stat_slouch.time(1) stat_slouch.time(end)];
cfg.title = num2str(toi(t));
cfg.comment = 'no';
ft_topoplotTFR(cfg,stat_slouch)


subplot(3,3,8:9)
cfg = [];
cfg.channel = soi_GLM_slouch;
cfg.parameter = 'stat';
cfg.title = '';
cfg.zlim = [-3 3];
ft_singleplotTFR(cfg,stat_slouch)
ylabel('frequency (Hz)')
xlabel('time (s)')
xticks(-1:0.5:0.5)
yticks([-4:4:18])
colormap(cm)
title('')
cb = colorbar;
cb.Label.String = 'slouch regressor (T statistic)';
cb.Ticks = [-3 0 3];

print(fig,fullfile(plotpth,'GLM_main_results_rt_tot_slouch'),'-dpng')






fig = figure;

subplot(121)
cfg = [];
cfg.ylim = [10 10];
cfg.parameter = 'stat';
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.zlim = 'maxabs';
ft_topoplotTFR(cfg,stat_int)
cb = colorbar;
cb.Label.String = 'rt*t-o-t regressor (T statistic)';
subplot(122)
ft_topoplotTFR(cfg,stat_const)
cb = colorbar;
cb.Label.String = 'constant (T statistic)';
colormap(cm)
print(fig,fullfile('GLM_inter_constant_tot_z'),'-dpng')



%% interaction - n.s.

toi = [];
toi_idx = [];
% find significant time-points
for t = 1:length(stat_int.time)

    if ~isempty(find(stat_int.mask(:,foi,t)))
        toi_idx = [toi_idx,t];
        toi = [toi, stat_rt.time(t)];
    end
end


chan_oi = cell(1,length(toi));
for t = 1:length(toi)
    chan_oi{t} = stat_int.label(find(sum(stat_int.mask(:,foi,toi_idx(t)),2)));
end

% save the SOIs
soi_GLM_int = stat_int.label(find(sum(logical(sum(stat_int.mask(:,foi,toi_idx),2)),3) >= 10));


cfg = [];
cfg.ylim = [8 12];
cfg.parameter = 'stat';
cfg.zlim = 'maxabs'
cfg.highlight = 'on';
cfg.layout = 'neuromag306cmb_helmet.mat';

fig = figure;
for t = 1:length(toi_idx)
    subplot(4,5,t)
    cfg.xlim = [toi(t) toi(t)];
    cfg.highlightchannel = chan_oi{t};
    cfg.title = num2str(toi(t));
    cfg.comment = 'no';
    cfg.zlim =[-3 3];
    ft_topoplotTFR(cfg,stat_int)
    title(num2str(toi(t)));
    colormap(cm)
    colorbar
end
print(fig,fullfile('GLM_int_regr_time'),'-dpng')


