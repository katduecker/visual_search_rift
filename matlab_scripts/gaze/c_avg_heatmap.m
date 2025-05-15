%% Eye movement heatmap

clear all; close all; clc;
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid','clean data all trials');
mergepth = fullfile(pth,'results','meg','2 merged edf mat');

occupth = fullfile(pth,'results','eyelink','heatmap');

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

propixx_res      = [1920 1080];                                   % propixx resolution
scr.w            = 72;                                       % screen width in cm
scr.h            = 40.5;                                     % screen height in cm
scr.d            = 142;
scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
scr.scrdegrw     = asind(scr.w/scr.ch);
scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pixel

hist_bins_x = (linspace(0, propixx_res(1), propixx_res(1)/10) - propixx_res(1)/2)./scr.onedegrpix;
hist_bins_y = (linspace(0, propixx_res(2), propixx_res(2)/10) - propixx_res(2)/2)./scr.onedegrpix;


ocu_pre = zeros(length(subj),length(hist_bins_x)-1, length(hist_bins_y)-1);
ocu_search = zeros(length(subj),length(hist_bins_x)-1, length(hist_bins_y)-1);

for s = 1:length(subj)
    
    load(fullfile(occupth,subj{s},'eye_heatmap.mat'))
    
    
    ocu_pre(s,:,:) = squeeze(mean(pdf_eye_pre,1));
    ocu_search(s,:,:) = squeeze(mean(pdf_eye_search,1));
end

fig = figure('Position', [0 0 1920/2 1080]);

for s = 1:length(subj)
    subplot(8,4,s)
    
    imagesc(hist_bins_x, hist_bins_y, squeeze(ocu_pre(s,:,:)))
    colormap(cm(round(length(cm)/2):end,:))
    caxis([0, 0.01])
    colorbar
    hold on
    axis image
    rectangle('Position',[-1, -1, 2, 2])
    rectangle('Position',[-5, -5, 10, 10])
    
end

print(fig,fullfile(occupth, 'heatmaps_subj_pre'),'-dpng')
print(fig,fullfile(occupth, 'heatmaps_subj_pre'),'-dsvg')


fig = figure('Position', [0 0 1920/2 1080]);

for s = 1:length(subj)
    subplot(8,4,s)
    
    imagesc(hist_bins_x, hist_bins_y, squeeze(ocu_search(s,:,:)))
    colormap(cm(round(length(cm)/2):end,:))
    caxis([0, 0.01])
    colorbar
    hold on
    axis image
    rectangle('Position',[-1, -1, 2, 2])
    rectangle('Position',[-5, -5, 10, 10])
end

print(fig,fullfile(occupth, 'heatmaps_subj_search'),'-dpng')
print(fig,fullfile(occupth, 'heatmaps_subj_search'),'-dsvg')

fig = figure('Position', [0, 0, 1920/1.5,1080/2]);
subplot(1,2,1)
imagesc(hist_bins_x, hist_bins_y, squeeze(mean(ocu_pre)))
colormap(cm(round(length(cm)/2):end,:))
caxis([0, 0.01])

axis image
rectangle('Position',[-1, -1, 2, 2])
rectangle('Position',[-5, -5, 10, 10])
cb = colorbar;
cb.Ticks = [0, 0.01];
ylabel('height (degree)')
xlabel('width (degree)')
yticks(-5:5:5)
xticks(-10:10:10)

subplot(1,2,2)
imagesc(hist_bins_x, hist_bins_y, squeeze(mean(ocu_search)))
colormap(cm(round(length(cm)/2):end,:))
caxis([0, 0.01])
yticks(-5:5:5)
xticks(-10:10:10)
axis image
rectangle('Position',[-1, -1, 2, 2])
rectangle('Position',[-5, -5, 10, 10])
xlabel('width (degree)')
cb = colorbar;
cb.Ticks = [0, 0.01];
print(fig,fullfile(occupth, 'avg_heatmaps'),'-dpng')
print(fig,fullfile(occupth, 'avg_heatmaps'),'-dsvg')