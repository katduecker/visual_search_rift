%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c4. plot grandaverage alpha TFR (Fig 4 a)

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
outpth = fullfile(pth,'results','meg','6 Alpha','pow align iaf');
iafpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Arial')
alphafigpth = fullfile(pth,'results','meg','6 Alpha','fig');
behavpth = fullfile(pth,'results','behavior');
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;


addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',201);
cm = flipud(cm);

load(fullfile(outpth,'tfr_ga.mat'))
varnames = {'gui 16', 'ung 16','gui 32','ung 32'};

fig = figure('Position',[0 0 1920 1080/2.5]);
subplot(131)
imagesc(timevec,freq,squeeze(mean(tfr_ga,1)))
xticks([-1.5:0.5:0.5])
xlim([-1.5 0.5])
yticks(-4:4:16)

axis xy
xlabel('time (s)')
ylabel('Individual Alpha Frequency + (Hz)')
cb = colorbar;
caxis([-56 -52])
cb.Label.String = 'log(power)';
cb.FontSize = 12;
cb.FontName = 'Arial';
cb.Ticks = -56:-52;
box off
colormap(cm)


print(fig,fullfile(alphafigpth,'tfr_ga'),'-dpng','-r0')
print(fig,fullfile(alphafigpth,'tfr_ga'),'-dsvg','-r0')

