%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d. plot grandaverage RIFT beamformer result (Fig. 3 b)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Source analysis
% a: align digitized headshape to T1 -> realign Tq
% b: Forward model/Lead field
% c: DICS beamformer

clear all; close all; clc
%% settings & paths
pth = 'Z:\Visual Search RFT';                                          % server path

toi = [-2.5 2];
dtpth = fullfile(pth, 'data'); % data path
addpath(fullfile('Z:\','fieldtrip'))
ft_defaults;

% location where fieldtrip is installed
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template', 'sourcemodel');
% load 4 mmm sourcemodel
template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm'));
templatedir = fullfile(ftdir, 'template','anatomy');
templmri = ft_read_mri(fullfile(templatedir,'single_subj_T1.nii'));
outpth = fullfile(pth,'results','meg','7 Beamformer');
plotpth = fullfile(outpth,'fig');

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

%% RFT

% colormap
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);

% 1. RFT general

rft_subj = cell(1,length(subj));
for s=1:length(subj)

    load(fullfile(outpth,subj{s},'dics_filt_rft.mat'))

    cfg = [];
    cfg.parameter = 'coh';
    rft_subj{s} = ft_sourcegrandaverage(cfg,rft_dics{:});

end

cfg = [];
cfg.parameter = 'coh';
ga_rft = ft_sourcegrandaverage(cfg,rft_subj{:});

ga_rft = rmfield(ga_rft,'cfg');


clear rft_subj src_rft rft_dics

x = maxk(ga_rft.coh,floor(sum(ga_rft.inside)*0.01));
mask = ga_rft.coh > x(end);
ga_rft.mask = mask;

clear x mask

% interpolate
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'coh';
src_rft_int = ft_sourceinterpolate(cfg,ga_rft,templmri);

cfg.parameter = 'mask';
src_rft_int_mask = ft_sourceinterpolate(cfg,ga_rft,templmri);

src_rft_int.mask = src_rft_int_mask.mask;

clear src_rft_int_mask
src_rft_int.coordsys = 'mni';

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'coh';
cfg.funcolorlim = 'zeromax';
cfg.maskparameter = 'mask';
cfg.crosshair = 'no';
cfg.funcolormap = cm(round(length(cm)/2):end,:);
ft_sourceplot(cfg,src_rft_int);


savefig(fullfile(plotpth,'ga_rft'))

clear src_*
close all

fig = openfig(fullfile(plotpth,'ga_rft'));
print(fig,fullfile(plotpth,'ga_rft'),'-dpng','-r600')
print(fig,fullfile(plotpth,'ga_rft'),'-dsvg','-r600')
delete(fullfile(plotpth,'ga_rft.fig'))

close all
