%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e1. plot grandaverage T boost Dsupp

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked April 2024

%% Source analysis
% a: align digitized headshape to T1 -> realign Tq
% b: Forward model/Lead field
% c: DICS beamformer

clear all; close all; clc
%% settings & paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path

toi = [-2.5 2];
dtpth = fullfile(pth, 'data'); % data path
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

% get template source localized coherence
load(fullfile(outpth,subj{1},'dics_filt_rft.mat'))


%% RFT

% colormap
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);

% prepare contrast t boosting d suppression

rftT_subj = cell(1,length(subj));
rftD_subj = rftT_subj;
rftungT_subj = rftT_subj;
rftungD_subj = rftT_subj;

for s=1:length(subj)
    load(fullfile(outpth,subj{s},'dics_tboost_dsupp.mat'))
    
    % fill in target
    rftT_subj{s} = rft_dics{1};
    rftT_subj{s}.avg.coh = source_guiT;
    % fill in distractor
    rftD_subj{s} = rft_dics{1};
    rftD_subj{s}.avg.coh = source_guiD;
    % fill in unguided
    rftungT_subj{s} = rft_dics{1};
    rftungT_subj{s}.avg.coh = source_ungT;
    rftungD_subj{s} = rft_dics{1};
    rftungD_subj{s}.avg.coh = source_ungD;
end

%grandaverage to 

nsubj = length(subj);
% run statistics over subjects %
cfg=[];
cfg.dim         = template.sourcemodel.dim;
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.parameter   = 'coh';
cfg.correctm    = 'cluster';
cfg.numrandomization = 1000;
cfg.alpha       = 0.05; % note that this only implies single-sided testing
cfg.tail        = 1;

cfg.design(1,:) = [1:nsubj 1:nsubj];
cfg.design(2,:) = [ones(1,nsubj)*1 ones(1,nsubj)*2];
cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

statT = ft_sourcestatistics(cfg, rftT_subj{:}, rftungT_subj{:});

cfg.tail = -1;
statD = ft_sourcestatistics(cfg, rftD_subj{:}, rftungD_subj{:});

clear rftung_subj rftD_subj
%% plot
% interpolate
temp_subjT = rftT_subj{1};
temp_subjT.avg.coh = statT.stat;
temp_subjT.mask = statT.mask;
temp_subjD = temp_subjT;
temp_subjD.avg.coh = statD.stat;
temp_subjD.mask = statD.mask; 

clear statT statD

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'coh';
srcT = ft_sourceinterpolate(cfg,temp_subjT,templmri);
srcD = ft_sourceinterpolate(cfg,temp_subjD,templmri);

cfg.parameter = 'mask';
srcT_mask = ft_sourceinterpolate(cfg,temp_subjT,templmri);
srcD_mask = ft_sourceinterpolate(cfg,temp_subjD,templmri);

srcT.mask = srcT_mask.mask;
srcD.mask = srcD_mask.mask;

clear srcT_mask srcD_mask

srcT.coordsys = 'mni';
srcD.coordsys = 'mni';

fig = figure;
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'coh';
%cfg.funcolorlim = 'zeromax';
cfg.maskparameter = 'mask';
cfg.figure = 'gca';
cfg.crosshair = 'no';
%cfg.funcolormap = cm(round(length(cm)/2):end,:);
ft_sourceplot(cfg,srcT);
% print(fig,fullfile(plotpth,'disc_tvsung'),'-dpng','-r600')
% print(fig,fullfile(plotpth,'disc_tvsung'),'-dsvg','-r600')
% close all
% 
fig = figure;
%cfg.funcolorlim = 'zeromin';
%cfg.funcolormap = cm(1:round(length(cm)/2),:);
ft_sourceplot(cfg,srcD);
print(fig,fullfile(plotpth,'disc_dvsung'),'-dpng','-r600')
print(fig,fullfile(plotpth,'disc_dvsung'),'-dsvg','-r600')

