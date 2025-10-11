%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d. plot grandaverage alpha high-low beamformer

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


% colormap
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('seq','Greens',101);
cm(cm>1) =1;
cm(cm<0)=0;

%cm = flipud(cm);


%% Alpha high vs low bsl
load(fullfile(outpth,subj{1},'dics_filt_iaf_search.mat'))
iaf_subj = cell(1,length(subj));
for s=1:length(subj)

%     load(fullfile(outpth,subj{s},'dics_filt_iaf.mat'))
% 
     load(fullfile(outpth,subj{s},'dics_iaf_high_low.mat'),'iaf_high_dics','iaf_low_dics')
     
     x = iaf_high_dics./iaf_low_dics - 1;
    iaf_subj{s} = iaf_dics;
    iaf_subj{s}.avg.pow = x;
    
    iaf_subj{s}.freq = 0;
    
    clear iaf_high_dics iaf_low_dics x
end

cfg = [];
cfg.parameter = 'pow';
ga_iaf = ft_sourcegrandaverage(cfg,iaf_subj{:});

ga_iaf = rmfield(ga_iaf,'cfg');


x = maxk(ga_iaf.pow,floor(sum(ga_iaf.inside)*0.01));
mask = ga_iaf.pow > x(end);
ga_iaf.mask = mask;

clear x mask

% interpolate
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'pow';
src_rft_int = ft_sourceinterpolate(cfg,ga_iaf,templmri);

cfg.parameter = 'mask';
src_rft_int_mask = ft_sourceinterpolate(cfg,ga_iaf,templmri);

src_rft_int.mask = src_rft_int_mask.mask;

clear src_rft_int_mask
src_rft_int.coordsys = 'mni';

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'zeromax';
cfg.maskparameter = 'mask';
cfg.crosshair = 'no';
cfg.funcolormap = cm;
ft_sourceplot(cfg,src_rft_int);


savefig(fullfile(plotpth,'ga_iaf_high_vs_low'))

close all

fig = openfig(fullfile(plotpth,'ga_iaf_high_vs_low'));
print(fig,fullfile(plotpth,'ga_iaf_high_vs_low_iaf'),'-dpng','-r600')
print(fig,fullfile(plotpth,'ga_iaf_high_vs_low_iaf'),'-dsvg','-r600')
delete(fullfile(plotpth,'ga_iaf_high_vs_low.fig'))

%% Alpha high vs low search


load(fullfile(outpth,subj{1},'dics_filt_iaf_search.mat'))
iaf_subj = cell(1,length(subj));
for s=1:length(subj)

%     load(fullfile(outpth,subj{s},'dics_filt_iaf.mat'))
% 
     load(fullfile(outpth,subj{s},'dics_iaf_high_low_search.mat'),'iaf_high_dics','iaf_low_dics')
     
     x = iaf_high_dics./iaf_low_dics - 1;
    iaf_subj{s} = iaf_dics;
    iaf_subj{s}.avg.pow = x;
    
    iaf_subj{s}.freq = 0;
    
    clear iaf_high_dics iaf_low_dics x
end

cfg = [];
cfg.parameter = 'pow';
ga_iaf = ft_sourcegrandaverage(cfg,iaf_subj{:});

ga_iaf = rmfield(ga_iaf,'cfg');


x = maxk(ga_iaf.pow,floor(sum(ga_iaf.inside)*0.01));
mask = ga_iaf.pow > x(end);
ga_iaf.mask = mask;

clear x mask

% interpolate
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'pow';
src_rft_int = ft_sourceinterpolate(cfg,ga_iaf,templmri);

cfg.parameter = 'mask';
src_rft_int_mask = ft_sourceinterpolate(cfg,ga_iaf,templmri);

src_rft_int.mask = src_rft_int_mask.mask;

clear src_rft_int_mask
src_rft_int.coordsys = 'mni';

cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolorlim = 'zeromax';
cfg.maskparameter = 'mask';
cfg.crosshair = 'no';
cfg.funcolormap = cm;
ft_sourceplot(cfg,src_rft_int);


savefig(fullfile(plotpth,'ga_iaf_high_vs_low'))

close all

fig = openfig(fullfile(plotpth,'ga_iaf_high_vs_low'));
print(fig,fullfile(plotpth,'ga_iaf_high_vs_low_iaf_search'),'-dpng','-r600')
print(fig,fullfile(plotpth,'ga_iaf_high_vs_low_iaf_search'),'-dsvg','-r600')
delete(fullfile(plotpth,'ga_iaf_high_vs_low.fig'))

close all