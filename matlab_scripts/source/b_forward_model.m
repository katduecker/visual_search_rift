%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b. estimate leadfield using aligned image for each participant

% Input
% - s: subject index

% Output:
% - sourcemodel
% - headmodel
% - lead field

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Source analysis
% a: align digitized headshape to T1 -> realign Tq
% b: Forward model/Lead field
% c: DICS beamformer

function b_forward_model(s)

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
dtpth = fullfile(pth,'data');
mxpth = fullfile(pth,'results','meg','1 maxfilter','1 maxfilter');
cohpth = fullfile(pth,'results','meg','5 COH hilb','coh','sinusoid','conditions','alpha RFT');
outpth = fullfile(pth,'results','meg','7 Beamformer/');
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

% location where fieldtrip is installed
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template', 'sourcemodel');
% load 8 mmm sourcemodel
template = load(fullfile(templatedir, 'standard_sourcemodel3d8mm'));

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

% load T1 scan aligned with headshape/fiducials (digitized during MEG)
load(fullfile(outpth,subj{s},'t1_align.mat'))

% segment brain
cfg = [];
cfg.output = {'brain'};
mri_segm = ft_volumesegment(cfg,mri_realigned3);


% warp mri to mni
cfg = [];
cfg.warpmni    = 'yes';
cfg.template   = template.sourcemodel; % source grid
cfg.nonlinear  = 'yes';
cfg.mri        = mri_realigned3;   % aligned, resliced mri
sourcemdl      = ft_prepare_sourcemodel(cfg);


% Headmodel & leadfield for each part

% list recording files (typically 6 per subject)
d = dir(fullfile(mxpth,subj{s}));
d = {d.name};
files = d(cell2mat(cellfun(@(x) ~isempty(x),regexp(d,'_sss.fif'),'UniformOutput',false)));


% load sensors for each experiment file & compute headmodel & leadfield
% grad structure (sensor/coil positions)

% average over positions
grad = [];
for fl = 1:length(files)
    grad = [grad;ft_read_sens(fullfile(mxpth,subj{s},files{fl}))];
end
% 
% mean channel position
mGrad = grad(1);
% average grad structure
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad(g).chanpos;
mGrad.coilpos = mGrad.coilpos + grad(g).coilpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
mGrad.coilpos = mGrad.coilpos./length(grad);
clear grad


% headmodel
cfg        = [];
cfg.grad   = mGrad;
cfg.method = 'singleshell';
cfg.channel = 'MEGGRAD';
headmodel = ft_prepare_headmodel(cfg,mri_segm);
headmodel = ft_convert_units(headmodel,'cm');

% leadfield
cfg = [];
cfg.grad = mGrad;
cfg.grid = sourcemdl;
cfg.headmodel = headmodel;
cfg.channel = 'MEGGRAD';

leadf= ft_prepare_leadfield(cfg);

delete(fullfile(outpth,subj{s},'head_leadf_aligned.mat'))

save(fullfile(outpth,subj{s},'head_leadf_avg.mat'),'headmodel','leadf','mGrad','-v7.3')

