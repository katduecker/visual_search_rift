%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a. align T1 image with digitized headshape
% uses ICP algorithm but needs to be adjusted manually

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Source analysis
% a: align digitized headshape to T1 -> realign Tq
% b: Forward model/Lead field
% c: DICS beamformer



% list subjects with T1
clear all; close all; clc
pth = 'Z:\Visual Search RFT';
dtpth = fullfile(pth,'data');
cohpth = fullfile(pth,'results','meg','5 COH hilb','coh','sinusoid','conditions','alpha RFT');
outpth = fullfile(pth,'results','meg','7 Beamformer');
addpath('Z:\fieldtrip')
ft_defaults;

% template mri
[ftver, ftdir] = ft_version;
templatedir = fullfile(ftdir, 'template','anatomy');

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

% Align T1 and headshape

for s = 1:length(subj)
   
    if exist(fullfile(outpth,subj{s},'t1_align.mat'),"file")
        continue
    end

  % if subject got individual T1 read in
    if exist(fullfile(dtpth,subj{s},'mri'))
        d = dir(fullfile(dtpth,subj{s},'mri'));
        d = {d.name};
        t1_file = (cell2mat(cellfun(@(x) ~isempty(x),regexp(d,'T1'),'UniformOutput',false)) ...
            + cell2mat(cellfun(@(x) ~isempty(x),regexp(d,'nii'),'UniformOutput',false))) == 2;
        % unzip
        x = gunzip(fullfile(dtpth,subj{s},'mri',d{t1_file}));


        % read in
        mri = ft_read_mri(x{1});
    % else, use template
    else
        mri = ft_read_mri(fullfile(templatedir,'single_subj_T1.nii'));
    end

    % read headshape from MEG file
    headshape = ft_read_headshape(fullfile(dtpth,subj{s},'meg','part1.fif'), 'unit', 'cm');

    % reslice
    mri_reslice = ft_volumereslice([],mri);
    % convert T1 scan to neuromag coordinate system

    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'neuromag';
    mri_realigned = ft_volumerealign(cfg,mri_reslice);
    mri_realigned = ft_convert_units(mri_realigned,'cm');
    close all;

    cfg = [];
    cfg.method = 'headshape';
    cfg.headshape.headshape = headshape;
    cfg.headshape.interactive = 'no';
    cfg.headshape.icp = 'yes';
    cfg.coordsys = 'neuromag';
    cfg.paramter = 'anatomy';
    cfg.viewresult = 'yes';

    mri_realigned2 = ft_volumerealign(cfg,mri_realigned);

    cfg = [];
    cfg.method = 'headshape';
    cfg.headshape.interactive = 'yes';
    cfg.headshape.headshape = headshape;
    mri_realigned3 = ft_volumerealign(cfg,mri_realigned2);

    mri_realigned3.coordsys = 'neuromag';
    close all;

    mkdir(fullfile(outpth,subj{s}));
    save(fullfile(outpth,subj{s},'t1_align.mat'),'mri_realigned3')

    clear mri*
end


