%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 1: Guided Visual Search is associated with a feature-based priority
% map in early visual cortex

%% Single trial coherence
% Plot average of single trial coherence 

% (c), Katharina Duecker
% last edited, Nov-29-2024

clear all;
rmpath(genpath('/rds/projects/2018/jenseno-entrainment/fieldtrip'))
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
datpth = fullfile(pth,'data');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

cohpth = fullfile(pth,'results','meg','8 COH single trl');

%% Template struct
d = dir(fullfile(datpth,subj{1},'meg'));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);
fs = 1000;
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(datpth,subj{1},'meg',f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = {'MEGGRAD'};
% load in data for this part
data = ft_preprocessing(cfg);

ERP = ft_timelockanalysis([], data);
%% load coherence

RIFT_subj = cell(1,31);
for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},'coh_single_trial_mcohere.mat'))
    
    coh_T = ERP;
    coh_T.avg = mean(cohT_bsl);
    coh_T.var = std(cohT);
    coh_T.dof = ones(size(coh_T.avg));
        
    coh_D = ERP;
    coh_D.avg = mean(cohD_bsl);
    coh_D.var = std(cohD);
    coh_D.dof = ones(size(coh_D.avg));
    
    % change layout to please fieldtrip
    coh_T.avg = repmat(coh_T.avg,1,1,2);
    coh_T.time = [1, 2];
    coh_D.avg = repmat(coh_D.avg,1,1,2);
    coh_D.time = [1, 2];
    
    cfg = [];
    cfg.method = 'sum';
    coh_T = ft_combineplanar(cfg,coh_T);
    coh_D = ft_combineplanar(cfg,coh_D);
   
    rift = ft_timelockgrandaverage([],coh_T,coh_D);
    
    RIFT_subj{s} = rift;
end

RIFT_GA = ft_timelockgrandaverage([],RIFT_subj{:});

addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);
% 
fig = figure('Position',[0 0 1940/2 1080/3]);
cfg = [];
cfg.layout = 'neuromag306cmb_helmet.mat';
cfg.marker = 'off';
cfg.figure = 'gca';
cfg.comment = 'no';
cfg.xlim = [1,1];
ft_topoplotER(cfg,RIFT_GA)
%caxis([0, 0.014])
colormap(cm)
cb = colorbar;
% cb.Limits = [0, 0.014];
% cb.Ticks =[0, 0.014];
mkdir(fullfile(cohpth,'fig'))
print(fig,fullfile(cohpth,'fig','single_trial_coh_mcohere_planar'),'-dsvg')
