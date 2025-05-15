%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a2. plot coherence for each subject (-> Supplementary figure 1);
% grandaverage (Fig 3a)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023


%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)



%% settings
clear all; close all; clc; beep off;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'SNR');

cohfigpth = fullfile(pth,'results','meg','5 COH hilb', 'fig', '0 SNR');
mkdir(cohfigpth)
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
datpth = fullfile(pth,'data');             % max filter
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

rmpath(genpath('/rds/projects/2018/jenseno-entrainment/fieldtrip'))
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))
addpath(fullfile(pth,'matlab scripts/',"cbrewer/"))
cm = cbrewer('div','RdBu',101);
cm = flipud(cm);

ft_defaults;

toi = [-2.5 2];                             % start and end of trial in sec
avgtoi = 0.5;
fw = 2;                                     % bandwidth bp filter
fs = 1000;
foi = 50:75;
% list subjects
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

addpath(fullfile(pth, 'raacampbell-shadedErrorBar-19cf3fe/'))
%% load in example subject to get ft structure
d = dir(fullfile(datpth,subj{1},'meg'));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);
% trial structure to load in trl
load(fullfile(mergepth, subj{1},'trl_overlap_meg_el_rsp.mat'))

trlstruct{1} = [meginfo.alltrl_bl{1}(:,3)-fs*2.5,meginfo.alltrl_bl{1}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{1}),1)-2.5*fs];
trlstruct{1}(trlstruct{1}(:,1) <0,1) = 1;

cfg = [];
cfg.dataset = fullfile(datpth,subj{1},'meg',f{1});
cfg.preproc.detrend = 'yes';
cfg.trl = trlstruct{1}(1,:);
cfg.channel = {'MEG','MISC004','MISC005'};
% load in data for this part
data = ft_preprocessing(cfg);

% fourier transform
cfg            = [];
cfg.output     = 'pow';
cfg.method     = 'mtmconvol';
%cfg.pad        = 'nextpow2';
cfg.taper      = 'hanning';
cfg.toi        = -1.5:2;
cfg.foi        = foi;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*.1;
cfg.tapsmofrq  = 3;
cfg.keeptrials = 'no';
cfg.channel    = {'MEG', 'MISC004','MISC005'};
freq       = ft_freqanalysis(cfg, data);

cfg = [];
cfg.channel = {'MEGGRAD', 'MISC004','MISC005'};
freqgrad = ft_selectdata(cfg,freq);
cfg.channel = {'MEGMAG', 'MISC004','MISC005'};
freqmag        = ft_selectdata(cfg, data);
freqgrad.time  = toi(1):1/fs:toi(2);
freqgrad.freq  = foi;
freqmag.time  = toi(1):1/fs:toi(2);
freqmag.freq  = foi;
freqmag.dimord = freqgrad.dimord;
clear data

%% Grads

filename = 'SNR_MEGGRAD.mat';

% Topos
cfg = [];
cfg.marker = 'off';
cfg.zlim = 'zeromax';
cfg.comment = 'no';
cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
cfg.colorbar = 'yes';
cfg.ylim = [60 60];

coh_all60 = cell(1,length(subj));
coh_all67 = coh_all60;
for s = 1:length(subj)
    load(fullfile(cohpth,subj{s},filename))                    % load coherence

    coh_all60{s} = freqgrad;
    coh_all67{s} = freqgrad;
    coh_all60{s}.powspctrm = coh60;
    coh_all67{s}.powspctrm = coh67;
    clear coh60 coh67
end



% % This is Supplementary Fig. 1
% fig = figure('Position',[0 0 920 1080]);
% for s = 1:length(subj)
%         load(fullfile(soipth,subj{s},'soi_stat.mat'))
% 
%         coh_for_plot = freqgrad;
%     % paste into freq structure and plot
%     coh_for_plot.powspctrm = coh_all60{s}.*10^2;
% 
%     cfg.layout = 'neuromag306planar_helmet.mat';
%     cfg.xlim = [0.0 0.5];
%     cfg.highlight = 'on';
%     cfg.highlightsymbol = 'o';
%     cfg.highlightcolor = [5,113,176]./255;
%     cfg.highlightsize = 4;
%     cfg.highlightchannel = find(ismember(coh_for_plot.label,soigrad));
%     subplot(8,4,s)
%     ft_topoplotTFR(cfg,coh_for_plot)
%     cb = colorbar;
%     cb.Ticks = cb.Limits;
%     cb.TickLabels = strsplit(sprintf('%0.1f ',cb.Ticks));
% 
% end
% 
% print(fig,fullfile(cohfigpth,'topo_coh60_allsubj_grad_fwdth35'),'-dpng','-r300')
% print(fig,fullfile(cohfigpth,'topo_coh60_allsubj_grad_fwdth35'),'-dsvg','-r600','-painters')
% 
% close all

coh_soi60 = cell(1,length(subj));
coh_soi67 = cell(1, length(subj));
coh_spec60 = cell(1, length(subj));
coh_spec67 = cell(1, length(subj));

fig = figure('Position',[0 0 920 1080]);
for s = 1:length(subj)
    
    % load SOI
    load(fullfile(soipth,subj{s},'soi_stat.mat'))
    soimag = soi_stat(logical(cell2mat(cellfun(@(x) strcmp(x(end),'1'),soi_stat,'UniformOutput',false))));
    soigrad = soi_stat(~ismember(soi_stat,soimag));
    
    [~,startsamp] = min(abs(freqgrad.time - 0.2));
    [~,endsamp] = min(abs(freqgrad.time - .5));     % minimum RT is about 500ms

    
    % paste into freq structure and plot
    coh_soi60{s} = mean(coh_all60{s}.powspctrm(ismember(freqgrad.label,soigrad),:,:),1,'omitnan');
    coh_spec60{s} = mean(coh_soi60{s}(:,:,startsamp:endsamp),3,'omitnan');

%     subplot(9,4,s)
%     plot(freqgrad.freq,coh_spec60{s})
%     hold on
    
    % paste into freq structure and plot
    coh_soi67{s} = mean(coh_all67{s}.powspctrm(ismember(freqgrad.label,soigrad),:,:),1,'omitnan');
    coh_spec67{s} = mean(coh_soi67{s}(:,:,startsamp:endsamp),3,'omitnan');

%     subplot(9,4,s)
%     plot(freqgrad.freq,coh_spec67{s})
end

% print(fig,fullfile(cohfigpth,'coh_allsubj_grad_spect_fwdth35'),'-dpng','-r300')
% print(fig,fullfile(cohfigpth,'coh_allsubj_grad_spect_fwdth35'),'-dsvg','-r600')
% 


%% Grandavergae (Fig 3a)

GA60 = ft_freqgrandaverage([],coh_all60{:});
GA67 = ft_freqgrandaverage([], coh_all67{:});

ga_spec60 = vertcat(coh_spec60{:});
ga_spec67 = vertcat(coh_spec67{:});

ga_tfr60 = mean(vertcat(coh_soi60{:}));
ga_tfr67 = mean(vertcat(coh_soi67{:}));

ga_tfr = squeeze(mean(vertcat(coh_soi60{:}, coh_soi67{:})));

close all
fig = figure('Position', [0, 0, 1980, 1080/3]);
set(gca, 'FontName', 'Arial')
subplot(131)
cfg = [];
cfg.marker = 'off';
% cfg.highlight = 'labels';
% cfg.highlightsymbol = 'o';
% cfg.highlightsize = 4;
% cfg.highlightcolor = [0 0.5 0];
cfg.zlim = [0, 0.02];
cfg.comment = 'no';
cfg.colormap = cm(floor(length(cm)/2)+1:end,:);
cfg.ylim = [60 60];
cfg.xlim = [0 .5];
cfg.zlim = 'maxmin';
cfg.layout = 'neuromag306planar_helmet.mat';
ft_topoplotTFR(cfg,GA60)   

cb = colorbar;
cb.Ticks = [0, 0.01, 0.02];

ft_topoplotTFR(cfg, GA60)

subplot(132)
% shadedErrorBar(freqgrad.freq, ga_spec60, {@mean, @std}, 'lineprops', {'-', 'MarkerFaceColor', [0 0.447 0.741]})
% hold on
% shadedErrorBar(freqgrad.freq, ga_spec67, {@mean, @std}, 'lineprops', {'-', 'MarkerFaceColor', [0.85 0.325 0.098]})

plot(freqgrad.freq, mean(ga_spec60))
hold on
plot(freqgrad.freq, mean(ga_spec67))
ylabel('coherence MEG - diode')
xlabel('frequency (Hz)')
legend('60 Hz', '67 Hz')
ylim([0, 0.02])
yticks([0, 0.01, 0.02])

subplot(133)
imagesc(freqgrad.time, freqgrad.freq, squeeze(ga_tfr60))
axis xy
ylabel('frequency (Hz)')
xlabel('time (s)')
xlim([-0.2, 0.5])
caxis([0, 0.02])

cb = colorbar;
cb.Ticks = [0, 0.01, 0.02];


print(fig,fullfile(cohfigpth,'topo_GA_grad_60_fwdth5'),'-dpng','-r300')

print(fig,fullfile(cohfigpth,'topo_GA_grad_60_fwdth5'),'-dsvg', '-painters')

