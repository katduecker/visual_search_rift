%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b. Find alpha SOI and individual alpha frequency
% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs
% - s: subject index

% Output
% - soi_grad: gradiometers with high alpha power
% iaf_grad: identified IAF in gradiometers

% [c] K. Duecker, PhD candidate Neuronal Oscillations group
% last changed: 23/03/2022
% katharina.duecker@gmail.com

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

function b_find_iaf_soi(s)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

inpth = fullfile(pth,'results','meg','6 Alpha','pow');

outpth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;


load(fullfile(inpth,subj{s},'data_winl_5.mat'),'TFR_alpha_avg')

cfg = [];
cfg.method = 'sum';
TFR_alpha_avg = ft_combineplanar(cfg,TFR_alpha_avg);

cfg = [];
cfg.channel = 'MEGGRAD';
cfg.latency = [-1 0];
cfg.avgovertime = 'yes';
cfg.frequency = [4 14];
TFR_alpha_avg = ft_selectdata(cfg,TFR_alpha_avg);

for sens = 1:length(TFR_alpha_avg.label)
    [v{sens},p{sens}] = findpeaks(TFR_alpha_avg.powspctrm(sens,:));
end

f = zeros(1,length(TFR_alpha_avg.label));
powf = zeros(1,length(TFR_alpha_avg.label));
for sens = 1:length(TFR_alpha_avg.label)
    if ~isempty(v{sens})
        if length(p{sens})
            [~,x] = max(v{sens});
            f(sens) = TFR_alpha_avg.freq(p{sens}(x));
            powf(sens) = v{sens}(x);
        else
            f(sens) = TFR_alpha_avg.freq(p{sens});
            powf(sens) = v{sens};
        end
    end
end

plot(TFR_alpha_avg.freq,TFR_alpha_avg.powspctrm)
[~,soi_idx] = maxk(powf,4);

soi_grad_cmb = TFR_alpha_avg.label(soi_idx);

% not combined
soi_grad = {};
for sg = 1:length(soi_grad_cmb)
    soi_grad = [soi_grad,{soi_grad_cmb{sg}(1:7)},{['MEG',soi_grad_cmb{sg}(9:end)]}];
end

iaf_grad = mode(f(soi_idx));

mkdir(fullfile(outpth,subj{s}))

save(fullfile(outpth,subj{s},'iaf_soi_most_common.mat'),'soi_grad_cmb','soi_grad','iaf_grad')

