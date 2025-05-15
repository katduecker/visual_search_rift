%% simulation GLM
% simulate GLM with surrogate signal to find best settings for GLM

clear all; close all; clc;
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha');
soipth = fullfile(alphapth,'iaf_soi');

cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');
outpth = fullfile(alphapth,'GLM RT');

% load frequencies to which IAF vector should be aligned to
load(fullfile(pth,'matlab scripts','alpha','alpha_align_vec.mat'))

s = 1;
%% load data condition
d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];

load(fullfile(inpth,subj{s},file{1}),'behav_array')

% extract specs condition
load(fullfile(pth,'experiment','trigdef.mat'))

rt = [behav_array{:,3}]';

%% condition specifier: guided, set size 32, target present
condi_specs = {'ti','32t','tp'};

condi = zeros(length(behav_array),3);

% find trigger condition
% store: 1: guided; 1: set 32; 1: target present

for c = 1:length(condi_specs)
    trig = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x), strfind(trigdef(:,2),condi_specs{c}),'UniformOutput',false)),1));
       
    idx = ismember([behav_array{:,1}],trig);
    condi(:,c) = idx;
end

tot = [1:length(behav_array)]./length(behav_array);     % time on task

%% Create alpha distribution

num_it = 1000;
model_perf = zeros(num_it,2,3);
b = zeros(3,num_it);

for i = 1:num_it

% simulate rt
b0 = .5;
b1 = randn(1);
b2 = randn(1);
b(:,i) = [b0;b1;b2];

lin_pred = b0 + b1*tot' + b2*rt;

% log link function
k = exp(2*lin_pred);

% synthetic RT
theta = 3; % Scale parameter
iaf = .5e-13 + gamrnd(k, theta).*10^(-13);

% log link function

X = [ones(length(behav_array),1),tot',rt];

beta_pow = X\iaf;

X = [ones(length(behav_array),1),tot',rt];
beta_zpow = X\zscore(iaf);

X = [ones(length(behav_array),1),tot',rt];

beta_log = X\log(iaf);

model_perf(i,:,1) = b(2:3,i) - beta_pow(2:3);
model_perf(i,:,2) = b(2:3,i) - beta_zpow(2:3);
model_perf(i,:,3) = b(2:3,i) - beta_log(2:3);

end


fig = figure;
subplot(211)
histogram(rt, 'NumBins', 100)
title('RT example participant')
subplot(212)
histogram(iaf, 'NumBins', 100)
title('simulated alpha power')


figure;
boxplot(squeeze(model_perf(:,2,:)))
ylabel('diff(b - b_hat)')
xticklabels({'pow', 'z(pow)', 'log(pow)'})


%% Simulate based on alpha

rmpath(genpath('/rds/projects/2018/jenseno-entrainment'))
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha','pow');

load(fullfile(alphapth,subj{1},'data_fourier_winl_10.mat'),'TFR')

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);

% select time of interest
cfg = [];
cfg.latency = [-1.5 .5];
cfg.frequency = 10;
IAF = ft_selectdata(cfg,TFR);

four_pow = squeeze(IAF.powspctrm(:,1,1,1));
z_four = zscore(four_pow);

load(fullfile(alphapth,subj{1},'data_winl_10.mat'),'TFR')

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);

% select time of interest
cfg = [];
cfg.latency = [-1.5 .5];
cfg.frequency = 10;
IAF = ft_selectdata(cfg,TFR);

pow = squeeze(IAF.powspctrm(:,1,1,1)).*10^12;

figure;
subplot(311)
histogram(pow,'NumBins', 100)
title('pow')
subplot(312)
histogram(four_pow,'NumBins',100)
title('fourier')
subplot(313)
histogram(z_four,'NumBins',100)
title('zscore fourier')

close all
num_it = 1000;
model_perf = zeros(num_it,2,4);
b = zeros(3,num_it);

for i = 1:num_it

% simulate rt
b0 = .5;
b1 = randn(1);
b2 = randn(1);
b(:,i) = [b0;b1;b2];

lin_pred = b0 + b1*tot' + b2*(pow);

% log link function
k = exp(lin_pred);

% synthetic RT
theta = 2; % Scale parameter
rt = gamrnd(k, theta);

X = [ones(length(behav_array),1),tot',pow];

beta_pow = X\rt;

X = [ones(length(behav_array),1),tot',zscore(pow)];
beta_zpow = X\rt;

X = [ones(length(behav_array),1),tot',four_pow];

beta_fou = X\rt;

X = [ones(length(behav_array),1),tot',z_four];
beta_zfou = X\rt;

model_perf(i,:,1) = b(2:3,i) - beta_pow(2:3);
model_perf(i,:,2) = b(2:3,i) - beta_fou(2:3);
model_perf(i,:,3) = b(2:3,i) - beta_zpow(2:3);
model_perf(i,:,4) = b(2:3,i) - beta_zfou(2:3);

end

figure
scatter(ones(1,length(model_perf)).*0.5,squeeze(model_perf(:,2,1)),'filled')
hold on
scatter(ones(1,length(model_perf)).*75, squeeze(model_perf(:,2,2)),'filled')
scatter(ones(1,length(model_perf)), squeeze(model_perf(:,2,3)),'filled')
xlim([0.25 1.25])

figure;
subplot(121)
boxplot(squeeze(model_perf(:,2,1:2)))
ylabel('diff(b - b_hat)')
xticklabels({'pow', 'fourier'})
subplot(122)
boxplot(squeeze(model_perf(:,2,3:4)))
xticklabels({'z(pow)', 'z(fourier)'})

std(squeeze(model_perf(:,2,1:2)))
std(squeeze(model_perf(:,2,3:4)))

[h,p,ci,stats]=ttest(squeeze(model_perf(:,2,3))-squeeze(model_perf(:,2,4)))

[h,p,ci,stats]=ttest(squeeze(model_perf(:,2,1))-squeeze(model_perf(:,2,2)))
mean(squeeze(model_perf(:,2,1))-squeeze(model_perf(:,2,2)))
% assuming that RT depends on power, the fit is best when using z-scored
% fourier spectrum

%% RIFT alpha

cohpth = fullfile(pth,'results','meg','8 COH single trl');
cohsoipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');

load(fullfile(cohpth,subj{1},'coh_single_trial_pearson.mat'))
load(fullfile(cohsoipth,subj{s},'soi_stat.mat'))
RIFT = corrT.avg(:,1)
figure
histogram(RIFT,'NumBins',100)

%% Generate correlated data

b = zeros(3,num_it);
model_perf = zeros(num_it,2,4);
for i = 1:num_it
    b0 = .3;
    b1 = randn(1);
    b2  = randn(1);
    RIFTsim = RIFT.*pow*b2*10^13+RIFT.*tot.*b1+b0;
    b(:,i) = [b0;b1;b2];
    X = [ones(size(RIFTsim,1),1),tot,pow];
    beta = X\RIFTsim;
    model_perf(i,:,1) = b(2:3,i)*10^13 - beta(2:3);
    
    X = [ones(size(RIFTsim,1),1),tot,four_pow];
    beta = X\RIFTsim;
    model_perf(i,:,2) = b(2:3,i)*10^13 - beta(2:3);
    
    X = [ones(size(RIFTsim,1),1),tot,zscore(pow)];
    beta = X\RIFTsim;
    model_perf(i,:,3) = b(2:3,i) - beta(2:3);
    
    X = [ones(size(RIFTsim,1),1),tot,zscore(four_pow)];
    beta = X\RIFTsim;
    model_perf(i,:,4) = b(2:3,i) - beta(2:3);
end
histogram(RIFTsim,'NumBins',100)

figure
subplot(121)
boxplot(squeeze(model_perf(:,2,1:2)))
subplot(122)
boxplot(squeeze(model_perf(:,2,3:4)))

std(squeeze(model_perf(:,2,1:2)))
std(squeeze(model_perf(:,2,3:4)))

% for RIFT zscored power vs power doesn't matter!
[h,p,ci,stats]=ttest(squeeze(model_perf(:,2,3))-squeeze(model_perf(:,2,4)))