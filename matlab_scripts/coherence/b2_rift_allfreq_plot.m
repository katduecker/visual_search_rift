%% settings
clear all; close all; clc;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
outpth = fullfile(pth,'results','meg','5 COH hilb','coh','conditions');
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

condi_all = {{'ti','32t','6067'},{'ti','32t','6760'}};
fwdth = 5;
filttype = {'but','twopass'};
% select current condition

timevec = linspace(-0.5,0.5,1001);
freqvec = 40:80;

cohT6067 = zeros(length(subj),length(freqvec),length(timevec));
cohD6067 = zeros(length(subj),length(freqvec),length(timevec));

cohT6760 = cohT6067;
cohD6760 = cohT6760;

for s = 1:length(subj)
    
    condi = condi_all{1};
    load(fullfile(outpth,subj{s},strcat('coh_',strjoin(condi,'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'_allfreq.mat')))
    
    cohT6067(s,:,:) = squeeze(mean(coh.cohTgrad(1:end-2,:,2000:3000),1));
    cohD6067(s,:,:) = squeeze(mean(coh.cohDgrad(1:end-2,:,2000:3000),1));
    
    condi = condi_all{2};
    load(fullfile(outpth,subj{s},strcat('coh_',strjoin(condi,'_'),'_freqw_',num2str(fwdth),'_',strjoin(filttype,'_'),'_allfreq.mat')))
    
    cohT6760(s,:,:) = squeeze(mean(coh.cohTgrad(1:end-2,:,2000:3000),1));
    cohD6760(s,:,:) = squeeze(mean(coh.cohDgrad(1:end-2,:,2000:3000),1));
end

cohT6067 = squeeze(mean(cohT6067,1));
cohD6067 = squeeze(mean(cohD6067,1));
cohT6760 = squeeze(mean(cohT6760,1));
cohD6760 = squeeze(mean(cohD6760,1));


figure('Position',[0,0,1940 540])
subplot(141)
imagesc(timevec,freqvec,cohT6067)
caxis([0, 0.05])
colorbar
axis xy
subplot(142)
imagesc(timevec,freqvec,cohD6067)
caxis([0, 0.05])
colorbar
axis xy
subplot(143)
imagesc(timevec,freqvec,cohT6760)
caxis([0, 0.05])
colorbar
axis xy
subplot(144)
imagesc(timevec,freqvec,cohD6760)
caxis([0, 0.05])
colorbar
axis xy


f1 = freqvec == 60;
f2 = freqvec == 67;

figure;
plot(timevec,cohT6067(f1,:),'k')
hold on
plot(timevec,cohD6067(f2,:),'Color',[0.5,0.5,0.5])
plot(timevec,cohD6760(f1,:),'b')
plot(timevec,cohT6760(f2,:),'Color','c')
legend('T60','D67','D60','T67')

figure;
plot(freqvec,mean(cohT6067(:,700:end),2))
hold on
plot(freqvec,mean(cohT6760(:,700:end),2))

plot(foi,mean(cohD6067(:,700:end),2))

