%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a2: plot behavioral performance per condition (Fig. 2 c&d)
% also export .csv file for statistical analyses in R

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Behavioural analyses
% a: performance per condition
% b. performance for alpha high vs low

clear all; close all; clc; beep off;
addpath('Z:\fieldtrip')
ft_defaults;
addpath('Z:\Visual Search RFT\Violinplot-Matlab-master')
pth = 'Z:\Visual Search RFT';
set(0,'defaultAxesFontSize',12,'defaultAxesFontName','Arial')

mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphapth = fullfile(pth,'results','meg','6 Alpha');
alphapowpth = fullfile(alphapth,'pow');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions','alpha RFT');
behavpth = fullfile(pth,'results','behavior');
behavfigpth = fullfile(pth,'results','behavior','fig');
mkdir(behavfigpth)
condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};
condi_label = {{'gui', '16'}, {'ung', '16'},{'gui', '32'},{'ung', '32'}};
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;


load(fullfile(behavpth,'sign_detect_condi.mat'))

fig = figure('Position',[0 0 1920/1.5 1080/2.5]);

subplot(121)
violinplot([mean_rt(:,2),mean_rt(:,1),mean_rt(:,4),mean_rt(:,3)],cellfun(@strjoin,condi_label,'Uniformoutput',false),'ViolinColor',col_palette,'ShowMean',true);
ylabel('mean RT (s)')
yticks(0:1:2)
ylim([0 2])
ytickformat('%.1f')

subplot(122)
violinplot([d_prime(:,2),d_prime(:,1),d_prime(:,4),d_prime(:,3)],cellfun(@strjoin,condi_label,'Uniformoutput',false),'ViolinColor',col_palette,'ShowMean',true);
ylabel('d''')
ylim([0 6])
yticks(0:3:6)
ytickformat('%.1f')

% subplot(133)
% violinplot([C(:,1),C(:,2),C(:,3),C(:,4)],cellfun(@strjoin,condi_label,'Uniformoutput',false),'ViolinColor',col_palette,'ShowMean',true);
% ylabel('criterion')
% ylim([-.4 2])
% yticks(0:1:2)
% ytickformat('%.1f')

print(fig,fullfile(behavfigpth,'perform_condition'),'-dpng','-r600')
print(fig,fullfile(behavfigpth,'perform_condition'),'-dsvg','-r600')


%% Create table for statistics

% reaction time
T = reshape(mean_rt,[],1);
% d'
T_d = reshape(d_prime,[],1);
% C
T_C = reshape(C,[],1);

condi_label = {{'ung', '16'},{'gui', '16'},{'ung', '32'},{'gui', '32'}};

%sanity check
[isequal(T(1:length(mean_rt)),mean_rt(:,1)), isequal(T(length(mean_rt)+1:length(mean_rt)*2),mean_rt(:,2)), isequal(T(length(mean_rt)*2+1:length(mean_rt)*3),mean_rt(:,3)), isequal(T(length(mean_rt)*3+1:length(mean_rt)*4),mean_rt(:,4))]

% fifth column: condition
T = num2cell([repmat([1:length(mean_rt)]',4,1),T,T_d,T_C,zeros(size(T,1),1)]);

% condition 1: ungui 16
T(1:length(mean_rt),5) =repmat({1},length(mean_rt),1);

% condition 2: gui 16
T(length(mean_rt)+1:length(mean_rt)*2,5) =  repmat({2},length(mean_rt),1);

% condition 4: ungui 32
T(length(mean_rt)*2+1:length(mean_rt)*3,5) =  repmat({3},length(mean_rt),1);

% condition 4: gui 32
T(length(mean_rt)*3+1:length(mean_rt)*4,5) =  repmat({4},length(mean_rt),1);


T = cell2table(T,'VariableNames',{'id','RT','d','C','condi'});
writetable(T,fullfile(behavpth,'perf_condi.csv'),'Delimiter',',')