%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d2. Alpha power per condition/subject (Supplementary Fig. 5)
% prepare statistical test
% violin plots

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

%% Alpha analysis pipeline
% a. TFR of alpha power over all trials (keep trials) & average trials
% b. Identify IAF and SOI
% c. align TFR to IAF
% d. contrast conditions
% e. control analysis: compare alpha for fast vs slow

clear all; close all; clc; beep off;
condi_all = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},...
    {'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;

condi_label = {'ung, 16', 'gui, 16', 'ung, 32', 'gui,32'};

%% settings
pth = 'Z:\Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
alphapth = fullfile(pth,'results','meg','6 Alpha');
outpth = fullfile(alphapth,'pow align iaf');
soipth = fullfile(alphapth,'iaf_soi');
addpath('Z:\Visual Search RFT\Violinplot-Matlab-master')


load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

alpha_pow_bsl = zeros(length(subj),length(condi_all));
alpha_pow_search = zeros(length(subj),length(condi_all));
alpha_pow_time = zeros(length(subj),length(condi_all),56);


T_bsl = zeros(length(subj)*length(condi_label),3);
T_search = zeros(length(subj)*length(condi_all),4);

h = 0;
for s = 1:length(subj)

    for c = 1:length(condi_all)
        
        h = h+1;
        condname = strjoin(condi_all{c},'_');

        load(fullfile(outpth,subj{s},[condname,'.mat']))

        alpha_pow_bsl(s,c) = iaf_bsl;

        T_search(h,1) = s;

        % set size = 32?
        if ~isempty(intersect(c,[3,4,7,8]))
            T_search(h,2) = 1;
        end


        % unguided/guided
        if ~mod(c,2)
            T_search(h,3) = 1;
        end

        % target present?
        if ~isempty(intersect(c,5:8))
            T_search(h,4) = 1;
        end

        alpha_pow_search(s,c) = iaf_search;
        T_search(h,5) = iaf_search;

        alpha_pow_time(s,c,:) = (iaf_time - mean(iaf_time))./std(iaf_time);


        clear hit_rate rt iaf_*
    end
end

% average alpha power over target absent and present
alpha_pow_bsl = (alpha_pow_bsl(:,1:4)+alpha_pow_bsl(:,5:8))./2;


h = 0;
for s = 1:length(subj)

    for c = 1:length(condi_label)
        
        h = h+1;

        T_bsl(h,1) = s;

        % set size = 32?
        if ~isempty(intersect(c,[3,4]))
            T_bsl(h,2) = 1;
        end


        % unguided/guided
        if ~mod(c,2)
            T_bsl(h,3) = 1;
        end

        T_bsl(h,4) = alpha_pow_bsl(s,c);


    end
end

col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;


fig = figure('Position',[0 0 1920 1080/2.5]);
subplot(131)
% baseline
violinplot(log(alpha_pow_bsl),condi_label,'ViolinColor',col_palette,'ShowMean',true);

ylabel('log(power) at IAF')
yticks(-56:4:-48)
ylim([-56 -48])

subplot(132)
% search
violinplot(log(alpha_pow_search(:,1:4)),condi_label,'ViolinColor',col_palette,'ShowMean',true);
xlabel('Target absent')

ylabel('log(power) at IAF')
yticks(-56:4:-48)
ylim([-56 -48])


subplot(133)
% search
violinplot(log(alpha_pow_search(:,5:end)),condi_label,'ViolinColor',col_palette,'ShowMean',true);
xlabel('Target present')

ylabel('log(power) at IAF')
yticks(-56:4:-48)
ylim([-56 -48])

print(fig,fullfile(alphapth,'fig','alpha_condition_avg_time'),'-dsvg')
print(fig,fullfile(alphapth,'fig','alpha_condition_avg_time'),'-dpng')
close all


T_bsl = array2table(T_bsl,'VariableNames',{'id','set size','guided','iafpow'});
T_search = array2table(T_search,'VariableNames',{'id','set size','guided','tp','iafpow'});

writetable(T_bsl,fullfile(alphapth,'iaf_pow_bsl_condi.csv'))
writetable(T_search,fullfile(alphapth,'iaf_pow_search_condi.csv'))




%%% Alpha power over time

alpha_pow_time = (alpha_pow_time(:,1:4,:)+alpha_pow_time(:,5:end,:))./2;


fig = figure('Position',[0 0 1920/2 1080/2.5]);
avg_alpha_time = squeeze(mean(alpha_pow_time,1));
timevec = -1.75:0.05:1;
col_palette = [0,114,178;86,180,233;213,94,0;230,159,0]./255;

for c = 1:length(condi_label)
    plot(timevec,avg_alpha_time(c,:),'Color',col_palette(c,:),'LineWidth',1.5)
    hold on
end
xlim([-1.5 0.5])
xlabel('time (s)')
ylabel('alpha power ')
print(fig,fullfile(alphapth,'fig','alpha_condition_time'),'-dsvg')
print(fig,fullfile(alphapth,'fig','alpha_condition_time'),'-dpng')
