%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a1: behavioral performance per condition

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Behavioural analyses
% a: performance per condition
% b. performance for alpha high vs low

%% set up paths
clear all; close all; clc; beep off;
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/RT')
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';

mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphapth = fullfile(pth,'results','meg','6 Alpha');
alphapowpth = fullfile(alphapth,'pow');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions','alpha RFT');
behavpth = fullfile(pth,'results','behavior');

condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));


% load trigger specs
load(fullfile(pth, 'experiment','trigdef.mat'))
d_prime = zeros(length(subj),4);
C = zeros(length(subj),4);
mean_rt = zeros(length(subj),4);
min_rt = zeros(length(subj),4);
d_prime_data = zeros(length(subj),2,4);

load(fullfile(pth, 'experiment','trigdef.mat'))


for s = 1:length(subj)
    % loop over subjects
    % select behavior that went into TFR analysis
    load(fullfile(alphapowpth,subj{s},'data_winl_5.mat'),'perf_TFR')
    
    %% both correct present and absent are encoded as "hit" -> change correct target absent trials into "ca" (correct absent)
    ta_idx = cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),'ta'),'UniformOutput',false));
    
    % find "hit" target absent trials
    h_ta_idx = (ismember([perf_TFR{:,1}],[trigdef{ta_idx}])' + strcmp(perf_TFR(:,2),'h')) == 2;
    perf_TFR(h_ta_idx,2) = {'ca'};
    
    save(fullfile(alphapowpth,subj{s},'data_winl_5.mat'),'perf_TFR','-append')
   
    %% Select trials
    % load trial indices that are included in alpha-RFT analyses (randomly
    % selected to ensure equal number of trials)
    
    for c = 1:4
        
         c_idx = (cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{1}),'UniformOutput',false)) + ...
                    cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),condi{c}{2}),'UniformOutput',false))) == 2;

        cur_trig = vertcat(trigdef{c_idx,1});
        
        trl_idx = ismember(vertcat(perf_TFR{:,1}),cur_trig);

        % select trials
        perf_trl = perf_TFR(trl_idx,:);
        

        %% Signal detection measures
        
        % proportion hits
        n_h = sum(strcmp(perf_trl(:,2),'h'))/(sum(strcmp(perf_trl(:,2),'h'))+sum(strcmp(perf_trl(:,2),'m')));
       
        % proportion false alarm
        n_fa = sum(strcmp(perf_trl(:,2),'fa'))/(sum(strcmp(perf_trl(:,2),'ca'))+sum(strcmp(perf_trl(:,2),'fa')));
        
        % adjust extreme values https://stats.stackexchange.com/questions/134779/d-prime-with-100-hit-rate-probability-and-0-false-alarm-probability
        if n_fa == 0
            n = (sum(strcmp(perf_trl(:,2),'ca'))+sum(strcmp(perf_trl(:,2),'fa')));
            n_fa = 0.5/n;
        end
        if n_h == 1
            n = (sum(strcmp(perf_trl(:,2),'h'))+sum(strcmp(perf_trl(:,2),'m')));
            n_h = (n-0.5)/n;
        end
        
        % z transform all scores
        
        z_nh = norminv(n_h,0,1);
        z_nfa = norminv(n_fa,0,1);
        
        d_prime_data(s,1,c) = n_h;
        d_prime_data(s,2,c) = n_fa;
        
        % sensitivity
        d_prime(s,c) = z_nh-z_nfa;
        
        % criterion: https://www.birmingham.ac.uk/Documents/college-les/psych/vision-laboratory/sdtintro.pdf
        % position shift between false alarm and hits
        C(s,c) = - (z_nh+z_nfa)/2;
        
        %% reaction time
        
        h_idx = strcmp(perf_trl(:,2),'h');
        mean_rt(s,c) = mean([perf_trl{h_idx,3}]);
        min_rt(s,c) = min([perf_trl{h_idx,3}]);
        clear perf_trl
    end
end

save(fullfile(behavpth,'sign_detect_condi.mat'),'d_prime','C','mean_rt','min_rt')

% transform d prime data into table and write to excel

T_dp = num2cell([squeeze(d_prime_data(:,:,1));squeeze(d_prime_data(:,:,2));squeeze(d_prime_data(:,:,3));squeeze(d_prime_data(:,:,4))]);
T_dp = [num2cell(repmat([1:length(subj)]',4,1)), T_dp,[repmat({strjoin(condi{1},'_')},length(subj),1);repmat({strjoin(condi{2},'_')},length(subj),1);...
    repmat({strjoin(condi{3},'_')},length(subj),1);repmat({strjoin(condi{4},'_')},length(subj),1)]];

T_dp = cell2table(T_dp,'VariableNames',{'subj','p_hits','p_miss','condi'});
writetable(T_dp,fullfile(behavpth,'data_dprime_condi.csv'),'Delimiter',',')