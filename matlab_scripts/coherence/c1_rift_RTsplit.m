%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c1. coherence per condition (guided, unguided; set size) split for fast
% and slow trials

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Input:
% - s: subject index
% - c_idx: condition index
% - fwdth: bandwith of bpfilter applied before Hilbert transform
% (recommended: 5 Hz)
% - filttype: {'but','twopass'}
% - data_trim: discard trials with RT+- 3*std(RT) (0 or 1)
% - split_ta_tp: split based on target absent/present? (0 or 1)


% Output: struct containing coherence, PSD, CSD for gradiometers and
% magnetometers for fast and slow trials in respective condition

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)

function c1_rift_RTsplit(s,c_idx,fwdth,filttype,split_ta_tp)



condi_all = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

condi = condi_all{c_idx};                   % select current condition

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid', 'no noise diode');
outpth = fullfile(pth,'results','meg','5 COH hilb','coh','RT');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
alphasoipth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');
soipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

mkdir(fullfile(outpth, subj{s}))

%% load data

d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];


% coherence separately for Target absent/present?
if split_ta_tp
    % if yes, separate trials for ta/tp
    ta_tp = {'ta','tp'};
    
    for ti = 1:length(ta_tp)
        
        % find relevant data files
        condi_files = zeros(length(files),1);
        for c = 1:length(condi)
            
            condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c}),'UniformOutput',false))';
            
        end
        
        condi_files = condi_files == length(condi);
        condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,[condi{2}(end-2:end-1),ta_tp{ti}]),'UniformOutput',false))';
        
        
        files6067 =  condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'6067'),'UniformOutput',false))';
        files6760 =  condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'6760'),'UniformOutput',false))';
        
        files6067 = files6067 == 3;
        files6760 = files6760 == 3;
        
        
        % load files 6067
        files6067 = files(files6067);
        
        load(fullfile(inpth,subj{s},files6067{1}));
        trl_idx = find(trlcur);
        
        % performance in first trials
        trl_perf6067 = perf_cur;
        
        % load data
        data_load = data_trig;
        % data_load = rmfield(data_load,'sampleinfo');
        % load & append data
        for f = 2:length(files6067)
            clear data_trig trlcur perf_cur
            load(fullfile(inpth,subj{s},files6067{f}));
            % append data
            data_load = ft_appenddata([],data_load,data_trig);
            %append performance
            trl_perf6067 = [trl_perf6067;perf_cur];
        end
        
        data6067 = data_load;
        
        clear data_load
        
        
        
        % load files 6760
        files6760 = files(files6760);
        
        load(fullfile(inpth,subj{s},files6760{1}));
        trl_idx = find(trlcur);
        
        % performance in first trials
        trl_perf6760 = perf_cur;
        
        % load data
        data_load = data_trig;
        % data_load = rmfield(data_load,'sampleinfo');
        % load & append data
        for f = 2:length(files6760)
            clear data_trig trlcur perf_cur
            load(fullfile(inpth,subj{s},files6760{f}));
            % append data
            data_load = ft_appenddata([],data_load,data_trig);
            %append performance
            trl_perf6760 = [trl_perf6760;perf_cur];
        end
        
        data6760 = data_load;
        
        clear data_load
        
        
        
        %% Median split fast slow
        
        
        rt_all = cell2mat([trl_perf6067(:,3);trl_perf6760(:,3)]);
        
        % trim data
        m = mean(rt_all);
        std_rt = std(rt_all);
        
%         trl_idx6067 = logical(([trl_perf6067{:,3}]< m-3*std_rt) + ([trl_perf6067{:,3}] > m+3*std_rt));
%         trl_idx6760 = logical(([trl_perf6760{:,3}]< m-3*std_rt) + ([trl_perf6760{:,3}] > m+3*std_rt));
%         
%         cfg = [];
%         cfg.trials = ~trl_idx6067;
%         data6067 = ft_selectdata(cfg,data6067);
%         cfg.trials = ~trl_idx6760;
%         data6760 = ft_selectdata(cfg,data6760);
        
        
%         trl_perf6067(trl_idx6067,:) = [];
%         trl_perf6760(trl_idx6760,:) = [];
%         
%         
%         rt_all = cell2mat([trl_perf6067(:,3);trl_perf6760(:,3)]);
        
        mRT = median(rt_all);
        
        hrate_all = [trl_perf6067(:,2);trl_perf6760(:,2)];
        
        trl_fast6067 = [trl_perf6067{:,3}] < mRT;
        trl_fast6760 = [trl_perf6760{:,3}] < mRT;
        
        trl_slow6067 = [trl_perf6067{:,3}] > mRT;
        trl_slow6760 = [trl_perf6760{:,3}] > mRT;
        
        % proportion ta/tp trials
        load(fullfile(pth, 'experiment','trigdef.mat'))
        
        perc_corr_fast = sum(strcmp(hrate_all(rt_all<mRT),'h'))/sum(rt_all<mRT);
        perc_corr_slow = sum(strcmp(hrate_all(rt_all>mRT),'h'))/sum(rt_all>mRT);
        
        
        % select trials
        cfg = [];
        cfg.trials = find(trl_fast6067);
        data6067_fast = ft_selectdata(cfg,data6067);
        cfg.trials = find(trl_fast6760);
        data6760_fast = ft_selectdata(cfg,data6760);
        cfg.trials = find(trl_slow6067);
        data6067_slow = ft_selectdata(cfg,data6067);
        cfg.trials = find(trl_slow6760);
        data6760_slow = ft_selectdata(cfg,data6760);
        
        clear data
        
        
        
        %% Coherence
        
        % coh in 6067 for T and D
        [coh6067.fast.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_fast,'diode T', {'MEGGRAD'},60, fwdth,filttype);
        [coh6067.fast.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_fast,'diode D', {'MEGGRAD'}, 67, fwdth,filttype);
        [coh6067.slow.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_slow,'diode T', {'MEGGRAD'},60, fwdth,filttype);
        [coh6067.slow.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_slow,'diode D', {'MEGGRAD'}, 67, fwdth,filttype);
        
        
        [coh6760.fast.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_fast,'diode T', {'MEGGRAD'},67, fwdth,filttype);
        [coh6760.fast.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_fast,'diode D',{'MEGGRAD'} , 60, fwdth,filttype);
        [coh6760.slow.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_slow,'diode T', {'MEGGRAD'},67, fwdth,filttype);
        [coh6760.slow.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_slow,'diode D', {'MEGGRAD'}, 60, fwdth,filttype);
        
        
        % baseline correction
        
        % fast
        bsl = mean(coh6067.fast.cohTgrad(:,:,1000:2500),3);
        coh6067.fast.bslcor.cohTgrad = coh6067.fast.cohTgrad - bsl;
        
        bsl = mean(coh6067.fast.cohDgrad(:,:,1000:2500),3);
        coh6067.fast.bslcor.cohDgrad = coh6067.fast.cohDgrad - bsl;
        
        bsl = mean(coh6760.fast.cohTgrad(:,:,1000:2500),3);
        coh6760.fast.bslcor.cohTgrad = coh6760.fast.cohTgrad - bsl;
        
        bsl = mean(coh6760.fast.cohDgrad(:,:,1000:2500),3);
        coh6760.fast.bslcor.cohDgrad = coh6760.fast.cohDgrad - bsl;
        
        % slow
        bsl = mean(coh6067.slow.cohTgrad(:,:,1000:2500),3);
        coh6067.slow.bslcor.cohTgrad = coh6067.slow.cohTgrad - bsl;
        
        bsl = mean(coh6067.slow.cohDgrad(:,:,1000:2500),3);
        coh6067.slow.bslcor.cohDgrad = coh6067.slow.cohDgrad - bsl;
        
        bsl = mean(coh6760.slow.cohTgrad(:,:,1000:2500),3);
        coh6760.slow.bslcor.cohTgrad = coh6760.slow.cohTgrad - bsl;
        
        bsl = mean(coh6760.slow.cohDgrad(:,:,1000:2500),3);
        coh6760.slow.bslcor.cohDgrad = coh6760.slow.cohDgrad - bsl;
        
        
      
        condname = strjoin(condi,'_');
        save(fullfile(outpth,subj{s},[condname(1:end-1),ta_tp{ti},'_RTsplit.mat']),'coh6067','coh6760','perc_corr_slow','perc_corr_fast')
        
    end
    
% if coherence calculated over both Target absent, target present, still
% split into both and recombine
else
    
    % if yes, separate trials for ta/tp
    ta_tp = {'ta','tp'};
    
    for ti = 1:length(ta_tp)
        
        % find relevant data files
        condi_files = zeros(length(files),1);
        for c = 1:length(condi)
            
            condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c}),'UniformOutput',false))';
            
        end
        
        condi_files = condi_files == length(condi);
        condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,[condi{2}(end-2:end-1),ta_tp{ti}]),'UniformOutput',false))';
        
        
        files6067 =  condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'6067'),'UniformOutput',false))';
        files6760 =  condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,'6760'),'UniformOutput',false))';
        
        files6067 = files6067 == 3;
        files6760 = files6760 == 3;
        
        
        % load files 6067
        files6067 = files(files6067);
        
        load(fullfile(inpth,subj{s},files6067{1}));
        trl_idx = find(trlcur);
        
        % performance in first trials
        trl_perf6067 = perf_cur;
        
        % load data
        data_load = data_trig;
        % data_load = rmfield(data_load,'sampleinfo');
        % load & append data
        for f = 2:length(files6067)
            clear data_trig trlcur perf_cur
            load(fullfile(inpth,subj{s},files6067{f}));
            % append data
            data_load = ft_appenddata([],data_load,data_trig);
            %append performance
            trl_perf6067 = [trl_perf6067;perf_cur];
        end
        
        data6067 = data_load;
        
        clear data_load
        
        
        
        % load files 6760
        files6760 = files(files6760);
        
        load(fullfile(inpth,subj{s},files6760{1}));
        trl_idx = find(trlcur);
        
        % performance in first trials
        trl_perf6760 = perf_cur;
        
        % load data
        data_load = data_trig;
        % data_load = rmfield(data_load,'sampleinfo');
        % load & append data
        for f = 2:length(files6760)
            clear data_trig trlcur perf_cur
            load(fullfile(inpth,subj{s},files6760{f}));
            % append data
            data_load = ft_appenddata([],data_load,data_trig);
            %append performance
            trl_perf6760 = [trl_perf6760;perf_cur];
        end
        
        data6760 = data_load;
        
        clear data_load
        
        
        
        %% Median split fast slow
        
        
        rt_all = cell2mat([trl_perf6067(:,3);trl_perf6760(:,3)]);
        
        % trim data
        m = mean(rt_all);
        std_rt = std(rt_all);
        
%         trl_idx6067 = logical(([trl_perf6067{:,3}]< m-3*std_rt) + ([trl_perf6067{:,3}] > m+3*std_rt));
%         trl_idx6760 = logical(([trl_perf6760{:,3}]< m-3*std_rt) + ([trl_perf6760{:,3}] > m+3*std_rt));
%         
%         cfg = [];
%         cfg.trials = ~trl_idx6067;
%         data6067 = ft_selectdata(cfg,data6067);
%         cfg.trials = ~trl_idx6760;
%         data6760 = ft_selectdata(cfg,data6760);
%         
%         
%         trl_perf6067(trl_idx6067,:) = [];
%         trl_perf6760(trl_idx6760,:) = [];
%         
%         
%         rt_all = cell2mat([trl_perf6067(:,3);trl_perf6760(:,3)]);
        
        mRT = median(rt_all);
        
        hrate_all = [trl_perf6067(:,2);trl_perf6760(:,2)];
        
        trl_fast6067 = [trl_perf6067{:,3}] < mRT;
        trl_fast6760 = [trl_perf6760{:,3}] < mRT;
        
        trl_slow6067 = [trl_perf6067{:,3}] > mRT;
        trl_slow6760 = [trl_perf6760{:,3}] > mRT;
        
        % proportion ta/tp trials
        load(fullfile(pth, 'experiment','trigdef.mat'))
        
        perc_corr_fast = sum(strcmp(hrate_all(rt_all<mRT),'h'))/sum(rt_all<mRT);
        perc_corr_slow = sum(strcmp(hrate_all(rt_all>mRT),'h'))/sum(rt_all>mRT);
        
        
        % select trials
        cfg = [];
        cfg.trials = find(trl_fast6067);
        data6067_fast{ti} = ft_selectdata(cfg,data6067);
        cfg.trials = find(trl_fast6760);
        data6760_fast{ti} = ft_selectdata(cfg,data6760);
        cfg.trials = find(trl_slow6067);
        data6067_slow{ti} = ft_selectdata(cfg,data6067);
        cfg.trials = find(trl_slow6760);
        data6760_slow{ti} = ft_selectdata(cfg,data6760);
        
    end
    
    data6067_fast = ft_appenddata(cfg,data6067_fast{:});
    data6067_slow = ft_appenddata(cfg,data6067_slow{:});
    
    data6760_fast = ft_appenddata(cfg,data6760_fast{:});
    data6760_slow = ft_appenddata(cfg,data6760_slow{:});
    
    %% Coherence
    
    % coh in 6067 for T and D
    [coh6067.fast.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_fast,'diode T', {'MEGGRAD'},60, fwdth,filttype);
    [coh6067.fast.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_fast,'diode D', {'MEGGRAD'}, 67, fwdth,filttype);
    [coh6067.slow.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_slow,'diode T', {'MEGGRAD'},60, fwdth,filttype);
    [coh6067.slow.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6067_slow,'diode D', {'MEGGRAD'}, 67, fwdth,filttype);
    
    
    [coh6760.fast.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_fast,'diode T', {'MEGGRAD'},67, fwdth,filttype);
    [coh6760.fast.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_fast,'diode D',{'MEGGRAD'} , 60, fwdth,filttype);
    [coh6760.slow.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_slow,'diode T', {'MEGGRAD'},67, fwdth,filttype);
    [coh6760.slow.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data6760_slow,'diode D', {'MEGGRAD'}, 60, fwdth,filttype);
    
    
    % baseline correction
    
    % fast
    bsl = mean(coh6067.fast.cohTgrad(:,:,1500:2500),3);
    coh6067.fast.bslcor.cohTgrad = coh6067.fast.cohTgrad - bsl;
    
    bsl = mean(coh6067.fast.cohDgrad(:,:,1500:2500),3);
    coh6067.fast.bslcor.cohDgrad = coh6067.fast.cohDgrad - bsl;
    
    bsl = mean(coh6760.fast.cohTgrad(:,:,1500:2500),3);
    coh6760.fast.bslcor.cohTgrad = coh6760.fast.cohTgrad - bsl;
    
    bsl = mean(coh6760.fast.cohDgrad(:,:,1500:2500),3);
    coh6760.fast.bslcor.cohDgrad = coh6760.fast.cohDgrad - bsl;
    
    % slow
    bsl = mean(coh6067.slow.cohTgrad(:,:,1500:2500),3);
    coh6067.slow.bslcor.cohTgrad = coh6067.slow.cohTgrad - bsl;
    
    bsl = mean(coh6067.slow.cohDgrad(:,:,1500:2500),3);
    coh6067.slow.bslcor.cohDgrad = coh6067.slow.cohDgrad - bsl;
    
    bsl = mean(coh6760.slow.cohTgrad(:,:,1500:2500),3);
    coh6760.slow.bslcor.cohTgrad = coh6760.slow.cohTgrad - bsl;
    
    bsl = mean(coh6760.slow.cohDgrad(:,:,1500:2500),3);
    coh6760.slow.bslcor.cohDgrad = coh6760.slow.cohDgrad - bsl;
    
    condname = strjoin(condi,'_');
    save(fullfile(outpth,subj{s},[condname,'_RTsplit.mat']),'coh6067','coh6760','perc_corr_slow')
end

end
