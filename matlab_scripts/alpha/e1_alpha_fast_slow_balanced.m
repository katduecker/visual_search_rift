%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e1. calculate TFRs for fast vs slow trials (Fig. 4 d)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 May 2026

% Inputs
% - s: subject index
% - c_idx: condition index
% - split_ta_tp (bool): split into target absent/present?

% Output
% - soi_grad: gradiometers with high alpha power
% iaf_grad: identified IAF in gradiometers


function e1_alpha_fast_slow_balanced(s)

num_bins = 4;
winl = 0.5;
condi_all = {{'ni','16ta'},{'ti','16ta'}, {'ni','32ta'},{'ti','32ta'},{'ni','16tp'},{'ti','16tp'}, {'ni','32tp'},{'ti','32tp'}};


%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
outpth = fullfile(pth,'results','meg','6 Alpha','rt_balanced');

glmpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');
addpath('/rds/projects/j/jenseno-visual-search-rft/visual_search_rift/fieldtrip')

load(fullfile(pth,'matlab_scripts/',"preproc_meg",'idx_subjoi.mat'));

load(fullfile(glmpth, 'glm_RT_soi_iaf_subj'))

soi = subj_soi{s};
time_oi = glm_time_sig;

%% load data

d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];


freqvec = 4:1/winl:30;
spec_fast = zeros(length(condi_all), length(freqvec));
spec_slow = spec_fast;

for c_idx = 1:length(condi_all)
    condi = condi_all{c_idx};
    % find relevant data files
    condi_files = zeros(length(files),1);
    for c = 1:length(condi)
        condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c}),'UniformOutput',false))';
    end
    condi_files = condi_files == length(condi);   

    c_files = files(condi_files);

    load(fullfile(inpth,subj{s},c_files{1}));
    trl_idx = find(trlcur);

    % performance in first trials
    trl_perf = perf_cur;

    % load data
    data_load = data_trig;
    % data_load = rmfield(data_load,'sampleinfo');
    % load & append data
    for f = 2:length(c_files)
        clear data_trig trlcur perf_cur
        load(fullfile(inpth,subj{s},c_files{f}));
        % append data
        data_load = ft_appenddata([],data_load,data_trig);
        %append performance
        trl_perf = [trl_perf;perf_cur];
    end

    data = data_load;

    clear data_load

    %% Balanced split
    rt_all = [trl_perf{:,3}];

    lblock = floor(length(rt_all)/num_bins);

    % sort according to t-o-t
    [~,I] = sort([trl_perf{:,end}]);

    % sort RT
    rt_all = rt_all(I);

    h = 1;

    rt_block = zeros(num_bins,lblock);
    idx = zeros(num_bins,lblock);
    for b = 1:num_bins

        rt_block(b,:) = rt_all(h:h+lblock-1);

        idx(b,:) = I(h:h+lblock-1);

        h = h+lblock;

    end

    idx_fast = [];
    idx_slow = [];
    for b = 1:num_bins
        % get the top and bottom 10%
        [~, p] = mink(rt_block(b,:), round(size(rt_block,2)*0.1));
        idx_fast = [idx_fast,p];
        [~, p] = maxk(rt_block(b,:), round(size(rt_block,2)*0.1));
        idx_slow = [idx_slow, p];
    end


    % select trials
    cfg = [];
    cfg.trials = idx_fast(:);
    data_fast = ft_selectdata(cfg,data);

    cfg.trials = idx_slow(:);
    data_slow = ft_selectdata(cfg,data);

    clear data


    %% Power 

    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.channels = {'MEGGRAD'};
    cfg.taper = 'hanning';
    cfg.foi = freqvec;
    cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
    cfg.toi = -1.75:0.05:1;
    cfg.keeptrials = 'no';
    TFRfast = ft_freqanalysis(cfg,data_fast);
    TFRslow = ft_freqanalysis(cfg,data_slow);
    
    TFRfast = ft_combineplanar([], TFRfast);
    TFRslow = ft_combineplanar([], TFRslow);

    cfg = [];
    cfg.avgovertime = 'yes';
    cfg.latency = time_oi;
    cfg.channel = soi;
    cfg.avgoverchan = 'yes';
    IAFfast = ft_selectdata(cfg,TFRfast);
    IAFslow = ft_selectdata(cfg,TFRslow);

    spec_slow(c_idx,:) = IAFslow.powspctrm./max(IAFslow.powspctrm);
    spec_fast(c_idx,:) = IAFfast.powspctrm./max(IAFslow.powspctrm);
    

clear IAFfast IAFslow idx_fast idx_slow data_fast data_slow I rt_block idx trl_perf
end
save(fullfile(outpth,subj{s},'spectra_RTsplit.mat'),'spec_fast', 'spec_slow', 'freqvec')


end
