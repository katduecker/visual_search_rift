%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% g. split data according to conditions

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Input
% - s: subject id
% - noise: amplitude of noise added to the diode (diode will be replaced w/
% perfect sine wave + noise (noise is important for coherence later)

% Output
% - mat files containing trials + performance for each condition

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions


function g_split_trl(s,noise)

%% settings
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path
ldiode = {'MISC004','MISC005'}; % diode label
rft_freq = [60,67];             % tagging frequencies
phshft = [0,0.65];              % phase shift per tagging frequency
fs = 1000;                      % sampling rate
toi = [-2.5,4];                 % extract data from when to when? (relative to search display onset)

dtpth = fullfile(pth, 'data'); % data path
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))
ft_defaults;

dtpth = fullfile(pth,'data');                                       % raw data
maxfpth = fullfile(pth,'results','meg', '1 maxfilter','1 maxfilter');             % max filter
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');   % scripts
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
icapth = fullfile(pth,'results','meg', '3 ICA', '1 all subj');

if noise
    respth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid');
else
    respth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid','no noise diode');
end
mkdir(respth)
%% load in data and discard artefactual trials
load(fullfile(scriptpth,'idx_subjoi.mat'))                          % usable subjects index

% if exist(fullfile(respth,subj{s}))
%     error('subject folder already exists')
% end

% load in trial structure
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))

% list maxfiltered data
d = dir(fullfile(maxfpth,subj{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);

% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    hdr = ft_read_header(fullfile(maxfpth,subj{s},f{p}));
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)+fs*toi(1),meginfo.alltrl_bl{p}(:,3)+toi(2)*fs,zeros(length(meginfo.alltrl_bl{p}),1)+toi(1)*fs];
    
    % if end of file is reached, don't load in full 4 s
    x = find(trlstruct{p}(:,2) > hdr.nSamples);
    if x
        trlstruct{p}(x,2) = hdr.nSamples;
    end
    clear hdr
end


% load data
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(maxfpth,subj{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    cfg.channel = 'MEG';
    % load in data for this part
    dtprt{p} = ft_preprocessing(cfg);
    grad{p} = dtprt{p}.grad;
    % diodes
    cfg.channel = {'MISC004', 'MISC005'};
    diodes_trl{p}= ft_preprocessing(cfg);
end

data = ft_appenddata([],dtprt{:});
data.hdr = dtprt{1}.hdr;
diodes = ft_appenddata([], diodes_trl{:});


%% head position relative to z axis
% parietal sensors as reference
pariet_sens = {'MEG0712','MEG0711','MEG0713','MEG0722','MEG0721','MEG0723','MEG0742','MEG0741','MEG0743','MEG0732','MEG0731','MEG0733','MEG1832','MEG1831','MEG1833','MEG2242','MEG2241','MEG2243','MEG2022','MEG2021','MEG2023'};

pari_idx = ismember(grad{1}.label,pariet_sens);

dist_z = cell(1,length(f));

for k = 1:length(dist_z)
    dist_z{k} = repmat(mean(grad{1}.chanpos(pari_idx,3)-grad{k}.chanpos(pari_idx,3)),length(dtprt{k}.trial),1);
end

dist_z_array = vertcat(dist_z{:});

clear dtprt diodes_trl
%%

% discard strange trials
cfg = [];
cfg.trials =find(meginfo.keeptrl_all);                   % select clean trials
data = ft_selectdata(cfg,data);
diodes = ft_selectdata(cfg,diodes);

% sampleinfo
trl_end_samp = cell2mat(cellfun(@length,data.trial,'UniformOutput',false)).*[1:length(data.trial)];
data.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];
diodes.sampleinfo = [[1;1+trl_end_samp(1:end-1)'],trl_end_samp'];

%% load in ICA and reject comps

% load ICA data and components
load(fullfile(icapth, [subj{s},'_ica.mat']))

cfg = [];
cfg.component = badcomps; % to be removed component(s)
dataclean = ft_rejectcomponent(cfg, dataICA, data);
clear dataICA data

% mean channel position
mGrad = grad{1};
% average grad structure
for g = 2:length(grad)
mGrad.chanpos = mGrad.chanpos + grad{g}.chanpos;
end
mGrad.chanpos = mGrad.chanpos./length(grad);
clear grad

% concatenate cleaned data and diodes
data = ft_appenddata([],dataclean,diodes);
data.grad = mGrad;

%% reject trials containing rt's outside reasonable range

% RT
rt_array = [rspinfo.trl{:,3}];
rt_array = rt_array(rspinfo.keeptrl_rsp);

% behavior (trigger, hit/miss/fa, RT);
behav_array = rspinfo.trl;%(rspinfo.keeptrl_rsp,:);

% add time on task
behav_array = [behav_array,num2cell([1:length(rspinfo.keeptrl_rsp)]'./length(rspinfo.keeptrl_rsp))];

behav_array = behav_array(rspinfo.keeptrl_rsp,:);

 % reject trials with RT < 200 ms or RT > 4 s
rej_trl = logical(double(rt_array< 0.2) + double(rt_array == 4));   % this didn't work in first run, check
meginfo.keeptrl_all(find(rej_trl)) = 0;

% reject all trials based on meginfo.keeptrl_all
alltrl_list = meginfo.alltrl_list(logical(meginfo.keeptrl_all),:);
behav_array = behav_array(logical(meginfo.keeptrl_all),:);

% add head position z-axis
dist_z_array = dist_z_array(logical(meginfo.keeptrl_all),:);

%% save all trials (no split into conditions to keep trial order
mkdir(fullfile(respth,'clean data all trials',subj{s}))
save(fullfile(respth,'clean data all trials',subj{s},'data.mat'),'data','behav_array','dist_z_array','-v7.3')

save(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'),'meginfo','-append')
% 
%% Split into trials based on condition

load(fullfile(pth,'experiment','trigdef.mat'))

for trg_idx = 5:length(trigdef)

    % find trials
    trlcur = ismember(alltrl_list(:,1),trigdef{trg_idx});

    % select data
    cfg = [];
    cfg.trials = trlcur;
    data_cur = ft_selectdata(cfg, data);

    % select performance
    perf_cur = behav_array(trlcur,:);
    
    % head position
    head_cur = dist_z_array(trlcur,:);

    % re-name diode based on target and distractor
    % select T: Target color: end-5: channel index (1: MISC004, 2: MISC005)
    cfg = [];
    chanT = ldiode{str2num(trigdef{trg_idx,2}(end-5))};
    cfg.channel = chanT;
    data_T = ft_selectdata(cfg,data_cur);
    data_T.label = {'diode T'};
    
    % replace with perfect sinusoid
    freqT = str2num(trigdef{trg_idx,2}(end-3:end-2));           % target frequency
    data_T = kd_replace_diode_sinu(data_T,1,freqT,phshft(rft_freq == freqT),-2.5,fs,noise);
    
    % select distractor diode (the channel that was not selected above)
    cfg = [];
    cfg.channel = ldiode{~strcmp(ldiode,chanT)};
    data_D = ft_selectdata(cfg,data_cur);
    data_D.label = {'diode D'};
    
    % replace with perfect sinusoid
    freqD = str2num(trigdef{trg_idx,2}(end-1:end));
    data_D = kd_replace_diode_sinu(data_D,1,freqD,phshft(rft_freq == freqD),-2.5,fs,noise);
    
    % only select MEG
    cfg= [];
    cfg.channel = 'MEG';
    data_cur = ft_selectdata(cfg,data_cur);


    % append with target and distractor diode data
    data_trig = ft_appenddata([],data_cur,data_T,data_D);
    data_trig.fsample = 1000;
    
    mkdir(fullfile(respth,subj{s}))
    
    % save: data current condition, trial index, performance
    save(fullfile(respth,subj{s},['data_',trigdef{trg_idx,2},'.mat']),'data_trig','trlcur','perf_cur','head_cur')

    clear data_D data_T data_cur data_trig
end
