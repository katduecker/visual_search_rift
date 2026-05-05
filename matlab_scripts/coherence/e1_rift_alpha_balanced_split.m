%% Investigating Guided Search using Rapid Invisible Frequency Tagging
% Paper 2: Alpha inhibition

%% Balanced median split RIFT ~ alpha

% (c), Katharina Duecker
% last edited, Oct-05-2025

function e1_rift_alpha_balanced_split(s, time_init)

% inputs
% s: subject index
% n_blocks: number of blocks experiment should be split into to account for
% t-o-t

fwdth = 3.5;
filttype={'firws', 'twopass'};


num_bins = 4;
winl = 0.5;
% don't split into target present/absent because RT is not part of the
% split!
condi_all = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};


rmpath(genpath('/rds/projects/2018/jenseno-entrainment/fieldtrip'))

pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';
addpath(fullfile(pth, 'fieldtrip'))
ft_defaults;
addpath(fullfile(pth,'matlab scripts','alpha'))


inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid');
soipth = fullfile(pth,'results','meg','6 Alpha','iaf_soi');

cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','balanced_split');
cohsoipth = fullfile(pth,'results','meg','5 COH hilb', 'soi','sinusoid');

glmpth = fullfile(pth,'results','meg','9 GLM', 'glm_spec');

% load participant id's
load(fullfile(pth,'matlab_scripts',"preproc_meg",'idx_subjoi.mat'));

% make output directory
mkdir(fullfile(cohpth,subj{s}));

% load SOI, IAFs, and signfiicant time 
load(fullfile(glmpth, 'glm_RT_soi_iaf_subj'), 'subj_soi', 'iaf_glm', 'glm_time_sig')

chan_subj = subj_soi{s};
iaf = iaf_glm(s);
% if pre-defined time_oi, take that; otherwise use GLM time points
time_oi = zeros(1,2);
time_oi(isnan(time_init)) = glm_time_sig(isnan(time_init));
time_oi(~isnan(time_init)) = time_init(~isnan(time_init));

load(fullfile(cohsoipth, subj{s}, 'soi_stat.mat'), 'soigrad')

%% load data 
%- this is a bit complicated but the split up data has information about which diode contains the T and D signal!

d = dir(fullfile(inpth,subj{s}));
files = {d.name};
files(1:2) = [];

load(fullfile(inpth,subj{s},files{1}))

data = data_trig;
behav_array = perf_cur;

for f = 2:length(files)
    clear data_trig trlcur perf_cur
    load(fullfile(inpth,subj{s},files{f}));
    % append data
    data = ft_appenddata([],data,data_trig);
    %append performance
    behav_array = [behav_array;perf_cur];
end

% sort according to t-o-t
[tot,I] = sort([behav_array{:,end}]);

behav_array = behav_array(I,:);

data.trial = data.trial(I);
data.time = data.time(I);

if isfield(data,'sampleinfo')
 data = rmfield(data,'sampleinfo');
end

%% find trigger current condition
trig_load = load(fullfile(pth,'experiment','trigdef.mat'));
trigdef = trig_load.trigdef;
    
%% TFR

cfg = [];
cfg.method = 'mtmconvol';
soi = {};

% splot sensors OI into separate gradiometers
for c = 1:length(chan_subj)
    chan = {chan_subj{c}(1:7)};
    soi = [soi;chan];
    chan = {['MEG',chan_subj{c}(9:end)]};
    soi = [soi;chan];
end

cfg.channel = soi;
cfg.foi = iaf;

cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:0.5;
cfg.keeptrials = 'yes';
cfg.taper = 'hanning';
TFR = ft_freqanalysis(cfg,data);

% combine planar
cfg = [];
cfg.method = 'sum';
TFR = ft_combineplanar(cfg,TFR);


% IAF
cfg = [];
cfg.latency = time_oi;
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
IAF = ft_selectdata(cfg,TFR);

clear TFR
%% Perform balanced split

lblock = floor(length(behav_array)/num_bins);

% find Target and Distractor frequencies
% trig 6067
trig6067 = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),'6067'),'UniformOutput',false)),1));
% trig6760
trig6760 = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x),regexp(trigdef(:,2),'6760'),'UniformOutput',false)),1));

cohT_high = zeros(length(condi_all),6501);
cohD_high = zeros(length(condi_all),6501);
cohT_low = zeros(length(condi_all),6501);
cohD_low = zeros(length(condi_all),6501);

coh6067_high = cell(1,length(condi_all));
coh6760_high = cell(1,length(condi_all));
coh6067_low = cell(1,length(condi_all));
coh6760_low = cell(1,length(condi_all));

hits_high = [];
hits_low = [];

trigcomp = cell(1,length(condi_all));
for c_idx = 1:length(condi_all)
    trigcomp{c_idx} = trigdef;
end

% number of trials per block
n_trials = zeros(length(condi_all),num_bins);
% loop over conditions
for c_idx = 1:length(condi_all)
    
    data_high6067 = cell(1,num_bins);
    data_low6067 = cell(1,num_bins);
    
    data_high6760 = cell(1,num_bins);
    data_low6760 = cell(1,num_bins);
    
    condi = condi_all{c_idx};

    
    % find trigger for current condition
    trig_idx = cell2mat(cellfun(@(x) ~isempty(x),regexp(trigcomp{c_idx}(:,2),strjoin(condi,'')),'UniformOutput',false));
    
    trigger = cell2mat(trigcomp{c_idx}(trig_idx,1));
    
    

    h = 1;
    for block_id = 1:num_bins
        
        % trials per block
        block_trials = behav_array(h:h+lblock-1,:);
        
        % which condition do trials belong to?
        condi_trials = cell2mat(block_trials(ismember([block_trials{:,1}],trigger),1));
        n_trials(c_idx,block_id) = length(condi_trials);
        
        block_hits = strcmp(block_trials(:,2),'h');
        if length(condi_trials)>1
            
            % alpha
            block_alpha = IAF.powspctrm(h:h+lblock-1);
            condi_alpha = block_alpha(ismember([block_trials{:,1}],trigger));
            
            % data
            cfg = [];
            cfg.trials = h:h+lblock-1;
            data_block = ft_selectdata(cfg,data);
            
            cfg.trials = ismember([block_trials{:,1}],trigger);
            data_condi = ft_selectdata(cfg,data_block);
            
            hits_high = [hits_high;block_hits(condi_alpha > median(condi_alpha))];
            hits_low = [hits_low;block_hits(condi_alpha < median(condi_alpha))];
            
            % split based on condition and whether it's Target 60 and
            % Distractor 67 or the other way around
            cfg.trials = (condi_alpha > median(condi_alpha)) + (ismember(condi_trials,trig6067)) == 2;
            
            if ~isempty(find(cfg.trials))
                data_high6067{block_id} = ft_selectdata(cfg,data_condi);
            end
            
            cfg.trials = (condi_alpha < median(condi_alpha))+ (ismember(condi_trials,trig6067)) == 2;
            
            if ~isempty(find(cfg.trials))
                
                data_low6067{block_id} = ft_selectdata(cfg,data_condi);
            end
            
            
            cfg.trials = (condi_alpha > median(condi_alpha)) + (ismember(condi_trials,trig6760)) == 2;
            
            if ~isempty(find(cfg.trials))
                
                data_high6760{block_id} = ft_selectdata(cfg,data_condi);
            end
            
            cfg.trials = (condi_alpha < median(condi_alpha))+ (ismember(condi_trials,trig6760)) == 2;
            
            if ~isempty(find(cfg.trials))
                
                data_low6760{block_id} = ft_selectdata(cfg,data_condi);
            end
        end
        
        h = h+lblock;

        
    end
    
    data_high6067 = ft_appenddata([],data_high6067{cell2mat(cellfun(@(x) ~isempty(x),data_high6067,'UniformOutput',false))});
    data_high6760 = ft_appenddata([],data_high6760{cell2mat(cellfun(@(x) ~isempty(x),data_high6760,'UniformOutput',false))});
    
    
    data_low6067 = ft_appenddata([],data_low6067{cell2mat(cellfun(@(x) ~isempty(x),data_low6067,'UniformOutput',false))});
    data_low6760 = ft_appenddata([],data_low6760{cell2mat(cellfun(@(x) ~isempty(x),data_low6760,'UniformOutput',false))});
    
    
 
    [coh6067_high{c_idx}.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data_high6067,'diode T', soigrad,60, fwdth,filttype);
    [coh6067_high{c_idx}.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data_high6067,'diode D', soigrad,67, fwdth,filttype);
    [coh6760_high{c_idx}.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data_high6760,'diode T', soigrad,67, fwdth,filttype);
    [coh6760_high{c_idx}.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data_high6760,'diode D', soigrad,60, fwdth,filttype);


    [coh6067_low{c_idx}.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data_low6067,'diode T', soigrad,60, fwdth,filttype);
    [coh6067_low{c_idx}.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data_low6067,'diode D', soigrad,67, fwdth,filttype);
    [coh6760_low{c_idx}.cohTgrad, ~, ~, ~] = kd_coh_hilb_fun(data_low6760,'diode T', soigrad,67, fwdth,filttype);
    [coh6760_low{c_idx}.cohDgrad, ~, ~, ~] = kd_coh_hilb_fun(data_low6760,'diode D', soigrad,60, fwdth,filttype);
    
    
    % baseline correct
    bsl = mean(coh6067_high{c_idx}.cohTgrad(:,:,1000:2500),3);
    coh6067_high{c_idx}.bslcor.cohTgrad = coh6067_high{c_idx}.cohTgrad - bsl;
    bsl = mean(coh6067_high{c_idx}.cohDgrad(:,:,1000:2500),3);
    coh6067_high{c_idx}.bslcor.cohDgrad = coh6067_high{c_idx}.cohDgrad - bsl;
    
    bsl = mean(coh6760_high{c_idx}.cohTgrad(:,:,1000:2500),3);
    coh6760_high{c_idx}.bslcor.cohTgrad = coh6760_high{c_idx}.cohTgrad - bsl;
    bsl = mean(coh6760_high{c_idx}.cohDgrad(:,:,1000:2500),3);
    coh6760_high{c_idx}.bslcor.cohDgrad = coh6760_high{c_idx}.cohDgrad - bsl;
    
    bsl = mean(coh6067_low{c_idx}.cohTgrad(:,:,1000:2500),3);
    coh6067_low{c_idx}.bslcor.cohTgrad = coh6067_low{c_idx}.cohTgrad - bsl;
    bsl = mean(coh6067_low{c_idx}.cohDgrad(:,:,1000:2500),3);
    coh6067_low{c_idx}.bslcor.cohDgrad = coh6067_low{c_idx}.cohDgrad - bsl;
    
    bsl = mean(coh6760_low{c_idx}.cohTgrad(:,:,1000:2500),3);
    coh6760_low{c_idx}.bslcor.cohTgrad = coh6760_low{c_idx}.cohTgrad - bsl;
    bsl = mean(coh6760_low{c_idx}.cohDgrad(:,:,1000:2500),3);
    coh6760_low{c_idx}.bslcor.cohDgrad = coh6760_low{c_idx}.cohDgrad - bsl;
    
    cohT_high(c_idx,:) = (squeeze(mean(coh6067_high{c_idx}.cohTgrad(1:end-2,1,:),1)) + squeeze(mean(coh6760_high{c_idx}.cohTgrad(1:end-2,1,:),1)))./2;
    cohD_high(c_idx,:) = (squeeze(mean(coh6067_high{c_idx}.cohDgrad(1:end-2,1,:),1)) + squeeze(mean(coh6760_high{c_idx}.cohDgrad(1:end-2,1,:),1)))./2;
    
    cohT_low(c_idx,:) = (squeeze(mean(coh6067_low{c_idx}.cohTgrad(1:end-2,1,:),1)) + squeeze(mean(coh6760_low{c_idx}.cohTgrad(1:end-2,1,:),1)))./2;
    cohD_low(c_idx,:) = (squeeze(mean(coh6067_low{c_idx}.cohDgrad(1:end-2,1,:),1)) + squeeze(mean(coh6760_low{c_idx}.cohDgrad(1:end-2,1,:),1)))./2;

end

toi_split = time_oi;
toi_split = strjoin(arrayfun(@(x) num2str(x),toi_split.*1000,'UniformOutput',false),'_to_');
save(fullfile(cohpth,subj{s},['balanced_split_',toi_split]),'cohT_high','cohD_high','cohT_low','cohD_low','n_trials', 'hits_high', 'hits_low')




