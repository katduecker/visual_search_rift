%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a1. Define trial structure

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% a. Define trial structure

% outputs: trl_overlap_meg_el_rsp.mat

% - meginfo: structure containing 
%   - alltrl_bl/alltrl_list: info on MEG trials, indices in blocks and lists
%   - keeptrl_all:           index of trials that should be rejected from edf, matfile and MEG
%   - meg_rt:                reaction time read out by MEG triggers
%   - trigtrl:               triggers read out for each trial (order should be trial type, 2, number >5000 = button press)

% - elinfo: structure containing
%   - eltrl (trial information): 
%            - trial type (trigger); 
%            - time points start - onset display - end of trial (trigger 4)
%            - samples corrsponding to time points above
%   - keeptrl_rsp: logical array, trials to be kept based on response file -> use before keeptrl_all!
%   - keeptrl_all: logical array, which trials to be kept of all (after rsp have been selected)

% - rspinfo: structure containing
%   - trl: cell containing trial type (trigger), 'h/m/fa', RT as read out by stim computer
%   - keeptrl_rsp (same as elinfo): logical array, trials to be kept based on response file -> use before meginfo.keeptrl_all!
%   - keeptrl_all (same as elinfo): logical array, which trials to be kept of all (after rsp have been selected)

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement -> just sanity check, kept trials w/ saccades
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RIFT response
% g. Split trials into conditions

clear all; close all; clc; beep off

%% PATHS
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                 % main path
mtpth = fullfile(pth,'matlab scripts');                       % matlab scripts path


edfconvpath = fullfile(pth,'edf-converter-master/');          % path to edf converter
addpath(edfconvpath)
command = 'edf2asc.exe';                                      % eyelink function edf2asc                                         

dtpth = fullfile(pth,'data');                                 % data path
maxfpth = fullfile(pth, 'results','meg','1 maxfilter');
trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');  % results path
mkdir(trl_merge_pth)

addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))                              % fieldtrip  
ft_defaults;

% list subjects
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds


% document: subject id, n trial meg, el adjusted, response adjusted
if exist(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))             % subject documentation exists already
    load(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))             % load
    if size(mergesubj,1) < length(subjfolds)
        mergesubj = [mergesubj,cell(length(subjfolds)-size(mergesubj,1),4)];
    end
else
    mergesubj = cell(length(subjfolds),4);                                               % else create empty cell
end

propixx_res = [1920 1080];                                   % propixx resolution

% load(fullfile(pth,'experiment','button_val.mat'),'buttons_right');
% buttons_right = buttons_right(1:3);


for s = 1:length(subjfolds)
    % current subject
    disp(['subj ', subjfolds{s}])
    mergesubj{s,1} = subjfolds{s};
    trl_merge_pth = fullfile(pth, 'results','meg','2 merged edf mat');
    
    %% Prep el file

    % convert edf file to mat struct using edf2mat
    % list files ending in edf.
    d = dir(fullfile(dtpth,subjfolds{s}));
    f = {d.name};
    idx = cellfun(@(x) regexp(x,'.edf'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    edf_files = f(idxx);
    file_out = 'el_struct.mat';

    % if el structure exists,load in
    if exist(fullfile(dtpth,subjfolds{s},file_out))

        load(fullfile(dtpth,subjfolds{s},file_out))
        el_fsample = 1000;

        % if el structure does not exist, load in edf file using edf2mat
        % and then convert matlab object to struct using kd_edf2matstruct
    else
        % if there is one edf file (experiment not interrupted)
        if length(edf_files) == 1
            el = kd_edf2mat2struct(edfconvpath, fullfile(dtpth,subjfolds{s}), edf_files{1},file_out);
            el_fsample = 1000;

            % if there is two
        elseif length(edf_files) == 2
            el = kd_edf2mat2struct_2files(edfconvpath, fullfile(dtpth,subjfolds{s}), edf_files,file_out);
            el_fsample = 1000;

            % error if there is none or >2
        else
            mergesubj{s,3} = 'unusual number of el files';
        end
    end
    

    %% Prep MEG 
    % list maxfiltered fif files
    mgfls = dir(fullfile(maxfpth,subjfolds{s}));
    f = {mgfls.name};
    
    % find fif file
    idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    f = f(idxx);
    
    % load events separated in parts
    for fl = 1:length(f)
        events{fl} = ft_read_event(fullfile(maxfpth,subjfolds{s},f{fl}));
    end
    
    load(fullfile(pth,'experiment','trigdef.mat'))
    
    %  List trial type fif
    alltrl_bl = {};
    trigtrl = {};

    % trl: type, start sample, end sample
    for fl = 1:length(events)
        allval = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).value];   % these are all values, extract button press
        allsmp = [events{fl}(find(strcmp('STI101',{events{fl}.type}))).sample];

        find(allval==1)
        % get trigger indices
        trig_bsl = find((allval >= 5)+(allval<=trigdef{end,1}) == 2);

        % if the following trigger is not 2, reject
        trl_idx = allval(trig_bsl+1) == 2;
        trig_bsl = trig_bsl(trl_idx);

        % store samples
        trl = [allval(trig_bsl)',allsmp(trig_bsl)' allsmp(trig_bsl+1)' allsmp(trig_bsl+2)'];
        % store trigger values
        trigtrl{fl} = [allval(trig_bsl)',allval(trig_bsl+1)',allval(trig_bsl+2)'];
        alltrl_bl{fl} = trl;
    end
    
    % if there is not enough identified trials, chances are the triggers
    % were messed up by P pressing button for too long
    alltrl_list = vertcat(alltrl_bl{:});
    trigtrl_list = vertcat(trigtrl{:});
  

    if length(trigtrl_list) < 800
        mergesubj{s,3} = 'trigger messed up';
        continue
    elseif length(trigtrl_list) > 960
        error('too many trials')
    end
    meg_rt = (alltrl_list(:,4) - alltrl_list(:,3))/1000;
    
    %% List trial type EL
    
    % extract el triggers
    idx = cellfun(@(x) regexp(x,'trig:'),el.Events.Messages.info,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    eltrials = el.Events.Messages.info(idxx);
    elsamples = el.Events.Messages.time(idxx);
    
    % extract trigger values
    extrtrig = cell2mat(cellfun(@str2num,cellfun(@(x) x(7:end),eltrials,'UniformOutput',false),'UniformOutput',false));
    
    % if practice trials are included (more than 24 1's = start of block) ->
    % only select experimental trials (last 24)
    begexp = find(extrtrig == 1);
    begexp = begexp(end-23:end);
    extrtrig = extrtrig(begexp(1):end);
    elsamples = elsamples(begexp(1):end);
    eltrials = eltrials(begexp(1):end);
    
    % eltrl: trigger, sample begin, sample end
    eltrl = [];
    trig_bsl = 1;
    while trig_bsl+2 <= length(extrtrig)
        
        % a: beginning of sample, a+1: onset search display, a+2: trial end (4
        if extrtrig(trig_bsl) >= 5 && extrtrig(trig_bsl) <= trigdef{end,1} && extrtrig(trig_bsl+1) == 2 && extrtrig(trig_bsl+2) == 4
            eltrl = [eltrl;extrtrig(trig_bsl) elsamples(trig_bsl) elsamples(trig_bsl+1) elsamples(trig_bsl+2)];
        end
        trig_bsl = trig_bsl + 1;
    end
    
    % reaction time according to eye tracker
    eltrl = [eltrl,zeros(size(eltrl,1),3)];
    elrt = zeros(size(eltrl,1),1);
    
    % extract samples corresponding to time
    fixdot = propixx_res./2; % coordinates fixation do
    
    % resample by hand if sampling rate higher than 1000 has been by accident, this shouldn't happen)
    bs = find(el.Samples.time == eltrl(1,2));
    el_fsample = 1000*length(bs);
    if el_fsample > 1000
        sr = el_fsample/1000;
        el_resamp = [];
        el_resamp.Samples.time = el.Samples.time(1:sr:length(el.Samples.time));
        el_resamp.Samples.posX = el.Samples.posX(1:sr:length(el.Samples.time));
        el_resamp.Samples.posY = el.Samples.posY(1:sr:length(el.Samples.time));
        el_resamp.Events.Esacc.start = el.Events.Esacc.start;
        el_resamp.Events.Ssacc.time = el.Events.Ssacc.time;
        el_resamp.Events.Esacc.duration  = el.Events.Esacc.duration;
        el_resamp.Events.Sblink.time = el.Events.Sblink.time;
        el_resamp.Events.Messages.info = el.Events.Messages.info;
        el_resamp.Events.Messages.time = el.Events.Messages.time;
        
        el = el_resamp;
        save(fullfile(dtpth,subjfolds{s},'el_struct.mat'),'el','-v7.3')
        el_fsample = 1000;
    end
    
    % store eye movement in matrix
    elxy = zeros(size(eltrl,1),5.5.*el_fsample,2);
    for e = 1:size(eltrl,1)
        % start sample: index of time
        bs = find(el.Samples.time == eltrl(e,2));   % first baseline sample
        vs = find(el.Samples.time == eltrl(e,3));   % first search sample
        es = find(el.Samples.time == eltrl(e,4));   % end sample

        % store sample info for trial
        eltrl(e,end-2:end) = [bs vs es];
        elxy(e,1:es(end)-bs+1,1) = el.Samples.posX(bs:es(end))';    % movement in x direction
        elxy(e,1:es(end)-bs+1,2) = el.Samples.posY(bs:es(end))';    % movement in y direction
        
        % check rt
        elrt(e) = (eltrl(e,4)-eltrl(e,3))/el_fsample;
    end
    
    %% Trial type mat file (responses stored during experiment)
    
    % list mat files
    d = dir(fullfile(dtpth,subjfolds{s}));
    
    f = {d.name};
    idx = cellfun(@(x) regexp(x,'mat'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    % delete practice mat file
    f = f(idxx);
    if length(f) > 1
        for h = 1:length(f)
            if strfind(f{h},'pract')
                delete(fullfile(dtpth,subjfolds{s},f{h}))
            end
        end
    end
    
    d = dir(fullfile(dtpth,subjfolds{s}));

    f = {d.name};
    idx = cellfun(@(x) regexp(x,'mat'),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    f = f(idxx);
    
    % find response files
    idx = cellfun(@(x) regexp(x,subjfolds{s}(end-3:end)),f,'UniformOutput',false);
    idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
    
    f = f(idxx);
    
    
    % if there are 2 response files, concatenate
    if length(f) == 1
        load(fullfile(dtpth,subjfolds{s},f{1}))
    else
        for i = 1:length(f)
            subj_all{i} = load(fullfile(dtpth,subjfolds{s},f{i}));
        end

        % set "subj" to first file and concatenate second file
        subj = subj_all{1}.subj;
        
        subj.exp.trials{1} = [];
        subj.exp.trials{2} = [];
        subj.srchdsp = {};
        subj.rspns = {};

        % 
        for i = 1:length(f)
            subj.exp.trials{1} = [subj.exp.trials{1};subj_all{i}.subj.exp.trials{1}];       % trial type
            subj.exp.trials{2} = [subj.exp.trials{2};subj_all{i}.subj.exp.trials{2}];       % baseline (always 1.5, this is included for historic reasons)

            subj.srchdsp = [subj.srchdsp;subj_all{i}.subj.srchdsp];                         % search display layout
            subj.rspns = [subj.rspns,subj_all{i}.subj.rspns];                               % responses
        end
        
        % number of blocks = 23
        subj.exp.trials{1} = subj.exp.trials{1}(end-23:end,:);
        subj.exp.trials{2} = subj.exp.trials{2}(end-23:end,:);
        subj.rspns = subj.rspns(end-23:end);
        subj.srchdsp = subj.srchdsp(end-23:end,:);
    end
    
    % find trial trigger, response and RT
    
    trltype = subj.exp.trials{1};
    
    rsptrig = zeros(size(trltype));
    rt = zeros(size(trltype));
    % loop over block
    for b = 1:size(trltype,1)
        % loop over trials
        for t = 1:size(trltype,2)
            curspex = cellfun(@num2str,trltype{b,t},'UniformOutput',false);
            curspex = [curspex{:}];
            
            % decode trial type string into trigger
            rsptrig(b,t) = trigdef{ismember(trigdef(:,2),curspex),1};

            % store reaction time
            if ~isempty(subj.rspns{b}{3,t})
                rt(b,t) = subj.rspns{b}{3,t};
            else 
                rt(b,t) = 4;      % if no response, RT = 4
            end
        end
    end

    rsptrig = rsptrig';   
    rsptrigc = reshape(rsptrig,[],1);
    
    % create response cell
    rsp = cell(size(rsptrigc,1),3);

    % convert to cell
    rsp(:,1) = arrayfun(@(x) x*1,rsptrigc,'UniformOutput',false);
    
    % helper cell
    h = cell(prod(size(rsptrig)),2);
    % fill in RT and performance
    c = 0;
    for b = 1:size(subj.rspns,2)
        for t = 1:size(subj.rspns{b},2)
            c = c + 1;

            % 1st: hit/miss/fa - 2nd: RT
            h(c,:) = {subj.rspns{b}{end,t}, subj.rspns{b}{3,t}};
            
            % if no RT = 4
            if isempty(subj.rspns{b}{3,t})
                h(c,:) = {subj.rspns{b}{end,t}, 4};
            end
        end
    end
       
    rsp(:,2:end) = h;
    
    meg_rt = floor(meg_rt.*100)./100;
    
    %% check and align response triggers and MEG triggers

    if length(rsp) > length(trigtrl_list)
        wm = 1;
        wr = 1;
        rej_rsp = [];

        % trigger stored in RSP
        rsp_trig = [rsp{:,1}]';
        % if response and MEG trigger are not the same -> delete those that
        % don't match
        while wm <= length(trigtrl_list) && wr < length(rsptrigc)

            % if response trigger and MEG trigger are not the sampe
            if rsp_trig(wr) ~= trigtrl_list(wm,1)
                % find the first one that is the same
                x = find(rsp_trig(wr:end)==trigtrl_list(wm,1),1);
                rej_rsp = [rej_rsp,wr:wr+x-2];
                wr = wr + x -1;
            else
                wr = wr +1;
                wm = wm +1;
            end

        end
    elseif length(rsp) < length(trigtrl_list)
        error('MEG triggers weird')
    else
        rej_rsp = [];
    end
    
    % convert rej_rsp to keep_rsp (from rejection to the ones to be kept)
    keep_rsp = ~ismember(1:size(rsp,1),rej_rsp);

    %% Trial rejection
    
    % RT read out by stim computer
    rt_array = [rsp{:,3}]';
    rt_array(rej_rsp) = [];
    
   
    % difference stim RT and MEG RT
    rt_diff = rt_array-meg_rt;

    % if difference larger than expected 14 ms > reject (something wonky
    % with triggers)
    idx_bigdiff = abs(rt_diff) > 0.014;
    rej_trl = logical(idx_bigdiff);
    
    keep_trl = ~rej_trl;

    %% save info in elinfo, meginfo, rspinfo
    
    mergesubj{s,2} = length(alltrl_list);
    
    % store all el info in one struct each per MEG, edf, response
    elinfo.eltrl = eltrl;                   % trial info
    elinfo.move = elxy;                     % movement over time
    elinfo.keeptrl_all = keep_trl;            % reject these for all
    elinfo.keep_rsp = keep_rsp;               % reject specifically for eyelink and rsp
    meginfo.keeptrl_all = keep_trl;           % reject for all
    meginfo.meg_rt = meg_rt;
    meginfo.trigtrl = trigtrl_list;
    meginfo.alltrl_list  = alltrl_list;
    meginfo.alltrl_bl = alltrl_bl;
    elinfo.fsample = el_fsample;
    
    % rsp info as struct
    rspinfo.trl = rsp;                      % trial info (trigger, performance, rt
    rspinfo.keeptrl_all = keep_trl;       % rejected for all
    rspinfo.keeptrl_rsp = keep_rsp;       % reject only for response

    mkdir(fullfile(pth,'results','meg', '2 merged edf mat', subjfolds{s}))


    mergesubj{s,3} = 'adjusted';
    mergesubj{s,4} = 'adjusted';

    mkdir(fullfile(trl_merge_pth, subjfolds{s}))
    save(fullfile(trl_merge_pth, subjfolds{s},...
        'trl_overlap_meg_el_rsp.mat'),'rspinfo','elinfo','meginfo','-v7.3')
    
   
    
    save(fullfile(trl_merge_pth,'docu_merge.mat'),'mergesubj')

    
end

save(fullfile(trl_merge_pth,'docu_merge.mat'),'mergesubj')
