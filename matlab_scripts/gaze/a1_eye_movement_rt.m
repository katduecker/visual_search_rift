%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b1. same as a1, but for fast vs slow trials

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Eye movement analysis
% a: ocular artefacts for alpha high low
% b: ocular artefacts for fast vs slow trials

function a1_eye_movement_rt(s)

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
pth = '/rds/projects/j/jenseno-visual-search-rft/visual_search_rift';

inpth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid');
mergepth = fullfile(pth,'results','meg','2 merged edf mat');
dtpth = fullfile(pth,'data');
cohpth = fullfile(pth,'results','meg','5 COH hilb', 'coh','sinusoid','conditions','alpha RFT');
occupth = fullfile(pth,'results','eyelink');
mkdir(occupth)
addpath(genpath(fullfile(pth,'matlab_scripts')))
condi = {{'ni','16t'},{'ti','16t'}, {'ni','32t'},{'ti','32t'}};

load(fullfile(pth,'matlab scripts/',"preproc_meg/",'idx_subjoi.mat'));

folds = subj;

% saccades matrix
sac_subj_condi_bslsearch_fs = zeros(length(condi),2,2);
% blink matrix
bl_subj_condi_bslsearch_fs = zeros(length(condi),2,2);

% eye bias matrix
bias_subj_condi_fs = zeros(length(condi),2);



%% load eye movement data
d = dir(fullfile(inpth,folds{s}));
d = {d.name};
files = d(strncmp(d,'dat',3));

load(fullfile(mergepth,folds{s},"trl_overlap_meg_el_rsp.mat"))

eye_move = elinfo.move(rspinfo.keeptrl_rsp,:,:);
eye_eltrl= elinfo.eltrl(rspinfo.keeptrl_rsp,:);

eye_move = eye_move(meginfo.keeptrl_all,:,:);
eye_eltrl = eye_eltrl(meginfo.keeptrl_all,:);


%% load response file with info on search display
d = dir(fullfile(dtpth,folds{s}));
d = {d.name};

% find matfile (there are some irregularities in how it's named, but it
% always contains subject code & .mat)
mat_file = d{(strncmp(d,folds{s}(end-3:end),4) + cell2mat(cellfun(@(x) ~isempty(x), regexp(d,'.mat'),'UniformOutput',false))) == 2};
load(fullfile(dtpth,folds{s},mat_file))


% bring trial defintions & search displays in correct order

trl_def = {};
trl_search_dsp = {};

% loop over blocks
for bl_idx = 1:size(subj.exp.trials{1},1)

    % store current trial definition
    cur_trl = cell(1,size(subj.exp.trials{1},2));
    for trl_idx = 1:size(subj.exp.trials{1},2)
        cur_trl{trl_idx} = strjoin(cellfun(@num2str, subj.exp.trials{1}{bl_idx,trl_idx}, 'UniformOutput',false),'');
    end

    trl_def = [trl_def;cur_trl'];
    % store search displays in this block
    trl_search_dsp = [trl_search_dsp;subj.srchdsp(bl_idx,:)'];
end

% decode trial info into trigger
% are triggers the same as in rspinfo?
load(fullfile(pth,'experiment','trigdef.mat'))

% find triggers related to trial definitions
trl_trig = zeros(size(trl_def));

for tg = 5:size(trigdef,1)

    idx = ismember(trl_def,trigdef{tg,2});

    trl_trig(idx) = trigdef{tg,1};
end

% same trial order as in rspinfo?
if ~isequal([rspinfo.trl{:,1}]',trl_trig)
    error('error in reshaping experiment file!')
end

% select trials to be kept
trl_trig = trl_trig(rspinfo.keeptrl_rsp);
trl_search_dsp = trl_search_dsp(rspinfo.keeptrl_rsp);

trl_trig = trl_trig(meginfo.keeptrl_all);
trl_search_dsp = trl_search_dsp(meginfo.keeptrl_all);

trl_trig = num2cell(trl_trig);

trl_search_dsp = [trl_trig,trl_search_dsp];

clear trl_def trl_trig



%% get saccade and blink timing
load(fullfile(dtpth,folds{s},'el_struct.mat'))

sac_time = el.Events.Ssacc.time;
blink_time = el.Events.Sblink.time;

clear el



%% list data files

d = dir(fullfile(inpth,folds{s}));
files = {d.name};
files(1:2) = [];


for c_idx = 1:length(condi)

    %% sort trials according to trial indices of the pre-processed data
    % select current trials

    ta_tp = {'ta','tp'};


    

    for ti = 1:length(ta_tp)
        
        % find relevant data files
        condi_files = zeros(length(files),1);
        for c = 1:length(condi{c_idx})
            
            condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,condi{c_idx}{c}),'UniformOutput',false))';
            
        end
        
        condi_files = condi_files + cell2mat(cellfun(@(x) ~isempty(x),regexp(files,[condi{c_idx}{2}(end-2:end-1),ta_tp{ti}]),'UniformOutput',false))';
        condi_files = condi_files == (length(condi{c_idx})+1);
        
        % load files 6067
        c_files = files(condi_files);
        
        load(fullfile(inpth,folds{s},c_files{1}));
        trl_idx = find(trlcur);

        perf_TFR_coh = perf_cur;

        eye_move_condi = eye_move(trlcur,:,:);
        eye_trl_condi = eye_eltrl(trlcur,:,:);
        srch_dsp_condi = trl_search_dsp(trlcur,:);

        % load & append performance and eye movement
        for f = 2:length(c_files)
            clear trlcur perf_cur
            load(fullfile(inpth,folds{s},c_files{f}), 'perf_cur','trlcur');

            perf_TFR_coh = [perf_TFR_coh;perf_cur];

            eye_move_condi = [eye_move_condi;eye_move(trlcur,:,:)];
            eye_trl_condi = [eye_trl_condi;eye_eltrl(trlcur,:)]; % eltrl contains samples and time points for baseline, search onset, search
            srch_dsp_condi = [srch_dsp_condi;trl_search_dsp(trlcur,:)];

        end

        % compare if data triggers and eye data triggers are the same
        if ~isequal([perf_TFR_coh{:,1}]',eye_trl_condi(:,1),[srch_dsp_condi{:,1}]')
            error('MEG and eyelink trials are not the same!')
        end


       

        %% select condition
        load(fullfile(pth, 'experiment','trigdef.mat'))
        rt = [perf_TFR_coh{:,3}];
        m_rt = median(rt);

        trl_fast = find(rt<m_rt);
        trl_slow = find(rt>m_rt);

        if length(trl_fast)>length(trl_slow)
            idx = randperm(length(trl_fast));
            trl_fast = sort(trl_fast(1:length(trl_slow)));
        elseif length(trl_fast)<length(trl_slow)
            idx = randperm(length(trl_slow));
            trl_slow = sort(trl_fast(1:length(trl_fast)));
        end

        % eye movement data and trial info for high vs low (first 500 ms)
        eye_move_fast = eye_move_condi(trl_fast,1500:2000-1,:);
        eye_move_slow = eye_move_condi(trl_slow,1500:2000-1,:);

        eye_trl_fast= eye_trl_condi(trl_fast,:);
        eye_trl_slow = eye_trl_condi(trl_slow,:);

        %% select search display for current condition

        srch_dsp_fast = srch_dsp_condi(trl_fast,:);
        srch_dsp_slow = srch_dsp_condi(trl_slow,:);

        %% find saccades that happened during trial
        
        % store number of saccades and blinks
        s_fast = zeros(size(eye_trl_fast,1),2);
        s_slow = zeros(size(eye_trl_fast,1),2);
        
        bl_fast = s_fast;
        bl_slow = s_fast;
        for t = 1:size(eye_trl_fast,1)
            
            % saccades during baseline
            s_fast(t,1) = sum((eye_trl_fast(t,2) < sac_time) + (eye_trl_fast(t,3) > sac_time) == 2);
            s_slow(t,1) = sum((eye_trl_slow(t,2) < sac_time) + (eye_trl_slow(t,3) > sac_time) == 2);
            
            % saccades during search (first 500 ms
            s_fast(t,2) = sum((eye_trl_fast(t,3) < sac_time) + (eye_trl_fast(t,3)+500 > sac_time) == 2);
            s_slow(t,2) = sum((eye_trl_slow(t,3) < sac_time) + (eye_trl_slow(t,3)+500 > sac_time) == 2);


            % blinks during baseline
            bl_fast(t,1) = sum((eye_trl_fast(t,2) < blink_time) + (eye_trl_fast(t,3) > blink_time) == 2);
            bl_slow(t,1) = sum((eye_trl_slow(t,2) < blink_time) + (eye_trl_slow(t,3) > blink_time) == 2);

            % saccades during search
            bl_fast(t,2) = sum((eye_trl_fast(t,3) < blink_time) + (eye_trl_fast(t,3)+500 > blink_time) == 2);
            bl_slow(t,2) = sum((eye_trl_slow(t,3) < blink_time) + (eye_trl_slow(t,3)+500 > blink_time) == 2);


        end


        %% Gaze bias towards Target colour


        % bin eye movement data
        % average location in 50 ms bins
        numbin = size(eye_move_slow,2)/50;


        eye_move_fast_bin = zeros(length(trl_fast),numbin,size(eye_move_fast,3));
        eye_move_slow_bin = eye_move_fast_bin;
        sb = 1;     % start of bin
        for b = 1:numbin
            eye_move_fast_bin(:,b,:) = mean(eye_move_fast(:,sb:sb+50-1,:),2);
            eye_move_slow_bin(:,b,:) = mean(eye_move_slow(:,sb:sb+50-1,:),2);
            sb = sb+50;
        end

        for t = 1:size(eye_trl_fast,1)



            %% Alpha high
            search_coord = srch_dsp_fast{t,2}.cxy;

            % x-distance to each stimulus (x) at each time bin (y)
            x_dist = squeeze(eye_move_fast_bin(t,:,1)) - search_coord(1,:)';
            % y distance
            y_dist = squeeze(eye_move_fast_bin(t,:,2)) - search_coord(2,:)';

            % euclidean distance
            dist = sqrt(x_dist.^2 + y_dist.^2);

            [min_dist, p_dist] = min(dist,[],1);


            % how many stimuli where looked at and how often?
            [count_time, which_stim] = groupcounts(p_dist');

            % stimulus that eye was closest to
            [~, p] = max(count_time);
            stim_idx = which_stim(p);

            % percentage of time spent close to a stimulus in T color
            % multiply logical array (Target yes no) with number of bins
            % for each target


            stim_id_fast(ti,t) = sum(strcmp(srch_dsp_fast{t,2}.stimidentity(which_stim),'lt').*(count_time'./sum(count_time)));


            %% Alpha low
            search_coord = srch_dsp_slow{t,2}.cxy;

            % x-distance to each stimulus (x) at each time bin (y)
            x_dist = squeeze(eye_move_slow_bin(t,:,1)) - search_coord(1,:)';
            % y distance
            y_dist = squeeze(eye_move_slow_bin(t,:,2)) - search_coord(2,:)';

            % euclidean distance
            dist = sqrt(x_dist.^2 + y_dist.^2);

            [min_dist, p_dist] = min(dist,[],1);


            % how many stimuli where looked at and how often?
            [count_time, which_stim] = groupcounts(p_dist');

            % stimulus that eye was closest to
            [~, p] = max(count_time);
            stim_idx = which_stim(p);

            % percentage of time spent close to a stimulus in T color
            % multiply logical array (Target yes no) with number of bins
            % for each target


            stim_id_slow(ti,t) = sum(strcmp(srch_dsp_slow{t,2}.stimidentity(which_stim),'lt').*(count_time'./sum(count_time)));

        end


    s_fast_ti(ti,1) = mean(s_fast(:,1));
    s_fast_ti(ti,2) = mean(s_fast(:,2));
    s_slow_ti(ti,1) = mean(s_slow(:,1));
    s_slow_ti(ti,2) = mean(s_slow(:,2));
    
    bl_fast_ti(ti,1) = mean(bl_fast(:,1));
    bl_fast_ti(ti,2) = mean(bl_fast(:,2));
    bl_slow_ti(ti,1) = mean(bl_slow(:,1));
    bl_slow_ti(ti,2) = mean(bl_slow(:,2));
    end

    % fast baseline
    sac_subj_condi_bslsearch_fs(c_idx,1,1) = mean(s_fast_ti(:,1));
    bl_subj_condi_bslsearch_fs(c_idx,1,1) = mean(bl_fast_ti(:,1));
    % fast search
    sac_subj_condi_bslsearch_fs(c_idx,2,1) = mean(s_fast_ti(:,2));
    bl_subj_condi_bslsearch_fs(c_idx,2,1) = mean(bl_fast_ti(:,2));

    % slow baseline
    sac_subj_condi_bslsearch_fs(c_idx,1,2) = mean(s_slow_ti(:,1));
    bl_subj_condi_bslsearch_fs(c_idx,1,2) = mean(bl_slow_ti(:,1));
    % slow search
    sac_subj_condi_bslsearch_fs(c_idx,2,2) = mean(s_slow_ti(:,2));
    bl_subj_condi_bslsearch_fs(c_idx,2,2) = mean(bl_slow_ti(:,2));


    % high alpha
    bias_subj_condi_fs(c_idx,1) = mean(stim_id_fast(:));

    % low alpha
    bias_subj_condi_fs(c_idx,2) = mean(stim_id_slow(:));

    clear s_high s_low bl_high bl_low



end

mkdir(fullfile(occupth,folds{s}))
save(fullfile(occupth,folds{s},'occu_fast_slow.mat'),'sac_subj_condi_bslsearch_fs','bl_subj_condi_bslsearch_fs','bias_subj_condi_fs')