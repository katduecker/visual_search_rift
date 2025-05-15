%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% c. Define delays in photodiode and reject trials with longer or shorter
% delays than usual (typically first trial in block; <17 ms or > 21ms)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Input
% s: subject index

% Output
% - sets meginfo.keeptrl_all = 0 for strange trials (more of a sanity
% check, shouldn't happen)

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement -> just sanity check, kept trials w/ saccades
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RIFT response
% g. Split trials into conditions

function c_delay_diode(s)
%% Paths
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';                                          % server path
ldiode = {'MISC004','MISC005'}; % diode label
dtpth = fullfile(pth, 'data'); % data path
addpath(fullfile('/rds/projects/j/jenseno-visual-search-rft/','fieldtrip'))
ft_defaults;
fs = 1000;
dtpth = fullfile(pth,'data');                                       % raw data
maxfpth = fullfile(pth,'results','meg', '1 maxfilter');             % max filter
scriptpth = fullfile(pth,'matlab scripts/','preprocessing MEG/');   % scripts
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure
icapth = fullfile(pth,'results','meg', '3 ICA', '2 subjoi');
respth = fullfile(pth,'results','meg','4 split conditions');

% list subj
d = dir(dtpth);
folds = {d.name};
subj = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds
% load in trial structure
load(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'))


% trial structure to load in trl
for p = 1:length(meginfo.alltrl_bl)
    trlstruct{p} = [meginfo.alltrl_bl{p}(:,3)-fs*2.5,meginfo.alltrl_bl{p}(:,3)+2*fs,zeros(length(meginfo.alltrl_bl{p}),1)-2.5*fs];
    trlstruct{p}(trlstruct{p}(:,1) <0,1) = 1;
end


% list maxfiltered data
d = dir(fullfile(maxfpth,subj{s}));
f = {d.name};
% find fif file
idx = cellfun(@(x) regexp(x,'fif'),f,'UniformOutput',false);
idxx = cell2mat(cellfun(@(x) ~isempty(x),idx,'UniformOutput',false));
f = f(idxx);


% load diodes
cfg = [];
for p = 1:length(f)
    cfg.dataset = fullfile(maxfpth,subj{s},f{p});
    cfg.preproc.detrend = 'yes';
    cfg.trl = trlstruct{p};
    % diodes
    cfg.channel = {'MISC004', 'MISC005'};
    diodes_trl{p}= ft_preprocessing(cfg);
end
diodes = ft_appenddata([], diodes_trl{:});

alltrl_list = meginfo.alltrl_list(meginfo.keeptrl_all);

%% Delays
% find first gradient that is not 0 -> delay
cfg = [];
cfg.latency = [0 0.1];
cfg.channel = {'MISC004', 'MISC005'};
misc = ft_selectdata(cfg,diodes);

% concatenate into a matrix
misc_mtrx = vertcat(misc.trial{:});

% matrices all zeros to fill in delays (gradient ~=0)
diode_delay_misc4 = zeros(size(misc_mtrx,1),1);
diode_delay_misc5 = zeros(size(misc_mtrx,1),1);

for m = 1:2:size(misc_mtrx,1)
    % misc 4
    % gradient
    misc_grad = round(gradient(misc_mtrx(m,:)),3);
    try
        % delay
        diode_delay_misc4(m) = find(misc_grad > 0,1);
    catch ME
        
    end
    
    % misc 5
    % gradient
    misc_grad = round(gradient(misc_mtrx(m+1,:)),3);
    try
        % delay
        diode_delay_misc5(m) = find(misc_grad > 0,1);
    catch ME
        
    end
    
end

diode_delay_misc4 = diode_delay_misc4(1:2:size(misc_mtrx,1));
diode_delay_misc5 = diode_delay_misc5(1:2:size(misc_mtrx,1));

diode_delay_off = logical((diode_delay_misc4==0)+(diode_delay_misc4<17)+(diode_delay_misc4>21)+(diode_delay_misc5==0)+(diode_delay_misc5<17)+(diode_delay_misc5>21));

meginfo.keeptrl_all(find(diode_delay_off)) = 0;
meginfo.keeptrl_all = logical(meginfo.keeptrl_all);
save(fullfile(mergepth, subj{s},'trl_overlap_meg_el_rsp.mat'),'-append','meginfo')


end


