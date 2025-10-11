%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% f2. list participants with a significant tagging response
% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Output
% - usable_idx: index of subject with soi
% - num_soi: number of identified soi
% - subj: subject id

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions


clear all; close all; clc
pth = 'Z:\Visual Search RFT';
% addpath('/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT/matlab scripts/edf-converter-master'); %edf2mat converter
% rawdata 
megpth = fullfile(pth,'data');
% maxfiltered data
dtpth = fullfile(pth,'results','meg', '1 maxfilter', '1 maxfilter');
% ICA comps
icapth = fullfile(pth,'results','meg', '3 ICA');
% coherence path
cohpth = fullfile(pth,'results','meg', '5 COH hilb','soi','sinusoid');
mergepth = fullfile(pth,'results','meg', '2 merged edf mat');       % path containing trial structure

% list subj
d = dir(dtpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
fs = 1000;
clear d folds

load(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge.mat'))

% empty_idx = logical(cell2mat(cellfun(@isempty,mergesubj(:,2),'UniformOutput',false)));
% mergesubj = mergesubj(~empty_idx,:)
mergesubj = [mergesubj,cell(length(subjfolds),1)];
num_soi = zeros(length(subjfolds),1);
for s = 1:length(subjfolds)

   
    if exist(fullfile(cohpth,subjfolds{s},'soi_stat_consv.mat'),'file')
        load(fullfile(cohpth,subjfolds{s},'soi_stat_consv.mat'))
        num_soi(s) = length(soi_stat);
        load(fullfile(mergepth, subjfolds{s},'trl_overlap_meg_el_rsp.mat'),'meginfo')



        if isempty(soi_stat)
            mergesubj{s,5} = 'no soi';
        else
            soimag = soi_stat(logical(cell2mat(cellfun(@(x) strcmp(x(end),'1'),soi_stat,'UniformOutput',false))));
            soigrad = soi_stat(~ismember(soi_stat,soimag));

            %save(fullfile(cohpth,subjfolds{s},'soi_stat_not_align.mat'),'-append','soimag','soigrad')

            if isempty(soigrad)
                mergesubj{s,5} = 'no grad';
            else
                mergesubj{s,5} = 'soi identified';
            end

            save(fullfile(cohpth,subjfolds{s},'soi_stat_consv.mat'),'-append','soigrad','soimag')
        end
        clear soi_stat
    else
        mergesubj{s,5} = 'some issue';
    end

end

save(fullfile(pth,'results','meg', '2 merged edf mat','docu_merge_soi_coh.mat'),'mergesubj')

numberusable = sum(strcmp(mergesubj(:,5),'soi identified'));

mergesubj{19,5} = 'diode fell off';
% list all p's with raw data
d = dir(megpth);
folds = {d.name};
subjfolds = folds(strncmp(folds,'202',3));
[~,usable_idx] = intersect(subjfolds,mergesubj(strcmp(mergesubj(:,5),'soi identified'),1));

subj = subjfolds(usable_idx);

save(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi_consv.mat'),'usable_idx','num_soi','subj', 'noise')

