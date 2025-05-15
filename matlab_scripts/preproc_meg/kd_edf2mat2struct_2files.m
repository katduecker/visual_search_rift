%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function for a.: convert edf file to matlab structure (for 2 edf
% files)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023


% Input
% - edfconvpath: path where edf2mat is stored
% - subjpth: subject data path
% - file_in: file name .edf file
% - file_out: output filename

% Output
% - soi_stat: sensors with significant tagging response

%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions     

function el = kd_edf2mat2struct_2files(edfconvpath, subjpth, pth_in, edf_files, file_out)

addpath(edfconvpath)
for i = 1:length(edf_files)
    obj{i} = Edf2Mat(fullfile(pth_in,edf_files{i}));
end

% delete converted file
try delete(fullfile(subjpth,'el_edf2mat.mat'))
catch ME
end
    
% concatenate relevant el fields
el = [];
el.timeline = [obj{1}.timeline;obj{2}.timeline];
el.normalizedTimeline = [obj{1}.normalizedTimeline;obj{2}.normalizedTimeline];
el.Events.Messages.info = [obj{1}.Events.Messages.info,obj{2}.Events.Messages.info];
el.Events.Messages.time = [obj{1}.Events.Messages.time,obj{1}.Events.Messages.time(end)+obj{2}.Events.Messages.time];
el.Samples.time = [obj{1}.Samples.time;obj{1}.Samples.time(end)+obj{2}.Samples.time];
el.Samples.posX = [obj{1}.Samples.posX; obj{2}.Samples.posX];
el.Samples.posY = [obj{1}.Samples.posX; obj{2}.Samples.posY];
el.Events.Esacc.start = [obj{1}.Events.Esacc.start, obj{2}.Events.Esacc.start];
el.Events.Ssacc.time = [obj{1}.Events.Ssacc.time,obj{1}.Events.Ssacc.time(end)+obj{2}.Events.Ssacc.time];
el.Events.Esacc.duration = [obj{1}.Events.Esacc.duration,obj{2}.Events.Esacc.duration];
el.Events.Sblink.time = [obj{1}.Events.Sblink.time,obj{2}.Events.Sblink.time];

file_out = fullfile(subjpth,file_out);

if ~exist(file_out)
    save(file_out, "el",'-v7.3')
end
