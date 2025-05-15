%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function for a.: convert edf file to matlab structure 

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

function el = kd_edf2mat2struct(edfconvpath, subjpth, file_in,file_out)

addpath(edfconvpath)
obj = Edf2Mat(fullfile(subjpth,file_in));

% delete converted file
try delete(fullfile(subjpth,'el_edf2mat.mat'))
catch ME
end
    
props = properties(obj);
for p = 1:numel(props)
    try
    el.(props{p})=obj.(props{p});
    catch ME
    end
end

file_out = fullfile(subjpth,file_out);

if ~exist(file_out)
    save(file_out, "el",'-v7.3')
end