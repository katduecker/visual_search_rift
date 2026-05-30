%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function kindly provided by Dr. Oscar Ferrante

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023



function [eye_dir] = of_eyelinkSetup(exp_dir, edf_file)

override = 0;

%add eyelink script folder (should be in main experiment folder)
addpath([exp_dir filesep 'Eyelink']);

%make directory if it doesn't already exist (local computer)
eye_dir = [exp_dir filesep 'Eyelink' filesep ];
if ~exist(eye_dir, 'dir')
    mkdir(eye_dir);
end

%check whether files already exist for this subject/session
if exist([exp_dir filesep 'Eyelink' filesep 'Data' filesep  edf_file '.edf'],'file')>0
    cont = input('Warning! Eyelink file will be overwritten, do you want to continue? (y/n) ','s');
    if cont == 'n'
        error('Session aborted')
    end
end

end
