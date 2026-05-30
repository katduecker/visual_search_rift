%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function kindly provided by Dr. Oscar Ferrante

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023


function of_eyelinkStop(el)

% STEP 7
% finish up: stop recording eye-movements,
% close graphics window, close data file and shut down tracker
Eyelink('StopRecording');
Eyelink('CloseFile');
% download data file

fprintf('Receiving data file ''%s''\n', el.edf_file);
%status=Eyelink('ReceiveFile');
status=Eyelink('ReceiveFile',el.edf_file,el.eye_dir,1); %transfer file to experiment directory
if status > 0
    fprintf('ReceiveFile status %d\n', status);
end
if 2==exist(el.edf_file, 'file')
    fprintf('Data file ''%s'' can be found in ''%s''\n', el.edf_file, el.eye_dir);
end

% Shutdown Eyelink:
Eyelink('Shutdown');
% Close window:
% sca;
% Restore keyboard output to Matlab:
ListenChar(0);

end


