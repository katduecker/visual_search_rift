%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function kindly provided by Dr. Oscar Ferrante

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023



function [defaults] = of_eyelinkStart(edf_file,ptb, win, rect,scr)
% Open screen for calibration, calibrate and start recording

    % STEP 1
    % Open a graphics window on the main screen
    % using the PsychToolbox's Screen function.    
    % use the shrunk version of the window
    %window=Screen('OpenWindow', cfg.screenNumber, [] ,cfg.el_rect);
%     win = Screen('OpenWindow', screenid);

    % STEP 2
    % Provide Eyelink with details about the graphics environment
    % and perform some initializations. The information is returned
    % in a structure that also contains useful defaults
    % and control codes (e.g. tracker state bit and Eyelink key values).
    % Psychtoolbox defaults function
    defaults = EyelinkInitDefaults(win);
    defaults.window = win;
    % Disable key output to Matlab window:
    %ListenChar(2);

    % STEP 3
    % Initialization of the connection with the Eyelink Gazetracker.
    % exit program if this fails.
    if ~EyelinkInit
        fprintf('Eyelink Init aborted.\n');         
        Eyelink('Shutdown');%shutdown Eyelink:
        sca;% close window:
        ListenChar(0);%restore keyboard output to Matlab:
        return;
    else
        disp('Eyelink initizalized')
    end

    % open file to record data to
    disp('Opening EDF file'); 
    status = Eyelink('Openfile',[edf_file '.edf'] );

    if ~status
        disp('EDF file opened on Eyelink computer')
    else
        error(['Could not open EDF file on Eyelink computer, error: ' int2str(status)])
    end

    % set custom parameters
    disp('Setting parameters')
    defaults = of_eyelinkParams(defaults,ptb, rect,scr);

    % STEP 4 
    % Calibrate the eye tracker
    disp('Starting calibration')
    EyelinkDoTrackerSetup(defaults);

    % do a final check of calibration using driftcorrection
    %     EyelinkDoDriftCorrection(el);

    % STEP 5
    % start recording eye position
    disp('Start recording')
    %Screen('Close', win);
    Eyelink('StartRecording');
    % record a few samples before we actually start displaying
    WaitSecs(0.1);
    % mark zero-plot time in data file
    disp('Sending message')
    Eyelink('Message', 'SYNCTIME');

%     sca
    ListenChar(0);

end
