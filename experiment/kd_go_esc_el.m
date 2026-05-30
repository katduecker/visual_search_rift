%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function:
% Wait for keys: go, escape, re-do calibration?
% can re-do calibration

% is called in x_demo6 if one of the keys is pressed

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023



function go = kd_go_esc_el(go,ptb,el)
% Input:
% - go = 0, start value of go -> if go = 1 loop will end
% - ptb: psychtoolbox settings
% - el: eyelink settings
KbQueueStart();             % Start listening
while ~go
    [prs, frstprs] = KbQueueCheck(); %check response, return whether pressed, and first press timestamp
    if prs
        % Note that two keys may have been pressed
        kycd = find(frstprs);
        % escape
        if kycd == ptb.keyb.escape
            sca;
            error('Experiment terminated by user.')
            if el.on
                Eyelink('Message', 'ABORTED');
            end
        
        % go on
        elseif kycd == ptb.keyb.start || logical(sum(ptb.rsp_pr == kycd))
            go = 1;
            
        % calibration?
        elseif kycd == ptb.keyb.el_key && el.on
            
            % change propixx to normal mode
            Datapixx('SetPropixxDlpSequenceProgram',0); % 2 for 480, 5 for 1440 Hz, 0 for normal
            Datapixx('RegWrRd');
            
            % Re-do calibration
            EyelinkDoTrackerSetup(el.defaults);
            
            % re-start recording
            disp('Re-start recording')
            Eyelink('StartRecording');
            
            % record a few samples before we actually start displaying
            WaitSecs(0.1);
            
            % mark zero-plot time in data file
            disp('Sending message')
            Eyelink('Message', 'Re-start recording');
            
            Datapixx('SetPropixxDlpSequenceProgram',2); % 2 for 480, 5 for 1440 Hz, 0 for normal
            Datapixx('RegWrRd');
            go = 1;
        end
    end
end