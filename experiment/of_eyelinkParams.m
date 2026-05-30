% eyelinkParams function
function [defaults] = of_eyelinkParams(defaults,ptb, rect,scr)
%Custom parameters for eyelink

    el = defaults;

    el.eye_used                = 'LEFT_EYE';
    el.calibrationtargetsize   = 1;
    el.calibrationtargetwidth  = 0.5;
    el.calibrationtargetcolour = 255;
    el.targetbeep              = 0;
    el.feedbackbeep            = 0;
    el.displayCalResults       = 1;
    el.msgfontcolour           = 255;
    % 100 -> 10*10 dispay spread (confirm measurements)
    el.eyeimgesize             = 100/(scr.scrdegrh*scr.scrdegrw)*100;  % percentage of screen (??)
    el.backgroundcolour        = ptb.black;

    disp('Updating Parameters')
    EyelinkUpdateDefaults(el);

    defaults = el;

    % make sure we're still connected.
    if Eyelink('IsConnected')~=1
        warning('eyelink is not connected! restart the tracker');
        Eyelink('Shutdown'); %shutdown Eyelink:
        el.online = 0;
        return;
    end

    % enable drift correction
      Eyelink('Command', 'driftcorrect_cr_disable = NO');

    % make sure that we get gaze data from the Eyelink
    Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,INPUT'); 

    % This Command is crucial to map the gaze positions from the tracker to
    % screen pixel positions to determine fixation
    Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld',  rect(1)-rect(1), rect(2)-rect(2), rect(3)-rect(1), rect(4)-rect(2));
    Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld',         rect(1)-rect(1), rect(2)-rect(2), rect(3)-rect(1), rect(4)-rect(2));

    % Use Psychophysical setting
    Eyelink('Command', 'recording_parse_type = GAZE');
    Eyelink('Command', 'saccade_velocity_threshold = 22');
    Eyelink('Command', 'saccade_acceleration_threshold = 3800');
    Eyelink('Command', 'saccade_motion_threshold = 0.0');
    Eyelink('Command', 'saccade_pursuit_fixup = 60');
    Eyelink('Command', 'fixation_update_interval = 0');

    % Other tracker configurations

    % these might crash:
    Eyelink('Command', 'heuristic_filter = 0');
    Eyelink('Command', 'pupil_size_diameter = YES');

    % use 9 point calibration (Default)
    Eyelink('Command', 'calibration_type = HV9');
    %Eyelink('Command', 'calibration_type = HV13');

    Eyelink('Command', 'generate_default_targets = YES');
    Eyelink('Command', 'enable_automatic_calibration = YES');
    Eyelink('Command', 'automatic_calibration_pacing = 1000');
    Eyelink('Command', 'binocular_enabled = NO');
    Eyelink('Command', 'use_ellipse_fitter = NO');
    Eyelink('Command', 'sample_rate = 1000');
    %Eyelink('Command', 'elcl_tt_power = %d', 2); % illumination, 1 = 100%, 2 = 75%, 3 = 50%

    switch el.eye_used
        case 'RIGHT_EYE'
            Eyelink('Command', 'file_event_filter = RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
            Eyelink('Command', 'link_event_filter = RIGHT,FIXATION,FIXUPDATE,SACCADE,BLINK,MESSAGE,INPUT');
        case  'LEFT_EYE'
            Eyelink('Command', 'file_event_filter = LEFT,FIXATION,SACCADE,BLINK,MESSAGE,INPUT');
            Eyelink('Command', 'link_event_filter = LEFT,FIXATION,FIXUPDATE,SACCADE,BLINK,MESSAGE,INPUT');
    end

    Eyelink('Command', 'file_sample_data  = GAZE,GAZERES,HREF,PUPIL,AREA,STATUS,INPUT');

end