function [] = bbox_test()
%
% Birimgham button box test
%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

% Clear the workspace and the screen
close all;
clear all;
sca;
fclose('all');
Screen('Preference', 'SkipSyncTests', 1); %must be 0 during experiment
AssertOpenGL;
KbName('UnifyKeyNames'); % for easy use of Keyboard keys
PsychDefaultSetup(2);    % call some default settings for setting up Psychtoolbox



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initalise Labjack / buttonbox

% Make the UD .NET assembly visible in MATLAB.
bBox.ljasm = NET.addAssembly('LJUDDotNet');
bBox.ljudObj = LabJack.LabJackUD.LJUD;

% Open the first found LabJack U3.
[bBox.ljerror, bBox.ljhandle] = bBox.ljudObj.OpenLabJackS('LJ_dtU3', 'LJ_ctUSB', '0', true, 0);

% Start by using the pin_configuration_reset IOType so that all pin
% assignments are in the factory default condition.
bBox.ljudObj.ePutS(bBox.ljhandle, 'LJ_ioPIN_CONFIGURATION_RESET', 0, 0, 0);

cfg.bBox=bBox;
%cfg.bBox.ljhandle=bBox.ljhandle;
cfg.bBox.ActiveKeys=16:17; %buttons to check

disp(['Checking button box keys ' int2str(min(cfg.bBox.ActiveKeys)) ' to ' int2str(max(cfg.bBox.ActiveKeys))]); 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial KeyBoard settings
KbName('UnifyKeyNames');
escKey = KbName('ESCAPE');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOOP
GetResp=1;
maxtime=600; 
t_start=GetSecs;

disp('Listening to button box in a continuous loop, press ESC to quit')

while GetResp %turns at eye movement
    
    [keyIsDown, t_keypress, keyCode] = BbCheck(cfg);
    [KbkeyIsDown, ~, KbkeyCode] = KbCheck;
    
    if KbkeyIsDown && KbkeyCode(escKey)
        disp('End of button box test')
        GetResp=0;
    end
    
    if keyIsDown
        disp(['Button ' int2str(keyCode) ' pressed at time ' num2str(t_keypress-t_start,2)]);
    else if GetSecs-t_start<maxtime %check whether we are still within allowed response time
            WaitSecs(0.001); %let's wait a whole millisecond and try again
        else
            GetResp=0;
        end
    end
end
        

disp('Waiting for button box key press (there is no escape...)')

[t_keypress, keyCode] = BbWait(cfg);
disp(['Button ' int2str(keyCode) ' pressed at time ' num2str(t_keypress-t_start,2)]);

disp('Great, done.')

        
    
%subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Button box resonse (blocking)
function [t,resp] = BbWait(cfg)

% loop till broken from
while true
    for i=1:length(cfg.bBox.ActiveKeys) %cycle through active keys/buttons
        
        cur=cfg.bBox.ActiveKeys(i);
        
        % get current state
        [~,cs]  = cfg.bBox.ljudObj.eDI(cfg.bBox.ljhandle, cur, 1); % check button
        
        % get time
        t       = GetSecs();
        
        % check if voltage exceeds threshold
        if cs < 1
            resp = cur;
            return
        end
    end
end
end

%Button box resonse (non-blocking)
function [KeyIsDown,t,resp] = BbCheck(cfg)

%check all active keys
for i=1:length(cfg.bBox.ActiveKeys) %cycle through active keys/buttons
    
    cur=cfg.bBox.ActiveKeys(i);
    
    % get current state
    [~,cs]  = cfg.bBox.ljudObj.eDI(cfg.bBox.ljhandle, cur, 1); % check button
    
    % check if voltage exceeds threshold
    if cs < 1
        KeyIsDown=1;
        t = GetSecs();% get time
        resp = cur;
        break;
    else
        KeyIsDown=0;
        t=[];
        resp=[];
    end
end
end

end