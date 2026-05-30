%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% b. set-up PTB (window, screen, screen size)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023


%% PTB & screen settings

function [ptb, scr, exp] = b_setup_ptb(ptb,scr,exp,inlab,demo_ver)

% input: 
% - ptb: psychtoolbox settings
% - scr: screen settings (see main experiment script)
% - inlab: are we in the lab?
% - demo: showing demo version? (or quadrants?)
% output:
% - ptb: psychtoolbox settings
% - scr: screen settings

% start propixx in normal mode
if inlab
    Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram',0); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
end

% open window
ptb.screens      = Screen('Screens');
ptb.screenNumber = max(ptb.screens);

ptb.black        = BlackIndex(ptb.screenNumber);
ptb.white        = WhiteIndex(ptb.screenNumber);
ptb.colo         = ptb.colo.*ptb.white;
% get screen resolution - also convert to degree!
if inlab
    [ptb.window, ptb.windowRect] = PsychImaging('OpenWindow', ptb.screenNumber, ptb.black);
else
    [ptb.window, ptb.windowRect] = PsychImaging('OpenWindow', ptb.screenNumber, ptb.black, [0 0 1600 1000]);
end

Screen('Flip', ptb.window);                                                 % Flip to clear
scr.res = Screen('Resolution', ptb.screenNumber);

% get ifi
scr.ifi          = Screen('GetFlipInterval', ptb.window);                   % Measure refresh rate of the monitor
%ptb.trigfr       = round(ptb.trigtime/ptb.ifi);                             % trigger timing in framjes

% priority function
topPriorityLevel = MaxPriority(ptb.window);

% set blend function for the screen 
Screen('BlendFunction', ptb.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('TextFont', ptb.window, 'Arial');
Screen('TextSize', ptb.window, scr.fsizeinstr);

% size of screen
[scr.winW, scr.winH] = Screen('WindowSize', ptb.window);

% Get the centre coordinate of the window
[scr.xC, scr.yC]     = RectCenter(ptb.windowRect);

% for propixx: compute quadrants & centers of quadrants
[scr.quadrCenter, scr.quadrSize] = kdcompQuadr(scr.winW, scr.winH, [], []);

% stim size in cm and pixel (from degree)
ch = sqrt(scr.d^2+scr.h^2);      % hypothenuse (height screen)
cw = sqrt(scr.d^2+scr.w^2);      % hypothenuse (width screen)
if ~inlab
    scr.res.width  = scr.winW;
    scr.res.height = scr.winH;
end

if ~demo_ver
    scr.res.width = scr.res.width/2;
    scr.res.height = scr.res.height/2;
end
% sin(alpha) = a/c (Winkel = Gegenkathete/Hyopthenuse; GK = size in cm)
scr.scrdegrh    = asind(scr.h/ch);                          % screen height in degree
scr.scrdegrw    = asind(scr.w/cw);                          % screen width
scr.onedegrcm   = scr.h/scr.scrdegrh;                       % one degree in cm
scr.onedegrpix  = round(scr.res.height/scr.scrdegrh);       % one degree in pix
scr.stimcm      = scr.stimdegr .* scr.onedegrcm;            % stim size in cm
scr.stimpix     = scr.stimdegr .* scr.onedegrpix;           % stim size in pix
scr.searchcm    = scr.searchdegr .* scr.onedegrcm;          % search disp in cm
scr.searchpix   = scr.searchdegr .* scr.onedegrpix;         % search disp in pix

% if ~inlab
%     shrinkfact = [scr.winW scr.winH]./[scr.res.width, scr.res.height];             % shrink by this much if not tested in lab)
%     scr.stimpix = round(scr.stimpix .* shrinkfact(2));
%     scr.searchpix = round(scr.searchpix .* shrinkfact);
% end

%% Replicate flicker in upper and lower corner

% size photodiode flicker in pix
exp.rft.phtsz = 30;

% flickering dot in lower right corner 
for s = 1:size(scr.quadrSize,1)
    exp.rft.locpht_b(s,:) = [scr.quadrSize(s,3)-exp.rft.phtsz, scr.quadrSize(s,4) - exp.rft.phtsz,scr.quadrSize(s,3),scr.quadrSize(s,4)];
end

% flickering dot in upper right corner 
for s = 1:size(scr.quadrSize,1)
    exp.rft.locpht_t(s,:) = [scr.quadrSize(s,3)-exp.rft.phtsz, scr.quadrSize(s,2),scr.quadrSize(s,3),scr.quadrSize(s,2)+exp.rft.phtsz];
end
