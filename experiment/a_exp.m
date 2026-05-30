%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% a. experiment file (start from here)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

PsychDefaultSetup(2);
AssertOpenGL;
sca; % close screens
beep off;
close all;
clearvars;


%% a: Settings

% experimental settings
inlab            = 1;                           % are we in the lab?
demo_ver         = 0;                           % demo version?
el.on            = 0;
exp.tag          = 'a';                         % what are we tagging:
% 'a': all stimuli;
% 't': target only;
% 'd': distractor only
hdcurs           = 1;                           % hide cursor?
exp.fxdt         = 1;                           % present fixation dot throughout trial?

exp.condinstr    = {'ni','ti'};                 % 'ni' : no instruction; 'ti': target colour instruction
exp.rft.freq     = [60 67];
exp.block        = 1;                           % block design? 1 or o
exp.prct.on      = 0;                           % practice trials?

% keep /number of freq for now (important for balancing rft freq)
exp.ntr          = 60/length(exp.rft.freq);     % number of trials per condition, per target and per abs/pres, has to be divisible by 2
%exp.ntr          = 20/length(exp.rft.freq);
exp.activekey    = 1;                          % active keys only?
exp.stim         = {'t'};                      % stimuli
exp.maxltrl      = 4;
exp.jit          = 0;                       % jittered baseline?
exp.ppxset       = 2;
% screen settings
scr.w            = 72;             % screen width in cm
scr.h            = 40.5;              % screen height in cm
scr.d            = 142;             % distance from screen in cm

scr.fsizeinstr   = 60;
scr.stimdegr     = [.6 .3];         % stim size in degree
scr.searchdegr   = [12 12];         % search display [w h] in degree
scr.sets         = [32, 16];    % set sizes

% colour settings (full saturation)
%ptb.colo         = [255 0 229; 0 255 255] ./ 255;
ptb.colo         = [255 230 0; 0 255 255] ./ 255;
%ptb.colo         = [0 255 255; 255 0 229] ./ 255;
ptb.colw         = {'yellow','teal'};

% trial settings
exp.trttl = exp.ntr*length(exp.stim)*2*length(exp.condinstr)*size(ptb.colo,1)*length(scr.sets)*length(exp.rft.freq);             % total number of trials
% x = divisors(exp.trttl/length(exp.condinstr));
% exp.blckl = max(x(x <= (exp.trttl/length(exp.condinstr)/8)));
exp.blckl = 40;
clear x
% subj.domeye = input('subj dominant eye (RIGHT_EYE/LEFT_EYE)?', 's');



%% Logfile, info on participant
% add relevant paths
% paths
if inlab
  mpth           = 'C:\Users\MEG-STIM-01\Paradigms\dueckerk\VS_rft';
else
  mpth           = 'C:\Users\katha\Desktop\PhD\2_alpha_vs';
end
addpath(fullfile(mpth,'experiment','MEG triggering'))
cdpth          = fullfile(mpth,'experiment');
rsltpth        = fullfile(mpth,'results');
bttnpth        = fullfile(rsltpth,'responses');
addpath(cdpth)
if ~exist(bttnpth)
  mkdir(bttnpth)
end


subj.dt     = datetime;
subj.code   = input('subj code:','s');

while exist(fullfile(bttnpth,[subj.code,'.mat']))
  disp('subject code already exists, please define a different one.')
  subj.code = input('subj code:','s');
end

subj.age    = input('subj age:');
subj.gender = input('subj sex:','s');
subj.hdsc = input('subj handedness score:');
subj.hddecile = input('subj handedness decile:','s');
subj.augmhdsc = input('subj augmented handedness score:');


%% set up PTB
if ~demo_ver
  scr.fsizeinstr = round(scr.fsizeinstr/4);
end
% psychtoolbox set up
[ptb, scr, exp] = b_setup_ptb(ptb,scr,exp,inlab,demo_ver);


%% c: Eye tracker init

if el.on
  el.eyelink_key = KbName('E');  %key used to toggle eyelink validation/calculation on or off during experiment.
  el.continue_key = KbName('C'); %skip eye check part in start/end box, inorder to continue exp
  el.edf_file = [subj.code,'.edf'];
  [el.eye_dir] = of_eyelinkSetup(rsltpth, subj.code);
end

%% d: Buttons
% initialise keys
ptb = tg_init_key(ptb,exp.activekey,inlab);


%% e: Conditions
[exp.trials, condef] = e_create_blocks(ptb,scr,exp,0);

% triggers
if inlab
  ptb.PortAddress = hex2dec('BFF8');                  % check this port address!
  ptb.ioObjTrig   = io64;
  ptb.status      = io64(ptb.ioObjTrig);              % initialize the interface to the inpoutx64 system drive
end


ptb.trigtime     = 0.02;                             % send each trigger for 20 ms
ptb.trigfr       = round(ptb.trigtime/scr.ifi);

% trigger 1:4: 1: experiment start, 2: display onset, 4: end of trial, 3:
% ???
ptb.trigdef = kd_trigdef(1:4,condef);

% practice trials
if exp.prct.on
  exp.prct.nblocks = 4;
  exp.prct.ntrials = 16;
  exp.prct.trials = e2_practice(ptb,scr,exp);
  subj.exp.prct.trials = exp.prct.trials;
end

subj.exp.trials = exp.trials;

exp.done = 'Fantastic job, you are done! Thank you for participating in this study. \n Please wait for the experimenter.';
exp.feedback = {'Perfect score!', 'That''s awesome!', 'Nice job!', 'Don''t worry too much! Next one will be better!','Please let the experimenter know if something is wrong.'};
exp.fbnum = [1,.9,.75,.65];
exp.break  = 'Please take a break to rest your eyes.';
exp.startexp = 'Well done! \n Are you ready to start the experiment?';

% blink break smiley
exp.blinkbreak = 3;                                  % blink break
exp.smiley.img = imread(fullfile(cdpth, 'smiley.bmp'));
exp.smiley.img = uint8(exp.smiley.img).*255;
exp.smiley.img(:,:,2) = exp.smiley.img.*255;
exp.smiley.tx = Screen('MakeTexture', ptb.window, exp.smiley.img);
[~, ~, exp.smiley.dstrct] = kdcompQuadr(scr.winW, scr.winH, scr.onedegrpix*2, scr.onedegrpix*2);
%% f: Search display

[exp.l,exp.t,scr] = f1_stim_grid(scr);

%% g: RFT

% generate RFT signals
phshft = [0, 0.65];             % phase shift in cycles
for f = 1:length(exp.rft.freq)
  exp.rft.sig{f} = kd_rft(exp.maxltrl, scr.res.hz, exp.rft.freq(f), exp.ppxset,.5,.5,phshft(f));
end

% add third non flicker frequency
exp.rft.sig{f+1} = ones(size(exp.rft.sig{1})).*0.5;
%% update empty subject file to store trial settings, responses
% store search displays in subj file
subj.srchdsp = cell(size(exp.trials{1}));
% store response

%% Demo

% start eyelink
if el.on
  [el.defaults] = of_eyelinkStart(subj.code,ptb, ptb.window, ptb.windowRect,scr);
  Eyelink('Message', 'practice_start');
  
end
% start propixx in high freq mode
if inlab
  % Datapixx('Open');
  Datapixx('SetPropixxDlpSequenceProgram',exp.ppxset); % 2 for 480, 5 for 1440 Hz, 0 for normal
  Datapixx('RegWrRd');
  
  %   % psychtoolbox set up
  %   [ptb, scr, exp] = b_setup_ptb(ptb,scr,exp,inlab,demo_ver);
end

% send trigger zero to reset
kd_trig(ptb,inlab,0)

% hide cursor
if hdcurs
  HideCursor;
end


% %  practice trials
if exp.prct.on                          % practice trials?
  g_draw_instruct(ptb,exp,scr,el)
  
  go = 0;
  go = kd_go_esc_el(go,ptb,el);
  % HPI coil instruction
  for q = 1:size(scr.quadrCenter,1)
    DrawFormattedText(ptb.window, 'Please find a comfortable position for your head and sit still. \n Your head position will now be recorded.' , 'center', 'center', ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
  end
  ptb.vbl = Screen('Flip', ptb.window);
  if el.on
    Eyelink('Message', 'HPI coil measurement');
  end
  
  go = 0;
  go = kd_go_esc_el(go,ptb,el);
  
  % get ready
  g2_get_ready(ptb,exp,scr,el,1)
  [ptb,subj] = x_demo6_blockrft_practice(ptb,scr,exp,subj,el,inlab,demo_ver);
  
  save(fullfile(bttnpth,[subj.code,'pract.mat']),'subj')
  if el.on
    Eyelink('Message', 'END OF SESSION');
    of_eyelinkStop(el);
  end
  sca;
else
  
  % HPI coil instruction
  for q = 1:size(scr.quadrCenter,1)
    DrawFormattedText(ptb.window, 'Please find a comfortable position for your head and sit still. \n Your head position will now be recorded.' , 'center', 'center', ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
  end
  ptb.vbl = Screen('Flip', ptb.window);
  if el.on
    Eyelink('Message', 'HPI coil measurement');
  end
  
  go = 0;
  go = kd_go_esc_el(go,ptb,el);
  
  % % eyelink calibration again
  % if el.on
  % %    change propixx to normal mode
  %     Datapixx('SetPropixxDlpSequenceProgram',0); % 2 for 480, 5 for 1440 Hz, 0 for normal
  %     Datapixx('RegWrRd');
  %
  %     %Re-do calibration
  %     EyelinkDoTrackerSetup(el.defaults);
  %
  %     %re-start recording
  %     disp('Re-start recording')
  %     Eyelink('StartRecording');
  %
  %     %record a few samples before we actually start displaying
  %     WaitSecs(0.1);
  %
  %     %mark zero-plot time in data file
  %     disp('Sending message')
  %     Eyelink('Message', 'Re-start recording');
  %
  %     Datapixx('SetPropixxDlpSequenceProgram',exp.ppxset); % 2 for 480, 5 for 1440 Hz, 0 for normal
  %     Datapixx('RegWrRd');
  %
  % end
  
  % get ready
  g2_get_ready(ptb,exp,scr,el,1)
  
  % experiment
  [ptb,subj] = x_demo6_blockrft(ptb,scr,exp,subj,el,inlab,demo_ver);
  sca;
  ShowCursor;
  save(fullfile(bttnpth,[subj.code,'.mat']),'subj')
  
  if inlab
    % Datapixx('Open');
    Datapixx('SetPropixxDlpSequenceProgram',0); % 2 for 480, 5 for 1440 Hz, 0 for normal
    Datapixx('RegWrRd');
    
    %   % psychtoolbox set up
    %   [ptb, scr, exp] = b_setup_ptb(ptb,scr,exp,inlab,demo_ver);
  end
  % save(fullfile(bttnpth,['pseudotrig_test.mat']),'pseudotrig')
  
  % check pseudotrigger
  % for f = 5:36
  % checktrig(f) = numel(find(pseudotrig == f));
  % end
  % sum(checktrig) % ok!
  ShowCursor;
  %
  
  % stop eyelink
  if el.on
    Eyelink('Message', 'END OF SESSION');
    of_eyelinkStop(el);
  end
end

% x_demo3(ptb,scr,exp)
% HideCursor;
% [ptb,subj,info] = x_demo4_rft(ptb,scr,exp,subj,el,inlab);
% % 
% subj = xx_checkdisps(ptb,scr,exp,subj,demo_ver)
% % check search display size
% % spread in x direction
% for d = 1:length(subj.srchdsp)
%     spreadx = [min(subj.srchdsp{d}.cxy(1,:) ./ scr.onedegrpix), ...
%         max(subj.srchdsp{d}.cxy(1,:) ./ scr.onedegrpix)];
%     spready = [min(subj.srchdsp{d}.cxy(2,:) ./ scr.onedegrpix), ...
%         max(subj.srchdsp{d}.cxy(2,:) ./ scr.onedegrpix)];
%     
%     spread(d,:) = [diff(spreadx), diff(spready)];
% end
% % 
% 

% %ShowCursor;
% if inlab
%     % stop eyelink
%     if el.on
%         Eyelink('Message', 'end of experiment');
%         el_Stop(el);
%     end
%     
%     % store subj logfile
%    
% end


%x_demo3(ptb,scr,exp);

% check number of trials - ok!
% for t = 1:size(exp.trials{1},2)
%     ct(t) = 0;
%     for et = 1:size(exp.trials{1},2)
%         if isequal(exp.trials{2}(t,:),exp.trials{2}(et,:))
%             ct(t) = ct(t) + 1;
%         end
%     end
% end
% 
% find(ct ~= exp.ntr)