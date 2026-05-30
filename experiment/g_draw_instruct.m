%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% g. draw instructions

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function g_draw_instruct(ptb,exp,scr,el)


% instructions

% texts etc
instruct1 = 'Your task will be to indicate whether or not a T is presented among L''s. Note that the all stimuli will be presented in different orientations. \n For "present", press the "index-finger-button", for "absent", press the "middle-finger-button". \n Press the "index-finger-button" to continue.';
instruct2 = 'At the beginning of each block, a T will be presented, either in a colour or in white. \n A coloured T reveals the target''s colour for the following block.';
instruct3 = 'While this information is 100% valid, remember that the T will not appear in every trial.';
instruct4 = 'A white T indicates that the T can be shown in any colour throughout the block.';
instruct5 = 'Each trial begins with the presentation of a dot at the centre of the screen. \n Please fixate on the dot and only move your eyes when it disappears.';
if exp.fxdt
  instruct5 = 'Before and during the presentation of the T''s and L''s, a white fixation dot will be presented on the screen. \n Please fixate on the dot the entire time and try to find the T without moving your eyes.';

end
instruct6a = 'Please try to avoid eye blinks as much as you can, at least while you''re searching for the T.';
instruct6b = 'This screen will be presented every 8-10 trials for 3 seconds. \n Please utilise this time to rest your eyes.';
instruct6c = 'To help you remember the colour instruction, \n the smiley face will presented in the same colour as the T presented at the beginning of the block.';
instruct7 = 'Please solve the task as fast as you can, while being as accurate as possible! \n Are you ready to start with a few practice trials?';

% 1st instruction screen
for q = 1:size(scr.quadrCenter,1)
    DrawFormattedText(ptb.window, instruct1, 'center', 'center', ptb.white,[],[],[],2,[],scr.quadrSize(q,:)); 
end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct1');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);

% 2nd instruction: colour cue
for q = 1:size(scr.quadrCenter,1)
    % instruction above cue T
    DrawFormattedText(ptb.window, instruct2, 'center', scr.quadrCenter(q,2)-scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
    % cue T example
    Screen('FillRect', ptb.window, ptb.colo(1,:)'.*0.5, [repmat(repmat(scr.quadrCenter(q,:),1,2),2,1)+exp.t.coord{1}]')
    
    % instruction below cue t
    DrawFormattedText(ptb.window, instruct3, 'center', scr.quadrCenter(q,2)+scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:)); 

end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct colour cue');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);

% 3rd instruction screen: white cue
for q = 1:size(scr.quadrCenter,1)
    % instruction above cue T
    DrawFormattedText(ptb.window, instruct4, 'center', scr.quadrCenter(q,2)-scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
    % cue T example
    Screen('FillRect', ptb.window, ptb.white, [repmat(repmat(scr.quadrCenter(q,:),1,2),2,1)+exp.t.coord{1}]')
    
end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct white cue');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);

% 4th instruction: fixation
for q = 1:size(scr.quadrCenter,1)
    % instruction above cue T
    DrawFormattedText(ptb.window, instruct5, 'center', scr.quadrCenter(q,2)-scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
    % fixation dot
    Screen('DrawDots', ptb.window, [0 0], 10, ptb.white, scr.quadrCenter(q,:),1);
    
end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct fixation');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);

% instruction 4: blink break
for q = 1:size(scr.quadrCenter,1)
    DrawFormattedText(ptb.window, instruct6a, 'center', scr.quadrCenter(q,2)-scr.onedegrpix*3, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
       
    DrawFormattedText(ptb.window, 'Blink break!', 'center', scr.quadrCenter(q,2)-scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
    DrawFormattedText(ptb.window, instruct6b, 'center', scr.quadrCenter(q,2)+scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));

    Screen('DrawTexture', ptb.window,exp.smiley.tx, [], exp.smiley.dstrct(q,:)); % smiley face
    

end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct blink break');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);

% instruction 5: blink break + colour
for q = 1:size(scr.quadrCenter,1)
    DrawFormattedText(ptb.window, instruct6c, 'center', scr.quadrCenter(q,2)-scr.onedegrpix*3, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
       
    Screen('DrawTexture', ptb.window,exp.smiley.tx, [], exp.smiley.dstrct(q,:)); % smiley face
    
    Screen('DrawTexture', ptb.window,exp.smiley.tx, [], exp.smiley.dstrct(q,:),[],[],[],ptb.colo(1,:)'.*0.5); % smiley face


end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct blink break');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);


% instruction 5
for q = 1:size(scr.quadrCenter,1)
    DrawFormattedText(ptb.window, instruct7, 'center', 'center', ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
end
ptb.vbl = Screen('Flip', ptb.window);
if el.on
    Eyelink('Message', 'Instruct go');
end

go = 0;
go = kd_go_esc_el(go,ptb,el);