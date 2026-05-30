%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% pracicte experiment loop

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function [ptb,subj] = x_demo6_blockrft_practice(ptb,scr,exp,subj,el,inlab,demo_ver)


subj.prct.srchdsp = cell(size(exp.prct.trials{1}));
subj.prct.rspns = cell(1,size(exp.prct.trials{1},1));

for b = 1:size(exp.prct.trials{4},2)
    
    subj.prct.rspns{b} = cell(5,size(exp.prct.trials{1},2));
    % Cue current block
    if exp.ppxset == 5
        for q = 1:size(scr.quadrCenter,1)
            Screen('FillRect', ptb.window, ptb.white, [repmat(repmat(scr.quadrCenter(q,:),1,2),2,1)+exp.t.coord{1}]')
        end
    else
        switch exp.prct.trials{4}{b}
            case 'ti'
                for q = 1:size(scr.quadrCenter,1)
                    Screen('FillRect', ptb.window, ptb.colo(exp.prct.trials{5}(b),:)'.*0.5, [repmat(repmat(scr.quadrCenter(q,:),1,2),2,1)+exp.t.coord{1}]')
                end
            case 'ni'
                for q = 1:size(scr.quadrCenter,1)
                    Screen('FillRect', ptb.window, ptb.white, [repmat(repmat(scr.quadrCenter(q,:),1,2),2,1)+exp.t.coord{1}]')
                end
        end
    end
    ptb.vbl = Screen('Flip', ptb.window);
    
    % instruct trigger
    kd_trig(el,ptb,inlab,1)
    WaitSecs(ptb.trigtime)
    % reset
    kd_trig(el,ptb,inlab,0)
    WaitSecs(1-ptb.trigtime)
    
    % loop over trials within block
    for t = 1:size(exp.prct.trials{1}(b,:),2)
        
        % jittered baseline
        % fixation dot
        for q = 1:size(scr.quadrCenter,1)
            Screen('DrawDots', ptb.window, [0 0], 10, ptb.white, scr.quadrCenter(q,:),1);
        end
        ptb.vbl = Screen('Flip', ptb.window);
        
        % specs if the current trial
        curspex = exp.prct.trials{2}(exp.prct.trials{1}(b,t),:);
        % current trigger
        % string helper to find trigger
        x = cellfun(@num2str,curspex,'UniformOutput',false);
        curstr = strcat(x{:});
        % compare cells to find trigger
        ctrig = ptb.trigdef{find(strcmp(ptb.trigdef(:,2),curstr)),1};
        clear x curstr
        kd_trig(el,ptb,inlab,ctrig)
        
        WaitSecs(ptb.trigtime)
        % reset
        kd_trig(el,ptb,inlab,0)
        
        % rest of jittered baseline
        % prepare search display while baseline is on screen
        
        % how long does it take to make display?
        startdisp= GetSecs;
        % search display
        srchdsp = f2_search_displ(scr,curspex{2},exp.t,exp.l,curspex{3},ptb.colo(curspex{4},:),ptb.colo(curspex{5},:),demo_ver);
        
        subj.prct.srchdsp{b,t} = srchdsp;

        enddisp = GetSecs;
        
        % wait rest of time to present jittered baseline
        WaitSecs(exp.prct.trials{3}(b,t)-ptb.trigtime-(enddisp-startdisp))
        
        % present search display
        starttr = GetSecs;
        kd_trig(el,ptb,inlab,2)
        
        KbQueueStart();             % Start listening
        resp = 0;
        x = 0;
        while resp == 0
            x = x + 1;              % current frame
            for q = 1:size(scr.quadrCenter,1)
                
                for s = 1:size(srchdsp.set,2)
                    % select rft signal based on target/distr colour
                    if srchdsp.colid(s) == 1 || srchdsp.colid(s) == 0
                        f = find(exp.rft.freq == curspex{6});
                    elseif srchdsp.colid(s) == 2
                        f = find(exp.rft.freq == curspex{7});
                    end
                    
                    if exp.ppxset == 5
                        Screen('FillRect', ptb.window, ptb.white* exp.rft.sig{f}(q,:,x), srchdsp.sethf{q}{s});
                    elseif exp.ppxset == 2
                        Screen('FillRect', ptb.window, srchdsp.col(:,s).* exp.rft.sig{f}(q,x), srchdsp.sethf{q}{s});
                    end
                end
            end
            % flickering dots for photodiodes in upper and lower right corner
            for q = 1:size(scr.quadrCenter,1)
                if curspex{4} == 1
                    fyell = find(exp.rft.freq == curspex{6});
                    fteal = find(exp.rft.freq == curspex{7});
    
                elseif curspex{4} == 2
                    fyell = find(exp.rft.freq == curspex{7});
                    fteal = find(exp.rft.freq == curspex{6});
                end
                % top: yellow
                
                if exp.ppxset == 5
                    Screen('FillOval', ptb.window, ptb.white.* exp.rft.sig{fyell}(q,:,x), exp.rft.locpht_t(q,:));
                    % bottom: teal
                    Screen('FillOval', ptb.window, ptb.white.* exp.rft.sig{fteal}(q,:,x), exp.rft.locpht_b(q,:));
                else
                    Screen('FillOval', ptb.window, ptb.colo(1,:).* exp.rft.sig{fyell}(q,x), exp.rft.locpht_t(q,:));
                    % bottom: teal
                    Screen('FillOval', ptb.window, ptb.colo(2,:).* exp.rft.sig{fteal}(q,x), exp.rft.locpht_b(q,:));
                end
            end
             % fixation dot (if no eye movement)
            if exp.fxdt
                for q = 1:size(scr.quadrCenter,1)
                    Screen('DrawDots', ptb.window, [0 0], 10, ptb.white, scr.quadrCenter(q,:),1);
                end
            end

            ptb.vbl = Screen('Flip', ptb.window);
            
            % reset trigger after specified num of frames
            if x == ptb.trigfr
                kd_trig(el,ptb,inlab,0)
            end
            
            %% Listen to buttons
            [prs, frstprs] = KbQueueCheck(); %check response, return whether pressed, and first press timestamp
            
            if prs
                % Note that two keys may have been pressed
                kycd = find(frstprs);
                if kycd == ptb.keyb.escape
                    sca;
                elseif kycd == ptb.keyb.el_key
                    Eyelink('Command', 'driftcorrect_cr_disable = NO');
                    
                    EyelinkDoDriftCorrection(el.defaults);
                    Eyelink('Message', 'drift_correction');
                else
                    % two or more buttons pressed
                    if length(kycd) > 1
                        [~,ind]= min(frstprs(kycd));
                        kycd = kycd(ind); % select first response
                    end
                    
                    % store response
                    % 1. trial type
                    subj.rspns{b}{1,t} = curspex{3};
                    % 2. set size 
                    subj.rspns{b}{2,1} = curspex{2};
                    % 2. time point after trial onset = reaction time
                    subj.rspns{b}{3,t} = frstprs(kycd) - starttr;
                    % 3. key code
                    subj.rspns{b}{4,t} = kycd;
                    % 4. key name
                    subj.rspns{b}{5,t} = KbName(kycd);
                    % 5. hit/miss/fa?
                    % target present
                    if strcmp(curspex{3}(2),'p')
                        % if subject indicated presence
                        if find(ptb.rsp_pr == subj.rspns{b}{4,t})
                            subj.rspns{b}{6,t} = 'h';
                        else
                            subj.rspns{b}{6,t} = 'm';
                        end
                        % target absent
                    elseif strcmp(curspex{3}(2),'a')
                        % if subject indicated absence
                        if find(ptb.rsp_ab == subj.rspns{b}{4,t})
                            subj.rspns{b}{6,t} = 'h';                    % correct absence
                        else
                            subj.rspns{b}{6,t} = 'fa';
                        end
                        
                    end
                    
                    % change resp to 1
                    resp = 1;
                end
            end
            if x == exp.maxltrl * scr.res.hz && resp == 0
                subj.rspns{b}{1,t} = curspex{3};
                crtime = GetSecs;
                subj.rspns{b}{2,t} = crtime - starttr;
                subj.rspns{b}{3,t} = [];
                subj.rspns{b}{4,t} = [];
                subj.rspns{b}{6,t} = 'm';
                resp = 1;
            end
           % WaitSecs(0.1)
        end
        KbQueueStop();
        % black screen
        ptb.vbl = Screen('Flip', ptb.window);
        % send end trigger
        kd_trig(el,ptb,inlab,4);
        WaitSecs(ptb.trigtime);
        % reset
        kd_trig(el,ptb,inlab,0);
        ptb.vbl = Screen('Flip', ptb.window + scr.ifi);
        WaitSecs(0.5);
        
        % blink break every 10 trials
        if ~mod(t,10) && t ~= size(exp.trials{1},2)
          for q = 1:size(scr.quadrCenter,1)
            DrawFormattedText(ptb.window, 'Blink break!', 'center', scr.quadrCenter(q,2)-scr.onedegrpix*2, ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
            
            switch exp.prct.trials{4}{b}
              case 'ti'
                Screen('DrawTexture', ptb.window,exp.smiley.tx, [], exp.smiley.dstrct(q,:),[],[],[],ptb.colo(exp.prct.trials{5}(b),:)'.*0.5); % smiley face
              case 'ni'
                Screen('DrawTexture', ptb.window,exp.smiley.tx, [], exp.smiley.dstrct(q,:)); % smiley face
            end
          end
          ptb.vbl = Screen('Flip', ptb.window);
          WaitSecs(exp.blinkbreak)
          
        end
    end
    
    % feedback: number of hits
    numhit = sum(strcmp(subj.rspns{b}(6,:),'h'));
    
    fbscreen = kd_gen_feedback(numhit,exp,1);

    if b ~= 4
        
        fbscreen = [fbscreen, '\n', exp.break, '\n Press the "index-finger-button" to continue whenever you''re ready..'];
    else
        
        fbscreen = [fbscreen, '\n', exp.startexp];
    end
    % break
    for q = 1:size(scr.quadrCenter,1)
        DrawFormattedText(ptb.window, fbscreen, 'center', 'center', ptb.white,[],[],[],2,[],scr.quadrSize(q,:)); % 1.5: spacing
        
    end
    ptb.vbl = Screen('Flip', ptb.window);
    
    KbQueueStart();
    go = 0;
    go = kd_go_esc_el(go,ptb,el);
    
    
    if b ~= size(exp.prct.trials{4},2)
        g2_get_ready(ptb,exp,scr,el,0)
    end
    KbQueueStop();
    
end

end