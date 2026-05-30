%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% g2. get ready screen with countdown 
% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function g2_get_ready(ptb,exp,scr,el,expbeg)

% - expbeg: start of experiment (1) or block (0)

instruct = 'Get ready!';

if expbeg
    instruct = [instruct, ' The experiment will begin in '];
else
    instruct = [instruct, ' The next round will start in '];
end

for x = 3:-1:1
    for q = 1:size(scr.quadrCenter,1)
        DrawFormattedText(ptb.window, [instruct,num2str(x)], 'center', 'center', ptb.white,[],[],[],2,[],scr.quadrSize(q,:));
    end
    ptb.vbl = Screen('Flip', ptb.window);
    WaitSecs(1);
end
