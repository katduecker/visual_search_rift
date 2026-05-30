%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function:
% generate feedback

% is called in x_demo6 after each block

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function fbscreen = kd_gen_feedback(numhit,exp,pract)

% - generates feedback based on number of hits

if pract
    ntrial = size(exp.prct.trials{1},2);
else
    ntrial = size(exp.trials{1},2);
end

fbscreen = ['Your score for this block is ', num2str(numhit), '/', num2str(ntrial),'!'];

if numhit/ntrial == exp.fbnum(1)
    fbscreen = [fbscreen, '\n',exp.feedback{1}];
elseif numhit/ntrial < exp.fbnum(end)
    fbscreen = [fbscreen, '\n',exp.feedback{end}];
else
    fbscreen = [fbscreen, '\n',exp.feedback{find(numhit/ntrial >= exp.fbnum,1)}];
end
