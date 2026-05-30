%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% initialize keys 
% [c]  Dr Tjerk Gutteling, 
% adaped by Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function ptb = tg_init_key(ptb,activeKeysOnly,inlab)

% INPUTS
% - ptb: psychtoolbox settings structure
% - activeKeysOnly: only listen to active keys? 0 or 1
% - inlab: are we in the lab? 0 or 1
% OUTPUT
% - ptb: adds subfield to ptb including all active keys

KbName('UnifyKeyNames');
% keyboard
ptb.keyb.escape = KbName('ESCAPE');
ptb.keyb.start  = KbName('LeftArrow');
ptb.keyb.space  = KbName('space');
ptb.rsp_pr        = [KbName('UpArrow')];    % keyboard (just demo checks) 'presence'
ptb.rsp_ab        = [KbName('RightArrow')];   % keyboard 'absence'

% Eyelink keys
ptb.keyb.el_key = KbName('E');              % key used to toggle eyelink validation/calculation on or off during experiment.
ptb.keyb.el_continue_key = KbName('C');     % skip eye check part in start/end box, inorder to continue exp

% NAtA
if inlab
    ptb.nata.pr_r     = KbName('7&');       % 'presence' right hand
    ptb.nata.ab_r     = KbName('8*');       % 'absence' right hand
    ptb.nata.pr_l     = KbName('4$');       % 'presence' left hand
    ptb.nata.ab_l     = KbName('3#');       % 'absence' left hand
    
    % concatenate to rsp matrices
    ptb.rsp_pr = [ptb.rsp_pr, ptb.nata.pr_r,ptb.nata.pr_l];
    ptb.rsp_ab = [ptb.rsp_ab, ptb.nata.ab_r,ptb.nata.ab_l];
end


device = 0; %default device

if activeKeysOnly % here we listen only to specific keys
    
    % keyboard
    ptb.keyb.active = [ptb.keyb.escape,ptb.keyb.start,ptb.keyb.space,ptb.rsp_pr,ptb.rsp_ab,ptb.keyb.el_key,ptb.keyb.el_continue_key];

    
    % create queue
    keylist=zeros(1,256); %Set all keys to zero (ignore)
    keylist(ptb.keyb.active)=1; %set active keys to 1 (listen)
    KbQueueCreate(device,keylist);%%Create queue, this is a time consuming operation (relatively), do while non-time critical
    disp('I will only listen to: '); disp(KbName(ptb.keyb.active));
else
    KbQueueCreate; %Default initialization, all keys listened to
end



end