%% Experiment: VS, alpha and RFT
%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e2. create some practice trials (this is done manually to ensure that the
% particiapants see a variation of all combinations)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Function
% Input
% - ptb:       psychtoolbox settings
% - scr:       structure including screen settings
% - condinstr: cell array with possible instructions (e.g. no instruction 'ni'...)
% - ntr:       number of trials per condition combination (number)
% - stim:      cell array stimuli, here {'t'}

% Output
% - randmat:   randomised output matrix with condition info for each trial
function trials = e2_practice(ptb,scr,exp)

% construct 64 trials in 4 blocks
trials{1} = reshape(1:4*16,[],4)';

% block color instruction

% target absence/presence
stimcat = repmat(exp.stim,1,2);
stimcat(1:length(stimcat)/2) = cellfun(@(x) strcat(x,'p'),stimcat(1:length(stimcat)/2),'Uniformoutput',false);          % present
stimcat(length(stimcat)/2+1:length(stimcat)) = cellfun(@(x) strcat(x,'a'),stimcat(length(stimcat)/2+1:length(stimcat)),'Uniformoutput',false);          % absent

stimcat = repmat(stimcat,1,prod(size(trials{1}))/length(stimcat));
stimcat = stimcat(randperm(length(stimcat),length(stimcat)));
stimcat = reshape(stimcat,size(trials{1},1),size(trials{1},2));

% set size
setsize([1,3],:) = repmat(scr.sets(1),2,size(stimcat,2));
setsize([2,4],:) = repmat(scr.sets(2),2,size(stimcat,2));

% colors
% target
colidx(1, :) = repmat(1,1,size(stimcat,2));
colidx(2,:) = randi([1,2],1,size(stimcat,2));
colidx(3,:) = randi([1,2],1,size(stimcat,2));
colidx(4, :) = repmat(2,1,size(stimcat,2));

% distractor
colidxd(1, :) = repmat(2,1,size(stimcat,2));
colidxd(2,find(colidx(2,:) == 1)) = 2;
colidxd(2,find(colidx(2,:) == 2)) = 1;
colidxd(3,find(colidx(3,:) == 2)) = 1;
colidxd(3,find(colidx(3,:) == 1)) = 2;
colidxd(4, :) = repmat(1,1,size(stimcat,2));

% freqs
freqst   = repmat(exp.rft.freq,4,size(trials{1},2)/length(exp.rft.freq));
freqsd   = repmat([exp.rft.freq(2), exp.rft.freq(1)],4,size(trials{1},2)/length(exp.rft.freq));



trials{4} = [{'ti'},{'ni'},{'ni'},{'ti'}];
trials{5} = [1,0,0,2];

allcond = [];
% make allcond matrix
for t = 1:size(trials{1},1)
    allcond = [allcond; repmat(trials{4}(t),size(trials{1},2),1),...
        num2cell(setsize(t,:))',stimcat(t,:)',num2cell(colidx(t,:))',num2cell(colidxd(t,:))',num2cell(freqst(t,:))',num2cell(freqsd(t,:))'];
end

trials{2} = allcond;
% jittered baseline?
bsl = ones(size(trials{2},1),1);
if exp.jit
    bsl = bsl + round(.5.*rand(size(trials{1})),3)+.5;
else
    bsl = bsl + 0.5;
end
% reshape baseline
trials{3} = reshape(bsl,size(trials{1},1),size(trials{1},2));