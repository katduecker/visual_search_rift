%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% e. create random blocks
% We want to keep the number of target absent and present the same per
% block
% balance which tagging frequency is applied to which colour
% balance number of guided vs unguided (here called 'ti' and 'ni' for
% historical reasons)
% balance target colour
% balance set size (constant within block)

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

%% Function
% Input
% - ptb:       psychtoolbox settings
% - scr:       structure including screen settings
% - condinstr: cell array with possible instructions (e.g. no instruction 'ni'...)
% - ntr:       number of trials per condition combination (number)
% - stim:      cell array stimuli, here {'t','l'}

% Output

function [trials, condef] = e_create_blocks(ptb,scr,exp,prct)

% INPUT
% - ptb: psychtoolbox settings
% - scr: screen settings
% - exp: experiment settings
% - prct: practice trials?

% OUTPUT
% - trials: trial definition
%       trials{1}: block x trial per block cell containing trial specs
%       trials{2}: baseline length
% - trigs: trigger values

%%
% 1. Define blocks: 'ti' vs 'ni' (instruction (=guided) vs no instruction (unguided))

% experiment
if ~prct
    nblcks = exp.trttl/exp.blckl;
    instr = repmat(exp.condinstr,1,nblcks/length(exp.condinstr))';
    instr = instr(randperm(length(instr),length(instr)));
    trials{1} = instr;
    
% practice
else
   % trials{1} = cell(exp.prct.nblocks,exp.prct.nblocks);
    instr = repmat(exp.condinstr,1,exp.prct.nblocks/length(exp.condinstr))';
    instr = instr(randperm(length(instr),length(instr)));
    trials{1} = instr;
end

% 2. balance set size over 'ti' and 'ni'
stsz = repmat(scr.sets,1,size(trials{1},1)/length(scr.sets)/length(exp.condinstr))';

stsz_all = [];
% add as second column
for j = 1:length(exp.condinstr)
    
    % find strings
    fstr = trials{1}(find(strcmp(trials{1}(:,1),exp.condinstr{j})),1);
    stsz = stsz(randperm(length(stsz),length(stsz)));
    stsz_all = [stsz_all;stsz];
    for r = 1:length(fstr)
        fstr{r} = [fstr(r),stsz(r)];
    end
    trials{1}(find(strcmp(trials{1}(:,1),exp.condinstr{j}))) = fstr;
end

if ~prct
    trials{1} = repmat(trials{1},1,exp.blckl);
else
    trials{1} = repmat(trials{1},1,exp.prct.ntrials);   
end

% 3. within block, balance 'ta' and 'tp' (target absent, present)
% add present and absent to string
pres = repmat(exp.stim,1,2);
pres(1:length(pres)/2) = cellfun(@(x) strcat(x,'p'),pres(1:length(pres)/2),'Uniformoutput',false);          % present
pres(length(pres)/2+1:length(pres)) = cellfun(@(x) strcat(x,'a'),pres(length(pres)/2+1:length(pres)),'Uniformoutput',false);          % absent

% repmat within block
pres_all = [repmat(pres(1),size(trials{1},1),size(trials{1},2)/length(pres)),repmat(pres(2),size(trials{1},1),size(trials{1},2)/length(pres))];
for p = 1:size(pres_all,1)
    pres_all(p,:) = pres_all(p,randperm(length(pres_all)));
end

% 4. 'ti': balance target color over blocks, same within trials & same
% number per set size

% merge ti and ni into 1 matrix
tcol = zeros(size(trials{1}));
dcol = tcol;
for s = 1:length(scr.sets)
    tcol_ti = repmat(1:2,1,round(size(trials{1},1)/size(ptb.colo,1)/2/length(scr.sets)));
    tcol_ti = tcol_ti(randperm(length(tcol_ti),length(tcol_ti)));
    dcol_ti(tcol_ti == 1) = 2;
    dcol_ti(tcol_ti == 2) = 1;
    
    tcol_ti = repmat(tcol_ti',1,size(trials{1},2));
    dcol_ti = repmat(dcol_ti',1,size(trials{1},2));

    % find instruction and set size (ugly work-around)
    x = cellfun(@(x) ismember(x{1,1},'ti'),trials{1}(:,1),'UniformOutput',false);
    x = cell2mat(x);
    x = x(:,1);
    y = cellfun(@(x) (x{2} == scr.sets(s)),trials{1}(:,1),'UniformOutput',false);
    y = cell2mat(y);
    tcol(find((x + y) == 2),:) = tcol_ti;
    dcol(find((x + y) == 2),:) = dcol_ti;
    
%     x = cellfun(@(x) ismember(x{1,1},'ni'),trials{1}(:,1),'UniformOutput',false);
%     x = cell2mat(x);
%     x = x(:,1);
    clear tcol_ti dcol_ti
end

% 5. 'ni': balance target color over within trials, but not the same per
% block
% balance for present/absent

niidx = find(strcmp(instr,'ni'));                                       % no instruction index
colidx = repmat(1:2,1,size(trials{1},2)/length(pres)/2);               % color index vector
for n = 1:length(niidx) 
    for p = 1:length(pres)        
        % present/absent index
        cpidx = find(strcmp(pres_all(niidx(n),:),pres{p}));
        %numel(cpidx)
        tcol(niidx(n),cpidx) = colidx(randperm(length(colidx)));        % randomize color vector 
        
        % distractor color
        dcol(niidx(n),tcol(niidx(n),:) == 1) = 2;
        dcol(niidx(n),tcol(niidx(n),:) == 2) = 1;
    end
    
end


% 6. within all blocks, balance tagging frequencies
% balance over colours & present/absent
freq_t = zeros(size(trials{1}));
freq_d = freq_t;
for b = 1:size(trials{1},1)
    
    for p = 1:length(pres)
        cpidx = find(strcmp(pres_all(b,:),pres{p}));
        if strcmp(instr{b},'ti')

            % target instruction: balance over target color = length of trial
            rftfreqs = repmat(exp.rft.freq,1,size(trials{1},2)/length(exp.rft.freq)/length(pres));
            freq_t(b,cpidx) = rftfreqs(randperm(length(rftfreqs)));
            
        elseif strcmp(instr{b},'ni')
            % no target instruction: balance over each color (half of the
            % trial)
            rftfreqs = repmat(exp.rft.freq,1,size(trials{1},2)/length(exp.rft.freq)/size(ptb.colo,1)/length(pres));
            
            freq_t(b,cpidx(find(tcol(b,cpidx) == 1))) = rftfreqs(randperm(length(rftfreqs)));
            freq_t(b,cpidx(find(tcol(b,cpidx) == 2))) = rftfreqs(randperm(length(rftfreqs)));
 
        end
        clear cpidx
    end
    freq_d(b,freq_t(b,:) == exp.rft.freq(1)) = exp.rft.freq(2);
    freq_d(b,freq_t(b,:) == exp.rft.freq(2)) = exp.rft.freq(1);
end

condasstring = {};          % safe conditions as string
% add specs to blocks
for b = 1:size(trials{1},1)
    for t = 1:size(trials{1},2)
        condasstring = [condasstring;{[trials{1}{b,t}{1},num2str(trials{1}{b,t}{2}),pres_all{b,t},num2str(tcol(b,t)),num2str(dcol(b,t)),num2str(freq_t(b,t)),num2str(freq_d(b,t))]}];
        trials{1}{b,t} = cat(2,trials{1}{b,t},pres_all{b,t},num2cell(tcol(b,t)),num2cell(dcol(b,t)),num2cell(freq_t(b,t)),num2cell(freq_d(b,t)));
    end
end

condef = unique(condasstring);         % unique conditions
%% Jittered baseline -> we used a fixed baseline of 1.5 s in experiment!

% jittered baseline?
bsl = ones(size(trials{1}));
% round to 3rd decimal place: ms
% jitter between .5 and 1
if exp.jit
    bsl = bsl + round(.5.*rand(size(trials{1})),3)+.5;
else
    bsl = bsl + 0.5;
end
trials{2} = bsl;


%% Trigger (only once)
% % 7. find unique trial definitions and generate trigger numbers
for b = 1:size(trials{1},1)
    for t = 1:size(trials{1},2)
        x = cellfun(@num2str,trials{1}{b,t},'UniformOutput',false);
        funiquetrials{b,t} = strcat(x{:});
    end
end
uniquetrls = unique(funiquetrials);

% same number of combinations? - yes!!
for u = 1:length(uniquetrls)
    numtrialtype(u) = numel(find(strncmp(funiquetrials,uniquetrls{u}(1:8),8)));
end
 

% trigs = [1:4,5:4+length(uniquetrls)];
% 
% trigdef = {'instruct','beginning of trial','break','end'}';
% trigdef = [trigdef;uniquetrls];
% 
% trigdef = [num2cell(trigs)',trigdef];
% 
% save('trigdef.mat','trigdef');
