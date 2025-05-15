%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% d: prep: alpha vs time-on-task!

% Input
% -s : subject index


% Output
% csv file per participant, containing RT and alpha power for each trial 

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 5 Aug 2023

%% Behavioural analyses
% a: performance per condition
% b. performance for alpha high vs low
% c: confirmatory analysis for median split sceptics: correlation alpha
% power ~ RT over trials
% d: does alpha power explain variance in RT and accuracy in addition to
% time-on-task?



function g_POW_GLM_condi(s)

onef = 0;

addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;

pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions','sinusoid','clean data all trials');
load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));
alphapth = fullfile(pth,'results','meg','6 Alpha');
soipth = fullfile(alphapth,'iaf_soi');

cohsoipth = fullfile(pth,'results','meg','5 COH hilb','soi','sinusoid');
outpth = fullfile(alphapth,'GLM of TFR int');
mkdir(outpth,subj{s})

load(fullfile(pth,'matlab scripts','coherence','occi_grad.mat'))
load(fullfile(pth,'matlab scripts','alpha','alpha_align_vec.mat'))

%% load data condition
d = dir(fullfile(inpth,subj{s}));
file = {d.name};
file(1:2) = [];

load(fullfile(inpth,subj{s},file{1}))

% extract specs condition
load(fullfile(pth,'experiment','trigdef.mat'))

condi_specs = {'ti','32t','tp'};

condi = zeros(length(behav_array),3);



% find trigger condition
% store: 1: guided; 1: set 32; 1: target present

for c = 1:length(condi_specs)
    trig = cell2mat(trigdef(cell2mat(cellfun(@(x) ~isempty(x), strfind(trigdef(:,2),condi_specs{c}),'UniformOutput',false)),1));
       
    idx = ismember([behav_array{:,1}],trig);
    condi(:,c) = idx;
end



% combination conditions
condi_comb = [[0,0,0];[0,0,1];[0,1,0];[0,1,1];[1,0,0];[1,0,1];[1,1,0];[1,1,1]];

% RT trimming

rt = [behav_array{:,3}]';

rej_trl = zeros(size(condi,1),1);

for c = 1:length(condi_comb)
    cur_cond = zeros(size(condi,1),1);

    for ti = 1:size(condi,1)
        cur_cond(ti) = isequal(condi(ti,:),condi_comb(c,:));
    end
        
    m = mean(rt(logical(cur_cond)));
    std_rt = std(rt(logical(cur_cond),1));
    
    rej_trl = rej_trl + ((cur_cond + (rt < m-3*std_rt) +  (rt > m+3*std_rt)) == 2);
end

% reject trials
condi = condi(~rej_trl,:);
tot = [1:length(behav_array)]./length(behav_array);

% demean tot and distance along z axis
tot = tot(~rej_trl)'-mean(tot(~rej_trl));

% change sign -> large positive should be large distance from initial
% position
dist_z_array = dist_z_array .* sign(dist_z_array);
dist_z = dist_z_array(~rej_trl) - mean(dist_z_array(~rej_trl));

% store rt
rt = rt(~rej_trl);

for c = 1:length(condi_comb)
    cur_cond = zeros(size(condi,1),1);

    for ti = 1:size(condi,1)
        cur_cond(ti) = isequal(condi(ti,:),condi_comb(c,:));
    end
    
    rt(logical(cur_cond)) = zscore(rt(logical(cur_cond)));
    
   % rt(logical(cur_cond)) = (rt(logical(cur_cond))-min(rt(logical(cur_cond))))./(max(rt(logical(cur_cond)))-min(rt(logical(cur_cond))));

end


% now add alpha power
load(fullfile(soipth,subj{s},'iaf_soi.mat'))
load(fullfile(cohsoipth,subj{s},'soi_stat.mat'))

% TFR
winl=0.5;
cfg = [];
cfg.method = 'mtmconvol';
%cfg.channel = soi_grad;
%cfg.channel = soi_occi;

cfg.channel = 'MEGGRAD';
cfg.taper = 'hanning';
cfg.foi = 4:1/winl:30;
cfg.t_ftimwin = ones(length(cfg.foi),1)*winl;
cfg.toi = -1.75:0.05:1;
cfg.keeptrials = 'yes';
cfg.trials = find(~rej_trl);

TFR_alpha = ft_freqanalysis(cfg,data);

cfg = [];
cfg.method = 'sum';
TFR_alpha = ft_combineplanar(cfg,TFR_alpha);


cfg = [];
cfg.latency = [-1.5 .5];
IAF = ft_selectdata(cfg,TFR_alpha);



%% 1/f correction

if onef
    %
    cfg = [];
    cfg.avgoverrpt = 'yes';
    cfg.avgovertime = 'yes';
    
    % separate search and baseline
    cfg.latency = [-1 -0.25];
    IAFbslavg = ft_selectdata(cfg,IAF);
    cfg.latency = [0 0.5];
    IAFstimavg = ft_selectdata(cfg,IAF);
    
    
    % baseline
    
    bslendtime = find(IAF.time == 0);
    
    logf = log10(IAF.freq);
    logpow = log10(IAF.powspctrm);
    logpowavg = squeeze(log10(IAFbslavg.powspctrm));
    
    for chan = 1:length(IAFstimavg.label)
        tmp_f = logf;                   % temporary freq vector that's changing
        tmp_pow_avg = logpowavg(chan,:);
        
        while true
            
            % linear fit
            b = [ones(size(tmp_f))' tmp_f'] \ tmp_pow_avg';
            linft = tmp_f*b(2) + b(1);
            
            firstpass = tmp_pow_avg - linft;
            mean_error = mean(firstpass(firstpass < 0));
            
            error_index = firstpass>abs(mean_error);
            
            if (sum(error_index)) == 0 || (numel(tmp_f) <= numel(logf)/2)
                
                linft = logf*b(2) + b(1);
                % repeat linear fit for each time point and trial
                linft_trl = repmat(reshape(linft,1,1,[],1),1,1,1,bslendtime);
                linft_trl = repmat(linft_trl,size(IAF.powspctrm,1),1,1,1);
                logpow(:,chan,:,1:bslendtime) = logpow(:,chan,:,1:bslendtime) - linft_trl;
                
                break
                
            else
                
                % remove frequencies and fit anew
                tmp_f = tmp_f(error_index ==0);
                tmp_pow_avg = tmp_pow_avg(error_index==0);
            end
        end
        
    end
    
    
    % stimulation
    
    logpowavg = squeeze(log10(IAFstimavg.powspctrm));
    
    for chan = 1:length(IAFstimavg.label)
        tmp_f = logf;                   % temporary freq vector that's changing
        tmp_pow_avg = logpowavg(chan,:);
        
        while true
            
            % linear fit
            b = [ones(size(tmp_f))' tmp_f'] \ tmp_pow_avg';
            linft = tmp_f*b(2) + b(1);
            
            firstpass = tmp_pow_avg - linft;
            mean_error = mean(firstpass(firstpass < 0));
            
            error_index = firstpass>abs(mean_error);
            
            if (sum(error_index)) == 0 || (numel(tmp_f) <= numel(logf)/2)
                
                linft = logf*b(2) + b(1);
                % repeat linear fit for each time point and trial
                linft_trl = repmat(reshape(linft,1,1,[],1),1,1,1,size(logpow,4)-bslendtime);
                linft_trl = repmat(linft_trl,size(IAF.powspctrm,1),1,1,1);
                logpow(:,chan,:,bslendtime+1:end) = logpow(:,chan,:,bslendtime+1:end) - linft_trl;
                
                break
                
            else
                
                % remove frequencies and fit anew
                tmp_f = tmp_f(error_index ==0);
                tmp_pow_avg = tmp_pow_avg(error_index==0);
            end
        end
        
    end
    
    % convert back
    IAF.powspctrm = 10.^logpow;
    
end

% Design matrix
regr_names = {'c','guided','set size','target present','rt','tot','slouch'};

X = ones(length(rt),7);

% guided
X(:,2) = condi(:,1).*2-1;
% set size
X(:,3) = condi(:,2).*2-1;
% target present
X(:,4) = condi(:,3).*2-1;

% rt
X(:,5) = rt;


% movement etc
X(:,6) = tot;
% slouch
X(:,7) = dist_z;

multi_coll = corr(X);

% model variance
modcov = inv(X'*X);

model_beta = zeros(size(X,2),length(IAF.label),length(IAF.freq),length(IAF.time));
model_T = model_beta;

for f = 1:length(IAF.freq)
    for ti = 1:length(IAF.time)
        for chan = 1:length(IAF.label)
           
        Y = zscore(IAF.powspctrm(:,chan,f,ti));

        beta = X\Y;
        n = size(Y,1);
        r = Y - (X*beta);                   % error
        SSE = sum(r.^2);                    % sum of squared errors
        MSE = SSE ./ (n-size(X,2));         % Mean sqaured error
        SE = diag(sqrt(MSE*modcov));
        T = beta ./ SE;

        model_beta(:,chan,f,ti) = beta;
        model_T(:,chan,f,ti) = T;

        end

    end
end

if onef
    save(fullfile(outpth,subj{s},'Power_GLM_onef_allset_condi.mat'),'model_beta','model_T','regr_names','multi_coll')
else
    save(fullfile(outpth,subj{s},'Power_GLM_allset_condi.mat'),'model_beta','model_T','regr_names','multi_coll')
end

