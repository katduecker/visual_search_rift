%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function: coherence with hilbert transform

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023

% Inputs: 
% - data: MEG data containing SOI and diodes
% - foi: frequencies of interest
% - frqwdth: width of passband (5 works well)
% - diodelabel: label of photodiode
% - filttype: cell array containing filter type and direction, e.g.
% {'but','twopass'}

% Output:
% coh_spect: coherence spectrum
% psd_meg: psd meg channels
% psd_misc: psd diode 
% csd_meg_misc: cross spectral density

%% Coherence
% a. Coherence per frequency (collapsed over trials)
% b. Coherence per condition for target and distractor colour -> priority
% map
% c. Coherence ~ reaction time: Coherence for fast vs. slow trials
% d. Coherence ~ alpha power: Coherence for alpha high vs low trials
% (median split)

function [coh_spct, psd_meg, psd_misc,csd_meg_misc] = kd_coh_hilb_fun(data, diodelabel,soi, foi, fwdth,filttype)

% loop over frequencies and do hilbert transform of data
for f = 1:length(foi)
    
    cfg = [];
    cfg.hilbert     = 'complex';
    cfg.keeptrials  = 'yes';
    % bp filter
    cfg.bpfilter = 'yes'; 
    cfg.bpfreq = [foi(f)-fwdth foi(f)+fwdth];
    cfg.bpfilttype = filttype{1};
    cfg.bpfiltdir = filttype{2};
    cfg.channel = {soi{:}, 'diode T', 'diode D'};
    data_foi = ft_preprocessing(cfg,data);
    
    % keep diode frequency
    cfg.channel = diodelabel;
    cfg.bpfreq = [foi(1)-fwdth foi(end)+fwdth];

    diode_foi = ft_preprocessing(cfg,data);
    
    data_foi = rmfield(data_foi, 'cfg');
    diode_foi = rmfield(diode_foi, 'cfg');
    
    meg_pad = zeros(length(data_foi.trial),length(data_foi.label),length(data_foi.time{1}));
    misc_pad = meg_pad;
    csd_pad = meg_pad;
   
    
    % loop over trials
    for t = 1:length(data.trial)
        
        % magnitude per trial
        
        % meg in soi
        meg  = data_foi.trial{t};
        % diode -> repmat to same size
        misc = repmat(diode_foi.trial{t},size(meg,1),1);
        
        %             % Mike X Cohen
        %             psd_meg_trl = abs(meg).^2;
        %             psd_misc_trl = abs(misc).^2;
        %
        %             csd_meg_misc_trl = abs(meg).*abs(misc).*exp(1i*(angle(meg)-angle(misc)));
        
        %% Yali
        % power spectral density in this trial
        psd_meg_trl = meg.*conj(meg);
        psd_misc_trl = misc.*conj(misc);
        
        % cross-spectral density in this trial
        csd_meg_misc_trl = meg .* conj(misc);
        
        % fill into trialxchannelxtime matrix
        meg_pad(t,:,1:length(psd_meg_trl)) = psd_meg_trl;
        misc_pad(t,:,1:length(psd_meg_trl)) = psd_misc_trl;
        csd_pad(t,:,1:length(psd_meg_trl)) = csd_meg_misc_trl;
        
    end
    
    mean_meg_pad = squeeze(mean(meg_pad,1,'omitnan'));
    mean_misc_pad = squeeze(mean(misc_pad,1,'omitnan'));
    mean_csd_pad = squeeze(mean(csd_pad,1,'omitnan'));
    
    % match fieldtrip format for plotting
    for c = 1:size(meg_pad,2)
        % psd: average over trials
        psd_meg(c,f,:) = mean_meg_pad(c,:);
        % average psd over trials -> copy such that size of psd of misc is the same
        % as psd of meg
        psd_misc(c,f,:) = mean_misc_pad(c,:);
        
        %csd by hand Yali
        csd_meg_misc(c,f,:) = abs(mean_csd_pad(c,:)).^2;

        % MXC
        %csd_meg_misc(c,f,:) = abs(mean_csd_pad(c,:));
    end
    clear data_foi
end
        
coh_spct = csd_meg_misc./(psd_meg.*psd_misc);           % coherence spectrum

