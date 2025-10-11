%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% replace diode signal with sinusoid, amplitde 0.5, mean 0.5; noise
% amplitude defined by input

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 2 Aug 2023


% Input
% - data: data structure (fieldtrip)
% - diode_idx: diode indices
% - f: frequencies
% - phshft: phase shift of sinusoid
% - start_bsl: time point of bsl start relative to display onset in s

% Output
% - data: data structure with diode replaced w/ perfect sine wave


%% Preprocessing
% a. Define trial structure 
% b. Semi-automatic artefact rejection
% c. Define delays in photodiode and reject strange trials 
% d. Identify eye movement
% e. Run ICA with maximum 68 components
% f. Find sensors with significant RFT response
% g. Split trials into conditions

% PhD project 2

% replace diode signal with sinusoid
% ranging from 0 to 1, amplitude 0.5, mean 0.5%% VS + RFT


% [c] Katharina Duecker

function data = kd_replace_diode_sinu(data,diode_idx,rft_freq,phshft,start_bsl,fs,noise)

for t = 1:length(data.trial)
    
    % this is redundant, as there is only 1 diode
    for di = 1:length(diode_idx)

        % find gradient diode
        grad_diode = round(gradient(data.trial{t}(diode_idx(di),abs(start_bsl*fs):end)),3);
        % find delay for each trial
        diode_delay = find(grad_diode>0,1);
        if phshft(di) > 0
            diode_delay = find(grad_diode>0,1)+1;
        end

        % create time vector
        %tvec = 0:1/fs:data.time{1}(end)-diode_delay/fs+1/fs;
        %rft_sig = [zeros(1,abs(start_bsl*fs)+diode_delay-1),sin(2*pi*rft_freq(di)*tvec+phshft(di)).*0.5+0.5];
        
        tvec = data.time{t}(1:end-diode_delay);
        rft_sig = [zeros(1,diode_delay),sin(2*pi*rft_freq(di)*tvec+phshft(di)).*0.5+0.5];
        
        data.trial{t}(diode_idx(di),:) = rft_sig + randn(1,length(data.trial{t})).*noise;         % replace diode signal
    end
end