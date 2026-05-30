%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function:
% generate RIFT signal at certainf requency, format output such that it
% works for current propixx setting

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023


function rft = kd_rft(lsig, scrhz, freq, ppxset,ampl,meanlum,phshft)

%% creates tagging signal
% Inputs:
% - lsig:     length of signal in seconds
% - scrhz:    screen refresh rate in Hz
% - freq:     tagging frequency
% - ppxset:   propixx set
% - ampl:     amplitude of signal
% - meanlum:  mean luminance
% - phshft: phase shift

% Output:
% - rft:      sinusoidal tagging signal

% scrhz = scrhz*12;
% t = linspace(0, lsig, scrhz*lsig);
% sig = ampl*sin(2 * pi * freq * t)+meanlum;
if ppxset == 5
  scrhz = scrhz*12;
  t = linspace(0, lsig, scrhz*lsig);
  sig = ampl*sin(2 * pi * freq * t)+meanlum;
  rft = reshape(sig, 4, 3, []);
  
elseif ppxset == 2
  scrhz = scrhz*4;
  t = linspace(0, lsig, scrhz*lsig);
  sig = ampl*sin(2 * pi * freq * t + phshft*pi)+meanlum;
  % work around: 480 -> amplitudes caving in, cut out perfect sine wave!
  % (generated with 1440 Hz)
  rft = reshape(sig, 4,[]);
else
  rft = [];
end


end

