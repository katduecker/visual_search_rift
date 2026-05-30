%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function:
% define triggers

% this is only called ones at the start of the experiment -> defines the
% first 4 triggers which are not condition specific

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023


function trigdef = kd_trifdef(trig,condef)

% conditions and trigger
trigs = trig(end)+1:size(condef)+trig(end);

trigdef = [num2cell(trigs)',condef];
