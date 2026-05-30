%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% helper function:
% send trigger

% this is only called ones at the start of the experiment -> defines the
% first 4 triggers which are not condition specific

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function kd_trig(el,ptb,inlab,trig)

if inlab
    io64(ptb.ioObjTrig,ptb.PortAddress,trig);   
    if trig ~= 0 && el.on
        Eyelink('Message', ['trig: ',num2str(trig)])
    end
end