function []=Triggering()


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel Port IO & triggers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%initialization
PortAddress = hex2dec('D050'); %check this port address!
ioObjTrig = io64;

% initialize the interface to the inpoutx64 system driver
status = io64(ioObjTrig);

%send 0 trigger (reset all pins)
io64(ioObjTrig,PortAddress,0); %trigger 0 (reset)

%you can now send triggers 1-255 using io64(ioObjTrig,PortAddress,->trigger number<-)
%make sure you reset to 0 shortly after sending it (e.g. 50 or 100ms later)

%% %%%%%%%%%%%%%%%%%%%%%%
% Start experiment here %
%%%%%%%%%%%%%%%%%%%%%% %%

%Example 1:

%start main loop....

%some experiment code..

%Now you can use triggering in your experiment
%example (sending trigger 16 for 100ms):

io64(cfg.ioObjTrig,cfg.PortAddress,16); %send trigger 16 (pin 5)
WaitSecs(0.1)
io64(cfg.ioObjTrig,cfg.PortAddress,0); %set back to zero 

%some more experiment code


%--------------------------

%Example 2:

%start main loop....
%some experiment code..

%alternatively trigger using a function defined below
cfg=sendTrigger(cfg,16); %send trigger 16

%some more experiment code

%we want to reset the trigger 100ms after the last trigger
%place this somewhere in your main loop
if cfg.triggerSent && GetSecs>(cfg.triggerTime+0.1)
    io64(cfg.ioObjTrig,cfg.PortAddress,0);
    cfg.triggerSent=0;
end

%end main loop

%function to send MEG and triggers, and log triggers in logfile
function [cfg]= sendTrigger(cfg,trig)
    %send trigger to MEG
    io64(cfg.ioObjTrig,cfg.PortAddress,trig);
    cfg.triggerTime=GetSecs;
    cfg.triggerSent=1;
    
    %here you can also log the trigger in a logfile, eyelink etc..    
end

end