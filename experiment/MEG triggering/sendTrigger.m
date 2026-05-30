function [cfg]= sendTrigger(cfg,trig)
    %send trigger to MEG
    io64(cfg.ioObjTrig,cfg.PortAddress,trig);
    cfg.triggerTime=GetSecs;
    cfg.triggerSent=1;
    
    %here you can also log the trigger in a logfile, eyelink etc..    
end