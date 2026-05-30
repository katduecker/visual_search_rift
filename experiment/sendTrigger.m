function [el] = sendTrigger(el,trig)
% send trigger to eyelink
if el.eyelink
    Eyelink('Message', ['Trigger_' int2str(trig)]);
end
end