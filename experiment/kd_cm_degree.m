function stimsizedegr = kd_cm_degree(scr,stimsizecm)

ch = sqrt(scr.d^2+scr.h^2);      % hypothenuse (height screen)
stimsizedegr = asind(stimsizecm/ch)