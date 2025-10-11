%% Minimum reaction time

function min_rt = kd_find_minrt(mergepth,folds)
min_rt = zeros(size(folds));
for s = 1:length(folds)
    load(fullfile(mergepth, folds{s},'trl_overlap_meg_el_rsp.mat'))
    
    trl = rspinfo.trl(rspinfo.keeptrl_rsp,:);
    RT_cond = [trl{meginfo.keeptrl_all,3}];
       
      
    % minimum RT -> 5% fastest trials
    min_rt_10 = mink(RT_cond,round(length(RT_cond)*0.05));
    
    % minimum RT -> fastest 1/3 of trials
    %min_rt_13 = mink(RT_cond,round(length(RT_cond)/3));
    min_rt(s) = mean(min_rt_10);
    clear RT_cond rspinfo trl trim_rt
end