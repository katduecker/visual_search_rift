%% Eye movement heatmap


function b1_eye_heatmap(s)

rmpath(genpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip/'))
addpath('/rds/projects/j/jenseno-visual-search-rft/fieldtrip')
ft_defaults;
pth = '/rds/projects/j/jenseno-visual-search-rft/Visual Search RFT';
inpth = fullfile(pth,'results','meg','4 split conditions', 'sinusoid','clean data all trials');
mergepth = fullfile(pth,'results','meg','2 merged edf mat');

occupth = fullfile(pth,'results','eyelink','heatmap');
addpath(genpath(fullfile(pth,'matlab scripts')))

load(fullfile(pth,'matlab scripts/',"preprocessing MEG/",'idx_subjoi.mat'));

mkdir(fullfile(occupth,subj{s}))


% screen settings

propixx_res      = [1920 1080];                                   % propixx resolution
scr.w            = 72;                                       % screen width in cm
scr.h            = 40.5;                                     % screen height in cm
scr.d            = 142;
scr.ch           = sqrt(scr.d^2+scr.h^2);                    % hypothenuse (height screen)
scr.scrdegrh     = asind(scr.h/scr.ch);                      % screen height in degree
scr.scrdegrw     = asind(scr.w/scr.ch);
scr.onedegrpix   = round(propixx_res(2)/scr.scrdegrh);       % one degree in number of pixel

% histogram edges

% start from center such that distance from center to edges is the same for
% both sides
% hist_bins_x = [sort(propixx_res(1)/2:-scr.onedegrpix/5:0), propixx_res(1)/2:scr.onedegrpix/5:propixx_res(1)];
% hist_bins_y = [sort(propixx_res(2)/2:-scr.onedegrpix/5:0), propixx_res(2)/2:scr.onedegrpix/5:propixx_res(2)];

hist_bins_x = linspace(0, propixx_res(1), propixx_res(1)/10);
hist_bins_y = linspace(0, propixx_res(2), propixx_res(2)/10);

%% load eyelink data
d = dir(fullfile(inpth,subj{s}));
d = {d.name};
files = d(strncmp(d,'dat',3));

load(fullfile(mergepth,subj{s},"trl_overlap_meg_el_rsp.mat"))


% create data strcuture

eye_data = struct();
eye_data.time = linspace(-.5,.5,1000);
eye_data.pos_x_y = elinfo.move(rspinfo.keeptrl_rsp,1000:2000,:);
eye_data.trial_info = elinfo.eltrl(rspinfo.keeptrl_rsp,:);

eye_data.pos_x_y = eye_data.pos_x_y(meginfo.keeptrl_all,:,:);
eye_data.trial_info = eye_data.trial_info(meginfo.keeptrl_all,:);


% eye data: 1500 samples baseline, then as many samples as RT (sampling f: 1000 Hz) confirm
% diff(eye_data.trial_info(:,2:4),[],2)
% 
% diff(eye_data.trial_info(:,5:end),[],2)

% eye movement pre-search
eye_data_pre = eye_data;
eye_data_pre.pos_x_y(:,1:500,:);


x_coord = eye_data_pre.pos_x_y(:,:,1);
y_coord = eye_data_pre.pos_x_y(:,:,2);

% this is the heatmap
for t = 1:size(x_coord,1)
   [pdf_eye_pre(t,:,:),Xedges,Yedges] = histcounts2(x_coord(t,:),y_coord(t,:),hist_bins_x,hist_bins_y,'Normalization','probability');
    
end

% repeated for during search
eye_data_search = eye_data;
eye_data_search.pos_x_y(:,500:end,:);


x_coord = eye_data_search.pos_x_y(:,:,1);
y_coord = eye_data_search.pos_x_y(:,:,2);

% this is the heatmap
for t = 1:size(x_coord,1)
   [pdf_eye_search(t,:,:),Xedges,Yedges] = histcounts2(x_coord(t,:),y_coord(t,:),hist_bins_x,hist_bins_y,'Normalization','probability');
end

subplot(121)
imagesc(Xedges,Yedges,squeeze(mean(pdf_eye_pre,1)))
subplot(122)
imagesc(Xedges,Yedges,squeeze(mean(pdf_eye_search,1)))

save(fullfile(occupth,subj{s},'eye_heatmap.mat'),'pdf_eye_pre','pdf_eye_search')
