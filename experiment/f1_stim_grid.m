%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% f1. define x & y coordinates to get L's and T's (in all orientations)
% define search grid positions for the different set sizes

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function [l,t,scr] = f1_stim_grid(scr)

%% Define item search grid and stimuli

% input
% - scr: screen settings

% output
% - l: coordinates stimulus l
% - t: coordin stim t
% - scr: add grid to screen settings

% 1. Stimuli
% hlf_stpx: use to centralise stim
hlf_stpx = scr.stimpix./2;

%% L
% normal L
x = [-hlf_stpx(1) scr.stimpix(2)+hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1) hlf_stpx(2)];
y = [-hlf_stpx(1)  -scr.stimpix(1)+hlf_stpx(2) scr.stimpix(2)-hlf_stpx(1) hlf_stpx(2)];
%l.y = [0 0 0 0];
l.coord{1} = [x;y] ;                 % normal L

% upside down
x = [-hlf_stpx(1) -scr.stimpix(2)-hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1) -hlf_stpx(2)];
y = [scr.stimpix(1)-hlf_stpx(1)  -hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1) scr.stimpix(1)-hlf_stpx(2)];
% hlf_stpx(2) = hlf_stpx(2)*(-1);
l.coord{2} = [x;y];%- hlf_stpx;       % upside down

% mirrored
x = [-hlf_stpx(1) scr.stimpix(2)-hlf_stpx(1)+hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1) hlf_stpx(1)+hlf_stpx(2)];
y = [scr.stimpix(1)-hlf_stpx(1)  -scr.stimpix(1)+hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1) hlf_stpx(1)+hlf_stpx(2)];
l.coord{3} = [x;y];                 % mirrored L

% mirrored upside down
x = [-hlf_stpx(1) scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1) -hlf_stpx(1)-hlf_stpx(2)];
y = [-hlf_stpx(1)  -hlf_stpx(1) scr.stimpix(2)-hlf_stpx(1) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2)];
l.coord{4} = [x;y];                 % mirrored L upside down


% T
% normal
x = [-hlf_stpx(1)-hlf_stpx(2) scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2) hlf_stpx(1)+hlf_stpx(2) -scr.stimpix(2)-hlf_stpx(2)];
y = [ceil((scr.stimpix(1)+scr.stimpix(2))/2)-scr.stimpix(2)/2-hlf_stpx(1)-hlf_stpx(2)  -hlf_stpx(1)-hlf_stpx(2) floor((scr.stimpix(1)+scr.stimpix(2))/2+scr.stimpix(2)/2)-hlf_stpx(1)-hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2)];
%t.y = [0 0 0 0];
t.coord{1} = [x;y];
% 90 degree
x = [-hlf_stpx(1)-hlf_stpx(2) -hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2) hlf_stpx(2)];
y = [scr.stimpix(1)-hlf_stpx(1)-hlf_stpx(2)  -scr.stimpix(2)-hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2) scr.stimpix(2)+hlf_stpx(2)];
%t.y = [0 0 0 0];
t.coord{2} = [x;y];
% 180 degree
x = [-hlf_stpx(1)-hlf_stpx(2) scr.stimpix(2)+hlf_stpx(2) scr.stimpix(1)+scr.stimpix(2)-hlf_stpx(1)-hlf_stpx(2) hlf_stpx(2)];
y = [ceil((scr.stimpix(1)+scr.stimpix(2))/2)-scr.stimpix(2)/2-hlf_stpx(1)-hlf_stpx(2)  -scr.stimpix(1)+hlf_stpx(2) floor((scr.stimpix(1)+scr.stimpix(2))/2+scr.stimpix(2)/2)-hlf_stpx(1)-hlf_stpx(2) hlf_stpx(2)];
t.coord{3} = [x;y];
% 270 degree
x = [-hlf_stpx(2) scr.stimpix(2)-hlf_stpx(2) scr.stimpix(1)-hlf_stpx(2) -hlf_stpx(2)];
y = [-scr.stimpix(2)-hlf_stpx(2)  -scr.stimpix(2)-hlf_stpx(2) -hlf_stpx(2) scr.stimpix(1)-hlf_stpx(2)];
t.coord{4} = [x;y];

%% Grid
% grid in degree distance from centre
% make grids for set size 32,24, 16

% size of display in degree
plsearch = prod(scr.searchdegr);

% scr.sets = set size
for st = 1:length(scr.sets)
    
    plgridpoint = plsearch/((scr.sets(st)));
    
    % outer edges of grid cells - half the length of the grid
    scr.grid.x{st} = [-scr.searchdegr(1)/2+round(sqrt(plgridpoint))/2:round(sqrt(plgridpoint)):scr.searchdegr(1)/2 - round(sqrt(plgridpoint))/2] ;
    scr.grid.y{st} = scr.grid.x{st};
    scr.grid.idco{st} = CombVec(scr.grid.x{st},scr.grid.y{st});
    scr.grid.xpix{st}  = scr.grid.x{st} .* scr.onedegrpix + scr.xC;
    scr.grid.ypix{st}  = scr.grid.y{st}.* scr.onedegrpix + scr.yC;
    scr.grid.co{st} = CombVec(scr.grid.xpix{st},scr.grid.ypix{st});
    
    % align to quadrant centre (grid for highfreqmode)
    for q = 1:size(scr.quadrCenter,1)
        scr.grid.highfreqx{st}{q}  = scr.grid.x{st} .* scr.onedegrpix + scr.quadrCenter(q,1);
        scr.grid.highfreqy{st}{q}  = scr.grid.y{st}.* scr.onedegrpix + scr.quadrCenter(q,2);
        scr.grid.highfreqco{st}{q} = CombVec(scr.grid.highfreqx{st}{q},scr.grid.highfreqy{st}{q});
    end
    
end
% permute
