%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% Set up presentation in high frequency projector
% compute quadrant centers for stimulus presentation and destination
% rectangles

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023


% Input arguments: 
% winW: window/screen width in pixels
% winH: window/screen height in pixels
% stimX: stimulus width/size in x direction in pixels
% stimY: stimulus height/size in y direction in pixels

function [quadrCenter, quadrSize, destinRect] = kdcompQuadr(winW, winH, stimX, stimY)

% cut screen in 4 quadrants
xQuad = winW/4;
yQuad = winH/4;

%% Quadrant centers
% top left quadrant
tl = [xQuad, yQuad];
% top right
tr = [3*xQuad, yQuad];
% bottom left
bl = [xQuad, 3*yQuad];
% bottom right
br = [3*xQuad, 3*yQuad];

quadrCenter    = [tl; tr; bl; br];

%compute quadrant size
% ecc: eccentricity of quadrants
ecc            = repmat([-tl,tl],4,1);
% repeat quadrCenter horizontally
repquadrCenter = repmat(quadrCenter,1,2);
quadrSize      = repquadrCenter + ecc;

%% Destination rectangles: 
if stimX
    % stimulus matrix
    stimSize = repmat([round(-stimX/2), round(-stimY/2), round(stimX/2), round(stimY/2)],4,1);
    % destination rectangle: center +- half the stimulus size
    destinRect = round(repquadrCenter + stimSize);
else
    destinRect = [];
    disp('Destination rectangle can not be computed. Not enough input arguments.')
end
end