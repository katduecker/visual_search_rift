%% Visual Search + RIFT
% Duecker, Shapiro, Hanslmayr, Wolfe, Pan, and Jensen

% f2. generate search display (coordinates for each stimulus, colour of
% each stimulus)
% saves coordinates for normal mode (can be compared to eyelink data) and
% high-frequency grid (one grid per quadrant)


% This function is called inside x_demo6 ... to generate as many search
% displays as there are trials before the start of the experiment

% [c] Katharina Duecker, katharina.duecker@gmail.com
% last changed/checked 3 Aug 2023

function [srchdsp] = f2_search_displ(scr,curss,t,l,curt,tcol,dcol,demo_ver)
%% Search grid and stimulus positions
% input
% - scr:        screen settings
% - rft:        frequency tagging settings
% - curss:      number of stimuli on current search display
% - t
% - l
% - curt:       current target presence or absence



% output:
% - search: structure containing current search settings
%       .set: rectangle coordinates (2 per stimulus)
%       .col: colour of rectangles  (+ luminance change by RFT)


%% Search display
% target present trials: 't', 'l'
% target absent trials: 'nt' (target: t but absent); 'nl' (target: l but
% absent)
% note: repmat: repeat in x and y direction (we have 4 coordinates
% [x,y,x,y] per rectangle, 2 rectangles per stimus

% normal mode grid
srchdsp.cxy = scr.grid.co{find(scr.sets == curss)};                      % copy positions


% high frequency grid
srchdsp.highfreqxy = scr.grid.highfreqco{find(scr.sets == curss)};
selor = randi(4,1,curss);                               % randomly select orientation
srchdsp.set = cell(1,curss);

% divide search display into quadrants
quad{1} = find(double(scr.grid.co{find(scr.sets == curss)}(1,:) < scr.xC) + double(scr.grid.co{find(scr.sets == curss)}(2,:) > scr.yC) == 2);
quad{2} = find(double(scr.grid.co{find(scr.sets == curss)}(1,:) > scr.xC) + double(scr.grid.co{find(scr.sets == curss)}(2,:) > scr.yC) == 2);
quad{3} = find(double(scr.grid.co{find(scr.sets == curss)}(1,:) > scr.xC) + double(scr.grid.co{find(scr.sets == curss)}(2,:) < scr.yC) == 2);
quad{4} = find(double(scr.grid.co{find(scr.sets == curss)}(1,:) < scr.xC) + double(scr.grid.co{find(scr.sets == curss)}(2,:) <scr.yC) == 2);


if curss == 24 && length([quad{:}]) > 24
    lq = cell2mat(cellfun(@length,quad,'UniformOutput',false));
    [v,p] = max(lq);
    delq = randperm(v,1);
    srchdsp.cxy(:,quad{p}(delq)) = [];
    for q = 1:size(scr.quadrCenter,1)
        srchdsp.highfreqxy{q}(:,quad{p}(delq)) = [];
    end
    quad{p}(delq) = [];
    
    lq = cell2mat(cellfun(@length,quad,'UniformOutput',false));
    
    mq = cell2mat(arrayfun(@(x) mod(x,2),lq,'UniformOutput',false));
    
end
% if curss > 16 % for largest set size: just schuffle
%     jitdegr = [scr.stimdegr(2)*.25 scr.stimdegr(2)]*100;
% elseif curss == 16            % smaller set size, dividable by 4 (i.e. 16, 8, 4)
%     jitdegr = [scr.stimdegr(2) scr.stimdegr(2)*2.5]*100;
% elseif curss < 16
%     jitdegr = [scr.stimdegr(2)/2 scr.stimdegr(2)*2]*100;
% end

srchdsp.colid = zeros(1,size(srchdsp.cxy,2));
srchdsp.col = zeros(3,curss);

if ~mod(curss,4) && curss > 4
    delcell = [];
    for q = 1:length(quad)
        if curss == 32 && ~isempty(quad{q})
            % delete cell
            delq = randperm(length(quad{q}),1);
            delcell = [delcell,quad{q}(delq)];
            quad{q}(delq) = [];
%         elseif curss == 24 && ~isempty(quad{q}) && length(quad{q})>6      
%             delq = randperm(length(quad{q}),length(quad{q})-6);
%             delcell = [delcell,quad{q}(delq)];
%             quad{q}(delq) = [];
        end
        
        col = [ones(1,length(quad{q})/2),ones(1,length(quad{q})/2)*2];
        col = col(randperm(length(col),length(col)));
        srchdsp.colid(quad{q}) = col;
        clear col
    end

    if curss == 32
        srchdsp.cxy(:,delcell) = [];
        for q = 1:size(scr.quadrCenter,1)
            srchdsp.highfreqxy{q}(:,delcell) = [];
        end
    end

    srchdsp.colid(srchdsp.colid == 0) = [];
    srchdsp.col(:,srchdsp.colid == 1) = repmat(tcol',1,curss/2);
    % position target
    allones = find(srchdsp.colid== 1);
    allones = allones(randperm(length(allones),length(allones)));
    srchdsp.colid(allones(1)) = 0;
    srchdsp.col(:,srchdsp.colid == 2) = repmat(dcol',1,curss/2);

    
elseif curss == 4
    % colour
    srchdsp.colid = [ones(1,2),ones(1,2)*2];
    srchdsp.colid = srchdsp.colid(randperm(4,4));
    srchdsp.col(:,srchdsp.colid == 1) = repmat(tcol',1,curss/2);
    srchdsp.col(:,srchdsp.colid == 2) = repmat(dcol',1,curss/2);

    allones = find(srchdsp.colid== 1);
    allones = allones(randperm(length(allones),length(allones)));
    srchdsp.colid(allones(1)) = 0;
    
    % divide curss by 4 = number of postions selected per quadrant
    selpos = [quad{1}(randperm(length(quad{1}),curss/4)),quad{2}(randperm(length(quad{2}),curss/4)),...
        quad{3}(randperm(length(quad{3}),curss/4)),quad{4}(randperm(length(quad{4}),curss/4))];
    
        % shuffle positions
        %selpos = selpos(randperm(length(selpos)));
        % select positions
        srchdsp.cxy = srchdsp.cxy(:,selpos);
%     else
%         srchdsp.cxy = srchdsp.cxy(:,randperm(length(srchdsp.cxy),length(srchdsp.cxy)));
%         %if too many of the same colour in one quadrant
elseif curss == 2                           % smalles set size = 2
        jitdegr = [scr.stimdegr(2)/2 scr.stimdegr(2)*1.25]*100;

    % colour
    srchdsp.colid = [0,2];
    srchdsp.colid = srchdsp.colid(randperm(2,2));
    srchdsp.col(:,srchdsp.colid == 0) = repmat(tcol',1,curss/2);
    srchdsp.col(:,srchdsp.colid == 2) = repmat(dcol',1,curss/2);
    % half display
    half1 = find(double(srchdsp.cxy(1,:) < scr.xC));
    half2 = find(double(srchdsp.cxy(1,:) > scr.xC));

    selpos = [half1(randperm(length(half1),1)),half2(randperm(length(half2),1))];
    
    % shuffle positions
    for i = 1:10
        selpos = selpos(randperm(length(selpos)));
    end
    % select positions
    srchdsp.cxy = srchdsp.cxy(:,selpos);
    
    % select positions in quadrants
    for q = 1:size(scr.quadrCenter,1)
        srchdsp.highfreqxy{q}(:,selpos);
    end
end



% number of stim that fit in between centre
% Jitter: make sure that there is always a minimum of 1.5 stim sizes between
% stimuli
stimexpand = abs(t.coord{1}(1,1)-t.coord{1}(1,3));

% distances grid center
distgridcent = mean(diff(scr.grid.xpix{find(scr.sets == curss)}));

% jitter:
distendgrid = distgridcent/2-stimexpand/2;              % distance to boarder of grid
maxjit = distendgrid - scr.stimpix(2);% maximum jitter -> up to half a stim size away from grid

% half of the stimuli in - half in + direction, normally distributed around
% maxjit/2
jitpix = randi([-10 10],1,curss)./10*maxjit;
jitpix = [jitpix(randperm(length(jitpix),length(jitpix)));jitpix(randperm(length(jitpix),length(jitpix)))];

if demo_ver
    srchdsp.cxy = srchdsp.cxy + round(jitpix);
else
    srchdsp.cxy = srchdsp.cxy + round(jitpix)*2;                   % positions + jitter (*2 for non flicker)
    
end

% jitter in quadrants
for q = 1:size(scr.quadrCenter,1)
    srchdsp.highfreqxy{q} = srchdsp.highfreqxy{q} + jitpix;
end
srchdsp.set = cell(1,curss);
srchdsp.stimidentity = cell(1,curss);

% quadrant sets
srchdsp.sethf = cell(1,4);
for q = 1:size(scr.quadrCenter,1)
    srchdsp.sethf{q} = cell(1,curss);
end

% place distractors first
if strncmp(curt,'t',1)
    
    for c = 1:curss
        if srchdsp.colid(c) ~= 0
            srchdsp.set{c} = [l.coord{selor(c)}+repmat(repmat(srchdsp.cxy(:,c)',1,2),2,1)]';
            if srchdsp.colid(c) == 1
                srchdsp.stimidentity{c} = 'lt';                    % identity
            elseif srchdsp.colid(c) == 2
                 srchdsp.stimidentity{c} = 'ld';                    % identity

            end

            for q = 1:size(scr.quadrCenter,1)
                srchdsp.sethf{q}{c} = [l.coord{selor(c)}+repmat(repmat(srchdsp.highfreqxy{q}(:,c)',1,2),2,1)]';
            end
        end % select current orientation l and place it on grid
    end
    
elseif strncmp(curt,'l',1)
    
    for c = 1:curss
        if srchdsp.colid(c) ~= 0
            srchdsp.set{c} = [t.coord{selor(c)}+repmat(repmat(srchdsp.cxy(:,c)',1,2),2,1)]';
            
            srchdsp.stimidentity{c} = 't';
            for q = 1:size(scr.quadrCenter,1)
                srchdsp.sethf{q}{c} = [t.coord{selor(c)}+repmat(repmat(srchdsp.highfreqxy{q}(:,c)',1,2),2,1)]';
            end
        end% select current orientation l and place it on grid
    end
    
end

% if it's a target present trial
c = find(srchdsp.colid == 0);
if strcmp(curt,'tp')
    srchdsp.set{c} = [t.coord{selor(1)}+repmat(repmat(srchdsp.cxy(:,c)',1,2),2,1)]';
    srchdsp.stimidentity{c} = 't';
    
    for q = 1:size(scr.quadrCenter,1)
        srchdsp.sethf{q}{c} = [t.coord{selor(1)}+repmat(repmat(srchdsp.highfreqxy{q}(:,c)',1,2),2,1)]';
    end% select current orientation t and place it on grid
elseif strcmp(curt,'lp')
    srchdsp.stimidentity{c} = 'l';

    srchdsp.set{c} = [l.coord{selor(1)}+repmat(repmat(srchdsp.cxy(:,c)',1,2),2,1)]';
    for q = 1:size(scr.quadrCenter,1)
        srchdsp.sethf{q}{c} = [l.coord{selor(1)}+repmat(repmat(srchdsp.highfreqxy{q}(:,c)',1,2),2,1)]';
    end
end% select current orientation t and place it on grid

% if it's a target absent trial: place distractor at "target position"
c = find(srchdsp.colid == 0);
if strcmp(curt,'ta')
    srchdsp.stimidentity{c} = 'lt';

    srchdsp.set{c} = [l.coord{selor(1)}+repmat(repmat(srchdsp.cxy(:,c)',1,2),2,1)]';
    for q = 1:size(scr.quadrCenter,1)
        srchdsp.sethf{q}{c} = [l.coord{selor(1)}+repmat(repmat(srchdsp.highfreqxy{q}(:,c)',1,2),2,1)]';
    end% select current orientation t and place it on grid
elseif strcmp(curt,'la')
     srchdsp.stimidentity{c} = 'tt';

    srchdsp.set{c} = [t.coord{selor(1)}+repmat(repmat(srchdsp.cxy(:,c)',1,2),2,1)]';
    for q = 1:size(scr.quadrCenter,1)
        srchdsp.sethf{q}{c} = [t.coord{selor(1)}+repmat(repmat(srchdsp.highfreqxy{q}(:,c)',1,2),2,1)]';
    end
end% select current orientation t and place it on grid



