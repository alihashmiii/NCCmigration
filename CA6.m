% based on Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/10/09
% Edited on 07/12/09
% Further developed by Linus Schumacher (LJS) from Oct 2012 onwards.
% Cells move through a 2D box, each cell is represented by an (x,y) coordinate in a width x height y box
% at each timestep the list of cells is worked through in a random order
% and each cell extends filopodia in a set number of random directions
% Domain growth in steps intervals can be added with param.growingDomain=1,
% Chemoattractant concentration can be solved (rather than fixed) with caSolve=1
% Cells move if cellsMove=1,
% Cells are inserted at x=0 if insertCells=1

% Cells have a fixed radius cellRadius and can't move through each other or out of the box
% Uses cell_movement5.m, chemotaxis_solve.m (with associated files),
% domain_growth.m, initiate_cells.m, plot_cells.m,
% requires a subfolder avi_mat to save results

function out = CA6(in,experiment)
global param cells % using global variables is much faster than saving & loading from disk -- LJS
param.experiment = experiment;

if isstruct(in)
    time = in.time;
else
    time = in;
end

tic
%% Model Type Inputs %%
param.growingDomain = 1;     % the domain grows
followerFraction = 1;        % proportion of cells that are followers (0<=follow_per<=1)
% this is only an estimated fraction. actual leader fraction currently
% turns out to be ca. (1 - followerFraction)/2 -- LJS
divide_cells = 0;       % the cells can divide - they divide more where there's more c'tant
conversionType = 0;       % type of conversion used: 0 is no conversion; 1 is time frustrated; 2 is proportion of better directions
numFilopodia = [2,2];  % the number of filopodia for lead cells and follower cells

%%% probably don't want to change these %%%
param.makeChemoattractant = 1;   % there is a chemoattranctant source term
param.zeroBC = 0;                % = 1: make the boundary conditions for the c'tant c(edge) = 0 (has smoothed initial conditions)
                            % else no flux boundary conditions
caSolve = 1;           % solve for the chemoattractant concentration
cellsMove = 1;             % the cells move
insertCells = 1;           % new cells are inserted at x=0

volumeExclusion = 1;    % 1 = cells can't overlap, 0 = they can -- LJS
standStill = 0; % 1 = cells don't move if they don't know where to go; 0 = cells move in a random direction if they don't know where to go

%% Outputs (videos and figures) %%
makeMovies = 0;
makeCaMovie = 0; % makes a movie of a surface plot of the chemo attractant concentration -- LJS
makeAllMovie = 0; % makes a movie of the cells with filopodia on top of a contourplot of the chemoattractant -- LJS
makeFrames = 0; % makes frames at 0, 12 and 24 hours (can be changed) of the cells on top of the ca -- LJS

%% General parameters %%
param.tstep = 1/60;                   % time step in hours
numTsteps = floor(time/param.tstep)+1;   % number of time steps
t_0 = 6; % time in hrs at which simulation starts
cellRadius = 7.5;              % radius in um (= 7.5um)
leadSpeed = 41.6;                     % speed of the leader cells in mu/h
followSpeed = 49.9;                 % speed of the follower cells in mu/h

param.domainHeight = 120;                   % maximum y value
filolength = cellRadius + 10*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
maxFilolength = 45; % maximum length of filopodium before follower dettaches from leader. default = filolength (non-extensible)

dist = [leadSpeed; followSpeed]*param.tstep;             % the distance moved in a timestep
sensingAccuracy = 0.01; % relative accuracy with which concentration can be measurem. dC/C has to be greater than this to be noticed. This is the baseline value for the starting concentration, scales with 1/sqrt(c) -- LJS

needNeighbours = 0; % cells only move (directed) if there are at least this many other cells within filolength -- LJS
% set direction of movement 'parallel' or 'toward' to that of cell being
contactGuidance = 'parallel'; 
%% experimental parameters %%
param.insert = 0;                     % signal that the chemoattractant has been inserted (for experiment 1)
switch param.experiment
    case 12 %VEGF transplant back half
        param.transplantTime = 12; % time at which CA-production will be locally increased
        param.transplantXLocation = 0; % left edge of square region in which CA-production will be increased
        param.secondaryChi = 5; % strength of increased CA production
    case 13 %VEGF transplant middle half
        param.transplantTime = 12; % time at which CA-production will be locally increased
        param.transplantXLocation = 70; % left edge of square region in which CA-production will be increased
        param.secondaryChi = 5; % strength of increased CA production
    case 11 %VEGF transplant back edge
        param.transplantTime = 12; % time at which CA-production will be locally increased
        param.transplantXLocation = 0; % left edge of square region in which CA-production will be increased
        param.secondaryChi = 5; % strength of increased CA production
    case 14 %increased VEGF at far (right-most) edge
        param.transplantTime = 12; % time at which CA-production will be locally increased
        param.transplantXLocation = 425; % left edge of square region in which CA-production will be increased
        param.secondaryChi = 5; % strength of increased CA production
    otherwise
        param.transplantTime = NaN; param.transplantXLocation = NaN; param.secondaryChi = NaN;
end
%% caSolve parameters %%
param.diffus = 0.1;%252e3;    % chemoattractant diffusivity (in (mu)^2/h), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
param.chi = 0.0001;                  % chemoattractant production term (usually 0.0001)
param.eatRate = 1000;                      % chemoattractant consumption rate
param.eatWidth = cellRadius;         % width of eating chemoattractant, equivalent to gaussian sigma

%% convert parameters
if isstruct(in)
    if ismember('conversionType',fields(in))
        conversionType = in.conversionType;
    end
end
if conversionType == 1
    param.numSteps = 5; % number of steps to not sucessfully find a direction, before changing roles (convert type 1)
    param.numDirections=NaN;
elseif conversionType == 2
    param.numSteps = numFilopodia(1); % number of directions to sample in (convert type 2) -- this is currently set in convert_cells.m
    param.numDirections = 1/param.numSteps; % fraction of directions needed to be better to maintain a leader profile (convert type 2) -- this is currently set in convert_cells.m
elseif conversionType == 4
    param.numSteps = [8, 8]; % timescale in minutes for switching [lead2follow, follow2lead]
    param.numDirections = NaN;
else
    param.numSteps=10;
    param.numDirections=NaN;
end
if isstruct(in)
    if ismember('numSteps',fields(in))
        param.numSteps = in.numSteps;
    end
    if ismember('numDirections',fields(in))
        param.numDirections = in.numDirections;
    end
end
%% insertCells parameters %%
insertTimeStep = 0.1;         % the time in hours between each new cell insertion
insertEverySteps = floor(insertTimeStep/param.tstep);    % how often are new cells inserted
insertNumCells = 1;                  % how many new cells are inserted at each timepoint
if insertCells==1
    numCellsInitial = 6;              % initial number of cells
else
    numCellsInitial = 1;
end
initYFrac = (param.domainHeight-2*cellRadius)/param.domainHeight; % fraction of y initiated with cells (so that they aren't too close to the top or bottom)
initXFrac = 0;                 % initial fraction of x with cells

%% adjust parameters if they have been provided in input %%
if isstruct(in)
    if ismember('leadSpeed',fields(in))
        leadSpeed = in.leadSpeed; % speed of the leader cells in mu/h
        dist = [leadSpeed; followSpeed]*param.tstep;             % the distance moved in a timestep
    end
    if ismember('followSpeed',fields(in))
        followSpeed = in.followSpeed; % speed of the follower cells in mu/h
        dist = [leadSpeed; followSpeed]*param.tstep;             % the distance moved in a timestep
    end
    if ismember('numFilopodia',fields(in))
        numFilopodia = in.numFilopodia; % the number of filopodia for lead cells and follower cells
    end
    if ismember('filolength',fields(in))
        filolength = in.filolength; % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
    end
    if ismember('maxFilolength',fields(in))
        maxFilolength = in.maxFilolength; % maximum length of filopodium before folloer dettaches from leader
    end
    if ismember('diffus',fields(in))
        param.diffus = in.diffus; % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
    end
    if ismember('chi',fields(in))
        param.chi = in.chi; % chemoattractant production term (usually 0.0001)
    end
    if ismember('eatRate',fields(in))
        param.eatRate = in.eatRate; % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
    end
    if ismember('eatWidth',fields(in))
        param.eatWidth = in.eatWidth; % width of eating chemoattractant, equivalent to gaussian sigma
    end
    if ismember('followerFraction',fields(in))
        followerFraction = in.followerFraction; % proportion of cells that are followers (0<=follow_per<=1)
    end
    if ismember('tstep',fields(in))
        param.tstep = in.tstep; % time step in hours
        numTsteps = floor(time/param.tstep)+1;   % number of time steps
        dist = [leadSpeed; followSpeed]*param.tstep;             % the distance moved in a timestep
    end
    if ismember('volumeExclusion',fields(in))
        volumeExclusion = in.volumeExclusion;
    end
    if ismember('standStill',fields(in))
        standStill = in.standStill;
    end
    if ismember('insertEverySteps',fields(in))
        insertEverySteps = in.insertEverySteps;
    end
    if ismember('sensingAccuracy',fields(in))
        sensingAccuracy = in.sensingAccuracy;
    end
    if ismember('needNeighbours',fields(in))
        needNeighbours = in.needNeighbours;
    end
    if ismember('makeMovies',fields(in))
        makeMovies = in.makeMovies;
    end
    if ismember('makeCaMovie',fields(in))
        makeCaMovie = in.makeCaMovie;
    end
    if ismember('makeAllMovie',fields(in))
        makeAllMovie = in.makeAllMovie;
    end
    if ismember('makeFrames',fields(in))
        makeFrames = in.makeFrames;
    end
    if ismember('numCellsInitial',fields(in))
        numCellsInitial = in.numCellsInitial;
    end
    if ismember('caSolve',fields(in))
        caSolve = in.caSolve;
    end
    if ismember('divide_cells',fields(in))
        divide_cells = in.divide_cells;
    end
    if ismember('insertCells',fields(in))
        insertCells = in.insertCells;
    end
    if ismember('t_0',fields(in))
        t_0 = in.t_0;
    end
    if ismember('domainHeight',fields(in))
        param.domainHeight = in.domainHeight;
    end
    if ismember('initYFrac',fields(in))
        initYFrac = in.initYFrac;
    end
    if ismember('contactGuidance',fields(in))
        contactGuidance = in.contactGuidance; 
    end
end

followStart = floor(18/param.tstep) - floor(18*followerFraction/param.tstep)+1 % the time step after which new cells will be followers, to aim for the desired fraction of followers at t = 18hours -- LJS

%% domain growth parameters %%
param.initialDomainLength= 300;                % initial width of the domain with growth (um)
domainLengths = ones(1,numTsteps).*param.initialDomainLength;  % initialise domain length vector

% These parameters are found from data from the Kulesa Lab, using least
% squares regression in make_domain_plot.m on data in lengths_data.m
param.Linf = 870;                            % end domain length
param.a = 0.0800;                            % "steepness" of the logistic domain growth -- LJS
param.t_s = 16; %parameter used in domain_growth, giving the point of highest growth
%%
presave_stuff % create results file with parameters to signal that this simulation is being worked on

%% initialise cells
if isstruct(in)&&ismember('cells',fields(in))
    cells = in.cells; % take cells passed in
    numCellsInitial = size(cells,2);
else % set up the initial cells so that they aren't too close to each other or the edge
    temp = initiate_cells(numCellsInitial,cellRadius,0,param.initialDomainLength,param.domainHeight,initXFrac,initYFrac,[],1);
    cells = temp.cells;
end

finalNumCells = (numCellsInitial + floor(numTsteps/insertEverySteps)*insertNumCells)*2;    % final number of cells expected (*2 for divisions and experimental insertions)

if followerFraction > 1
        cellsFollow = true(finalNumCells,1); % cells are followeres by default
    else
        cellsFollow = false(finalNumCells,1); % cells are leaders by default.
end    
if isstruct(in)&&ismember('cellsFollow',fields(in))
        cellsFollow(1:numCellsInitial) = in.cellsFollow(1:numCellsInitial); % take cell states passed in
end
% cells being inserted after a certain time-point will be set to followers.
% This is a better approximation of leader fraction than pre-setting based
% on expected total cell numbers, which are too high when many cells cannot
% be inserted due to jamming -- LJS

%% initialise vectors and time %%
t_save = t_0+(0:param.tstep:param.tstep*numTsteps);
xlat_save = cell(1,numTsteps); % spatial lattices (lat)
ylat_save = cell(1,numTsteps);
ca_save = cell(1,numTsteps); % chemoattractant (ca)
cells_save = cell(numTsteps,1);
filopodia_save = cell(numTsteps,1);
cellsFollow_save = cell(numTsteps,1);
attach = zeros(finalNumCells,1,'uint16'); % indices of which cell each cell is attached to
if isstruct(in)&&ismember('attach',fields(in))
    attach(1:numCellsInitial) = in.attach(1:numCellsInitial); % take cell attachments passed in
end
theta = NaN(finalNumCells,1); % cells' movement directions-- LJS
attach_save = cell(1,numTsteps);
if isstruct(in)&&ismember('ca_new',fields(in))
    param.ca_new = in.ca_new; % take CA field passed in
    param.insert = 1;
end
if isstruct(in)&&ismember('xlat_new',fields(in));
    xlat_new = in.xlat_new; % take lattice passed in - not sure how important this is
    param.insert = 1;
else
    xlat_new=[];
end
moved = false(numTsteps,finalNumCells);
happiness = NaN(numTsteps,finalNumCells); % for integrate-and-switch cell behaviour conversion
%% begin timesteps %%
for timeCtr=1:numTsteps
    %% after t=followStart then all subsequent cells are followers
    if timeCtr==followStart
        if param.experiment==35 % makes only half of trailers followers, half leaders
            cellsFollow(length(cells(1,:))+1:2:end) = 1;
        else
            cellsFollow(length(cells(1,:))+1:end) = 1;
        end
        disp('followers start here')
    end
    
    %% Cells insertions for experiments
    transplant_cells
    
    %% If we are inserting new cells, do so here %%
    if mod(timeCtr,insertEverySteps)==0
        fprintf(['t = ' num2str(t_save(timeCtr+1)) '\r'] )
        if (insertCells==1)&&((param.experiment==0)||(param.experiment>3)||(in.it==1)||(in.ablate_type~=2)||t_save(timeCtr)<in.ablate_time)
            temp = initiate_cells(insertNumCells,cellRadius,0,param.initialDomainLength,param.domainHeight,0,initYFrac,cells,volumeExclusion);
            cells = temp.cells;
        end
    end
    
    %% chemoattractant %%
    if caSolve==1
        % give parameters for the solver (depending on whether this is the first run or not)
        if timeCtr==1
            ind = int64(0); % starts integration at t=0
            iwk = zeros(580230,1,'int64'); % is used by the solver for outputting the efficiency of integration, check documentation at 5.4-4: http://www.nag.co.uk/numeric/MB/manual_21_1/pdf/D03/d03ra.pdf#lnk_leniwk -- LJS
            rwk = zeros(1880000,1); % it's unclear from NAG documentation what this parameter is used for, but it needs to be a double array of a certain size -- LJS
        elseif ((param.experiment==1)||(param.experiment==2))&&(in.it==2)&&(t_save(timeCtr)==in.changeTime)
            %% experiments 1 and 2: inserting chemoattractant % backward compatibility not test -- LJS
            insert_tissue
        else
            ind = int64(1); % continuing integration from the previous solution
        end
        % run the solver
        if ((param.experiment==1)||(param.experiment==2))&&(in.it==2)&&(t_save(timeCtr)==in.changeTime)
            temp = chemotaxis_solve(t_save(timeCtr),t_save(timeCtr+1),ind,iwk,rwk,param.initialDomainLength,param.domainHeight,xlat_new,length(ylat_new),param.insert);            
        elseif timeCtr>1
            temp = chemotaxis_solve(t_save(timeCtr),t_save(timeCtr+1),ind,iwk,rwk,param.initialDomainLength,param.domainHeight,xlat_new,length(ylat_save{timeCtr-1}),param.insert);
        else % initialise chemoattractant
            temp = chemotaxis_solve(t_save(timeCtr),t_save(timeCtr+1),ind,iwk,rwk,param.initialDomainLength,param.domainHeight,xlat_new,32,param.insert);
        end
        if temp.ifail~=0
            save(['results/' saveInfo '_' num2str(timeCtr) '_solverWarningLog.mat'],'temp')
        end
        % take output 
        xlat_save{timeCtr} = temp.xsave;
        ylat_save{timeCtr} = temp.ysave;
        iwk = temp.iwk;
        rwk = temp.rwk;
        
        ca_save{timeCtr} = temp.chemotaxis;
    else
        %% Fixed chemoattractant %%
        if param.growingDomain==1
            % Domain Growth Happens at every timestep
            % and starts 6 hours before migration -- LJS
            [~, domainLengths(timeCtr), ~] = domain_growth(cells(1,:),t_save(timeCtr),param.tstep,param.Linf,param.a,param.initialDomainLength,param.t_s);
        end
        
        xlat_save{timeCtr} = 0:domainLengths(timeCtr)/100:domainLengths(timeCtr);
        ylat_save{timeCtr} = 0:param.domainHeight/100:param.domainHeight;
        ca = zeros(length(xlat_save{timeCtr}),length(ylat_save{timeCtr}));
        for i=1:length(xlat_save{timeCtr})
            for j=1:length(ylat_save{timeCtr})
                ca(i,j) = xlat_save{timeCtr}(i)./domainLengths(timeCtr);
            end
        end
        ca_save{timeCtr} = ca;
    end
    %% divide cells %%
    if (divide_cells==1)
        temp = cells_divide(cellsFollow,cellRadius,domainLengths(timeCtr),0,param.domainHeight,ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},param.tstep);
        cellsFollow = temp.cellsFollow;
    end
    %% domain growth %%
    if param.growingDomain==1
        % Domain Growth Happens at every timestep
        % and starts 6 hours before migration -- LJS
        [cells(1,:), domainLengths(timeCtr), ~] = domain_growth(cells(1,:),t_save(timeCtr),param.tstep,param.Linf,param.a,param.initialDomainLength,param.t_s);
    end
    
    %% move cells %%
    if cellsMove==1
        if timeCtr==1
            temp = new_move_cells(cellsFollow,[],attach,theta,...
                ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},...
                cellRadius,filolength,maxFilolength,param.eatWidth,param.domainHeight,dist,domainLengths(timeCtr),numFilopodia,...
                volumeExclusion, standStill,sensingAccuracy,needNeighbours,contactGuidance);
        else
            temp = new_move_cells(cellsFollow,filopodia,attach,theta,...
                ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},...
                cellRadius,filolength,maxFilolength,param.eatWidth,param.domainHeight,dist,domainLengths(timeCtr),numFilopodia,...
                volumeExclusion, standStill,sensingAccuracy,needNeighbours,contactGuidance);
        end
        attach = temp.attach;
        cellsFollow = temp.cellsFollow;
        filopodia = temp.filopodia;
        theta = temp.theta;
        moved(timeCtr,:) = [temp.moved, false(1,length(moved(1,:))-length(temp.moved))]; % with padding for not-yet-existing cells -- LJS
        if conversionType==4
            if timeCtr ==1
                happiness(timeCtr,1:length(cells(1,:))) = ~cellsFollow(1:length(cells(1,:))); % leaders start at happiness 1, followers at zero (their respective switching thresholds, i.e. max /min) -- LJS
                if isstruct(in)&&ismember('happiness',fields(in))
                    happiness(1,1:numCellsInitial) = in.happiness; % take cell happiness passed in
                end
            else
                newCellIdcs = isnan(happiness(timeCtr - 1,1:length(cells(1,:)))); % newly existing cells have previous happiness nan -- LJS
                happiness(timeCtr,newCellIdcs) = ~cellsFollow(newCellIdcs); % leaders start at happiness 1, followers at zero (their respective switching thresholds, i.e. max /min) -- LJS
                happiness(timeCtr,temp.sensed&~newCellIdcs) = min(1,... % the maximum happiness is 1 -- LJS
                    happiness(timeCtr-1,temp.sensed&~newCellIdcs) + param.tstep*60/param.numSteps(2)); % cells that sensed CA become happier, with maximum 1 -- LJS
                happiness(timeCtr,~temp.sensed&~newCellIdcs) = max(0,... % minimum happiness is 0 -- LJS
                    happiness(timeCtr-1,~temp.sensed&~newCellIdcs) - param.tstep*60/param.numSteps(1)); % cells that haven't sensed CA become sadder, with minimum 0 -- LJS
                % new cells change from current happiness, as previous is
                % nan
                happiness(timeCtr,temp.sensed&newCellIdcs) = min(1,... % the maximum happiness is 1 -- LJS
                    happiness(timeCtr,temp.sensed&newCellIdcs)+ param.tstep*60/param.numSteps(2)); % cells that sensed CA become happier, with maximum 1 -- LJS
                happiness(timeCtr,~temp.sensed&newCellIdcs) = max(0,... % minimum happiness is 0 -- LJS
                    happiness(timeCtr,~temp.sensed&newCellIdcs) - param.tstep*60/param.numSteps(1)); % cells that haven't sensed CA become sadder, with minimum 0 -- LJS
            end
        end
        attach_save{timeCtr} = attach(1:size(cells,2));
        cellsFollow_save{timeCtr} = cellsFollow(1:size(cells,2));
        filopodia_save{timeCtr} = filopodia;
    end
    cells_save{timeCtr}=cells;
    
    %% cells can convert from leaders <-> followers
    if (conversionType~=0)
        out = convert_cells(cellsFollow,timeCtr,cells_save,filolength,moved,happiness,ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},...
            param.eatWidth,conversionType,numFilopodia);
        cellsFollow = out.cellsFollow;
        moved = out.moved;
        cellsFollow_save{timeCtr} = cellsFollow(1:size(cells,2));
    end
end

toc
numCellsFinal = size(cells,2);
disp(['Number of cells: ' num2str(numCellsFinal)])
%% Save stuff
save_stuff

%% make movies %%
% make_figure

if makeMovies==1    
    caCmap = load('cmap_blue2cyan.txt');
    %%% make frames %%%
    if makeFrames==1
        make_frames
        disp('made frames')
%         open(['avi_mat/frames/frames3',saveInfo,'.fig'])
    end
    
    %%% make camovie.avi %%%
    if makeCaMovie==1
        make_ca_movie
    end
    
    %%% make cells+ca movie (allmovie.avi)%%%
    if makeAllMovie==1
        make_all_movie_hidden
    end
    
    if makeFrames==1
        open(['avi_mat/frames/',saveInfo,'.fig'])
    end
end

delete('avi_mat/*.mat');
