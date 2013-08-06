% Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/10/09
% Edited on 07/12/09
% Further developed by Linus Schumacher (LJS) from Oct 2012 onwards.
% Cells move through a 2D box, each cell is represented by an (x,y) coordinate in a width x height y box
% at each timestep the list of cells is worked through in a random order
% and each cell extends filopodia in a set number of random directions
% Domain growth in steps intervals can be added with growingDomain=1,
% Chemoattractant concentration can be solved (rather than fixed) with ca_solve=1
% Cells move if cells_move=1,
% Cells are inserted at x=0 if insert_cells=1

% Cells have a fixed radius cellRadius and can't move through each other or out of the box
% Uses cell_movement5.m, chemotaxis_solve.m (with associated files),
% domain_growth.m, initiate_cells.m, plot_cells.m,
% requires a subfolder avi_mat to save results

function out = CA6(in,experiment)
global param cells % using global variables is much faster than saving & loading from disk -- LJS

if isstruct(in)
    time = in.time;
else
    time = in;
end

tic
%% Model Type Inputs %%
growingDomain = 1;     % the domain grows
followerFraction = 0.7;        % proportion of cells that are followers (0<=follow_per<=1)
divide_cells = 0;       % the cells can divide - they divide more where there's more c'tant
convert_type = 0;       % type of conversion used: 0 is no conversion; 1 is time frustrated; 2 is proportion of better directions
numFilopodia = [6,2];  % the number of filopodia for lead cells and follower cells

%%% probably don't want to change these %%%
makeChemoattractant = 1;   % there is a chemoattranctant source term
zeroBC = 0;                % = 1: make the boundary conditions for the c'tant c(edge) = 0 (has smoothed initial conditions)
                            % else no flux boundary conditions
ca_solve = 1;           % solve for the chemoattractant concentration
cells_move = 1;             % the cells move
insert_cells = 1;           % new cells are inserted at x=0

%% Outputs (videos and figures) %%
movies = 0;
ca_movie = 0; % makes a movie of a surface plot of the chemo attractant concentration -- LJS
all_movie = 0; % makes a movie of the cells with filopodia on top of a contourplot of the chemoattractant -- LJS
frames = 0; % makes frames at 0, 12 and 24 hours (can be changed) of the cells on top of the ca -- LJS

%% General parameters %%
tstep = 5/2/60;                   % time step in hours
numTsteps = floor(time/tstep)+1;   % number of time steps
cellRadius = 7.5;              % radius in um (= 7.5um)
leadSpeed = 41.6;                     % speed of the leader cells in mu/h
followSpeed = 49.9;                 % speed of the follower cells in mu/h

domainHeight = 120;                   % maximum y value
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
inity_frac = (domainHeight-2*cellRadius)/domainHeight; % fraction of y initiated with cells (so that they aren't too close to the top or bottom)
dist = [leadSpeed; followSpeed]*tstep;             % the distance moved in a timestep

%% experimental parameters %%
insert = 0;                     % signal that the chemoattractant has been inserted (for experiment 1)

%% ca_solve parameters %%
diffus = 1;%252e3;    % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
chi = 0.0001;                  % chemoattractant production term (usually 0.0001)
eatRate = 1;                      % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
eatWidth = cellRadius;         % width of eating chemoattractant, equivalent to gaussian sigma

%% adjust parameters if they have been provided in input %%
if isstruct(in)
    if ismember('leadSpeed',fields(in))
        leadSpeed = in.leadSpeed; % speed of the leader cells in mu/h
    end
    if ismember('followSpeed',fields(in))
        followSpeed = in.followSpeed; % speed of the follower cells in mu/h
    end
    if ismember('numFilopodia',fields(in))
        numFilopodia = in.numFilopodia; % the number of filopodia for lead cells and follower cells
    end
    if ismember('filolength',fields(in))
        filolength = in.filolength; % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
    end
    if ismember('diffus',fields(in))
        diffus = in.diffus; % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
    end
    if ismember('chi',fields(in))
        chi = in.chi; % chemoattractant production term (usually 0.0001)
    end
    if ismember('eatRate',fields(in))
        eatRate = in.eatRate; % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
    end
    if ismember('eatWidth',fields(in))
        eatWidth = in.eatWidth; % width of eating chemoattractant, equivalent to gaussian sigma
    end
    if ismember('followerFraction',fields(in))
        followerFraction = in.followerFraction; % proportion of cells that are followers (0<=follow_per<=1)
    end
end
%% convert parameters
if convert_type == 1
    num_steps = 10; % number of steps to not sucessfully find a direction, before changing roles (convert type 1)
    num_directions=[];
elseif convert_type == 2
    num_steps = numFilopodia(1); % number of directions to sample in (convert type 2)
    num_directions = 1/num_steps; % fraction of directions needed to be better to maintain a leader profile (convert type 2)
else
    num_steps=[];
    num_directions=[];
end
    
%% insert_cells parameters %%
insert_time_step = 0.1;         % the time in hours between each new cell insertion
insert_step = floor(insert_time_step/tstep);    % how often are new cells inserted
num_cells = 1;                  % how many new cells are inserted at each timepoint
if insert_cells==1
    n = 6;              % initial number of cells
else
    n = 1;
end
initx_frac = 0;                 % initial fraction of x with cells
followStart = floor(18/tstep) - floor(18*followerFraction/tstep)+1 % the time step after which new cells will be followers, to aim for the desired fraction of followers at t = 18hours -- LJS

%% domain growth parameters %%
initialDomainLength= 300;                % initial width of the domain with growth (um)
domainLengths = ones(1,numTsteps).*initialDomainLength;  % initialise domain length vector

% These parameters are found from data from the Kulesa Lab, using least
% squares regression in make_domain_plot.m on data in lengths_data.m
Linf = 870;                            % end domain length
a = 0.0800;                            % "steepness" of the logistic domain growth -- LJS
t_start = -16; %parameter used in domain growth
param = [Linf, a, diffus, eatWidth, growingDomain, initialDomainLength, makeChemoattractant,...
    chi, domainHeight, zeroBC, insert, tstep, t_start, eatRate, num_steps, num_directions];
save avi_mat/param param

%% set up the initial cells so that they aren't too close to each other or
%% the edge %%
temp = initiate_cells(n,cellRadius,0,initialDomainLength,domainHeight,initx_frac,inity_frac,[]);
cells = temp.cells;
end_num_cells = (n + floor(numTsteps/insert_step)*num_cells)*2;    % final number of cells expected (*2 for divisions and experimental insertions)
cellsFollow = false(end_num_cells,1); % cells are leaders by default. For fixed fractions of followers, all
% cells being inserted after a certain time-point will be set to followers.
% This is a better approximation of leader fraction than pre-setting based
% on expected total cell numbers, which are too high when many cells cannot
% be inserted due to jamming -- LJS

%% initialise vectors and time %%
t_save = 0:tstep:tstep*numTsteps;
xlat_save = cell(1,numTsteps); % spatial lattices (lat)
ylat_save = cell(1,numTsteps);
ca_save = cell(1,numTsteps); % chemoattractant (ca)
cells_save = cell(numTsteps,1);
filopodia_save = cell(numTsteps,1);
cellsFollow_save = cell(numTsteps,1);
attach = zeros(end_num_cells,1,'uint16'); % indices of which cell each cell is attached to
theta = NaN(end_num_cells,1); % cells' movement directions-- LJS
attach_save = cell(1,numTsteps);
xlat_new=[];
barrier = zeros(numTsteps,1); % x-position of barrier for experiments 4&5, may be obsolete -- LJS
moved = false(numTsteps,end_num_cells);
num_better_foll_save = []; %% these may be obsolete -- LJS
num_foll_save = 0; %% these may be obsolete -- LJS
num_better_lead_save = []; %% these may be obsolete -- LJS
num_lead_save = 0; %% these may be obsolete -- LJS
%% begin timesteps %%
for timeCtr=1:numTsteps
    %% after t=followStart then all subsequent cells are followers
    if timeCtr==followStart
        cellsFollow(length(cells(1,:))+1:end) = 1;
        disp('followers start here')
    end
    
    %% Cells and barrier insertions for experiments
    insert_cells_and_barrier
    
    %% If we are inserting new cells, do so here %%
    if mod(timeCtr,insert_step)==0
        fprintf(['t = ' mat2str(t_save(timeCtr+1)) '\r'] )
        if (insert_cells==1)&&((experiment==0)||(experiment>3)||(in.it==1)||(in.ablate_type~=2)||t_save(timeCtr)<in.ablate_time)
            temp = initiate_cells(num_cells,cellRadius,0,initialDomainLength,domainHeight,0,inity_frac,cells);
            cells = temp.cells;
        end
    end
    
    %% chemoattractant %%
    if ca_solve==1
        cells_in = cells;
        % give parameters for the solver (depending on whether this is the first run or not)
        if t_save(timeCtr)==0
            ind = int64(0); % starts integration at t=0
            iwk = zeros(580230,1,'int64'); % is used by the solver for outputting the efficiency of integration, check documentation at 5.4-4: http://www.nag.co.uk/numeric/MB/manual_21_1/pdf/D03/d03ra.pdf#lnk_leniwk -- LJS
            rwk = zeros(1880000,1); % it's unclear from NAG documentation what this parameter is used for, but it needs to be a double array of a certain size -- LJS
        elseif ((experiment==1)||(experiment==2))&&(in.it==2)&&(t_save(timeCtr)==in.changeTime)
            %% experiments 1 and 2: inserting chemoattractant
            insert_tissue
        else
            ind = int64(1); % continuing integration from the previous solution
        end
        % run the solver
        if ((experiment==1)||(experiment==2))&&(in.it==2)&&(t_save(timeCtr)==in.changeTime)
            temp = chemotaxis_solve(t_save(timeCtr),t_save(timeCtr+1),ind,iwk,rwk,cells_in,initialDomainLength,domainHeight,length(xlat_new),length(ylat_new),insert);
        elseif timeCtr>1
            temp = chemotaxis_solve(t_save(timeCtr),t_save(timeCtr+1),ind,iwk,rwk,cells_in,initialDomainLength,domainHeight,length(xlat_new),length(ylat_save{timeCtr-1}),insert);
        else
            temp = chemotaxis_solve(t_save(timeCtr),t_save(timeCtr+1),ind,iwk,rwk,cells_in,initialDomainLength,domainHeight,length(xlat_new),50,insert);
        end

        % take output 
        xlat_save{timeCtr} = temp.xsave;
        ylat_save{timeCtr} = temp.ysave;
        iwk = temp.iwk;
        rwk = temp.rwk;
        
        ca_save{timeCtr} = temp.chemotaxis;
    else
        %% Fixed chemoattractant %%
        if growingDomain==1
            % Domain Growth Happens at every timestep
            if ((experiment==4)||(experiment==5))
                temp = domain_growth(cells(1,:),t_save(timeCtr),tstep,Linf,a,initialDomainLength,barrier(timeCtr),t_start);
                barrier(timeCtr+1)=temp.barrier;
            else
                temp = domain_growth(cells(1,:),t_save(timeCtr),tstep,Linf,a,initialDomainLength,[],t_start);
            end
            domainLengths(timeCtr) = temp.domainLength;
        end
        
        xlat_save{timeCtr} = 0:domainLengths(timeCtr)/100:domainLengths(timeCtr);
        ylat_save{timeCtr} = 0:domainHeight/100:domainHeight;
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
        temp = cells_divide(cells,cellsFollow,cellRadius,domainLengths(timeCtr),0,domainHeight,ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},tstep);
        cells = temp.cells;
        cellsFollow = temp.cellsFollow;
    end
    %% domain growth %%
    if growingDomain==1
        % Domain Growth Happens at every timestep
        if ((experiment==4)||(experiment==5))
            temp = domain_growth(cells(1,:),t_save(timeCtr),tstep,Linf,a,initialDomainLength,barrier(timeCtr),t_start);
            barrier(timeCtr+1)=temp.barrier;
        else
            temp = domain_growth(cells(1,:),t_save(timeCtr),tstep,Linf,a,initialDomainLength,[],t_start);
        end
        domainLengths(timeCtr) = temp.domainLength;
        cells(1,:) = temp.cells_next;
        cells(2,:) = cells(2,:);
    end
    
    %% move cells %%
    if cells_move==1
        if timeCtr==1
            temp = new_move_cells(cells,cellsFollow,[],attach,theta,...
                ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},...
                cellRadius,filolength,eatWidth,domainHeight,dist,domainLengths(timeCtr),barrier(timeCtr),experiment,t_save(timeCtr),in,numFilopodia);
        else
            temp = new_move_cells(cells,cellsFollow,filopodia,attach,theta,...
                ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},...
                cellRadius,filolength,eatWidth,domainHeight,dist,domainLengths(timeCtr),barrier(timeCtr),experiment,t_save(timeCtr),in,numFilopodia);
        end
        attach = temp.attach;
        cellsFollow = temp.cellsFollow;
        filopodia = temp.filopodia;
        theta = temp.theta;
        cells = temp.cells;
        moved(timeCtr,:) = [temp.moved, false(1,length(moved(1,:))-length(temp.moved))];
        
        attach_save{timeCtr} = attach;
        cellsFollow_save{timeCtr} = cellsFollow;
        filopodia_save{timeCtr} = filopodia;
    end
    cells_save{timeCtr}=cells;
    
    %% cells can convert from leaders <-> followers
    if (convert_type~=0)&&((experiment==0)||experiment==3||(in.it==1)||(t_save(timeCtr)==in.changeTime))
        out = convert_cells(cells,cellsFollow,attach_save,timeCtr,cells_save,filolength,moved,ca_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},...
            eatWidth,filopodia,convert_type,param,num_better_foll_save,num_foll_save,num_better_lead_save,num_lead_save);
        cellsFollow = out.cellsFollow;
        moved = out.moved;
        num_better_foll_save = out.num_better_foll_save;
        num_foll_save = out.num_foll_save;
        num_better_lead_save = out.num_better_lead_save;
        num_lead_save = out.num_lead_save;
        cellsFollow_save{timeCtr} = cellsFollow;
    end
    

    if (cells~=cells_save{timeCtr})
        size(cells)
        size(cells_save)
        disp('cells error')
        pause
    end
    if (timeCtr>=followStart)
        if (attach~=attach_save{timeCtr})
            disp('attach error')
            pause
        end
        if (cellsFollow~=cellsFollow_save{timeCtr})
            disp('follow error')
            pause
        end
        if (filopodia~=filopodia_save{timeCtr})
            disp('filopodia error')
            pause
        end
    end
end
%%
% hist(num_better_foll_save/num_foll_save)
% hist(num_better_lead_save/num_lead_save)
out.num_better_foll_save = num_better_foll_save;
out.num_foll_save = num_foll_save;
out.num_better_lead_save = num_better_lead_save;
out.num_lead_save = num_lead_save;
toc
% whitebg('white')
disp(['Number of cells: ' num2str(size(cells,2))])
%% Save stuff
save_stuff

%% make movies %%
% make_figure

if movies==1    
    %%% make frames %%%
    if frames==1
        make_frames
        disp('made frames')
    end
    
    %%% make camovie.avi %%%
    if ca_movie==1
        make_ca_movie
    end
    
    %%% make cells+ca movie (allmovie.avi)%%%
    if all_movie==1
        make_all_movie
    end
    
    close all
    open(['avi_mat/frames/frames3',saveInfo,'.fig'])
end

delete('avi_mat/*.mat');
