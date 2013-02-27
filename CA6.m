% Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/10/09
% Edited on 07/12/09
% Cells move through a 2D box, each cell is represented by an (x,y) coordinate in a width x height box
% at each timestep the list of cells is worked through in a random order
% and each cell extends filopodia in random directions until a better
% concentration is found or max_fil is reached
% Domain growth in steps intervals can be added with growing_domain=1,
% Chemoattractant concentration can be solved (rather than fixed) with ca_solve=1
% Cells move if cells_move=1,
% Cells are inserted at x=0 if insert_cells=1

% Cells have a fixed radius cell_radius and can't move through each other or out of the box
% Uses cell_movement5.m, chemotaxis_solve.m (with associated files),
% domain_growth.m, initiate_cells.m, plot_cells.m,
% requires a subfolder avi_mat to save results

function out = CA6(in,experiment)

if (experiment~=0)&&(experiment~=6)
    time = in.time;
else
    time = in;
end

tic
%% Model Type Inputs %%
growing_domain = 1;     % the domain grows
follow_perc = 7/8;        % proportion of cells that are followers (0<=follow_per<=1)
divide_cells = 0;       % the cells can divide - they divide more where there's more c'tant
convert_type = 0;       % type of conversion used: 0 is no conversion; 1 is time frustrated; 2 is proportion of better directions
metropolis = 0;         % cells sometimes move even if it is unfavourable
num_filopodia = [6,2];  % the number of filopodia for lead cells and follower cells
satiate = 0;            % number of timesteps after which the cells are satiated

%%% probably don't want to change these %%%
make_chemoattractant = 1;   % there is a chemoattranctant source term
zero_bc = 0;                % = 1: make the boundary conditions for the c'tant c(edge) = 0 (has smoothed initial conditions)
                            % else no flux boundary conditions
if experiment==6
    ca_solve = 0;           % don't bother solving for the chemoattractant concentration
else
    ca_solve = 1;           % solve for the chemoattractant concentration
end
% ca_solve = 0;
cells_move = 1;             % the cells move
insert_cells = 1;           % new cells are inserted at x=0

%% Outputs (videos and figures) %%
movies = 1;
ca_movie = 0; % makes a movie of a surface plot of the chemo attractant concentration -- LJS
all_movie = 1; % makes a movie of the cells with filopodia on top of a contourplot of the chemoattractant -- LJS
frames = 1; % makes frames at 0, 12 and 24 hours (can be changed) of the cells on top of the ca -- LJS

%% General parameters %%
tstep = 0.05;                   % time step in hours
tsteps = floor(time/tstep)+1;   % number of time steps
cell_radius = 7.5;              % radius in um (= 7.5um)
speed = 45;                     % speed of the leader cells (=2*45um/hr because half the time they look the wrong way and don't move)
height = 120;                   % maximum y value
filolength = cell_radius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS)
% filolength = 50;   % for trouble-shooting, this is what Louise used -- LJS
inity_perc = (height-2*cell_radius)/height; % percentage of y initiated with cells (so that they aren't too close to the top or bottom)
dist = speed*tstep;             % the distance moved in a timestep

%% diffusion parameters
if experiment==6
    dist = 1;                   % the distance moved in a timestep
    diffusion = speed / dist * tstep;   % rate of movement given speed, distance moved and timestep
else
    diffusion = 0;
end

%% experimental parameters %%
insert = 0;                     % signal that the chemoattractant has been inserted (for experiment 1)

%% ca_solve parameters %%
diffus = 0.02888;               % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s, for membrane bound VEGF unknown/near zero -- LJS
d = 30;                         % width of sensing chemoattractant
chi = 0.0001;                   % chemoattractant production term (usually 0.0001)
mult = 1;                      % lambda, usually 0.045

%% cells_move parameters %%
e = d;                          % width of eating chemoattractant (smaller e has larger radius)

%% convert parameters
if convert_type == 1
    num_steps = 10; % number of steps to not sucessfully find a direction, before changing roles (convert type 1)
    num_directions=[];
elseif convert_type == 2
    num_steps = 16; % number of directions to sample in (convert type 2)
    num_directions = 0.25; % fraction of directions needed to be better to maintain a leader profile (convert type 2)
else
    num_steps=[];
    num_directions=[];
end
    
%% insert_cells parameters %%
insert_time_step = 0.2;         % the time in hours between each new cell insertion
insert_step = floor(insert_time_step/tstep);    % how often are new cells inserted
num_cells = 2;                  % how many new cells are inserted at each timepoint
if insert_cells==1
    n = num_cells;              % initial number of cells
else
    n = 10;
end
initx_perc = 0;                 % initial percentage of x with cells
% follow_start = floor(time/tstep) - floor(time*follow_perc/tstep)+1
follow_start = floor(24/tstep) - floor(24*follow_perc/tstep)+1

%% domain growth parameters %%
width = 300;                % initial width of the domain with growth (um)
domain_length=ones(1,tsteps).*width;  % initialise domain length vector

% These parameters are found from data from the Kulesa Lab, using least
% squares regression in make_domain_plot.m on data in lengths_data.m
Linf = 870;                            % end domain length
% Linf = 830;
a = 0.0800;                            % how fast the domain grows (?)
t_start = -16;
% t_start = -11;
param = [Linf, a, diffus,e, growing_domain,width,make_chemoattractant,chi,height,zero_bc,insert,tstep,t_start,d,mult,num_steps,num_directions];
save avi_mat/param param

%% set up the initial cells so that they aren't too close to each other or
%% the edge %%
temp = initiate_cells(n,cell_radius,0,width,height,initx_perc,inity_perc,[]);
cells = temp.cells;
end_num_cells = (n + floor(tsteps/insert_step)*num_cells)*2;    % final number of cells expected (*2 for divisions and experimental insertions)
cells_follow = [zeros(end_num_cells -floor(56* follow_perc),1); ones(floor(end_num_cells * follow_perc),1)];
% so the first cells out are leaders and the subsequent ones are followers

%% initialise vectors and time %%
t_save = 0:tstep:tstep*tsteps;
xlat_save = cell(1,tsteps);
ylat_save = cell(1,tsteps);
ca_save = cell(1,tsteps);
cells_save = cell(tsteps,1);
filopodia_save = cell(tsteps,1);
cells_follow_save = cell(tsteps,1);
attach = zeros(end_num_cells,1);
attach_save = cell(1,tsteps);
xlat_new=[];
barrier = zeros(tsteps,1);
moved = zeros(tsteps,end_num_cells);
num_better_foll_save = [];
num_foll_save = 0;
num_better_lead_save = [];
num_lead_save = 0;
time_in_domain = zeros(1,length(cells));
%% begin timesteps %%
for k=1:tsteps
    
    if k==floor(tsteps/2)
        multiplier = 0;
        saved.multip=1;
        save avi_mat/saved saved
        save avi_mat/multiplier multiplier
    end
    %% after t=follow_start then all subsequent cells are followers
    if k==follow_start
        cells_follow(length(cells(1,:))+1:end) = ones(size(cells_follow(length(cells(1,:))+1:end)));
        disp('followers start here')
    end
    
    %% Cells and barrier insertions for experiments
    insert_cells_and_barrier
    
    %% If we are inserting new cells, do so here %%
    if mod(k,insert_step)==0
        fprintf(['t = ' mat2str(t_save(k+1)) '\r'] )
        if (insert_cells==1)&&((experiment==0)||(experiment>3)||(in.it==1)||(in.ablate_type~=2)||t_save(k)<in.ablate_time)
            temp = initiate_cells(num_cells,cell_radius,0,width,height,0,inity_perc,cells);
            cells = temp.cells;
        end
    end
    
    %% chemoattractant %%
    if ca_solve==1
        %% cells stop consuming c'tant after satiate timesteps
        if satiate~=0
            time_in_domain(end+1:length(cells)) = 1;
            time_in_domain = time_in_domain + ones(1,length(time_in_domain));
            cells_in = cells(:,time_in_domain<satiate);
        else
            cells_in = cells;
        end
        % give parameters for the solver (depending on whether this is the first run or not)
        if t_save(k)==0
            if isunix==1
                ind = int64(0); % beginning at t=0
                iwk = zeros(580230,1,'int64');
            else
                ind = int32(0);
                iwk = zeros(580230,1,'int32');
            end
            rwk = zeros(1880000,1);
        elseif ((experiment==1)||(experiment==2))&&(in.it==2)&&(t_save(k)==in.change_time)
            %% experiments 1 and 2: inserting chemoattractant
            insert_tissue
        else
            if isunix==1
                ind = int64(1); % continuing from the previous solution
            else
                ind = int32(1); % continuing from the previous solution
            end
        end
        % run the solver
        if ((experiment==1)||(experiment==2))&&(in.it==2)&&(t_save(k)==in.change_time)
            temp = chemotaxis_solve(t_save(k),t_save(k+1),ind,iwk,rwk,cells_in,width,height,length(xlat_new),length(ylat_new),insert);
        elseif k>1
            temp = chemotaxis_solve(t_save(k),t_save(k+1),ind,iwk,rwk,cells_in,width,height,length(xlat_new),length(ylat_save{k-1}),insert);
        else
            temp = chemotaxis_solve(t_save(k),t_save(k+1),ind,iwk,rwk,cells_in,width,height,length(xlat_new),50,insert);
        end

        % take output 
        xlat_save{k} = temp.xsave;
        ylat_save{k} = temp.ysave;
        iwk = temp.iwk;
        rwk = temp.rwk;
        
        ca_save{k} = temp.chemotaxis;
    else
        %% Fixed chemoattractant %%
        if growing_domain==1
            % Domain Growth Happens at every timestep
            if ((experiment==4)||(experiment==5))
                temp = domain_growth(cells(1,:),t_save(k),tstep,Linf,a,width,barrier(k),t_start);
                barrier(k+1)=temp.barrier;
            else
                temp = domain_growth(cells(1,:),t_save(k),tstep,Linf,a,width,[],t_start);
            end
            domain_length(k) = temp.domain_length;
        end
        
        xlat_save{k} = 0:domain_length(k)/100:domain_length(k);
        ylat_save{k} = 0:height/100:height;
        ca = zeros(length(xlat_save{k}),length(ylat_save{k}));
        for i=1:length(xlat_save{k})
            for j=1:length(ylat_save{k})
                ca(i,j) = xlat_save{k}(i)./domain_length(k);
            end
        end
        ca_save{k} = ca;
    end
    %% divide cells %%
    if (divide_cells==1)
        temp = cells_divide(cells,cells_follow,cell_radius,domain_length(k),0,height,ca_save{k},xlat_save{k},ylat_save{k},tstep);
        cells = temp.cells;
        cells_follow = temp.cells_follow;
    end
    %% domain growth %%
    if growing_domain==1
        % Domain Growth Happens at every timestep
        if ((experiment==4)||(experiment==5))
            temp = domain_growth(cells(1,:),t_save(k),tstep,Linf,a,width,barrier(k),t_start);
            barrier(k+1)=temp.barrier;
        else
            temp = domain_growth(cells(1,:),t_save(k),tstep,Linf,a,width,[],t_start);
        end
        domain_length(k) = temp.domain_length;
        cells(1,:) = temp.cells_next;
        cells(2,:) = cells(2,:);
    end
    
    %% move cells %%
    if cells_move==1
        if k==1
            temp = new_move_cells(cells,cells_follow,[],attach,[],...
                ca_save{k},xlat_save{k},ylat_save{k},...
                cell_radius,filolength,d,height,dist,domain_length(k),barrier(k),experiment,t_save(k),in,metropolis,num_filopodia, diffusion);
        else
            temp = new_move_cells(cells,cells_follow,filopodia,attach,theta,...
                ca_save{k},xlat_save{k},ylat_save{k},...
                cell_radius,filolength,d,height,dist,domain_length(k),barrier(k),experiment,t_save(k),in,metropolis,num_filopodia, diffusion);
        end
        attach = temp.attach;
        cells_follow = temp.cells_follow;
        filopodia = temp.filopodia;
        theta = temp.theta;
        cells = temp.cells;
        moved(k,:) = [temp.moved, zeros(1,length(moved(1,:))-length(temp.moved))];
        
        attach_save{k} = attach;
        cells_follow_save{k} = cells_follow;
        filopodia_save{k} = filopodia;
    end
    cells_save{k}=cells;
    
    %% cells can convert from leaders <-> followers
    if (convert_type~=0)&&((experiment==0)||experiment==3||(in.it==1)||(t_save(k)==in.change_time))
        out = convert_cells(cells,cells_follow,attach_save,k,cells_save,filolength,moved,ca_save{k},xlat_save{k},ylat_save{k},...
            d,filopodia,convert_type,param,num_better_foll_save,num_foll_save,num_better_lead_save,num_lead_save);
        cells_follow = out.cells_follow;
        moved = out.moved;
        num_better_foll_save = out.num_better_foll_save;
        num_foll_save = out.num_foll_save;
        num_better_lead_save = out.num_better_lead_save;
        num_lead_save = out.num_lead_save;
        %cells_follow = convert_cells(cells,cells_follow,attach_save,xlat_save{k},ylat_save{k},ca_save{k});
        cells_follow_save{k} = cells_follow;
    end
    

    if (cells~=cells_save{k})
        size(cells)
        size(cells_save)
        disp('cells error')
        pause
    end
    if (k>=follow_start)
        if (attach~=attach_save{k})
            disp('attach error')
            pause
        end
        if (cells_follow~=cells_follow_save{k})
            disp('follow error')
            pause
        end
        if (filopodia~=filopodia_save{k})
            disp('filopodia error')
            pause
        end
    end
end
%%
hist(num_better_foll_save/num_foll_save)
hist(num_better_lead_save/num_lead_save)
out.num_better_foll_save = num_better_foll_save;
out.num_foll_save = num_foll_save;
out.num_better_lead_save = num_better_lead_save;
out.num_lead_save = num_lead_save;
toc
whitebg('white')

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
    open(['avi_mat/frames/frames3',save_info,'.fig'])
end
%% calculate average directionality %%%
% for i=1:n
%     total = sum(sqrt((cells_save{2:end}(1,i)-cells_save{1:end-1}(1,i)).^2+(cells_save{2:end}(2,i)-cells_save{1:end-1}(2,i)).^2));
%     straight = sqrt((cells_save{end}(1,i)-cells_save{1}(1,i))^2 + (cells_save{end}(2,i)-cells_save{1}(2,i))^2);
%     directionality(i) = straight/total;
% end
% directionality
% average_direc = mean(directionality)
% average_speed =
% mean((cells_save(1,1:n,end)-cells_save(1,1:n,1))./t_save(end))

delete('avi_mat/cells.mat','avi_mat/param.mat','avi_mat/plotsol.mat','avi_mat/xsave.mat','avi_mat/ysave.mat')
