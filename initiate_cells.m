% Louise Dyson D.Phil project CA program, 05/11/09
% Initiate n cells with radius cellRadius, so that they don't intersect
% cells are initiated with x and y coordinates in cells(1,:,1) and
% cells(2,:,1), resp.
% cells is given m+1 timesteps to come

function out = initiate_cells(numCells2Insert,cellRadius,followerFraction,initialDomainLength,domainHeight,initXFrac,initYFrac,cells_in,volumeExclusion)
%% set up the initial cells so that they aren't too close to each other or the edge %%%
cells = [cells_in NaN(2,numCells2Insert)];
cellsFollow = false(length(cells(1,:,1)),1);
[~,j] = size(cells_in);
if j==0 % if we're inserting the first set of cells, position them uniformly, to have repeatable initial conditions
    deltaY = domainHeight*initYFrac/(numCells2Insert + 1); % spacing between uniformly spread cells
    cells(:,1:numCells2Insert,1) = [cellRadius(ones(1,numCells2Insert));...
        max(cellRadius,(1 - initYFrac)/2*domainHeight + deltaY):deltaY:...
        min(domainHeight - cellRadius,domainHeight*(1 + initYFrac)/2 - deltaY)];
end
%% If there are follower cells, initialise them here %%%
% I'm fairly sure this part is only relevant if we want to seed more than
% just the very left edge with cells at the start - leave this here for now
% as it might be important for the old transplantation experiments
if followerFraction~=0
    %% Iterate until the cells vector is full to the correct amount %%%
    while isnan(cells(1,floor(end*followerFraction),1))
        %% postulate coordinates for a new cell %%%
        if initXFrac==0
            cellx=cellRadius;
        else
            cellx = rand().*initXFrac.*followerFraction.*initialDomainLength;
        end
        celly = (rand().*initYFrac+(0.5-initYFrac/2)).*domainHeight;
        %% Check that there is no overlap with existing cells %%%
        if (j>0)%&&(volumeExclusion==1)
            %%% (x_cells,y_cells) are the pre-existing cell coordinates %%%
            x_cells = cells(1,1:j,1);
            y_cells = cells(2,1:j,1);
            if (sqrt(min((cellx - x_cells).^2+(celly-y_cells).^2))>2*cellRadius)...
                    &&(cellx>=cellRadius)&&(cellx<initialDomainLength-cellRadius)&&(celly>cellRadius)&&(celly<domainHeight-cellRadius)
                cells(:,j+1,1) = [cellx;celly];
                j = j+1;
            end
        elseif (j==0)%||(volumeExclusion==0) % if it's the first cell or if we don't have volume exclusion
            cells(:,j+1,1) = [cellx;celly];
            j = j+1;
        end
    end
    cellsFollow(1:floor(end*followerFraction))=ones(floor(length(cells(1,:,1))*followerFraction),1);
end

%% Initialise all other cells here
if isnan(cells(1,end,1))
    % (x_cells,y_cells) are the pre-existing cell coordinates
    x_cells = cells(1,1:j,1);
    y_cells = cells(2,1:j,1);
    %% try to insert all the cells that need to be inserted -- LJS
    % generate a number of candidate y-values -- LJS
    numSamples = 30;
    yValues = ((linspace(0,1,numSamples).*initYFrac + (0.5 - initYFrac/2))*domainHeight);
    % try putting cells at these y-values in a random order -- LJS
    yIndcs = randperm(length(yValues));
    for yCtr = yIndcs
        celly = yValues(yCtr);
        %% postulate x-coordinates for a new cell %%%
        if initXFrac==0
            cellx = cellRadius;
        else % this next bit is probably only needed for odd cases where we want to insert cells beyond the left edge of the domain, and if we want to use it it should be changed to not randomly sample x-coordinates
            cellx = (rand().*initXFrac.*(1 - followerFraction) + initXFrac*followerFraction).*initialDomainLength;
        end
        %% If there is no overlap, take that coordinate
        if (sqrt(min((cellx - x_cells).^2+(celly-y_cells).^2))>2*cellRadius)...
                &&(cellx>=cellRadius)&&(cellx<initialDomainLength-cellRadius)...
                &&(celly>cellRadius)&&(celly<domainHeight-cellRadius)
            % if it fits, insert cell
            cells(:,j+1,1) = [cellx; celly];
            j = j+1;
            % if this was the last cell, break out of the loop, else keep
            % going through the other yValues
            if ~isnan(cells(1,end))
                break
            end
            % update coordinates of cells that might be in the way
            x_cells = cells(1,1:j,1);
            y_cells = cells(2,1:j,1);
        end
    end
end
%%
if length(cells(1,any(cells,1)==0))==1
    disp([mat2str(length(cells(1,any(cells,1)==0))),' cell could not be initialised'])
elseif length(cells(1,any(cells,1)==0))>1
    disp([mat2str(length(cells(1,any(cells,1)==0))),' cells could not be initialised'])
end
cells(:,any(cells,1)==0)=[];

out.cells = cells;
out.cellsFollow = cellsFollow;
if isempty(cells_in)==1
    disp('cells initiated')
end