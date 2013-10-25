% Louise Dyson D.Phil project CA program, 05/11/09
% Initiate n cells with radius cellRadius, so that they don't intersect
% cells are initiated with x and y coordinates in cells(1,:,1) and
% cells(2,:,1), resp.
% cells is given m+1 timesteps to come

function out = initiate_cells(n,cellRadius,followerFraction,initialDomainLength,domainHeight,initx_frac,inity_frac,cells_in,volumeExclusion)
%% set up the initial cells so that they aren't too close to each other or the edge %%%
cells = [cells_in NaN*ones(2,n)];
cellsFollow = false(length(cells(1,:,1)),1);
[~,j] = size(cells_in);

%% If there are follower cells, initialise them here %%%
if followerFraction~=0
    %% Iterate until the cells vector is full to the correct amount %%%
    while isnan(cells(1,floor(end*followerFraction),1))
        %% postulate coordinates for a new cell %%%
        if initx_frac==0
            cellx=cellRadius;
        else
            cellx = rand().*initx_frac.*followerFraction.*initialDomainLength;
        end
        celly = (rand().*inity_frac+(0.5-inity_frac/2)).*domainHeight;
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
%% Iterate until the cells vector is full %%%
it = 0;
maxit = 50;
while isnan(cells(1,end,1))&&(it<maxit)
    %% postulate coordinates for a new cell %%%
    if initx_frac==0
        cellx = cellRadius;
    else
        cellx = (rand().*initx_frac.*(1-followerFraction)+initx_frac*followerFraction).*initialDomainLength;
    end
    celly = (rand().*inity_frac+(0.5-inity_frac/2))*domainHeight;

    %% Check that there is no overlap with existing cells %%%
    if (j>0)%&&(volumeExclusion==1)
        %%% (x_cells,y_cells) are the pre-existing cell coordinates %%%
        x_cells = cells(1,1:j,1);
        y_cells = cells(2,1:j,1);
        %% If there is no overlap, take that coordinate else increase 'it' %%
        if (sqrt(min((cellx - x_cells).^2+(celly-y_cells).^2))>2*cellRadius)...
                &&(cellx>=cellRadius)&&(cellx<initialDomainLength-cellRadius)&&(celly>cellRadius)&&(celly<domainHeight-cellRadius)
            cells(:,j+1,1) = [cellx;celly];
            j = j+1;
        else
            it = it+1;
        end
    elseif (j==0)%||(volumeExclusion==0) % if it's the first cell or if we don't have volume exclusion
        cells(:,j+1,1) = [cellx;celly];
        j = j+1;
    end
end
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