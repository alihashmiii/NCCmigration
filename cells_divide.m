% Louise Dyson D.Phil project second CA program, 14/02/10
% takes in cell (x_cell,y_cell) positions and follower information (cellsFollow) 
% and outputs new cell positions including new cells due to division

function out = cells_divide(cells,cellsFollow,cellRadius,domainLength,bottom,top,ca,ca_x,ca_y,tstep)

%% find the starting number of cells %%
n = length(cells(1,:));

%% find the chemoattractant concentration at each cell %%
chemo = find_ca(cells,ca_x,ca_y,ca);%,cellRadius);

%% a proportion of the cells divide %%
rate_div = 0.1/1.05;
if length(chemo)<3
    cells_to_divide = find(floor(2.*rand(1,length(chemo)))==1);
%     cells_to_divide = find(floor(2.*rand(1,length(cells)))==1); % half
%     the cells divide
else
    cells_should_divide = find((chemo-min(chemo))/(max(chemo)-min(chemo))>0.9); % find the positions in 'cells' of the the cells to divide
%    cells_should_not_divide = find((chemo-min(chemo))/(max(chemo)-min(chemo))<=0.8);
    cells_to_divide = cells_should_divide(floor(rate_div*tstep/(1-rate_div*tstep).*rand(1,length(cells_should_divide)))==1); % cells_should_not_divide(floor(1.00*rand(1,length(cells_should_not_divide)))==1)];
end

%% expand the cells and cellsFollow matrices %%
cells = [cells zeros(2,length(cells_to_divide))];
cellsFollow = [cellsFollow; false(length(cells_to_divide),1)];

%% iterate through the cells that are dividing %%
max_it = 10;
for i=1:length(cells_to_divide)
    
    r = cells_to_divide(i); % position in 'cells' of the cell that is dividing
    
    %% try different directions of division until an empty area is found %%
    success=0;
    it = 0;
    while (success==0)&&(it<max_it)
        theta = rand()*2*pi; % direction that the new cell will appear in
        
        %% find the putative new coordinates
        new_x = cells(1,r) + cos(theta)*2*cellRadius;
        new_y = cells(2,r) + sin(theta)*2*cellRadius;
        
        %% is there a cell or an edge in the new area?
        if ((sqrt(min((new_x- cells(1,:)).^2+(new_y-cells(2,:)).^2))>2*cellRadius)...
                &&(new_x>cellRadius)&&(new_x<domainLength-cellRadius)&&(new_y>cellRadius+bottom)&&(new_y<top-cellRadius))
            cells(1,n+i) = new_x;
            cells(2,n+i) = new_y;
            cellsFollow(n+i) = cellsFollow(r);
            success = 1;
        else
            it = it+1;
        end
    end
    if it==max_it
        cells(:,end) = [];
        cellsFollow(end)=[];
        n = n-1;
    else
        disp('cell divides')
    end
end

%% set output %%
out.cells = cells;
out.cellsFollow = cellsFollow;