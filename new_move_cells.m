function out = new_move_cells(cellsFollow,filopodia,attach,theta,...
    ca_save,xlat,ylat,cellRadius, filolength, maxFilolength, eatWidth, ...
    domainHeight, dist, domainLength, numFilopodia,volumeExclusion, ...
    standStill, sensingAccuracy, needNeighbours, contactGuidance, currentTime, dan)
global cells param
%% iterate through the cell movement in a random order %%%
cell_order = randperm(length(cells(1,:)));
moved = false(1,length(cells(1,:)));
sensed = false(1,length(cells(1,:)));
for i =1:length(cell_order)
    move = 0;
    cellIdx = cell_order(i);  % look at the ith cell
    other_cells = cells(:,(1:end)~=cellIdx);
    
    %% calculate neighbours within reach (work in progress) -- LJS
    distance = sqrt((cells(1,cellIdx) - other_cells(1,:)).^2 + (cells(2,cellIdx) - other_cells(2,:)).^2);
    numberOfNeighbours = nnz(distance <= filolength);
    
    %% Decide whether to try to move
    if cellsFollow(cellIdx)~=1
        %% if it's a leader
        [filopodia(cellIdx,1:numFilopodia(1),:),move,theta(cellIdx),~] = cell_movement5((rand(1,numFilopodia(1))*2 - 1)*pi,cells(1,cellIdx),cells(2,cellIdx),ca_save,xlat,ylat,...
            eatWidth,filolength,numFilopodia(1),[],sensingAccuracy);
        sensed(cellIdx) = move; % -- LJS
    else
        %% if it's a follower -- LJS
        if attach(cellIdx)~=0
            %% if it's a chained follower
            if (cells(1,attach(cellIdx)) - cells(1,cellIdx))^2 + (cells(2,attach(cellIdx)) - cells(2,cellIdx))^2 < (maxFilolength + cellRadius)^2
                %% if it can reach the cell ahead
                % set (first) filopodium position to closet point on membrane of cell being followed -- LJS
                phi = atan2((cells(2,attach(cellIdx)) - cells(2,cellIdx)),(cells(1,attach(cellIdx)) - cells(1,cellIdx))); % the angle towards the cell being followed -- LJS
                filopodia(cellIdx,1,1) = cells(1,attach(cellIdx)) - cellRadius*cos(phi);
                filopodia(cellIdx,1,2) = cells(2,attach(cellIdx)) - cellRadius*sin(phi);
                if strcmp(contactGuidance,'parallel')
                    % set direction of movement parallel to that of cell being
                    % followed -- LJS
                    theta(cellIdx) = theta(attach(cellIdx));
                elseif strcmp(contactGuidance,'toward')
                    % set direction of movement towards the cell being followed -- LJS
                    theta(cellIdx) = phi;
                end
                move = 1;
            else
                %% if the cell ahead is too far away, then dettach the chain
                % set (first) filopodial position and movement angle in the direction of cell
                % centre of lost cell -- LJS
                theta(cellIdx) = atan2((cells(2,attach(cellIdx)) - cells(2,cellIdx)),(cells(1,attach(cellIdx)) - cells(1,cellIdx)));
                filopodia(cellIdx,1,:) = cells(:,cellIdx)' + filolength.*[cos(theta(cellIdx)) sin(theta(cellIdx))];
                attach = dettach(cellIdx,attach);
            end
            % set any other filopodia in random direction (for now, until I
            % have a better idea) -- LJS
            if numFilopodia(2) > 1
                phi = (rand(1,numFilopodia(2)-1)*2 - 1)*pi;
                filopodia(cellIdx,2:numFilopodia(2),1) = cells(1,cellIdx) + filolength.*cos(phi);
                filopodia(cellIdx,2:numFilopodia(2),2) = cells(2,cellIdx) + filolength.*sin(phi);
            end
        end
        if attach(cellIdx)==0
            %% if it's an unchained follower (this uncludes previously chained, now dettached followers from the previous if-statement -- LJS)
            % look for other cells
            [foundCellidx,filopodia] = cell_movement5_follow((rand(1,numFilopodia(2))*2 - 1)*pi,cellIdx,cells(1,:),cells(2,:),cellRadius,...
                filolength,filopodia);
            if isempty(foundCellidx)~=1
                %% if another cell was found then find the head of that chain
                head = foundCellidx;
                while attach(head)~=0
                    head = attach(head);
                end
                if cellsFollow(head)==0    %% if the head is a leader, then attach and move
                    attach(cellIdx) = foundCellidx;
                    % set filopodium position to closet point on membrane of cell being followed -- LJS
                    phi = atan2((cells(2,attach(cellIdx)) - cells(2,cellIdx)),(cells(1,attach(cellIdx)) - cells(1,cellIdx))); % the angle towards the cell being followed -- LJS
                    filopodia(cellIdx,1,1) = cells(1,attach(cellIdx)) - cellRadius*cos(phi);
                    filopodia(cellIdx,1,2) = cells(2,attach(cellIdx)) - cellRadius*sin(phi);
                    if strcmp(contactGuidance,'parallel')
                        % set direction of movement parallel to that of cell being
                        % followed -- LJS
                        theta(cellIdx) = theta(attach(cellIdx));
                    elseif strcmp(contactGuidance,'toward')
                        % set direction of movement towards the cell being followed -- LJS
                        theta(cellIdx) = phi;
                    end
                    move = 1;
                end
            end
        end
        % check if favourable chemoattractant gradients can be sensed, for
        % phenotype conversion -- LJS
        [~,sensed(cellIdx),~,~] = cell_movement5([],cells(1,cellIdx),cells(2,cellIdx),ca_save,xlat,ylat,...
            eatWidth,filolength,numFilopodia(2),squeeze(filopodia(cellIdx,:,:)),sensingAccuracy);
    end
    if numberOfNeighbours < needNeighbours % check if a cell should wait around for others
        move = 0;
    end
    % end
    %
    %% Try to move
    % for i =1:length(cell_order)
    %     cellidx = cell_order(i);  % look at the ith cell
    %     other_cells = cells(:,(1:end)~=cellidx);
    
    if (standStill==0)&&(move==0) % if standStill = 0, cells move in a random direction
        theta(cellIdx) = (rand()*2 - 1)*pi; % pick a random direction for movement
    end
    if (move==1)||((standStill==0)&&(move==0))
        if move==1, moved(cellIdx)=1; end
        if ((param.experiment==40)||(param.experiment==41)||(param.experiment==42)...
                ||(param.experiment==43))&&(cells(1,cellIdx)<=1/3*domainLength) % move at reduced speed
            switch param.experiment
                case 40 % slow down is constant over time
                    slowDown = 1;
                case 41 % slow down is diluted with tissue growth
                    slowDown = param.initialDomainLength/domainLength;
                case 42 % slow down has it's own simple dynamics, first increases, then decreases with time
                    tPeakSlowdown = 12;
                    minSlowdown = 0.5;
                    slowDown = max(minSlowdown,...
                        (tPeakSlowdown - abs(currentTime - tPeakSlowdown))/tPeakSlowdown);
                case 43 % slow down is proportional to DAN conc in new loc (0 or 1)
                    try_x = cells(1,cellIdx) + cos(theta(cellIdx))*dist(1);
                    try_y = cells(2,cellIdx) + sin(theta(cellIdx))*dist(1);
                    xrange = xlat(xlat>=0&xlat<=max(xlat)/3);
                    % xsave and ysave are the coordinates of the chemoattractant lattice
                    % points. dan has the same y-coordinates, but only makes up a 1/3 of the x range
                    slowDown = mean(mean(dan(find_occupancy(xrange,ylat,try_x,try_y,cellRadius)))); % take the mean twice in case the cell sits on multiple lattice points
                otherwise
                    slowDown = 0;
            end
            new_x = cells(1,cellIdx) + cos(theta(cellIdx))*...
                (dist(1) - (dist(1) - dist(3))*slowDown);
            new_y = cells(2,cellIdx) + sin(theta(cellIdx))*...
                (dist(1) - (dist(1) - dist(3))*slowDown);
        else
            if (cellsFollow(cellIdx)==1) %if it's a follower
                new_x = cells(1,cellIdx) + cos(theta(cellIdx))*dist(2);
                new_y = cells(2,cellIdx) + sin(theta(cellIdx))*dist(2);
            else
                new_x = cells(1,cellIdx) + cos(theta(cellIdx))*dist(1);
                new_y = cells(2,cellIdx) + sin(theta(cellIdx))*dist(1);
            end
        end
        
        if volumeExclusion==1 %% if there is no cell or edge in the way, then move
            diff = [new_x-other_cells(1,:); new_y-other_cells(2,:)];
            if (length(cell_order)==1)||(min(vnorm(diff))>2*cellRadius)
                if (new_x>cellRadius)&&(new_x<domainLength-cellRadius)&&(new_y>cellRadius)&&(new_y<domainHeight-cellRadius)&&(new_x>min(xlat)+cellRadius) % this last condition may be important in tissue transplantation experiments?
                    cells(1,cellIdx) = new_x;
                    cells(2,cellIdx) = new_y;
                end
            end
        elseif volumeExclusion==0 % move regardless of other cell's position
            if (new_x>cellRadius)&&(new_x<domainLength-cellRadius)&&(new_y>cellRadius)&&(new_y<domainHeight-cellRadius)&&(new_x>min(xlat)+cellRadius) % this last condition may be important in tissue transplantation experiments?
                cells(1,cellIdx) = new_x;
                cells(2,cellIdx) = new_y;
            end
        end
    end
    if volumeExclusion==1
        if min(sqrt((cells(1,cellIdx)-cells(1,(1:end)~=cellIdx)).^2 + (cells(2,cellIdx)-cells(2,(1:end)~=cellIdx)).^2))<=(2*cellRadius)
            disp('error - cells moving through each other!')
            pause
        end
    end
end

%% save stuff
out.attach = attach;
out.cellsFollow = cellsFollow;
out.filopodia = filopodia;
out.theta = theta;
out.moved = moved;
out.sensed = sensed;