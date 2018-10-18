function out = move_cells_cont_states(param,cells,followerness,filopodia,attach,theta,...
    ca_save,xlat,ylat,cellRadius, filolength, maxFilolength, eatWidth, ...
    domainHeight, dist, domainLength, numFilopodia,volumeExclusion, ...
    standStill, sensingAccuracy, needNeighbours, contactGuidance, currentTime, dan)
%% iterate through the cell movement in a random order %%%
cell_order = randperm(length(cells(1,:)));
moved = false(1,length(cells(1,:)));
leaderness = 1-followerness;
for i =1:length(cell_order)
    moveDirected = 0;
    cellIdx = cell_order(i);  % look at the ith cell
    other_cells = cells(:,(1:end)~=cellIdx);
    
    %% calculate neighbours within reach (work in progress) -- LJS
    distance = sqrt((cells(1,cellIdx) - other_cells(1,:)).^2 + (cells(2,cellIdx) - other_cells(2,:)).^2);
    numberOfNeighbours = nnz(distance <= filolength);
    
    %% check attachment and compute directions of filopodia
    % this also computes contact guidance direction
    
%     % stay attached to cell with probability ~ its leaderness
%     if attach(cellIdx)~=0
%         p = rand();
%         regularizer = 0;
%         if p>(regularizer + leaderness(attach(cellIdx)))/(regularizer + 1)
%             attach(cellIdx) = 0; % dettach
%         end
%     end
    if attach(cellIdx)~=0
        %% if it's already in filopodial contact with other cell
        if (cells(1,attach(cellIdx)) - cells(1,cellIdx))^2 + (cells(2,attach(cellIdx)) - cells(2,cellIdx))^2 < (maxFilolength + cellRadius)^2
            %% if it can reach the cell ahead
            % set (first) filopodium position to closet point on membrane of cell being followed -- LJS
            phi = atan2((cells(2,attach(cellIdx)) - cells(2,cellIdx)),(cells(1,attach(cellIdx)) - cells(1,cellIdx))); % the angle towards the cell being followed -- LJS
            filopodia(cellIdx,1,1) = cells(1,attach(cellIdx)) - cellRadius*cos(phi);
            filopodia(cellIdx,1,2) = cells(2,attach(cellIdx)) - cellRadius*sin(phi);
            if strcmp(contactGuidance,'parallel')
                % set direction of movement parallel to that of cell being
                % followed -- LJS
                thetaContactGuidance = theta(attach(cellIdx));
            elseif strcmp(contactGuidance,'toward')
                % set direction of movement towards the cell being followed -- LJS
                thetaContactGuidance = phi;
            end
            moveDirected = 1;
        else
            %% if the cell ahead is too far away, then dettach
            % set movement angle in the direction of cell centre of lost cell -- LJS
            thetaContactGuidance = atan2((cells(2,attach(cellIdx)) - cells(2,cellIdx)),(cells(1,attach(cellIdx)) - cells(1,cellIdx)));
            %             attach = dettach(cellIdx,attach);
            attach(cellIdx) = 0;
            % set dettached (first) filopodium in random direction
            phi = (rand(1,1)*2 - 1)*pi;
            filopodia(cellIdx,1,1) = cells(1,cellIdx) + filolength.*cos(phi);
            filopodia(cellIdx,1,2) = cells(2,cellIdx) + filolength.*sin(phi);
        end
        % set any other filopodia in random direction -- LJS
        if numFilopodia > 1
            phi = (rand(1,numFilopodia-1)*2 - 1)*pi;
            filopodia(cellIdx,2:numFilopodia,1) = cells(1,cellIdx) + filolength.*cos(phi);
            filopodia(cellIdx,2:numFilopodia,2) = cells(2,cellIdx) + filolength.*sin(phi);
        end
    else % sample numFilopodia random directions
        phi = (rand(1,numFilopodia)*2 - 1)*pi;
        filopodia(cellIdx,:,1) = cells(1,cellIdx) + filolength.*cos(phi);
        filopodia(cellIdx,:,2) = cells(2,cellIdx) + filolength.*sin(phi);
    end
    
    %% Compute movement direction based on gradient signalling
    
    % sample gradient signal in directions of filopodia
    [~,moveCtaxis,thetaChemotaxis,~,deltaC] = ...
        sense_gradients([],cells(1,cellIdx),cells(2,cellIdx),ca_save,xlat,ylat,...
        eatWidth,filolength,numFilopodia,squeeze(filopodia(cellIdx,:,:)),sensingAccuracy);
    if moveCtaxis, moveDirected=1; end
    leaderness(cellIdx) = max(0,deltaC); % should be between 0 and 1
    
    %% check if extended filopodia have touched another cell
    % this also computes contact guidance direction
    if attach(cellIdx)==0
        % if cell is not in contact with another (this uncludes previously chained, now dettached followers from the previous if-statement -- LJS)
        % look for other cells
        [foundCellidx,filopodia] = sense_contacts([],cellIdx,...
            cells(1,:),cells(2,:),cellRadius,filolength,filopodia);
        if isempty(foundCellidx)~=1
            % if another cell was found then attach
            attach(cellIdx) = foundCellidx;
            % set (first) filopodium position to closet point on membrane of cell being followed -- LJS
            phi = atan2((cells(2,attach(cellIdx)) - cells(2,cellIdx)),...
                (cells(1,attach(cellIdx)) - cells(1,cellIdx))); % the angle towards the cell being followed -- LJS
            filopodia(cellIdx,1,1) = cells(1,attach(cellIdx)) - cellRadius*cos(phi);
            filopodia(cellIdx,1,2) = cells(2,attach(cellIdx)) - cellRadius*sin(phi);
            if strcmp(contactGuidance,'parallel')
                % set direction of movement parallel to that of cell being
                % followed -- LJS
                thetaContactGuidance = theta(attach(cellIdx));
            elseif strcmp(contactGuidance,'toward')
                % set direction of movement towards the cell being followed -- LJS
                thetaContactGuidance = phi;
            end
            moveDirected = 1;
        else % no cell found for contact guidance
            thetaContactGuidance = (rand()*2 - 1)*pi;
        end
    end
    
    if numberOfNeighbours < needNeighbours % check if a cell should wait around for others
        moveDirected = 0;
    end
    
    %% Try to move
    
    if (standStill==0)&&(moveDirected==0) % if standStill = 0, cells move in a random direction
        theta(cellIdx) = (rand()*2 - 1)*pi; % pick a random direction for movement
    else
%         % combine directional information from gradient and contact guidance
%         xc = cos(thetaContactGuidance);
%         yc = sin(thetaContactGuidance);
%         xg = cos(thetaChemotaxis);
%         yg = sin(thetaChemotaxis);
%         xcombined = xc*(1 - leaderness(cellIdx)) + xg*leaderness(cellIdx);
%         ycombined = yc*(1 - leaderness(cellIdx)) + yg*leaderness(cellIdx);
%         theta(cellIdx) = atan2(ycombined,xcombined);
        
                % randomly choose between gradient and contact guidance
                p = rand();
                if p<=leaderness(cellIdx)
                    theta(cellIdx) = thetaChemotaxis;
                else
                    theta(cellIdx) = thetaContactGuidance;
                end
    end
    
    if (moveDirected==1)||((standStill==0)&&(moveDirected==0))
        if moveDirected==1, moved(cellIdx)=1; end
        if ((param.experiment==40)||(param.experiment==41)||(param.experiment==42)...
                ||(param.experiment==43)||(param.experiment==44))&&(cells(1,cellIdx)<=1/3*domainLength) % move at reduced speed
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
                case 44 % slow down has intrinsic dynamics, but is also broken down by cells
                    tPeakSlowdown = 12;
                    minSlowdown = 0.5;
                    slowDownDynamics = max(minSlowdown,...
                        (tPeakSlowdown - abs(currentTime - tPeakSlowdown))/tPeakSlowdown);
                    try_x = cells(1,cellIdx) + cos(theta(cellIdx))*dist(1);
                    try_y = cells(2,cellIdx) + sin(theta(cellIdx))*dist(1);
                    xrange = xlat(xlat>=0&xlat<=max(xlat)/3);
                    % xsave and ysave are the coordinates of the chemoattractant lattice
                    % points. dan has the same y-coordinates, but only makes up a 1/3 of the x range
                    slowDown = slowDownDynamics.*mean(mean(dan(find_occupancy(xrange,ylat,try_x,try_y,cellRadius)))); % take the mean twice in case the cell sits on multiple lattice points
                otherwise
                    slowDown = 0;
            end
            new_x = cells(1,cellIdx) + cos(theta(cellIdx))*...
                (dist(1) - (dist(1) - dist(3))*slowDown);
            new_y = cells(2,cellIdx) + sin(theta(cellIdx))*...
                (dist(1) - (dist(1) - dist(3))*slowDown);
        else
            dist_combined = leaderness(cellIdx)*dist(1) + (1 - leaderness(cellIdx))*dist(2);
            new_x = cells(1,cellIdx) + cos(theta(cellIdx))*dist_combined;
            new_y = cells(2,cellIdx) + sin(theta(cellIdx))*dist_combined;
        end
        % turn setting new position respecting boundary conditions into
        % function
        if volumeExclusion==1&&length(cell_order)>1 %% if there is no cell or edge in the way, then move
            diff = [new_x-other_cells(1,:); new_y-other_cells(2,:)];
            if min(vnorm(diff))>2*cellRadius
                if new_x>cellRadius&&new_x<domainLength-cellRadius&&new_y>cellRadius&&new_y<domainHeight-cellRadius&&new_x>min(xlat)+cellRadius % this last condition may be important in tissue transplantation experiments?
                    cells(1,cellIdx) = new_x;
                    cells(2,cellIdx) = new_y;
                elseif new_y<cellRadius||new_y>domainHeight-cellRadius
                    theta(cellIdx) = -theta(cellIdx); % reflective boundary
                    if new_y>domainHeight-cellRadius
                        new_y = 2*(domainHeight-cellRadius) - new_y;
                    elseif new_y<cellRadius
                        new_y = 2*cellRadius - new_y;
                    end
                    diff = [new_x-other_cells(1,:); new_y-other_cells(2,:)];
                    if min(vnorm(diff))>2*cellRadius
                        cells(1,cellIdx) = new_x;
                        cells(2,cellIdx) = new_y;
                    end
                end
            end
        elseif volumeExclusion==0||length(cell_order)==1 % move regardless of other cell's position
            if (new_x>cellRadius)&&(new_x<domainLength-cellRadius)&&(new_y>cellRadius)&&(new_y<domainHeight-cellRadius)&&(new_x>min(xlat)+cellRadius) % this last condition may be important in tissue transplantation experiments?
                cells(1,cellIdx) = new_x;
                cells(2,cellIdx) = new_y;
            end
        end
    end
    if volumeExclusion==1
        if min(sqrt((cells(1,cellIdx)-cells(1,(1:end)~=cellIdx)).^2 + (cells(2,cellIdx)-cells(2,(1:end)~=cellIdx)).^2))<=(2*cellRadius)
            warning('error - cells moving through each other!')
        end
    end
    % check max filo length
    assert((filopodia(cellIdx,1,1) - cells(1,cellIdx))^2 + (filopodia(cellIdx,1,2) - cells(2,cellIdx))^2 <= (maxFilolength + cellRadius)^2)
end


%% save stuff
out.cells = cells;
out.attach = attach;
out.filopodia = filopodia;
out.theta = theta;
out.moved = moved;
out.leaderness = leaderness;
