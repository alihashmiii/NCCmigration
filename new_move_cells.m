function out = new_move_cells(cells,cellsFollow,filopodia,attach,theta,...
                            ca_save,xlat,ylat,...
                            cellRadius, filolength, eatWidth, domainHeight, dist, domainLength, barrier, experiment, t_save, in, numFilopodia)
%% iterate through the cell movement in a random order %%%
cell_order = randperm(length(cells(1,:)));
moved = zeros(1,length(cells(1,:)));
for i =1:length(cell_order)
    move = 0;
    cellidx = cell_order(i);  % look at the ith cell
    other_cells = cells(:,(1:end)~=cellidx);
    
    %% Decide whether to try to move
    if cellsFollow(cellidx)~=1   
        %% if it's a leader
            [filopodia(cellidx,:,:),E,move,theta(cellidx),~] = cell_movement5((rand(1,numFilopodia(1))*2 - 1)*pi,cells(1,cellidx),cells(2,cellidx),ca_save,xlat,ylat,...
                eatWidth,filolength,numFilopodia(1),[]);        
    elseif attach(cellidx)~=0     
        %% if it's a chained follower that can reach the cell ahead
        if (cells(1,attach(cellidx)) - cells(1,cellidx))^2 + (cells(2,attach(cellidx)) - cells(2,cellidx))^2 < (filolength + cellRadius)^2
%             % set angle of movement towards cell being followed -- LJS
%             theta(cellidx) = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx)));
            % set angle of movement parallel to that of cell being followed
            theta(cellidx) = theta(attach(cellidx));
            % set (first) filopodium position to closet point on membrane of cell being followed -- LJS
            phi = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx))); % the angle towards the cell being followed -- LJS
            filopodia(cellidx,1,1) = cells(1,attach(cellidx)) - cellRadius*cos(phi);
            filopodia(cellidx,1,2) = cells(2,attach(cellidx)) - cellRadius*sin(phi);
        else
            %% if the cell ahead is too far, then dettach the chain
            % set (first) filopodial position and movement angle in the direction of cell
            % centre of lost cell -- LJS
            theta(cellidx) = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx)));
            filopodia(cellidx,1,:) = cells(:,cellidx)' + filolength.*[cos(theta(cellidx)) sin(theta(cellidx))];
            attach = dettach(cellidx,attach);
        end
        % set any other filopodia in random direction (for now, until I
        % have a better idea) -- LJS
        if numFilopodia(2) > 1
            phi = (rand(1,numFilopodia(2)-1) - 1)*pi;
            filopodia(cellidx,2:numFilopodia(2),1) = cells(1,cellidx) + filolength.*cos(phi);
            filopodia(cellidx,2:numFilopodia(2),2) = cells(2,cellidx) + filolength.*sin(phi);
        end
        move = 1;
    else
        %% if it's an unchained follower
        [foundCellidx,filopodia] = cell_movement5_follow((rand(1,numFilopodia(2))*2 - 1)*pi,cellidx,cells(1,:),cells(2,:),cellRadius,...
            filolength,filopodia,barrier,experiment);
        if isempty(foundCellidx)~=1
            %% if another cell was found then find the head of that chain
            head = foundCellidx;
            while attach(head)~=0
                head = attach(head);
            end
            if cellsFollow(head)==0    %% if the head is a leader, then attach and move
                attach(cellidx) = foundCellidx;
                % set angle of movement parallel to that of cell being
                % followed -- LJS
                theta(cellidx) = theta(attach(cellidx));
                % set filopodium position to closet point on membrane of cell being followed -- LJS
                phi = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx))); % the angle towards the cell being followed -- LJS
                filopodia(cellidx,1,1) = cells(1,attach(cellidx)) - cellRadius*cos(phi);
                filopodia(cellidx,1,2) = cells(2,attach(cellidx)) - cellRadius*sin(phi);
                move = 1;
            end
        end
    end
    
    %% Try to move
    if (cellsFollow(cellidx)==1)&&(move==0)
%         currently unchained followers do not move
    end
    if (move==1)
        moved(cellidx)=1;
        if (cellsFollow(cellidx)==1)&&(move==1) %if it's a follower
                new_x = cells(1,cellidx) + cos(theta(cellidx))*dist(2);
                new_y = cells(2,cellidx) + sin(theta(cellidx))*dist(2);
        else
            new_x = cells(1,cellidx) + cos(theta(cellidx))*dist(1);
            new_y = cells(2,cellidx) + sin(theta(cellidx))*dist(1);
        end
        
        %% if there is no cell or edge in the way, then move
        diff = [new_x-other_cells(1,:); new_y-other_cells(2,:)];
        if (length(cell_order)==1)||(min(vnorm(diff))>2*cellRadius)
            if (new_x>cellRadius)&&(new_x<domainLength-cellRadius)&&(new_y>cellRadius)&&(new_y<domainHeight-cellRadius)&&(new_x>min(xlat)+cellRadius)
                % if we have inserted a barrier and the cell trys
                % to move through it, stop them
                if ((experiment==4)||(experiment==5))&&(t_save>in.barrier_time)
                    previous_side = cells(1,cellidx) + (-1)^max((cells(1,cellidx)-barrier)/abs(cells(1,cellidx)-barrier),0)*cellRadius;
                    new_side = new_x+ (-1)^max((new_x-barrier)/abs(new_x-barrier),0)*cellRadius;
                    if (barrier-previous_side)*(barrier-new_side)>0
                        cells(1,cellidx) = new_x;
                        cells(2,cellidx) = new_y;
                    end
                else
                    cells(1,cellidx) = new_x;
                    cells(2,cellidx) = new_y;
                    if (move==0)&&(cellsFollow(cellidx)==1)
                        disp('follower moved anyway')
                    elseif move==0
                        disp('leader moved anyway')
                    end
                end
            end
        end
    end
    if min(sqrt((cells(1,cellidx)-cells(1,(1:end)~=cellidx)).^2 + (cells(2,cellidx)-cells(2,(1:end)~=cellidx)).^2))<=(2*cellRadius)
        disp('error - cells moving through each other!')
        pause
    end
end

%% save stuff
out.attach = attach;
out.cellsFollow = cellsFollow;
out.filopodia = filopodia;
out.theta = theta;
out.cells = cells;
out.moved = moved;