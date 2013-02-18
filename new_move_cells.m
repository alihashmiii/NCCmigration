function out = new_move_cells(cells,cells_follow,filopodia,attach,theta,...
                            ca_save,xlat,ylat,...
                            cell_radius, filolength, d, height, dist, domain_length, barrier, experiment, t_save, in, metropolis, num_filopodia, diffusion)
% dettached = 0;
%% iterate through the cell movement in a random order %%%
cell_order = randperm(length(cells(1,:)));
moved = zeros(1,length(cells(1,:)));
for i =1:length(cell_order)
    move = 0;
    cellidx = cell_order(i);  % look at the rth cell
    other_cells = cells(:,(1:end)~=cellidx);
    
    %% Decide whether to try to move
    if cells_follow(cellidx)~=1   
        %% if it's a leader
        if diffusion==0
            [filopodia(cellidx,:),E,move,theta(cellidx),~] = cell_movement5(rand(1,num_filopodia(1))*2*pi,cells(1,cellidx),cells(2,cellidx),ca_save,xlat,ylat,...
                d,filolength,num_filopodia(1),[]);
        elseif rand(1,1)<diffusion
            move = 1;
        else
            move = 0;
        end            
    elseif attach(cellidx)~=0     
        %% if it's a chained follower that can reach the cell ahead
        if (cells(1,attach(cellidx)) - cells(1,cellidx))^2 + (cells(2,attach(cellidx)) - cells(2,cellidx))^2 < (filolength + cell_radius)^2
%             % set angle of movement towards cell being followed -- LJS
%             theta(cellidx) = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx)));
            % set angle of movement parallel to that of cell being followed
            theta(cellidx) = theta(attach(cellidx));
            % set filopodium position to closet point on membrane of cell being followed -- LJS
            phi = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx))); % the angle towards the cell being followed -- LJS
            filopodia(cellidx,1) = cells(1,attach(cellidx)) - cell_radius*cos(phi);
            filopodia(cellidx,2) = cells(2,attach(cellidx)) - cell_radius*sin(phi);
            move = 1;
        else
            %% if the cell ahead is too far, then dettach the chain
            % set filopodial position and movement angle in the direction of cell
            % centre of lost cell -- LJS
            theta(cellidx) = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx)));
            filopodia(cellidx,:) = cells(:,cellidx)' + filolength.*[cos(theta(cellidx)) sin(theta(cellidx))];
            attach = dettach(cellidx,attach);
        end
        
    else
        %% if it's an unchained follower
        theta(cellidx) = rand()*2*pi;
        [c,filopodia(cellidx,:)] = cell_movement5_follow(theta(cellidx),cellidx,cells(1,:),cells(2,:),cell_radius,...
            filolength,filopodia,barrier,experiment);
        if isempty(c)~=1
            %% if another cell was found then find the head of that chain
            head = c;
            while attach(head)~=0
                head = attach(head);
            end
            if cells_follow(head)==0    %% if the head is a leader, then attach and move
                attach(cellidx) = c;
                % set angle of movement parallel to that of cell being followed
                theta(cellidx) = theta(attach(cellidx));
                % set filopodium position to closet point on membrane of cell being followed -- LJS
                phi = atan2((cells(2,attach(cellidx)) - cells(2,cellidx)),(cells(1,attach(cellidx)) - cells(1,cellidx))); % the angle towards the cell being followed -- LJS
                filopodia(cellidx,1) = cells(1,attach(cellidx)) - cell_radius*cos(phi);
                filopodia(cellidx,2) = cells(2,attach(cellidx)) - cell_radius*sin(phi);
                move = 1;
            end
        end
    end
    
    %% Try to move
    r1 = rand();
    T = 0.0001; %%% ???? -- LJS
    if (cells_follow(cellidx)==1)&&(move==0)
        x_fil = cells(1,cellidx)+ filolength*cos(theta(cellidx));
        y_fil = cells(2,cellidx)+ filolength*sin(theta(cellidx));
        E = min(sqrt((x_fil-other_cells(1,:)).^2+(y_fil-other_cells(2,:)).^2));
        T = 50;
%     elseif (cells_follow(r)==0)&&(move==0)
%         exp(-E/T)
    end
    if (move==1)%||((metropolis==1)&&(r1<exp(-E/T))&&(cells_follow(r)==1))
        moved(cellidx)=1;
    end
    if (move==1)||((metropolis==1)&&(r1<exp(-E/T)))%&&(cells_follow(r)==1))
        if (cells_follow(cellidx)==1)&&(move==1)
                new_x = cells(1,cellidx) + cos(theta(cellidx))*dist;
                new_y = cells(2,cellidx) + sin(theta(cellidx))*dist;
        else
            new_x = cells(1,cellidx) + cos(theta(cellidx))*dist;
            new_y = cells(2,cellidx) + sin(theta(cellidx))*dist;
        end
        
        %% if there is no cell or edge in the way, then move
        diff = [new_x-other_cells(1,:); new_y-other_cells(2,:)];
        if (length(cell_order)==1)||(min(vnorm(diff))>2*cell_radius)
            if (new_x>cell_radius)&&(new_x<domain_length-cell_radius)&&(new_y>cell_radius)&&(new_y<height-cell_radius)&&(new_x>min(xlat)+cell_radius)
                % if we have inserted a barrier and the cell trys
                % to move through it, stop them
                if ((experiment==4)||(experiment==5))&&(t_save>in.barrier_time)
                    previous_side = cells(1,cellidx) + (-1)^max((cells(1,cellidx)-barrier)/abs(cells(1,cellidx)-barrier),0)*cell_radius;
                    new_side = new_x+ (-1)^max((new_x-barrier)/abs(new_x-barrier),0)*cell_radius;
                    if (barrier-previous_side)*(barrier-new_side)>0
                        cells(1,cellidx) = new_x;
                        cells(2,cellidx) = new_y;
%                         moved(r) = 1;
                    end
                else
                    cells(1,cellidx) = new_x;
                    cells(2,cellidx) = new_y;
%                     moved(r) = 1;
                    if (move==0)&&(cells_follow(cellidx)==1)
                        disp('follower moved anyway')
                    elseif move==0
                        disp('leader moved anyway')
                    end
                end
            end
        end
    end
    if min(sqrt((cells(1,cellidx)-cells(1,(1:end)~=cellidx)).^2 + (cells(2,cellidx)-cells(2,(1:end)~=cellidx)).^2))<=(2*cell_radius)
        disp('error - cells moving through each other!')
        pause
    end
end

%% save stuff
out.attach = attach;
out.cells_follow = cells_follow;
out.filopodia = filopodia;
out.theta = theta;
out.cells = cells;
out.moved = moved;