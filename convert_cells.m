function out = convert_cells(cells,cellsFollow,attach_save,timeCtr,cells_save,filolength,moved,happiness,ca_save,xlat,ylat,eatWidth,filopodia,conversionType,param,...
    num_better_foll_save,num_foll_save,num_better_lead_save,num_lead_save,numFilopodia)

if conversionType==1
    %% using amount of time not found a gradient / another cell
    numSteps = param(15); % number of steps to not sucessfully find a direction, before changing roles (convert type 1)

    for cellCtr=1:length(cells(1,:))
        [~,numCellsInTime] = cellfun(@size,cells_save);
        if (timeCtr>numSteps)&&(numCellsInTime(timeCtr-numSteps)>=cellCtr)
            if ~any(moved(timeCtr-numSteps:timeCtr,cellCtr)) % if none moved -- LJS
                if cellsFollow(cellCtr)==1
                    disp('foll->lead')
                    cellsFollow(cellCtr)=0;
                else
                    disp('lead->foll')
                    cellsFollow(cellCtr)=1;
                end
                moved(timeCtr,cellCtr)=1;
            end
        end
    end
    
elseif conversionType==2
    %% Using the presence of a c'tant gradient (with integral measures of c'tant)
    for cellCtr=1:length(cells(1,:))
        if cellsFollow(cellCtr)==1 % if it's a follower
            numSteps = numFilopodia(2); % number of directions to sample in (convert type 2)
            numDirections = 1/numSteps; % fraction of directions that need to be better for phenotype switch -- LJS
        elseif cellsFollow(cellCtr)==0 % if it's a leader
            numSteps = numFilopodia(1); % number of directions to sample in (convert type 2)
            numDirections = 1/numSteps; % fraction of directions that need to be better for phenotype switch -- LJS
        end
        r = rand()*2*pi;
        theta = (2*pi/numSteps:2*pi/numSteps:2*pi) + r; % this seems to sample evenly distributed with only a random offset... is that what we want? -- LJS
        [~,~,~,~,num_better] = cell_movement5(theta,cells(1,cellCtr),cells(2,cellCtr),ca_save,xlat,ylat,eatWidth,filolength,numSteps,[]);
        
        if num_better>=(numDirections*numSteps)
            %             if (rand()<0.7)
            if cellsFollow(cellCtr)==1
                disp('follow -> lead')
                cellsFollow(cellCtr)=0;
            end
            %             else
            %                 cellsFollow(i) = 1;
            %             end
        elseif num_better<(numDirections*numSteps)
            if cellsFollow(cellCtr)==0
                disp('lead -> follow')
                cellsFollow(cellCtr)=1;
            end
            if cellsFollow(cellCtr)==1
                num_better_foll_save = [num_better_foll_save, num_better/numSteps];
                num_foll_save = num_foll_save +1;
            else
                num_better_lead_save = [num_better_lead_save, num_better/numSteps];
                num_lead_save = num_lead_save +1;
            end
            %         elseif rand()<0.5
            %             cellsFollow(i) = 1-cellsFollow(i);
        end
    end
elseif conversionType==3 %% Conversion type 3 was because Ruth kept asking if we couldn't just set some concentration of local chemoattractant at which cells would convert between the types. It doesn't really work, because the overall levels are being diluted as the domain expands, so you can't set just one threshold. (LJS: but what if hte consumption was high enough so that dilution wasn't much of an issue?
    %% using maximum chemoattractant gradient (with point measures of c'tant)
    dx = 10;
    temp_thet = 0:0.1*pi:2*pi;
    for cellCtr=1:length(cells(1,:))
        temp = [cells(1,cellCtr) + dx*cos(temp_thet); cells(2,cellCtr) + dx*sin(temp_thet)];
        temp = [cells(:,cellCtr) temp];
        
        %         chemo = find_ca(temp,xlat_save,ylat_save,ca_save);
        chemo = find_ca(temp,xlat,ylat,ca_save); % ?? -- LJS
        if chemo(1)==0
            grad = sign(chemo(2:end)-chemo(1)).*max(abs(chemo(2:end)-chemo(1)))/dx;
        else
            grad = sign(chemo(2:end)-chemo(1)).*max(abs(chemo(2:end)-chemo(1)))/(dx*chemo(1));
        end
        if cellsFollow(cellCtr)==1
            temp_grad(cellCtr,timeCtr) = grad(1);
            if (timeCtr>100)
                [~,tempj] = find(abs(temp_grad)==max(abs(temp_grad(cellCtr,timeCtr-100:timeCtr))));
                %                     temp = temp_grad(i,k-100:k);
                %                     temp = temp(temp~=0);
                follow_grad(cellCtr,timeCtr) = temp_grad(cellCtr,tempj(1));
            end
            %                 if (k>10)&&(follow_grad(i,k)>1)
            %                     disp('follow -> lead')
            %                     cellsFollow(i)=0;
            %                 end
        else
            temp_grad(cellCtr,timeCtr) = grad(1);
            if timeCtr>100
                [~,tempj] = find(abs(temp_grad)==max(abs(temp_grad(cellCtr,timeCtr-100:timeCtr))));
                %                     temp = temp_grad(i,k-100:k);
                %                     temp = temp(temp~=0);
                leader_grad(cellCtr,timeCtr) = temp_grad(cellCtr,tempj(1));
                %                     leader_grad(i,k) = sign(temp_grad(i,k-100:k))*max(abs(temp_grad(i,k-100:k)));
            end
            %                 if (k>10)&&(leader_grad(i,k)<0.001)
            %                     disp('lead -> follow')
            %                     cellsFollow(i)=1;
            %                 end
        end
        if grad>0.8
            cellsFollow(cellCtr) = 0;
        elseif grad<0.05
            cellsFollow(cellCtr) = 1;
        end
        %             pause
    end
elseif conversionType == 4
     %% integrate-and-switch, or time-frustration with hysteresis -- LJS
    numSteps = param(15); % number of steps between happy and unhappy (threshold to switch) -- LJS
    cells2switch = happiness(timeCtr,:)<=0; % unhappy cells to switch type
    if any(cells2switch)
        cellsFollow(cells2switch) = ~cellsFollow(cells2switch); % switch cell type
        happiness(timeCtr,cells2switch) = numSteps; % cells that switch become maximally happy (hysteresis) -- LJS
        moved(timeCtr,cells2switch) = 1;
        
        if any(cellsFollow(cells2switch)), disp([num2str(nnz(cellsFollow(cells2switch))) 'lead->foll']); end
        if any(~cellsFollow(cells2switch)), disp([num2str(nnz(~cellsFollow(cells2switch))) 'foll->lead']); end
    end
    
end
out.num_better_foll_save = num_better_foll_save;
out.num_foll_save = num_foll_save;
out.num_better_lead_save = num_better_lead_save;
out.num_lead_save = num_lead_save;
out.cellsFollow = cellsFollow;
out.moved = moved;
out.happiness = happiness;

