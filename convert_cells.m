function out = convert_cells(cells,cellsFollow,attach_save,timeCtr,cells_save,filolength,moved,ca_save,xlat,ylat,d,filopodia,convert_type,param,...
    num_better_foll_save,num_foll_save,num_better_lead_save,num_lead_save)

%% using quorum sensing
% R = 7.5;
% for i=1:length(cells(1,:))
%     % find out how many cells there are within filopodia reach
%     num_cells = sum(sqrt((cells(1,i)-cells(1,(1:end)~=i)).^2 + (cells(2,i)-cells(2,(1:end)~=i)).^2)<=filolength);
%
%     % if there are many, and it was a leader, become a follower
%     if (cellsFollow(i)==1) && (num_cells<3)
%         disp('foll->lead')
%         cellsFollow(i)=0;
%         num_cells
% %     elseif (cellsFollow(i)==0) && (num_cells>15)
% %         disp('lead->foll')
% %         cellsFollow(i)=1;
% %         num_cells
%     end
% end

if convert_type==1
    %% using amount of time not found a gradient / another cell
    num_steps = param(15); % number of steps to not sucessfully find a direction, before changing roles (convert type 1)
    [size_attach_save,~] = cellfun(@size,attach_save(1:timeCtr)); % don't consider timesteps that haven't happened yet
    change = find(size_attach_save~=size_attach_save(1),1,'first');
    if isempty(change)
        temp_attach_save = cell2mat(attach_save);
    else
        temp_attach_save = cell2mat(attach_save(change:end)); % if an experiment changed the size of attach_save, only consider after that
    end
    for cellCtr=1:length(cells(1,:))
        [~,num_cells_k] = cellfun(@size,cells_save);
        if (timeCtr>num_steps)&&(num_cells_k(timeCtr-num_steps)>=cellCtr)
            if ~any(moved(timeCtr-num_steps:timeCtr,cellCtr)) % if none moved -- LJS
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
    % out.cellsFollow = cellsFollow;
    % out.moved = moved;
    
elseif convert_type==2
    %% Using the presence of a c'tant gradient (with integral measures of c'tant)

    num_steps = param(15); % number of directions to sample in (convert type 2)
    m = param(16);
    r = rand()*2*pi;
    theta = (2*pi/num_steps:2*pi/num_steps:2*pi) + r;
    for cellCtr=1:length(cells(1,:))
        [~,~,~,~,num_better] = cell_movement5(theta,cells(1,cellCtr),cells(2,cellCtr),ca_save,xlat,ylat,d,filolength,num_steps,[]);

        if num_better>=(m*num_steps)
%             if (rand()<0.7)
                if cellsFollow(cellCtr)==1
                    disp('follow -> lead')
                    cellsFollow(cellCtr)=0;
                end
%             else
%                 cellsFollow(i) = 1;
%             end
        elseif num_better<(m*num_steps)
            if cellsFollow(cellCtr)==0
                disp('lead -> follow')
                cellsFollow(cellCtr)=1;
            end
        if cellsFollow(cellCtr)==1
            num_better_foll_save = [num_better_foll_save, num_better/num_steps];
            num_foll_save = num_foll_save +1;
        else
            num_better_lead_save = [num_better_lead_save, num_better/num_steps];
            num_lead_save = num_lead_save +1;
        end
%         elseif rand()<0.5
%             cellsFollow(i) = 1-cellsFollow(i);
        end
    end
elseif convert_type==3 %% Conversion type 3 was because Ruth kept asking if we couldn't just set some concentration of local chemoattractant at which cells would convert between the types. It doesn't really work, because the overall levels are being diluted as the domain expands, so you can't set just one threshold. (LJS: but what if hte consumption was high enough so that dilution wasn't much of an issue?
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
end
out.num_better_foll_save = num_better_foll_save;
out.num_foll_save = num_foll_save;
out.num_better_lead_save = num_better_lead_save;
out.num_lead_save = num_lead_save;
out.cellsFollow = cellsFollow;
out.moved = moved;
