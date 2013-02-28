function out = convert_cells(cells,cellsFollow,attach_save,k,cells_save,filolength,moved,ca_save,xlat,ylat,d,filopodia,convert_type,param,...
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
    n = param(15);
    [size_attach_save,~] = cellfun(@size,attach_save(1:k)); % don't consider timesteps that haven't happened yet
    change = find(size_attach_save~=size_attach_save(1),1,'first');
    if isempty(change)
        temp_attach_save = cell2mat(attach_save);
    else
        temp_attach_save = cell2mat(attach_save(change:end)); % if an experiment changed the size of attach_save, only consider after that
    end
    for i=1:length(cells(1,:))
        [~,num_cells_k] = cellfun(@size,cells_save);
        if (k>n)&&(num_cells_k(k-n)>=i)
            moved_in_n = sum(moved(k-n:k,i));
            if (moved_in_n==0)
                if cellsFollow(i)==1
                    disp('foll->lead')
                    cellsFollow(i)=0;
                else
                    disp('lead->foll')
                    cellsFollow(i)=1;
                end
                moved(k,i)=1;
            end
        end
    end
    % out.cellsFollow = cellsFollow;
    % out.moved = moved;
    
elseif convert_type==2
    %% Using the presence of a c'tant gradient (with integral measures of c'tant)

    n = param(15);
    m = param(16);
    r = rand()*2*pi;
    theta = (2*pi/n:2*pi/n:2*pi) + r;
    for i=1:length(cells(1,:))
        [~,~,~,~,num_better] = cell_movement5(theta,cells(1,i),cells(2,i),ca_save,xlat,ylat,d,filolength,n,[]);

        if num_better>=(m*n)
%             if (rand()<0.7)
                if cellsFollow(i)==1
                    disp('follow -> lead')
                    cellsFollow(i)=0;
                end
%             else
%                 cellsFollow(i) = 1;
%             end
        elseif num_better<(m*n)
            if cellsFollow(i)==0
                disp('lead -> follow')
                cellsFollow(i)=1;
            end
        if cellsFollow(i)==1
            num_better_foll_save = [num_better_foll_save, num_better/n];
            num_foll_save = num_foll_save +1;
        else
            num_better_lead_save = [num_better_lead_save, num_better/n];
            num_lead_save = num_lead_save +1;
        end
%         elseif rand()<0.5
%             cellsFollow(i) = 1-cellsFollow(i);
        end
    end
elseif convert_type==3
    %% using maximum chemoattractant gradient (with point measures of c'tant)
    dx = 10;
    temp_thet = 0:0.1*pi:2*pi;
    for i=1:length(cells(1,:))
        temp = [cells(1,i) + dx*cos(temp_thet); cells(2,i) + dx*sin(temp_thet)];
        temp = [cells(:,i) temp];
    
        chemo = find_ca(temp,xlat_save,ylat_save,ca_save);
        if chemo(1)==0
            grad = sign(chemo(2:end)-chemo(1)).*max(abs(chemo(2:end)-chemo(1)))/dx;
        else
            grad = sign(chemo(2:end)-chemo(1)).*max(abs(chemo(2:end)-chemo(1)))/(dx*chemo(1));
        end
        if cellsFollow(i)==1
            temp_grad(i,k) = grad(1);
            if (k>100)
                [~,tempj] = find(abs(temp_grad)==max(abs(temp_grad(i,k-100:k))));
                %                     temp = temp_grad(i,k-100:k);
                %                     temp = temp(temp~=0);
                follow_grad(i,k) = temp_grad(i,tempj(1));
            end
            %                 if (k>10)&&(follow_grad(i,k)>1)
            %                     disp('follow -> lead')
            %                     cellsFollow(i)=0;
            %                 end
        else
            temp_grad(i,k) = grad(1);
            if k>100
                [~,tempj] = find(abs(temp_grad)==max(abs(temp_grad(i,k-100:k))));
                %                     temp = temp_grad(i,k-100:k);
                %                     temp = temp(temp~=0);
                leader_grad(i,k) = temp_grad(i,tempj(1));
                %                     leader_grad(i,k) = sign(temp_grad(i,k-100:k))*max(abs(temp_grad(i,k-100:k)));
            end
            %                 if (k>10)&&(leader_grad(i,k)<0.001)
            %                     disp('lead -> follow')
            %                     cellsFollow(i)=1;
            %                 end
        end
                    if grad>0.8
                        cellsFollow(i) = 0;
                    elseif grad<0.05
                        cellsFollow(i) = 1;
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
