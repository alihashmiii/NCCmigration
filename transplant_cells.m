if ((param.experiment==1)||(param.experiment==2))&&(in.it==2)&&(t_save(timeCtr)==in.changeTime)
    
    fprintf(['experiment = ',mat2str(param.experiment),'\r'])
    %% Ablations
    % Ablation type 1 deletes the cells backwards or forwards from
    % a point
    if (param.experiment==1)&&(in.ablate_type==1)
        %% experiment 1
        exp.experiment = param.experiment;
        make_plot(cells, cellsFollow, xlat_save{timeCtr-1},ylat_save{timeCtr-1},t_save(timeCtr-1),ca_save{timeCtr-1},filopodia,attach_save{timeCtr-1},cellRadius,filolength,sensingAccuracy,0,caCmap,0,param,[])
        sort(cells(1,:))
        x_start = input('where should we ablate back from?: ');
        
        ablate_out = ablate(cells,attach,x_start,'back',cellsFollow);
        cells = ablate_out.cells;
        attach = ablate_out.attach;
        cellsFollow = ablate_out.cellsFollow;
        y_start = 0;

    elseif (param.experiment==2)
        %% experiment 2
        exp.experiment = param.experiment;
        disp('up to date')
        make_plot(cells, cellsFollow, xlat_save{timeCtr-1},ylat_save{timeCtr-1},t_save(timeCtr-1),ca_save{timeCtr-1},filopodia,attach_save{timeCtr-1},cellRadius,filolength,sensingAccuracy,0,caCmap,0,param,[])
        sort(cells(1,:))
        disp(['the inserted region is ',mat2str(in.extracted_x_end-in.extracted_x_start),'um long'])
        x_start = input('where should we insert the region (x_start) [front region]?: ');
        disp(['the inserted region is ',mat2str(in.extracted_y_end-in.extracted_y_start),'um high'])
        y_start = input('where should we insert the region (y_start): ');
        
        if (in.ablate_type==1)
            ablate_out = ablate(cells,attach,x_start,'forwards',cellsFollow);
            cells = ablate_out.cells;
            attach = ablate_out.attach;
            cellsFollow = ablate_out.cellsFollow;
            
            cells(:,cells(1,:)>x_start)=[];    % take out the front cells
        end
    end
    %% Insertions
    %% experiments 1 and 2

    exp.experiment = param.experiment;
    make_plot(cells, cellsFollow, xlat_save{timeCtr-1},ylat_save{timeCtr-1},t_save(timeCtr-1),ca_save{timeCtr-1},filopodia,...
        attach_save{timeCtr-1},cellRadius,filolength,sensingAccuracy,0,caCmap,0,param,[])
    sort(cells(1,:))
    disp(['the inserted region is ',mat2str(in.extracted_x_end-in.extracted_x_start),'um long'])
    if param.experiment==2
        insert_x = x_start;
    else
        insert_x = -100;
    end
    if param.experiment==1
        y_start = 0;
    end
    insert_y = y_start;
    
    cells_in = in.cells_in;
    cells_in(1,:) = cells_in(1,:) - in.extracted_x_start + insert_x;    % move the inserted cell region to start at insert_x
    cells_in(2,:) = cells_in(2,:) - in.extracted_y_start + insert_y;    % move the inserted cell region to start at insert_y
    
    filopodia_in = in.filopodia_in;
    filopodia_in(:,1) = filopodia_in(:,1) - in.extracted_x_start + insert_x;
    filopodia_in(:,2) = filopodia_in(:,2) - in.extracted_y_start + insert_y;
    
    cellsFollow(length(cells(1,:))+1:length(cells(1,:))+length(cells_in(1,:)))=in.cellsFollow;
    cells = [cells cells_in];
    filopodia = [filopodia; filopodia_in];
    for i=1:length(in.attach)
        if in.attach(i)==0
            attach(end+1) = in.attach(i);
        elseif in.attach(i)<length(in.cellsFollow)
            attach(end+1) = in.attach(i)+length(cellsFollow);
        else
            attach(end+1) = 0;
        end
    end
    disp('inserted cells')
elseif (param.experiment==3)&&(t_save(timeCtr)==in.ablate_time)&&(in.ablate_type==1)
    %% experiment 3
    % Ablation type 1 deletes the cells backwards or forwards from
    % a point
    ablate_out = ablate(cells,attach,in.ablate_start,'back',cellsFollow);
    cells = ablate_out.cells;
    attach = ablate_out.attach;
    cellsFollow = ablate_out.cellsFollow;
end
