
%%  (only needed for backwards compatibility)  save stuff for tissue transplant simulations %%
if (param.experiment==1)||(param.experiment==2)
    if in.it==1
        %% make a plot and allow user to select the region to extract
        make_plot(cells, cellsFollow, xlat_save{saveCtr-1},ylat_save{saveCtr-1},t_save(saveCtr-1),ca_save{saveCtr-1},filopodia,attach_save{saveCtr-1},cellRadius,filolength,sensingAccuracy,0,caCmap,0,param,[])
        sort(cells(1,:))
        
        if param.experiment==1
            temp = find(cellsFollow==1)';
            sort(cells(1,temp(temp<length(cells(1,:)))))
            x_start = input('where should the extracted region start (x_start) [front region], (x_end=x_start+100)?:  ');
            x_end = x_start + 100;
            y_start =  in.y_start;   % start of the region to be extracted _in um_
            y_end = in.y_end;       % end of the region to be extracted _in um_
        else
            temp = find(cellsFollow==1)';
            sort(cells(1,temp(temp<length(cells(1,:)))))
            x_start = input('where should the extracted region start (x_start) [back region]?:  ');
            x_end = input('where should the extracted region end (x_end) [back region]?:  ');
            y_start = input('where should the extracted region start (y_start)?:  ');% start of the region to be extracted _in um_
            y_end = input('where should the extracted region end (y_end)?:  ');% end of the region to be extracted _in um_
        end
        
        % interpolate the chemoattractant to give a good description of the region
        xlat_new_fine = (0:5:max(xlat_save{end}));
        if (param.experiment==1)||((y_start==0)&&(y_end==param.domainHeight))
            ca_new = interp2(xlat_save{end},ylat_save{end},ca_save{end}',xlat_new_fine,ylat_save{end});
            ca_new = ca_new';
            extracted_xlat = find((xlat_new_fine>=x_start)&(xlat_new_fine<=x_end));
            extracted_ylat = 1:length(ylat_save{end});
            extracted_x = xlat_new_fine(extracted_xlat);
            extracted_y = ylat_save{end};
        else
            ylat_new_fine = (0:5:param.domainHeight);
            ca_new = interp2(xlat_save{end},ylat_save{end},ca_save{end}',xlat_new_fine,ylat_new_fine');
            extracted_xlat = find((xlat_new_fine>=x_start)&(xlat_new_fine<=x_end));
            extracted_ylat = find((ylat_new_fine>=y_start)&(ylat_new_fine<=y_end));
            extracted_y = ylat_new_fine(extracted_ylat);
        end
        
        % save the cell positions and filopodia (cells that are across the cut are not transplanted)
        extracted_cell_ind = find((cells_save{end}(1,:)>=x_start+cellRadius)&(cells_save{end}(1,:)<=x_end-cellRadius)...
            &(cells_save{end}(2,:)>=y_start+cellRadius)&(cells_save{end}(2,:)<=y_end-cellRadius));
        
        out.cells_save = cells_save{end}(:,extracted_cell_ind);
        out.filopodia_save = filopodia_save{end}(extracted_cell_ind,:);
        out.cellsFollow_save = cellsFollow_save{end}(extracted_cell_ind);
        out.attach = attach(extracted_cell_ind) -(length(cells_save{end}(1,:)) - length(out.cells_save(1,:)));
        
        % save the interpolated xlat, ylat and chemoattractant
        out.xlat_save = xlat_new_fine(extracted_xlat);
        out.ylat_save = extracted_y;
        out.ca_save = ca_new(extracted_xlat,extracted_ylat);
        out.extracted_x_start = x_start;
        out.extracted_x_end = x_end;
        out.extracted_y_start = y_start;
        out.extracted_y_end = y_end;
        makeMovies=0;
    elseif in.it==2
        makeMovies=1;
    end
end
%% saveInfo included in the naming of files %%%
if isstruct(in)&&ismember('saveInfo',fields(in))
    saveInfo = in.saveInfo;
elseif ~exist('saveInfo','var')
    if conversionType~=0
        if conversionType==4
            saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_exp_',num2str(param.experiment),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),'_steps_',num2str(numSteps(1)),'_',num2str(numSteps(2)),...
                '_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
        else
            saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_exp_',num2str(param.experiment),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),'_steps_',num2str(numSteps),...
                '_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
        end
    else
        saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_exp_',num2str(param.experiment),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),...
            '_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
    end
end
%% save the simulation results %%
if (param.experiment==0||param.experiment==3||param.experiment>=10)||(((param.experiment==1)||(param.experiment==2))&&(in.it~=1))
    % downsample as appropriate
    if saveEvery>1
        t_save = t_save(1:saveEvery:end);
        xlat_save = xlat_save(1:saveEvery:end); % spatial lattices (lat)
        ylat_save = ylat_save(1:saveEvery:end);
        ca_save = ca_save(1:saveEvery:end); % chemoattractant (ca)
        if ~isempty(dan) % dan tunneling simulation
            dan_save = dan_save(1:saveEvery:end);
        end
        if exist('happiness','var')
            happiness = happiness(1:saveEvery:end);
        end
        cells_save = cells_save(1:saveEvery:end);
        filopodia_save = filopodia_save(1:saveEvery:end);
        cellsFollow_save = cellsFollow_save(1:saveEvery:end);
        attach_save = attach_save(1:saveEvery:end);
        moved = moved(1:saveEvery:end,:);
        domainLengths = domainLengths(1:saveEvery:end);
    end
    numSavepoints = ceil(numTsteps/saveEvery);
    out.saveEvery = saveEvery;
    % convert floats to single precision for saving, to reduces disk space
    % used
    out.t_save = t_save;
    for saveCtr=1:numSavepoints
        xlat_save{saveCtr} = single(xlat_save{saveCtr});
        ylat_save{saveCtr} = single(ylat_save{saveCtr});
        ca_save{saveCtr} = single(ca_save{saveCtr});
        cells_save{saveCtr} = single(cells_save{saveCtr});
        filopodia_save{saveCtr} = single(filopodia_save{saveCtr});
        if exist('happiness','var')
            happiness = single(happiness);
        end
        if ~all(islogical(cellsFollow_save{saveCtr}))
            cellsFollow_save{saveCtr} = single(cellsFollow_save{saveCtr});
        end
    end
    out.xlat_save = xlat_save;
    out.ylat_save = ylat_save;
    out.ca_save = ca_save;
    out.dan_save = dan_save;
    out.cells_save = cells_save;
    out.filopodia_save = filopodia_save;
    out.domainLengths = domainLengths;
    out.saveInfo = saveInfo;
    out.numTsteps = numTsteps;
    out.numSavepoints = numSavepoints;
    out.growingDomain = param.growingDomain;
    out.followerFraction = followerFraction;
    out.tstep = param.tstep;
    out.cellsFollow_save = cellsFollow_save;
    out.cellRadius = cellRadius;
    out.param.domainHeight = param.domainHeight;
    out.filolength = filolength;
    out.attach_save = attach_save;
    out.param_names1 = 'Linf, a, diffus,e, growingDomain,initialDomainLength,makeChemoattractant';
    out.param_names2 = 'chi,param.domainHeight,zeroBC,insert,tstep,t_s,d,eatRate';
    out.param = param; %param = [Linf, a, diffus, eatWidth, growingDomain, initialDomainLength, makeChemoattractant, chi, param.domainHeight, zeroBC, insert, tstep, t_s, eatRate, numSteps, numDirections];
    out.moved = moved(:,1:numCellsFinal);
    if exist('happiness','var')
        out.happiness = happiness(:,1:numCellsFinal);
    end
    out.leadSpeed = leadSpeed;
    out.followSpeed = followSpeed;
    out.numFilopodia = numFilopodia;
    out.sensingAccuracy = sensingAccuracy;
    out.followerFraction = followerFraction;
    out.divide_cells = divide_cells;
    out.param.experiment = param.experiment;
    [~, computerName] = system('hostname -s');
    computerName = computerName(1:end-1); % remove newline character
    save(['results/',saveInfo,'.mat'],'out')
    system(['rm -f results/',saveInfo,'_running_on_',computerName,'.mat'])
    fprintf(['saved results to results/',saveInfo,'.mat \n'])
    
end
