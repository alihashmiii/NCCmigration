
%% save stuff for tissue transplant experiments %%
if (experiment==1)||(experiment==2)
    if in.it==1
        %% make a plot and allow user to select the region to extract
        make_plot(cells, cellsFollow, xlat_save{k-1},ylat_save{k-1},ca_save{k-1},filopodia,attach_save{k-1},cellRadius,0,0)
        sort(cells(1,:))
        
        if experiment==1
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
        
        %% interpolate the chemoattractant to give a good description of
        %% the region
        xlat_new_fine = (0:5:max(xlat_save{end}));
        if (experiment==1)||((y_start==0)&&(y_end==domainHeight))
            ca_new = interp2(xlat_save{end},ylat_save{end},ca_save{end}',xlat_new_fine,ylat_save{end});
            ca_new = ca_new';
            extracted_xlat = find((xlat_new_fine>=x_start)&(xlat_new_fine<=x_end));
            extracted_ylat = 1:length(ylat_save{end});
            extracted_x = xlat_new_fine(extracted_xlat);
            extracted_y = ylat_save{end};
        else
            ylat_new_fine = (0:5:domainHeight);
            ca_new = interp2(xlat_save{end},ylat_save{end},ca_save{end}',xlat_new_fine,ylat_new_fine');
            extracted_xlat = find((xlat_new_fine>=x_start)&(xlat_new_fine<=x_end));
            extracted_ylat = find((ylat_new_fine>=y_start)&(ylat_new_fine<=y_end));
            extracted_y = ylat_new_fine(extracted_ylat);
        end
        
        %% save the relevant data
        %% save the cell positions and filopodia (cells that are across the cut are not transplanted)
        extracted_cell_ind = find((cells_save{end}(1,:)>=x_start+cellRadius)&(cells_save{end}(1,:)<=x_end-cellRadius)...
            &(cells_save{end}(2,:)>=y_start+cellRadius)&(cells_save{end}(2,:)<=y_end-cellRadius));
        
        out.cells_save = cells_save{end}(:,extracted_cell_ind);
        out.filopodia_save = filopodia_save{end}(extracted_cell_ind,:);
        out.cellsFollow = cellsFollow_save{end}(extracted_cell_ind);
        out.attach = attach(extracted_cell_ind) -(length(cells_save{end}(1,:)) - length(out.cells_save(1,:)));

        %% save the interpolated xlat, ylat and chemoattractant
        out.xlat_save = xlat_new_fine(extracted_xlat);
        out.ylat_save = extracted_y;
        out.ca_save = ca_new(extracted_xlat,extracted_ylat);
        out.extracted_x_start = x_start;
        out.extracted_x_end = x_end;
        out.extracted_y_start = y_start;
        out.extracted_y_end = y_end;
        movies=0;
    elseif in.it==2
        movies=1;
    end
end
%%% saveInfo included in the naming of files %%%
if isstruct(in)&&ismember('saveInfo',fields(in))
    saveInfo = in.saveInfo;
elseif ~exist('saveInfo','var')
    if conversionType~=0
        saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),'_steps_',num2str(num_steps),...
        '_numleadfil_',mat2str(numFilopodia(1)),'_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
    else
            saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),...
        '_numleadfil_',mat2str(numFilopodia(1)),'_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
    end
end
%% save the results %%
if (experiment==0||experiment==3)||(((experiment==1)||(experiment==2))&&(in.it~=1))
    % convert floats to single precision for saving, to reduces disk space
    % used
    out.t_save = t_save;
    for timeCtr=1:numTsteps
        xlat_save{timeCtr} = single(xlat_save{timeCtr});
        ylat_save{timeCtr} = single(ylat_save{timeCtr});
        ca_save{timeCtr} = single(ca_save{timeCtr});
        cells_save{timeCtr} = single(cells_save{timeCtr});
        filopodia_save{timeCtr} = single(filopodia_save{timeCtr});
    end
    out.xlat_save = xlat_save;
    out.ylat_save = ylat_save;
    out.ca_save = ca_save;
    out.cells_save = cells_save;
    out.filopodia_save = filopodia_save;
    out.domainLengths = domainLengths;
    out.saveInfo = saveInfo;
    out.numTsteps = numTsteps;
    out.growingDomain = growingDomain;
    out.followerFraction = followerFraction;
    out.tstep = tstep;
    out.cellsFollow = cellsFollow_save;
    out.cellRadius = cellRadius;
    out.domainHeight = domainHeight;
    out.filolength = filolength;
    out.attach_save = attach_save;
    out.param_names1 = 'Linf, a, diffus,e, growingDomain,initialDomainLength,makeChemoattractant';
    out.param_names2 = 'chi,domainHeight,zero_bc,insert,tstep,t_start,d,eatRate';
    out.param = param; %param = [Linf, a, diffus, eatWidth, growingDomain, initialDomainLength, makeChemoattractant, chi, domainHeight, zeroBC, insert, tstep, t_start, eatRate, num_steps, num_directions];
    out.moved = moved;
    out.happiness = happiness;
    
    % these are parameters we might sweep
    out.leadSpeed = leadSpeed;
    out.followSpeed = followSpeed;
    out.numFilopodia = numFilopodia;
    
    out.growingDomain = growingDomain;
    out.followerFraction = followerFraction;
    out.divide_cells = divide_cells;
    out.experiment = experiment;
    
    save(['/mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'.mat'],'out')
    delete(['/mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'_running.mat'])
    fprintf(['saved results to /mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'.mat \n'])

    
end
