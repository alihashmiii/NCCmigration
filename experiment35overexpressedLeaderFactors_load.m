% load
% L.J. Schumacher 06.04.14

close all
clear

time = 18;
numRepeats = 20;
maxRuns2plot = 20;
numParamCombinations = 27;
% to calculate the density profile of cells and chemoattractant along the x-direction
cellRadius = 7.5;              % radius in um (= 7.5um)
filolength = cellRadius + 9*2;
xBins = 0:50:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(numParamCombinations,numRepeats,3,length(xBins));
caDistribution = NaN(numParamCombinations,numRepeats,64);
xlat_save = NaN(64,1);
% preallocate variables for saving collated results
actualLeaderFraction = NaN(numParamCombinations,1);
sensingAccuracies = NaN(numParamCombinations,1);
neighboursNeeds = NaN(numParamCombinations,1);

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;
followFracValues = [7/8, 15/16, 1];
sensingAccuracyValues = [0.1, 0.01, 0.001];
needNeighboursValues = [0, 1, 2];
followFracColors = jet(length(followFracValues));

caCmap = load('cmap_blue2cyan.txt');

for sensingAccuracy = sensingAccuracyValues
    profilesFig = figure('Visible','off');
    profiles2getherFig = figure('Visible','off');
    neighbourRelationsFig = figure('Visible','off');
    for followFracCtr = 1:length(followFracValues)
        followerFraction = followFracValues(followFracCtr);
        for needNeighboursCtr = 1:length(needNeighboursValues)
            needNeighbours = needNeighboursValues(needNeighboursCtr);
            runsFig = figure('Visible','off');
            %% load and plot data for every run of this parameter combination
            for repCtr = 1:numRepeats
                loadInfo = ['experiment35/exp35_followFrac_' num2str(followerFraction,precision) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
                    '_Run_' num2str(repCtr)];
                try % sometimes we get corrupt files, which crashes the script
                    load(['results/' loadInfo '.mat'])
                catch
                    delete(['results/' loadInfo '.mat']) % delete the corrupt file
                    experiment35overexpressedLeaderFactors; % recreate the missing results file
                    load(['results/' loadInfo '.mat']) % load again
                end                        % make a plot of all repeats
                if repCtr <= maxRuns2plot
                    subplot(min(numRepeats,maxRuns2plot)/2 + 2,2,repCtr+2)
                    make_plot(out.cells_save{end},out.cellsFollow{end},out.xlat_save{end},out.ylat_save{end}, ...
                        out.ca_save{end},out.filopodia_save{end},out.numFilopodia,out.attach_save{end},out.cellRadius,filolength,sensingAccuracy,0,caCmap,1)
                    title([num2str(size(out.cells_save{end},2)) ' cells, ' num2str(min([size(out.cells_save{end},2) nnz(out.cellsFollow{end}==0)])) ' leaders.'])
                end
                % calculate migration profile
                numberOfCells = size(out.cells_save{end},2);
                cellDistributions(paramCtr,repCtr,1,:) = histc(out.cells_save{end}(1,out.cellsFollow{end}(1:numberOfCells)==0),xBins); % leaders
                cellDistributions(paramCtr,repCtr,2,:) = histc(out.cells_save{end}(1,(out.cellsFollow{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)~=0)),xBins); % followers, attached
                cellDistributions(paramCtr,repCtr,3,:) = histc(out.cells_save{end}(1,(out.cellsFollow{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)==0)),xBins); % followers, attached
                caDistribution(paramCtr,repCtr,:) = mean(out.ca_save{end},2);
                if paramCtr==1, xlat_save = out.xlat_save{end}; end % load the x-coordinated of the CA profile, only once as they're always the same
                % calculate neighbour relationships
                if paramCtr==1&&repCtr==1
                    neighbours = neighbourRelationships(out.cells_save{end});
                elseif repCtr==1
                    neighbours(paramCtr) = neighbourRelationships(out.cells_save{end});
                else
                    tempNeighbours = neighbourRelationships(out.cells_save{end});
                    neighbours(paramCtr).numbers = neighbours(paramCtr).numbers + tempNeighbours.numbers;
                    neighbours(paramCtr).distances = neighbours(paramCtr).distances + tempNeighbours.distances;
                    neighbours(paramCtr).areas = neighbours(paramCtr).areas + tempNeighbours.areas;
                end
            end
            % plot migration profile
            subplot(min(numRepeats,maxRuns2plot)/2 + 2,2,[1 2])
            plot_migration_profile
            xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
            %legend([H3(3);H3(2);H3(1)],'lead','follow','lost');
            
            % title has parameter values and actual leader fraction
            actualLeaderFraction(paramCtr) = sum(mean(squeeze(cellDistributions(paramCtr,:,1,:)))); % mean number of leader cells
            actualLeaderFraction(paramCtr) = actualLeaderFraction(paramCtr)/(actualLeaderFraction(paramCtr) + sum(sum(mean(squeeze(cellDistributions(paramCtr,:,2:3,:)))))); % divide by mean total number of cells
            title(['Exp3.5: leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', sensAcc=' num2str(sensingAccuracy) ', needNbrs=' num2str(needNeighbours)])
            
            % save plot
            filename = ['results/experiment35/figures/exp35_followFrac_' num2str(followerFraction,precision) ...
                '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
                '_allRuns'];
            pos = get(runsFig,'Position');
            pos(4) = 3/2*pos(3);% adjust height to 3/2 width
            set(runsFig,'PaperUnits','centimeters','Position',pos);
            exportfig(runsFig,[filename '.eps'],exportOptions);
            system(['epstopdf ' filename '.eps']);
            system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
            close(runsFig);
            
            %% plot summary migration profiles
            set(0,'CurrentFigure',profilesFig);
            subplot(length(followFracValues),length(needNeighboursValues),length(needNeighboursValues)*(followFracCtr - 1) + needNeighboursCtr)
            plot_migration_profile
            % xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
            % %                     legend([H3;H1;H2],'leaders','followers','chemoattractant');
            
            % title has parameter values and actual leader fraction
            if needNeighboursCtr==1&&followFracCtr==1
                title(['Exp3.5: sensAcc=' num2str(sensingAccuracy)])
            elseif followFracCtr==1
                title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision)])
            else
                title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ' needNbrs=' num2str(needNeighbours)])
            end
            
            set(0,'CurrentFigure',profiles2getherFig);
            subplot(length(needNeighboursValues),1,needNeighboursCtr)
            if followFracCtr==1, hold on, end
            plot(xBins,squeeze(mean(sum(cellDistributions(paramCtr,:,:,:),3),2)),'Color',followFracColors(followFracCtr,:));
            if followFracCtr==length(followFracValues)
                xlabel('x/\mum'), ylabel('N(cells)'), legend(num2str(actualLeaderFraction((length(needNeighboursValues):length(needNeighboursValues):length(needNeighboursValues)*length(followFracValues)) - length(needNeighboursValues)*length(followFracValues) + paramCtr),precision))
                ylim([0 10]), xlim([0 735]), set(gca,'YTick',[0 2 4 6 8 10]), grid on
            end
            
            %% plot neighbour relations
            set(0,'CurrentFigure',neighbourRelationsFig);
            subplot(length(needNeighboursValues),3,1 + (needNeighboursCtr - 1)*3)
            if followFracCtr ==1, hold on, end
            plot(1:length(neighbours(paramCtr).numbers),neighbours(paramCtr).numbers./numRepeats,'Color',followFracColors(followFracCtr,:))
            if followFracCtr==length(followFracValues)
                xlabel('#neighbours'), ylabel('N(cells)'), grid on
            end
            
            subplot(length(needNeighboursValues),3,2 + (needNeighboursCtr - 1)*3)
            if followFracCtr ==1, hold on, end
            plot(neighbours(paramCtr).distancesBinEdges, neighbours(paramCtr).distances./numRepeats,'Color',followFracColors(followFracCtr,:))
            if followFracCtr==length(followFracValues)
                title(['Exp3.5: sensAcc=' num2str(sensingAccuracy) ', needNbrs=' num2str(needNeighbours)])
                xlabel('distance/\mum'), ylabel('N(cells)')
                xlim([0 max(neighbours(paramCtr).distancesBinEdges)]), grid on
            end
            
            subplot(length(needNeighboursValues),3,3 + (needNeighboursCtr - 1)*3)
            if followFracCtr ==1, hold on, end
            plot(neighbours(paramCtr).areasBinEdges, neighbours(paramCtr).areas./numRepeats,'Color',followFracColors(followFracCtr,:))
            if followFracCtr==length(followFracValues)
                xlabel('area/\mum^2'), ylabel('N(cells)'), legend(num2str(actualLeaderFraction((length(needNeighboursValues):length(needNeighboursValues):length(needNeighboursValues)*length(followFracValues)) - length(needNeighboursValues)*length(followFracValues) + paramCtr),precision))
                xlim([0 max(neighbours(paramCtr).areasBinEdges)/2]), grid on
            end
            
            %% save summary of results
            sensingAccuracies(paramCtr) = sensingAccuracy;
            neighboursNeeds(paramCtr) = needNeighbours;
            paramCtr = paramCtr + 1;
        end
    end
    % for each sensing Accuracy make a plot with migration profiles for all set leader fractions and
    % needNeighbours - in total numParamCombinations plots with maxRuns2plot+1 subplots each
    pos = get(profilesFig,'Position');
    %pos(4) = 3/2*pos(3);% adjust height to 3/2 width
    set(profilesFig,'PaperUnits','centimeters','Position',pos);
    filename = ['results/experiment35/figures/exp35_sensingAcc_' num2str(sensingAccuracy) '_migrationProfiles'];
    exportfig(profilesFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
    close(profilesFig);
    
    pos = get(profiles2getherFig,'Position');
    pos(4) = 3/2*pos(3);% adjust height to 3/2 width
    set(profiles2getherFig,'PaperUnits','centimeters','Position',pos);
    filename = ['results/experiment35/figures/exp35_sensingAcc_' num2str(sensingAccuracy) '_migrationProfiles2gether'];
    exportfig(profiles2getherFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
    close(profiles2getherFig);
    
    pos = get(neighbourRelationsFig,'Position');
    pos(4) = 3/2*pos(3);% adjust height to 3/2 width
    set(neighbourRelationsFig,'PaperUnits','centimeters','Position',pos);
    filename = ['results/experiment35/figures/exp35_sensingAcc_' num2str(sensingAccuracy) '_neighbourRelations'];
    exportfig(neighbourRelationsFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
    close(neighbourRelationsFig);
end

save('results/experiment35/figures/experiment35collatedResults','xBins','cellDistributions','caDistribution','actualLeaderFraction','sensingAccuracies','neighboursNeeds')