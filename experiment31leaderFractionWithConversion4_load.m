% load
% L.J. Schumacher 28.10.13

close all
clear all

time = 18;
numRepeats = 20;
maxRuns2plot = 20;
numParamCombinations = 36;
% to calculate the density profile of cells and chemoattractant along the x-direction
cellRadius = 7.5;
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
xBins = 0:(cellRadius + filolength):24*(cellRadius + filolength); % bins for counting cell num vs. x profiles
cellDistributions = NaN(numParamCombinations,numRepeats,3,length(xBins));
caDistribution = NaN(numParamCombinations,numRepeats,64);
xlat_save = NaN(64,1);
% preallocate variables for saving collated results
actualLeaderFraction = NaN(numParamCombinations,1);
sensingAccuracies = NaN(numParamCombinations,1);
neighboursNeeds = NaN(numParamCombinations,1);

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;
conversionType = 4;
lead2followValues = [4 8];
follow2leadValues = [8 16 30];
follow2leadColors = jet(length(follow2leadValues));
caCmap = load('cmap_blue2cyan.txt');

for sensingAccuracy = [0.1, 0.01]
    for needNeighbours = [0, 1 ,2]
        profilesFig = figure('Visible','off');
        profiles2getherFig = figure('Visible','off');
        neighbourRelationsFig = figure('Visible','off');
            %% load and plot data for every run of this parameter combination
            for lead2followCtr = 1:length(lead2followValues)
                lead2follow = lead2followValues(lead2followCtr);
            for follow2leadCtr = 1:length(follow2leadValues)
                follow2lead = follow2leadValues(follow2leadCtr);
                numSteps = [lead2follow, follow2lead];
                runsFig = figure('Visible','off');
                for repCtr = 1:numRepeats
                    loadInfo = ['experiment31conversion4/exp31'...
                        '_conversion_' num2str(conversionType) '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours)...
                        '_Run_' num2str(repCtr)];
                    try % sometime we get corrupt files, which crashes the script
                        load(['results/' loadInfo '.mat'])
                    catch
                        delete(['results/' loadInfo '.mat']) % delete the corrupt file
                        experiment31leaderFractionWithConversion4; % recreate the missing results file
                        load(['results/' loadInfo '.mat']) % load again
                    end
                    % make a plot of all repeats
                    if repCtr <= maxRuns2plot
                        subplot(min(numRepeats,maxRuns2plot)/2 + 2,2,repCtr+2)
                        make_plot(out.cells_save{end},out.cellsFollow_save{end},out.xlat_save{end},out.ylat_save{end}, ...
                            out.ca_save{end},out.filopodia_save{end},out.numFilopodia,out.attach_save{end},out.cellRadius,filolength,sensingAccuracy,0,caCmap,1)
                        title([num2str(size(out.cells_save{end},2)) ' cells, ' num2str(min([size(out.cells_save{end},2) nnz(out.cellsFollow_save{end}==0)])) ' leaders.'])
                    end
                    % calculate migration profile
                    numberOfCells = size(out.cells_save{end},2);
                    cellDistributions(paramCtr,repCtr,1,:) = histc(out.cells_save{end}(1,out.cellsFollow_save{end}(1:numberOfCells)==0),xBins); % leaders
                    cellDistributions(paramCtr,repCtr,2,:) = histc(out.cells_save{end}(1,(out.cellsFollow_save{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)~=0)),xBins); % followers, attached
                    cellDistributions(paramCtr,repCtr,3,:) = histc(out.cells_save{end}(1,(out.cellsFollow_save{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)==0)),xBins); % followers, dettached
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
                title(['Exp3.1: leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', sensAcc=' num2str(sensingAccuracy) ', needNbrs=' num2str(needNeighbours) ', lead2follow=' num2str(lead2follow) ', follow2lead=' num2str(follow2lead) ])
                
                % save plot
                filename = ['results/experiment31conversion4/figures/exp31conv4' ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
                    '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
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
                subplot(length(follow2leadValues),length(lead2followValues),length(lead2followValues)*(follow2leadCtr - 1) + lead2followCtr)
                plot_migration_profile
                % xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
                % %                     legend([H3;H1;H2],'leaders','followers','chemoattractant');
                
                % title has parameter values and actual leader fraction
                if follow2leadCtr==1&&lead2followCtr==1
                    title(['Exp3.1: lead2follow=' num2str(lead2follow) ', sensAcc=' num2str(sensingAccuracy) ' needNbrs=' num2str(needNeighbours) ', leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ])
                elseif follow2leadCtr==1&&lead2followCtr==2
                    title(['lead2follow=' num2str(lead2follow) ', leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ])
                else
                    title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', follow2lead=' num2str(follow2lead)])
                end
                
                set(0,'CurrentFigure',profiles2getherFig);
                subplot(length(lead2followValues),1,lead2followCtr)
                if follow2leadCtr==1, hold on, end
                plot(xBins,squeeze(mean(sum(cellDistributions(paramCtr,:,:,:),3),2)),'Color',follow2leadColors(follow2leadCtr,:));
                if follow2leadCtr==length(follow2leadValues)
                    title(['Exp3.1: sensAcc=' num2str(sensingAccuracy) ', needNbrs=' num2str(needNeighbours) ', lead2follow=' num2str(lead2follow)])
                    xlabel('x/\mum'), ylabel('N(cells)'), legend(num2str(follow2leadValues'))
                    ylim([0 10]), xlim([0 800]), set(gca,'YTick',[0 2 4 6 8 10]), grid on
                end
                
                %% plot neighbour relations
                set(0,'CurrentFigure',neighbourRelationsFig);
                subplot(length(lead2followValues),3,1 + (lead2followCtr - 1)*3)
                if follow2leadCtr ==1, hold on, end
                plot(1:length(neighbours(paramCtr).numbers),neighbours(paramCtr).numbers./numRepeats,'Color',follow2leadColors(follow2leadCtr,:))
                if follow2leadCtr==length(follow2leadValues)
                    xlabel('#neighbours'), ylabel('N(cells)'), grid on
                end
                
                subplot(length(lead2followValues),3,2 + (lead2followCtr - 1)*3)
                if follow2leadCtr ==1, hold on, end
                plot(neighbours(paramCtr).distancesBinEdges, neighbours(paramCtr).distances./numRepeats,'Color',follow2leadColors(follow2leadCtr,:))
                if follow2leadCtr==length(follow2leadValues)
                    title(['Exp3.1: sensAcc=' num2str(sensingAccuracy) ', needNbrs=' num2str(needNeighbours) ', lead2follow=' num2str(lead2follow)])
                    xlabel('distance/\mum'), ylabel('N(cells)')
                    xlim([0 max(neighbours(paramCtr).distancesBinEdges)]), grid on
                end
                
                subplot(length(lead2followValues),3,3 + (lead2followCtr - 1)*3)
                if follow2leadCtr ==1, hold on, end
                plot(neighbours(paramCtr).areasBinEdges, neighbours(paramCtr).areas./numRepeats,'Color',follow2leadColors(follow2leadCtr,:))
                if follow2leadCtr==length(follow2leadValues)
                    xlabel('area/\mum^2'), ylabel('N(cells)'), legend(num2str(follow2leadValues'))
                    xlim([0 max(neighbours(paramCtr).areasBinEdges)/2]), grid on
                end
                
                %% save summary of results
                sensingAccuracies(paramCtr) = sensingAccuracy;
                neighboursNeeds(paramCtr) = needNeighbours;
                paramCtr = paramCtr + 1;
            end
            end
        % make a plot with migration profiles for all parameter combinations
        pos = get(profilesFig,'Position');
        pos(4) = 3/2*pos(3);% adjust height to 3/2 width
        set(profilesFig,'PaperUnits','centimeters','Position',pos);
        filename = ['results/experiment31conversion4/figures/exp31conv4_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) '_migrationProfiles'];
        exportfig(profilesFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
        close(profilesFig);
        
        pos = get(profiles2getherFig,'Position');
%         pos(4) = 3/2*pos(3);% adjust height to 3/2 width
        set(profiles2getherFig,'PaperUnits','centimeters','Position',pos);
        filename = ['results/experiment31conversion4/figures/exp31conv4_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) '_migrationProfiles2gether'];
        exportfig(profiles2getherFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
        
        pos = get(neighbourRelationsFig,'Position');
%         pos(4) = 3/2*pos(3);% adjust height to 3/2 width
        set(neighbourRelationsFig,'PaperUnits','centimeters','Position',pos);
        filename = ['results/experiment31conversion4/figures/exp31conv4_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) '_neighbourRelations'];
        exportfig(neighbourRelationsFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        system(['cp ' filename '.pdf results/PDFs/' filename '.pdf']); % copying finished plots to a place where Dropbox will sync them
        close(neighbourRelationsFig);
    end
end
save('results/experiment31conversion4/figures/experiment31conv4collatedResults','xBins','cellDistributions','xlat_save','caDistribution','actualLeaderFraction','sensingAccuracies','neighboursNeeds')