% plot migration profiles for simulated NRP1 knockdown.
% L.J. Schumacher 05.09.14

close all
clear all

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
defaultFollowValues = [2];
sensingAccuracyValues = [0.1, 0.01];
numParamCombinations = length(defaultFollowValues)*length(sensingAccuracyValues);
timePoints2plot = [12, 24];

xBins = 0:50:800; % bins for counting cell num vs. x profiles
% parameters for neighbourhood analysis
neighbourCutoff = 84.34;
cellRadius = 7.5;

for defaultFollow = defaultFollowValues
    for sensAccCtr = 1:length(sensingAccuracyValues)
        sensingAccuracy = sensingAccuracyValues(sensAccCtr);
        for lead2follow = [4 8]
            for follow2lead = [4 8]
                numSteps = [lead2follow, follow2lead];
                migrationProfilesFig = figure;
                hold on
                neighbourRelationsFig = figure;
                for ctr = 1:2
                    subplot(1,2,ctr), hold on
                end
                for perturbation = 0:1
                    % preallocate variables for saving collated results
                    cellDistributions = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                        numRepeats,3,length(xBins),length(timePoints2plot));
                    actualLeaderFraction = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                        numRepeats,length(timePoints2plot));
                    numCells = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                        numRepeats,length(timePoints2plot));
                    neighbourAreas = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                        numRepeats,21,length(timePoints2plot));
                    neighbourDistances = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                        numRepeats,length(2*cellRadius:cellRadius:neighbourCutoff),length(timePoints2plot));
                    
                    %% load data
                    numSteps = [lead2follow, follow2lead];
                    for repCtr = 1:numRepeats
                        if perturbation==0 % load control simulation
                            loadInfo = ['experiment31conversion4/exp31'...
                                '_conversion_4_defaultFollow_' num2str(defaultFollow) ...
                                '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                                '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                            try % sometime we get corrupt files, which crashes the script
                                load(['results/' loadInfo '.mat'])
                            catch
                                delete(['results/' loadInfo '.mat']) % delete the corrupt file
                                experiment31leaderFractionWithConversion4; % recreate the missing results file
                                load(['results/' loadInfo '.mat']) % load again
                            end
                        else % load transplant simulations
                            loadInfo = ['experiment37reducedChemotaxis/exp37' ...
                                '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                                '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                                '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                            try % sometime we get corrupt files, which crashes the script
                                load(['results/' loadInfo '.mat'])
                            catch
                                delete(['results/' loadInfo '.mat']) % delete the corrupt file
                                experiment37reducedChemotaxisWithConversion4; % recreate the missing results file
                                load(['results/' loadInfo '.mat']) % load again
                            end
                        end
                        for timeCtr = 1:length(timePoints2plot)
                            % load cell positions into variables
                            cells = out.cells_save{end}; % all cells
                            numberOfCells = size(cells,2);
                            followIdcs = out.cellsFollow_save{end}(1:numberOfCells);
                            attachIdcs = out.attach_save{end}(1:numberOfCells);
                            leaders = cells(:,followIdcs==0);
                            followers = cells(:,followIdcs==1&attachIdcs~=0);
                            losts = cells(:,followIdcs==1&attachIdcs==0);
                            
                            actualLeaderFraction(defaultFollow + 1,sensAccCtr,repCtr,timeCtr) =...
                                size(leaders,2)/numberOfCells;
                            numCells(defaultFollow + 1,sensAccCtr,repCtr,timeCtr) = numberOfCells;
                            
                            % calculate migration profile
                            cellDistributions(defaultFollow + 1,sensAccCtr,repCtr,1,:,timeCtr) =...
                                histc(leaders(1,:),xBins); % leaders
                            cellDistributions(defaultFollow + 1,sensAccCtr,repCtr,2,:,timeCtr) =...
                                histc(followers(1,:),xBins); % followers, attached
                            cellDistributions(defaultFollow + 1,sensAccCtr,repCtr,3,:,timeCtr) =...
                                histc(losts(1,:),xBins); % followers, attached
                            
                            % calculate neighbour relationships
                            neighbours = neighbourRelationships(out.cells_save{end},neighbourCutoff);
                            neighbourDistances(defaultFollow + 1,sensAccCtr,repCtr,:,timeCtr) =...
                                neighbours.distances./sum(neighbours.distances);
                            neighbourAreas(defaultFollow + 1,sensAccCtr,repCtr,:,timeCtr) =...
                                neighbours.areas./sum(neighbours.areas);
                        end
                    end
                    %% plot migration profile
                    for timeCtr = 1:length(timePoints2plot)
                        f_L = mean(actualLeaderFraction(defaultFollow + 1,sensAccCtr,:,timeCtr));
                        n_C = mean(numCells(defaultFollow + 1,sensAccCtr,:),timeCtr);
                        % plot migration profile
                        set(0,'CurrentFigure',migrationProfilesFig);
                        % plot leaders
                        plotHandle = plot(xBins,squeeze(mean(cellDistributions(defaultFollow + 1,sensAccCtr,:,1,:,timeCtr),3)),...
                            'LineWidth',2,'LineStyle','-');
                        % plot followers
                        plot(xBins,squeeze(mean(sum(cellDistributions(defaultFollow + 1,sensAccCtr,:,2:3,:,timeCtr),4),3)),...
                            'LineWidth',2,'LineStyle','--','Color',get(plotHandle,'Color'));
                        %% plot neighbour relations
                        set(0,'CurrentFigure',neighbourRelationsFig);
                        
                        subplot(1,2,1)
                        plot(neighbours.distancesBinEdges, squeeze(mean(neighbourDistances(defaultFollow + 1,sensAccCtr,:,:,timeCtr),3)))
                        
                        subplot(1,2,2)
                        plot(neighbours.areasBinEdges, squeeze(mean(neighbourAreas(defaultFollow + 1,sensAccCtr,:,:,timeCtr),3)))
                    end
                end
                set(0,'CurrentFigure',migrationProfilesFig);
                grid off
                set(gca,'GridLineStyle','-')
                legend('control (leaders)','control (followers)',...
                    'reduced chemotaxis (leaders)','reduced chemotaxis (followers)')
                xlabel('distance along stream (\mum)')
                ylabel('number of cell (per 50\mum)')
                ylim([0 16])
                
                set(0,'CurrentFigure',neighbourRelationsFig);
                subplot(1,2,1)
                xlabel('distance/\mum'), ylabel('P')
                xlim([2*cellRadius,neighbourCutoff])
                ylim([0 0.5])
                legend('control', 'reduced chemotaxis')
                subplot(1,2,2)
                xlabel('area/\mum^2'), ylabel('P')
                xlim([3*sqrt(3)/2*cellRadius.^2,1e4])
                ylim([0 0.5])
                xlim([0 5e3])
                legend('control', 'reduced chemotaxis')
                %% export figure
                exportOptions = struct('Format','eps2',...
                    'Width','15.0',...
                    'Color','rgb',...
                    'Resolution',300,...
                    'FontMode','fixed',...
                    'FontSize',10,...
                    'LineWidth',2);
                
                filename = ['results/experiment37reducedChemotaxis/figures/migrationProfiles_defaultFollow_' num2str(defaultFollow) '_sensAcc_' num2str(sensingAccuracy)];
                set(migrationProfilesFig,'PaperUnits','centimeters');
                exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
                system(['epstopdf ' filename '.eps']);
                
                exportOptions = struct('Format','eps2',...
                    'Width','15.0',...
                    'Height','10.0',...
                    'Color','rgb',...
                    'Resolution',300,...
                    'FontMode','fixed',...
                    'FontSize',10,...
                    'LineWidth',2);
                
                filename = ['results/experiment37reducedChemotaxis/figures/neighbourRelations_' num2str(defaultFollow) '_sensAcc_' num2str(sensingAccuracy)];
                set(neighbourRelationsFig,'PaperUnits','centimeters');
                exportfig(neighbourRelationsFig,[filename '.eps'],exportOptions);
                system(['epstopdf ' filename '.eps']);
            end
        end
    end
end