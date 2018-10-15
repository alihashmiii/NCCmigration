% plot migration profiles for simulated VEGF transplants.
% L.J. Schumacher 05.09.14

close all
clear all

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
defaultFollowValues = [0 1 2];
lead2follow = [8];
follow2lead = [8];
sensingAccuracyValues = [0.1, 0.01];
experiments = [0 12 11];
numParamCombinations = length(defaultFollowValues)*length(sensingAccuracyValues)...
    *length(experiments);

xBins = 0:50:800; % bins for counting cell num vs. x profiles
neighbourCutoff = 160;
cellRadius = 7.5;

for defaultFollow = defaultFollowValues
    for sensAccCtr = 1:length(sensingAccuracyValues)
        sensingAccuracy = sensingAccuracyValues(sensAccCtr);
        migrationProfilesFig = figure;
        hold on
        neighbourRelationsFig = figure;
        for ctr = 1:2
        subplot(1,2,ctr), hold on
        end
        for expCtr = 1:length(experiments)
            experiment = experiments(expCtr);
            % preallocate variables for saving collated results
            cellDistributions = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats,3,length(xBins));
            actualLeaderFraction = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats);
            numCells = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats);
%             neighbourNumbers = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
%                 numRepeats,10);
            neighbourAreas = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats,21);
            neighbourDistances = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats,length(2*cellRadius:cellRadius:neighbourCutoff));
            
            %% load data
                    numSteps = [lead2follow, follow2lead];
                    for repCtr = 1:numRepeats
                        if experiment==0 % load control simulation
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
                            loadInfo = ['experiment31transplants/exp' num2str(experiment) ...
                                '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                                '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                                '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                            try % sometime we get corrupt files, which crashes the script
                                load(['results/' loadInfo '.mat'])
                            catch
                                delete(['results/' loadInfo '.mat']) % delete the corrupt file
                                experiment31transplantsWithConversion4; % recreate the missing results file
                                load(['results/' loadInfo '.mat']) % load again
                            end
                        end
                        
                        % load cell positions into variables
                        cells = out.cells_save{end}; % all cells
                        numberOfCells = size(cells,2);
                        followIdcs = out.cellsFollow_save{end}(1:numberOfCells);
                        attachIdcs = out.attach_save{end}(1:numberOfCells);
                        leaders = cells(:,followIdcs==0);
                        followers = cells(:,followIdcs==1&attachIdcs~=0);
                        losts = cells(:,followIdcs==1&attachIdcs==0);
                        
                        actualLeaderFraction(defaultFollow + 1,sensAccCtr,repCtr) = size(leaders,2)/numberOfCells;
                        numCells(defaultFollow + 1,sensAccCtr,repCtr) = numberOfCells;
                        
                        % calculate migration profile
                        cellDistributions(defaultFollow + 1,sensAccCtr,repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
                        cellDistributions(defaultFollow + 1,sensAccCtr,repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
                        cellDistributions(defaultFollow + 1,sensAccCtr,repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached

                        % calculate neighbour relationships
                        neighbours = neighbourRelationships(out.cells_save{end},neighbourCutoff);
%                         neighbourNumbers(defaultFollow + 1,sensAccCtr,repCtr,:) = neighbours.numbers./sum(neighbours.numbers);
                        neighbourDistances(defaultFollow + 1,sensAccCtr,repCtr,:) = neighbours.distances./sum(neighbours.distances);
                        neighbourAreas(defaultFollow + 1,sensAccCtr,repCtr,:) = neighbours.areas./sum(neighbours.areas);
                    end
                    %% plot migration profile
                    f_L = mean(actualLeaderFraction(defaultFollow + 1,sensAccCtr,:));
                    n_C = mean(numCells(defaultFollow + 1,sensAccCtr,:));
                    % plot migration profile
                    set(0,'CurrentFigure',migrationProfilesFig);
                    % plot leaders
                    plotHandle = plot(xBins,squeeze(mean(cellDistributions(defaultFollow + 1,sensAccCtr,:,1,:),3)),...
                        'LineWidth',2,'LineStyle','-');
                    % plot followers
                    plot(xBins,squeeze(mean(sum(cellDistributions(defaultFollow + 1,sensAccCtr,:,2:3,:),4),3)),...
                        'LineWidth',2,'LineStyle','--','Color',get(plotHandle,'Color'));
                    %% plot neighbour relations
                    set(0,'CurrentFigure',neighbourRelationsFig);
%                     subplot(1,3,1)
%                     plot(1:length(neighbours.numbers), squeeze(mean(neighbourNumbers(defaultFollow + 1,sensAccCtr,:,:),3)))
                    
                    subplot(1,2,1)
                    plot(neighbours.distancesBinEdges, squeeze(mean(neighbourDistances(defaultFollow + 1,sensAccCtr,:,:),3)))
                    
                    subplot(1,2,2)
                    plot(neighbours.areasBinEdges, squeeze(mean(neighbourAreas(defaultFollow + 1,sensAccCtr,:,:),3)))
        end
        set(0,'CurrentFigure',migrationProfilesFig);
        grid off
        set(gca,'GridLineStyle','-')
        legend('control (leaders)','control (followers)','within (leaders)','within (followers)','adjacent (leaders)','adjacent (followers)')
        xlabel('distance along stream (\mum)')
        ylabel('number of cell (per 50\mum)')
              
        set(0,'CurrentFigure',neighbourRelationsFig);
%         subplot(1,3,1)
%         xlabel('#neighbours'), ylabel('P')
%         legend('control', 'within','adjacent')
        subplot(1,2,1)
        xlabel('distance/\mum'), ylabel('P')
        xlim([2*cellRadius,neighbourCutoff])
        legend('control', 'within','adjacent')
        ylim([0 0.5])
        subplot(1,2,2)
        xlabel('area/\mum^2'), ylabel('P')
        xlim([3*sqrt(3)/2*cellRadius.^2,1e4])
        xlim([0 6e3])
        ylim([0 0.7])
        legend('control', 'within','adjacent')
        %% export figure
        exportOptions = struct('Format','eps2',...
            'Width','9.0',...
            'Color','rgb',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',10,...
            'LineWidth',2);
        
        filename = ['manuscripts/VEGF/figures/Fig5I_defaultFollow_' num2str(defaultFollow) '_sensAcc_' num2str(sensingAccuracy)];
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
        
        filename = ['manuscripts/VEGF/figures/FigS5B_defaultFollow_' num2str(defaultFollow) '_sensAcc_' num2str(sensingAccuracy)];
        set(neighbourRelationsFig,'PaperUnits','centimeters');
        exportfig(neighbourRelationsFig,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        
    end
end