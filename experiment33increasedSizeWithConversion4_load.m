% load simulation results and save data for plotting
% L.J. Schumacher 05.09.14, 14.08.15

close all
clear all

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
defaultFollowValues = [1 2];
switchingTimes = [4 8];
sensingAccuracyValues = [0.1, 0.01];
experiments = [33 34];
insertStepsValues = [1 3; 12 24];
xBins = 0:50:800; % bins for counting cell num vs. x profiles

for defaultFollow = defaultFollowValues
    for sensAccCtr = 1:length(sensingAccuracyValues)
        sensingAccuracy = sensingAccuracyValues(sensAccCtr);
        for switchingTime = switchingTimes
            numSteps = [switchingTime, switchingTime];
            figure
            hold on
            for expCtr = 1:length(experiments)
                experiment = experiments(expCtr);
                for insertStepsCtr = 1:2
                    insertEverySteps = insertStepsValues(expCtr,insertStepsCtr);
                    %% load data
                    actualLeaderFraction = NaN(numRepeats,1);
                    numCells = NaN(numRepeats,1);
                    for repCtr = 1:numRepeats
                        loadInfo = ['experiment' num2str(experiment) '/exp' num2str(experiment) ...
                            '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                            '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                            '_sensingAcc_' num2str(sensingAccuracy) '_insertSteps_' num2str(insertEverySteps) ...
                            '_Run_' num2str(repCtr)];
                        try % sometime we get corrupt files, which crashes the script
                            load(['results/' loadInfo '.mat'])
                        catch
                            delete(['results/' loadInfo '.mat']) % delete the corrupt file
                            switch experiment
                                case 33
                                    experiment33increasedSizeWithConversion4; % recreate the missing results file
                                case 34
                                    experiment34decreasedSizeWithConversion4; % recreate the missing results file
                            end
                            load(['results/' loadInfo '.mat']) % load again
                        end
                        
                        % load cell positions into variables
                        cells = out.cells_save{end}; % all cells
                        numberOfCells = size(cells,2);
                        followIdcs = out.cellsFollow_save{end}(1:numberOfCells);
                        attachIdcs = out.attach_save{end}(1:numberOfCells);
                        leaders = cells(:,followIdcs==0);
                        followers = cells(:,followIdcs==1&attachIdcs~=0);
                        losts = cells(:,followIdcs==1&attachIdcs==0);
                        
                        actualLeaderFraction(repCtr) = size(leaders,2)/numberOfCells;
                        numCells(repCtr) = numberOfCells;
                    end
                    %% plot data
                    h = plot(numCells,actualLeaderFraction,'o');
                    plot(numCells,actualLeaderFraction.*numCells/100,'+','Color',h.Color);
                end
            end
            grid on
            box on
            set(gca,'GridLineStyle','-')
            xlabel('number of cells')
            ylabel('fraction (o) and number/100 (+) of leaders')
            %% export figure
            exportOptions = struct('Format','eps2',...
                'Width','17.0',...
                'Color','rgb',...
                'Resolution',300,...
                'FontMode','fixed',...
                'FontSize',10,...
                'LineWidth',2);
            
            filename = ['results/experiment33/figures/leaderFractionPopulationSize' ...
                '_defaultFollow_' num2str(defaultFollow) '_switchT_' num2str(switchingTime) ...
                '_sensingAcc_' num2str(sensingAccuracy)];
            pos = get(gcf,'Position');
            set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
            exportfig(gcf,[filename '.eps'],exportOptions);
            system(['epstopdf ' filename '.eps']);
        end
    end
end