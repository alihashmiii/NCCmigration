% plot migration profiles for simulated VEGF transplants.
% L.J. Schumacher 05.09.14

close all
clear all

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
defaultFollowValues = [0 1 2];
lead2follow = [4];
follow2lead = [4];
sensingAccuracyValues = [0.1, 0.01];
experiments = [0 12 11];
numParamCombinations = length(defaultFollowValues)*length(sensingAccuracyValues)...
    *length(experiments);

xBins = 0:50:800; % bins for counting cell num vs. x profiles

for defaultFollow = defaultFollowValues
    for sensAccCtr = 1:length(sensingAccuracyValues)
        sensingAccuracy = sensingAccuracyValues(sensAccCtr);
        figure
        hold all
        for expCtr = 1:length(experiments)
            experiment = experiments(expCtr);
            % preallocate variables for saving collated results
            cellDistributions = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats,3,length(xBins));
            actualLeaderFraction = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats);
            numCells = NaN(length(defaultFollowValues),length(sensingAccuracyValues),...
                numRepeats);
            
            
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
                    end
                    %% plot migration profile
                    f_L = mean(actualLeaderFraction(defaultFollow + 1,sensAccCtr,:));
                    n_C = mean(numCells(defaultFollow + 1,sensAccCtr,:));
                    % plot migration profile
                    plot(xBins,squeeze(mean(sum(cellDistributions(defaultFollow + 1,sensAccCtr,:,:,:),4),3)),...
                        'LineWidth',2);
        end
        grid on
        set(gca,'GridLineStyle','-')
        legend('control', 'within','adjacent')
        xlabel('distance along stream (\mum)')
        ylabel('number of cell (per 50\mum)')
        %% export figure
        exportOptions = struct('Format','eps2',...
            'Width','18.0',...
            'Color','rgb',...
            'Resolution',300,...
            'FontMode','fixed',...
            'FontSize',10,...
            'LineWidth',2);
        
        filename = ['manuscripts/VEGF/figures/Fig4K_defaultFollow_' num2str(defaultFollow) '_sensAcc_' num2str(sensingAccuracy)];
        pos = get(gcf,'Position');
        % pos(4) = 1/2*pos(3); % adjust height to fraction of width
        set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
        exportfig(gcf,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
        
    end
end