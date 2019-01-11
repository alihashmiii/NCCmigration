% plot migration outcomes for continuous vs discrete states

close all
clear all
addpath('../')
addpath('../simulationscripts')

time = 18;
numRepeats = 20;

% simulation parameters
sensingAccuracyValues = [0.1, 0.01 ];
experiments = {'control','choice'};

% auxiliary variables for plotting and loading
xBins = 0:50:1000; % bins for counting cell num vs. x profiles
plotBins = xBins(2:end) - mean(diff(xBins))/2;
cellRadius = 7.5;
time2plot = [24];
precision = 2; % significant figures for filenames and plot labels etc.
loadpath = '../results/';
for sensAccCtr = 1:length(sensingAccuracyValues)
    sensingAccuracy = sensingAccuracyValues(sensAccCtr);
    migrationProfilesFig = figure;
    hold on
    
    for expCtr = 1:length(experiments)
        experiment = experiments(expCtr);
        % preallocate variables for saving collated results
        cellDistributions = NaN(length(sensingAccuracyValues),...
            numRepeats,length(xBins)-1);
        numCells = NaN(length(sensingAccuracyValues),...
            numRepeats);
        
        %% load data
        for repCtr = 1:numRepeats
            if strcmp(experiment,'control')% load control simulations
                filename = ['experiment31conversion4/exp31'...
                    '_conversion_4_defaultFollow_2_numSteps_8_8' ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                load([loadpath filename '.mat'])
            else % load  simulation
                filename = ['experiment31contStates_choice/exp31'...
                    '_contStates_choice_sensingAcc_' num2str(sensingAccuracy)...
                    '_Run_' num2str(repCtr)];
                load([loadpath filename '.mat'])
            end
            
            % load cell positions into variables
            timeIdx = find(out.t_save >= time2plot,1,'first');
            cells = out.cells_save{timeIdx}; % all cells
            
            numCells(sensAccCtr,repCtr) = size(cells,2);
            
            % calculate migration profile
            cellDistributions(sensAccCtr,repCtr,:) = histcounts(cells(1,:),xBins); % leaders
            
        end
        %% plot migration profile
        % plot migration profile
        set(0,'CurrentFigure',migrationProfilesFig);
        plot(plotBins,squeeze(mean(cellDistributions(sensAccCtr,:,:),2)),...
            'LineWidth',2);
    end
    set(0,'CurrentFigure',migrationProfilesFig);
    grid off
    set(gca,'GridLineStyle','-')
    legend('discrete states','continuous states')
    xlabel('distance (\mum)')
    ylabel('number of cells (per 50\mum)')
    
    %% export figure
    exportOptions = struct('Format','eps2',...
        'Width','9.0',...
        'Color','rgb',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    
    filename = ['../manuscripts/JTB/figures/Fig2_contStates_choice_sensAcc_' num2str(100*sensingAccuracy)];
    set(migrationProfilesFig,'PaperUnits','centimeters');
    exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    
end
