% plot migration profiles for continuous vs discrete states

close all
clear all
addpath('../')
addpath('../simulationscripts')

time = 18;
numRepeatsNew = 40;
numRepeatsControl = 20;

% simulation parameters
sensingAccuracyValues = [0.1, 0.01];
experiments = {'control','choice','combination'};

% auxiliary variables for plotting and loading
xBins = 0:50:1100; % bins for counting cell num vs. x profiles
plotBins = xBins(2:end) - mean(diff(xBins))/2;
cellRadius = 7.5;
time2plot = [24];
precision = 2; % significant figures for filenames and plot labels etc.
loadpath = '../results/';

exportOptions = struct('Format','eps2',...
    'Width','9.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

for sensAccCtr = 1:length(sensingAccuracyValues)
    sensingAccuracy = sensingAccuracyValues(sensAccCtr);
    migrationProfilesFig = figure;
    hold on
%     cellNumFig = figure;
%     hold on
    contactTimeFig = figure;
    hold on
    for expCtr = 1:length(experiments)
        experiment = experiments{expCtr};
        if strcmp(experiment,'control')
            numRepeats = numRepeatsControl;
        else
            numRepeats = numRepeatsNew;
        end
        
        % preallocate variables for saving collated results
        cellDistributions = NaN(numRepeats,length(xBins)-1);
        numCells = NaN(numRepeats,1);
        contactTimes = cell(numRepeats,1);
        %% load data
        for repCtr = 1:numRepeats
            if strcmp(experiment,'control')% load control simulations
                filename = ['experiment31conversion4/exp31'...
                    '_conversion_4_defaultFollow_2_numSteps_8_8' ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                load([loadpath filename '.mat'])
            else % load  simulation
                filename = ['experiment31contStates_' experiment '/exp31'...
                    '_contStates_' experiment '_sensingAcc_' num2str(sensingAccuracy)...
                    '_Run_' num2str(repCtr)];
                load([loadpath filename '.mat'])
            end
            
            % load cell positions into variables
            timeIdx = find(out.t_save >= time2plot,1,'first');
            cells = out.cells_save{timeIdx}; % all cells
            
            numCells(repCtr) = size(cells,2);
            
            % calculate migration profile
            cellDistributions(repCtr,:) = histcounts(cells(1,:),xBins); % leaders
            
            % compute contact times
            contactTimes{repCtr} = computeContactTimes(out.attach_save(1:timeIdx),...
                ceil(out.numTsteps/out.numSavepoints));
        end
        %% plot migration profile
        % plot migration profile
        set(0,'CurrentFigure',migrationProfilesFig);
        plot(plotBins,squeeze(mean(cellDistributions,1)),...
            'LineWidth',2);
%         % plot cell number distribution
%         set(0,'CurrentFigure',cellNumFig);
%         histogram(numCells,'BinWidth',10,'Normalization','Probability')
        % plot contact time distribution
        set(0,'CurrentFigure',contactTimeFig);
        histogram(vertcat(contactTimes{:})/60,'Normalization','Probability',...
            'DisplayStyle','stairs','LineWidth',2)
    end
    %% format and export figures
    % migration profiles
    set(0,'CurrentFigure',migrationProfilesFig);
    box on
    legend({'discrete states','signal choice','signal combination'},'Location','North')
    xlabel('distance (\mum)')
    ylabel('number of cells (per 50\mum)')
    xlim([0 1000])
    
    filename = ['../manuscripts/JTB/figures/Fig2_contStates_combination_sensAcc_' num2str(100*sensingAccuracy)];
    set(migrationProfilesFig,'PaperUnits','centimeters');
    exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    
    % contact times
    set(0,'CurrentFigure',contactTimeFig);
    box on
%     legend({'discrete states','signal choice','signal combination'},'Location','NorthEast')
    xlabel('cell-cell contact time (hrs)')
    ylabel('relative frequency')
    xlim([0 time])
    
    filename = ['../manuscripts/JTB/figures/Fig2B_contStates_combination_sensAcc_' num2str(100*sensingAccuracy)];
    set(contactTimeFig,'PaperUnits','centimeters');
    exportfig(contactTimeFig,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
    
end
