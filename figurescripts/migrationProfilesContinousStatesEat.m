% plot migration outcome for continuous vs discrete states

close all
clear all
addpath('../')

time = 18;
numRepeats = 40;

% simulation parameters
sensingAccuracy = 0.1;
eatRates = [10, 25, 50, 75, 100, 125, 250, 500, 1000, 1500, 2000];
nVals = length(eatRates);
guidanceModes = {'choice','combination'};

% auxiliary variables for plotting and loading
cellRadius = 7.5;
time2plot = [24];
precision = 2; % significant figures for filenames and plot labels etc.
loadpath = '../results/';
plotSnapshot = true;

lineStyles = {'-','--'};
plotColors = lines(2);
migrationProfilesFig = figure;
hold on
for gdmCtr = 1:length(guidanceModes)
    guidanceMode = guidanceModes{gdmCtr};
    % preallocate variables for saving collated results
    numCells = NaN(length(eatRates),numRepeats);
    xMax = NaN(length(eatRates),numRepeats);
    hold on
    for eatCtr = 1:nVals
        eatRate = eatRates(eatCtr);
        
        %% load data
        for repCtr = 1:numRepeats
            filename = ['experiment31contStates_eat/exp31' ...
                '_contStates_' guidanceMode '_eat_' num2str(eatRate,precision) ...
                '_sensingAcc_' num2str(sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            load([loadpath filename '.mat'])
            
            % load cell positions into variables
            timeIdx = find(out.t_save >= time2plot,1,'first');
            cells = out.cells_save{timeIdx}; % all cells
            
            numCells(eatCtr,repCtr) = size(cells,2);
            xMax(eatCtr,repCtr) = max(cells(1,:));
        end
        % determine which snapshot to plot, if any
        if plotSnapshot&&gdmCtr==1&&eatRate>=40&&eatRate<=150
            % calculate rep that is closest to mean in relative terms
            numCellsRel = abs(numCells(eatCtr,:) - mean(numCells(eatCtr,:),2))./mean(numCells(eatCtr,:),2);
            xMaxRel = abs(xMax(eatCtr,:) - mean(xMax(eatCtr,:),2))./mean(xMax(eatCtr,:),2);
            [minDev, repToPlot] = min(numCellsRel + xMaxRel);
            % reload this rep
            if repToPlot~=repCtr
                filename = ['experiment31contStates_eat/exp31' ...
                    '_contStates_' guidanceMode '_eat_' num2str(eatRate,precision) ...
                    '_sensingAcc_' num2str(sensingAccuracy,precision) '_Run_' num2str(repToPlot)];
                load([loadpath filename '.mat'])
                timeIdx = find(out.t_save >= time2plot,1,'first');
            end
            % plot snapshots of simulations
            plot_snapshot(out,timeIdx)
        end
    end
    
    yyaxis left
    errorbar(eatRates,mean(numCells,2),std(numCells,0,2)/sqrt(numRepeats),...
        lineStyles{gdmCtr},'Color',plotColors(1,:),'LineWidth',2);
    yyaxis right
    errorbar(eatRates,mean(xMax,2),std(xMax,0,2)/sqrt(numRepeats),...
        lineStyles{gdmCtr},'Color',plotColors(2,:),'LineWidth',2);
end
box on
xlabel('chemoattractant consumption rate \lambda')
set(gca,'xscale','log')
yyaxis left
ylim([60 150])
ylabel('number of cells')
yyaxis right
ylim([550 1000]);
ylabel('max. dist. migrated (\mum)')
% add a box to highlight features
rectangle('Position',[40,650,125,200],'Curvature',0.25)
%% export figure
exportOptions = struct('Format','eps2',...
    'Width','9.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = ['../manuscripts/JTB/figures/Fig3_contStates_eat_'...
    'sensAcc_' num2str(100*sensingAccuracy)];
set(migrationProfilesFig,'PaperUnits','centimeters');
exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
