% plot migration profiles for continuous vs discrete states

close all
clear all
addpath('../')

time = 18;
numRepeats = 40;

% simulation parameters
sensingAccuracy = 0.1;
psaValues = [0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 1];
nVals = length(psaValues);
guidanceModes = {'choice','combination'};

% auxiliary variables for plotting and loading
% xBins = 0:50:1100; % bins for counting cell num vs. x profiles
% plotBins = xBins(2:end) - mean(diff(xBins))/2;
cellRadius = 7.5;
time2plot = [24];
precision = 2; % significant figures for filenames and plot labels etc.
loadpath = '../results/';
% load('~/Dropbox/Utilities/colormaps_ascii/increasing_warm/cmap_RdOrYl.txt')
% set(0,'defaultAxesColorOrder',cmap_RdOrYl(round(linspace(1,200,nVals)),:));
lineStyles = {'-','--'};
plotColors = lines(2);
migrationProfilesFig = figure;
hold on
for gdmCtr = 1:length(guidanceModes)
    guidanceMode = guidanceModes{gdmCtr};
    % preallocate variables for saving collated results
    numCells = NaN(length(psaValues),numRepeats);
    xMax = NaN(length(psaValues),numRepeats);
    hold on
    for psaCtr = 1:nVals
        psa = psaValues(psaCtr);
        
        %% load data
        for repCtr = 1:numRepeats
            filename = ['experiment31contStates_Psa/exp31' ...
                '_contStates_' guidanceMode '_Psa_' num2str(psa,precision) ...
                '_sensingAcc_' num2str(sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            load([loadpath filename '.mat'])
            
            % load cell positions into variables
            timeIdx = find(out.t_save >= time2plot,1,'first');
            cells = out.cells_save{timeIdx}; % all cells
            
            numCells(psaCtr,repCtr) = size(cells,2);
            xMax(psaCtr,repCtr) = max(cells(1,:));
            
        end
    end
%     yyaxis left
%     plot(psaValues,mean(numCells,2),...
%         lineStyles{gdmCtr},'Color',plotColors(1,:),'LineWidth',2);
%     yyaxis right
%     plot(psaValues,mean(xMax,2),...
%         lineStyles{gdmCtr},'Color',plotColors(2,:),'LineWidth',2);
%         add errorbars? looks too crowded...
        yyaxis left
    errorbar(1-psaValues,mean(numCells,2),std(numCells,0,2)/sqrt(numRepeats),...
        lineStyles{gdmCtr},'Color',plotColors(1,:),'LineWidth',2);
    yyaxis right
    errorbar(1-psaValues,mean(xMax,2),std(xMax,0,2)/sqrt(numRepeats),...
        lineStyles{gdmCtr},'Color',plotColors(2,:),'LineWidth',2);
end
box on
xlabel('detachment probability P_d')
% xticks(0.8:0.05:1)
yyaxis left
ylim([50 130])
ylabel('number of cells')
yyaxis right
ylim([600 1000]);
ylabel('max. dist. migrated (\mum)')

%% export figure
exportOptions = struct('Format','eps2',...
    'Width','9.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = ['../manuscripts/JTB/figures/Fig3_contStates_Psa_'...
    'sensAcc_' num2str(100*sensingAccuracy)];
set(migrationProfilesFig,'PaperUnits','centimeters');
exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
