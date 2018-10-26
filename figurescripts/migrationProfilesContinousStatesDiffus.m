% plot migration profiles for continuous vs discrete states

close all
clear all
addpath('../')
addpath('../simulationscripts')

time = 18;
numRepeats = 40;

% simulation parameters
sensingAccuracyUnscaled = 0.1;
diffusValues = [0.1, 1, 10, 1e2, 1e3, 1e4, 10^(4.5), 1e5];
nVals = length(diffusValues);
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
    numCells = NaN(length(diffusValues),numRepeats);
    xMax = NaN(length(diffusValues),numRepeats);
    hold on
    for diffusCtr = 1:nVals
        diffus = diffusValues(diffusCtr);
        sensingAccuracy = sensingAccuracyUnscaled*sqrt(0.1/diffus); % sensing accuracy scales with diffusivity
        
        %% load data
        for repCtr = 1:numRepeats
            filename = ['experiment31contStates_diffus/exp31' ...
                '_contStates_' guidanceMode '_D_' num2str(diffus,precision) ...
                '_sensingAcc_' num2str(sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            load([loadpath filename '.mat'])
            
            % load cell positions into variables
            timeIdx = find(out.t_save >= time2plot,1,'first');
            cells = out.cells_save{timeIdx}; % all cells
            
            numCells(diffusCtr,repCtr) = size(cells,2);
            xMax(diffusCtr,repCtr) = max(cells(1,:));
            
        end
    end
    yyaxis left
    plot(log10(diffusValues),mean(numCells,2),...
        lineStyles{gdmCtr},'Color',plotColors(1,:),'LineWidth',2);
    yyaxis right
    plot(log10(diffusValues),mean(xMax,2),...
        lineStyles{gdmCtr},'Color',plotColors(2,:),'LineWidth',2);
% %     add errorbars? looks too crowded...
%         yyaxis left
%     errorbar(log10(diffusValues),mean(numCells,2),std(numCells,0,2)/sqrt(numRepeats),...
%         lineStyles{gdmCtr},'Color',plotColors(1,:),'LineWidth',2);
%     yyaxis right
%     errorbar(log10(diffusValues),mean(xMax,2),std(xMax,0,2)/sqrt(numRepeats),...
%         lineStyles{gdmCtr},'Color',plotColors(2,:),'LineWidth',2);
    
end
box on
xlabel('log10(D)')
xlim([-1 5])
yyaxis left
ylim([60 130])
ylabel('number of cells')
yyaxis right
ylim([800 1000]);
ylabel('max. dist. migrated (\mum)')

%% export figure
exportOptions = struct('Format','eps2',...
    'Width','8.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = ['../manuscripts/JTB/figures/FigS2_contStates_diffus_'...
    'sensAcc_' num2str(100*sensingAccuracyUnscaled)];
set(migrationProfilesFig,'PaperUnits','centimeters');
exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
