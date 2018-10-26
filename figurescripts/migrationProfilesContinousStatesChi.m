% plot migration profiles for continuous vs discrete states

close all
clear all
addpath('../')

time = 18;
numRepeats = 40;

% simulation parameters
sensingAccuracy = 0.1;
chiValues = [1e-2, 1e-3, 1e-4, 1e-5, 10^(-5.5), 1e-6];
nVals = length(chiValues);
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
    numCells = NaN(length(chiValues),numRepeats);
    xMax = NaN(length(chiValues),numRepeats);
    hold on
    for chiCtr = 1:nVals
        chi = chiValues(chiCtr);
        
        %% load data
        for repCtr = 1:numRepeats
            filename = ['experiment31contStates_chi/exp31' ...
                '_contStates_' guidanceMode '_chi_' num2str(chi,precision) ...
                '_sensingAcc_' num2str(sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            load([loadpath filename '.mat'])
            
            % load cell positions into variables
            timeIdx = find(out.t_save >= time2plot,1,'first');
            cells = out.cells_save{timeIdx}; % all cells
            
            numCells(chiCtr,repCtr) = size(cells,2);
            xMax(chiCtr,repCtr) = max(cells(1,:));
            
        end
    end
    yyaxis left
    plot(log10(chiValues),mean(numCells,2),...
        lineStyles{gdmCtr},'Color',plotColors(1,:),'LineWidth',2);
    yyaxis right
    plot(log10(chiValues),mean(xMax,2),...
        lineStyles{gdmCtr},'Color',plotColors(2,:),'LineWidth',2);
    % add errorbars?
    
end
box on
xlabel('log10(\chi)')
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

filename = ['../manuscripts/JTB/figures/FigS2_contStates_chi_'...
    'sensAcc_' num2str(100*sensingAccuracy)];
set(migrationProfilesFig,'PaperUnits','centimeters');
exportfig(migrationProfilesFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
