% load & plot parameter sensitivity of number of cells and distance
% migrated
% L.J. Schumacher 17.02.2015

close all
clear
%%
time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;

conversionType = 0;
followerFractionValues = [0, 1];
sensingAccuracy = 0.1;
needNeighbours = 0;

dataSets = {'experiment31/';
    'experiment31Eat/';
    'experiment31D/';
    'experiment31Chi/';
    'experiment31halfFilo/';
    'experiment37reducedChemotaxis/';
    'experiment31/'};
labels = {'reference';'0.1\lambda';'1e6D';'1e4\chi';'0.5(l_{filo} - R)';'0.5n_{filo}';'0.01\Deltac/c'};
plotMarkers = {'o','^','>','x','<','v','+'};
plotColors = lines(length(dataSets));

numberOfCells = NaN(length(dataSets),length(followerFractionValues),numRepeats);
distanceMigrated = NaN(length(dataSets),length(followerFractionValues),numRepeats);
plotHandles = NaN(size(dataSets));

mainfigure = figure;
subplot(3,4,[1:3,5:7])
hold on
for dataSetCtr = length(dataSets):-1:1
    for followFracCtr = 1:length(followerFractionValues)
        followerFraction = followerFractionValues(followFracCtr);
        % file names differ for different data sets, depend on looped params
        fileNames = {['exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours)];
            ['exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy)];
            ['exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy/1e3)];
            ['exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy)];
            ['exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy)];
            ['exp37_conversion_' num2str(conversionType) '_followFrac_' num2str(followerFraction) ...
            '_sensingAcc_' num2str(sensingAccuracy)];
            ['exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_0.001_needNeighbours_' num2str(needNeighbours)]};
        % loop through realisations of particular parameter combination
        for repCtr = 1:numRepeats
            loadInfo = [dataSets{dataSetCtr} fileNames{dataSetCtr} '_Run_' num2str(repCtr)];
            load(['results/' loadInfo '.mat'])
            
            % load cell positions into variables
            cells = out.cells_save{end}; % all cells
            numberOfCells(dataSetCtr,followFracCtr,repCtr) = size(cells,2);
            distanceMigrated(dataSetCtr,followFracCtr,repCtr) = max(cells(1,:));
        end
        % scatter plot individual simulation outcomes
        if followerFraction==1
            plotHandles(dataSetCtr) = plot(squeeze(distanceMigrated(dataSetCtr,followFracCtr,:)),...
                squeeze(numberOfCells(dataSetCtr,followFracCtr,:)),plotMarkers{dataSetCtr},'Color',plotColors(dataSetCtr,:));
        else
            plot(squeeze(distanceMigrated(dataSetCtr,followFracCtr,:)),...
                squeeze(numberOfCells(dataSetCtr,followFracCtr,:)),...
                plotMarkers{dataSetCtr},'Color',plotColors(dataSetCtr,:)/2*(1 + followerFraction));
        end
    end
    lineHandle = plot(mean(distanceMigrated(dataSetCtr,:,:),3),mean(numberOfCells(dataSetCtr,:,:),3),...
        'LineWidth',2);
    % undocumented trick for shaded lines: http://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
    drawnow
    set(lineHandle.Edge, 'ColorBinding','interpolated','colordata',...
        uint8(255*[plotColors(dataSetCtr,:)'*(1 + followerFractionValues)/2; ones(size(followerFractionValues))]))
end
legend(plotHandles,labels,'Location','NorthWest')
xlabel('furthest distance migrated (\mum)')
ylabel('number of cells')
box on
%% plot individual parameter perturbations in separate panels
panelIndcs = [4,8,9:12];
for dataSetCtr = 1:length(dataSets)
    if dataSetCtr == 1;
        for panelCtr = 1:6
            subplot(3,4,panelIndcs(panelCtr))
            hold on
            for followFracCtr = 1:length(followerFractionValues)
                followerFraction = followerFractionValues(followFracCtr);
                % scatter plot individual simulation outcomes
                plot(squeeze(distanceMigrated(dataSetCtr,followFracCtr,:)),...
                    squeeze(numberOfCells(dataSetCtr,followFracCtr,:)),...
                    plotMarkers{dataSetCtr},'Color',plotColors(dataSetCtr,:)/2*(1 + followerFraction));
            end
            lineHandle = plot(mean(distanceMigrated(dataSetCtr,:,:),3),mean(numberOfCells(dataSetCtr,:,:),3),...
                'LineWidth',2);
            % undocumented trick for shaded lines: http://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
            drawnow
            set(lineHandle.Edge, 'ColorBinding','interpolated','colordata',...
                uint8(255*[plotColors(dataSetCtr,:)'*(1 + followerFractionValues)/2; ones(size(followerFractionValues))]))
%             xlabel('distance (\mum)')
%             ylabel('number of cells')
            ylim([20 160])
            set(gca,'ytick',[40:40:160])
            box on
        end
    else
        subplot(3,4,panelIndcs(dataSetCtr-1))
        for followFracCtr = 1:length(followerFractionValues)
            followerFraction = followerFractionValues(followFracCtr);
            % scatter plot individual simulation outcomes
            plot(squeeze(distanceMigrated(dataSetCtr,followFracCtr,:)),...
                squeeze(numberOfCells(dataSetCtr,followFracCtr,:)),...
                plotMarkers{dataSetCtr},'Color',plotColors(dataSetCtr,:)/2*(1 + followerFraction));
        end
        lineHandle = plot(mean(distanceMigrated(dataSetCtr,:,:),3),mean(numberOfCells(dataSetCtr,:,:),3),...
            'LineWidth',2);
        % undocumented trick for shaded lines: http://undocumentedmatlab.com/blog/plot-line-transparency-and-color-gradient
        drawnow
        set(lineHandle.Edge, 'ColorBinding','interpolated','colordata',...
            uint8(255*[plotColors(dataSetCtr,:)'*(1 + followerFractionValues)/2; ones(size(followerFractionValues))]))
    end
end
%% export figure
exportOptions = struct('Format','eps2',...
    'Width','18',...
    'Height','18',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = '../../Thesis/figures/parameterSensitivityNoPlasticity';
set(mainfigure,'PaperUnits','centimeters');
exportfig(mainfigure,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);