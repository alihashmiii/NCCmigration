% plot result of simulation of neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width
% has a zone of reduced cell speeds at the start (grows with tissue)
% DAN zone near domain entrance slows down cells, increases then decreases
% with time linearly

clear
close all
%%
timeToPlot = 24;
numReps = 10;

precision = 2;

diffusivities = [1];
slowSpeeds = [30];
insertEveryStepsValues = [6 10 15 30]; % corresponding to 6, 4 and 2 (attempted) cell insertions per hour at tstep = 1min
numInsertValues = length(insertEveryStepsValues);

fileName = 'exp42_slowEntryDynamicDAN';

plotColors = summer(numInsertValues);

exportOptions = struct('Format','eps2',...
    'Width','13.8',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

for cntGdn = {'parallel'}
    contactGuidance = char(cntGdn);
    for diffus = diffusivities
        sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        for slowSpeedCtr = 1:length(slowSpeeds)
            slowSpeed = slowSpeeds(slowSpeedCtr);
            danZoneFig = figure; hold on
            streamWidthFig = figure; hold on
            cellFractions = NaN(numInsertValues,numReps);
            streamWidths = NaN(numInsertValues,numReps);
            refCellFractions = NaN(numInsertValues,numReps);
            refStreamWidths = NaN(numInsertValues,numReps);
            for insertEveryStepsCtr = 1:numInsertValues
                insertEverySteps = insertEveryStepsValues(insertEveryStepsCtr);
                %% load reference distribution
                for repCtr = 1:numReps
                    if insertEverySteps==6
                        loadInfo = ['exp39_widerDomainStripe/exp39_widerDomainStripe' ...
                            '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                            '_speed_40_contactGuidance_' contactGuidance '_Run_' num2str(repCtr)];
                    else
                        loadInfo = ['exp39_widerDomainStripe/exp39_widerDomainStripe' ...
                            '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                            '_speed_40_insertEvry_' num2str(insertEverySteps) ...
                            '_contactGuidance_' contactGuidance '_Run_' num2str(repCtr)];
                    end
                    load(['results/' loadInfo '.mat'])
                    % find time index to plot
                    timeIdx = find(out.t_save>timeToPlot,1,'first');
                    % load cell positions into variables
                    cells = out.cells_save{timeIdx}; % all cells
                    % calculate cell fractions in/out of DAN zone vs time
                    refCellFractions(insertEveryStepsCtr,repCtr) = nnz(cells(1,:)<=max(out.xlat_save{timeIdx})/3)/size(cells,2);
                    % calc stream width
                    refStreamWidths(insertEveryStepsCtr,repCtr) = std(cells(2,cells(1,:)<=max(out.xlat_save{timeIdx})/3));
                end
                %% load DAN data
                for repCtr = 1:numReps
                    if insertEverySteps==6
                        loadInfo = [fileName '/' fileName ...
                            '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                            '_slowSpeed_' num2str(slowSpeed) ...
                            '_contactGuidance_' contactGuidance '_Run_' num2str(repCtr)];
                    else
                        loadInfo = [fileName '/' fileName ...
                            '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                            '_slowSpeed_' num2str(slowSpeed) '_insertEvry_' num2str(insertEverySteps)...
                            '_contactGuidance_' contactGuidance '_Run_' num2str(repCtr)];
                    end
                    load(['results/' loadInfo '.mat'])
                    % find time index to plot
                    timeIdx = find(out.t_save>timeToPlot,1,'first');
                    % load cell positions into variables
                    cells = out.cells_save{timeIdx}; % all cells
                    % calculate cell fractions in/out of DAN zone vs time
                    cellFractions(insertEveryStepsCtr,repCtr) = nnz(cells(1,:)<=max(out.xlat_save{timeIdx})/3)/size(cells,2);
                    % calc stream width
                    streamWidths(insertEveryStepsCtr,repCtr) = std(cells(2,cells(1,:)<=max(out.xlat_save{timeIdx})/3));
                end
            end
            %
%             errorbar(danZoneFig.Children,60./insertEveryStepsValues,mean(refCellFractions,2),std(cellFractions,0,2)/sqrt(numReps))
            errorbar(danZoneFig.Children,60./insertEveryStepsValues,mean(cellFractions,2),std(cellFractions,0,2)/sqrt(numReps))
            xlabel(danZoneFig.Children,'attempted cell insertions (per hour)')
            danZoneFig.Children.XDir = 'reverse';
            ylabel(danZoneFig.Children,'fraction of cells in DAN zone')
            xlim(danZoneFig.Children,[1 10])
            danZoneFig.Children.XGrid = 'on';
            danZoneFig.Children.YGrid = 'on';
            danZoneFig.Children.Box = 'on';
            danZoneFig.Children.YLim = [0.82, 0.94];
            danZoneFig.Children.YTickLabel = num2str(danZoneFig.Children.YTick','%1.2f');
%             hl = legend(danZoneFig.Children,{'no DAN';['v_{DAN} = ' num2str(slowSpeed)]},...
%                 'Location','NorthWest');
%
            errorbar(streamWidthFig.Children,60./insertEveryStepsValues,mean(refStreamWidths,2),std(refStreamWidths,0,2)/sqrt(numReps))
            errorbar(streamWidthFig.Children,60./insertEveryStepsValues,mean(streamWidths,2),std(streamWidths,0,2)/sqrt(numReps))
            xlabel(streamWidthFig.Children,'attempted cell insertions (per hour)')
            streamWidthFig.Children.XDir = 'reverse';
            ylabel(streamWidthFig.Children,'stream width in DAN zone')
            xlim(streamWidthFig.Children,[1 10])
            streamWidthFig.Children.XGrid = 'on';
            streamWidthFig.Children.YGrid = 'on';
            streamWidthFig.Children.Box = 'on';
            hl = legend(streamWidthFig.Children,{'no DAN';['v_{DAN} = ' num2str(slowSpeed)]},...
                'Location','NorthWest');
            %% save figure as .fig file
            savename = ['~/Dropbox/projects/cellMigration/DAN/figures/' fileName ...
                '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                '_slowSpeed_' num2str(slowSpeed) ...
                '_contactGuidance_' contactGuidance '_danZoneNumbers'];
            saveas(danZoneFig,[savename '.fig'])
            %
            set(danZoneFig,'PaperUnits','centimeters','color','none');
            exportfig(danZoneFig,[savename '.eps'],exportOptions);
            system(['epstopdf ' savename '.eps']);
            %
            savename = ['~/Dropbox/projects/cellMigration/DAN/figures/' fileName ...
                '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                '_slowSpeed_' num2str(slowSpeed) ...
                '_contactGuidance_' contactGuidance '_streamWidths'];
            saveas(streamWidthFig,[savename '.fig'])
            %
            set(streamWidthFig,'PaperUnits','centimeters','color','none');
            exportfig(streamWidthFig,[savename '.eps'],exportOptions);
            system(['epstopdf ' savename '.eps']);
        end
    end
end