% plot various migration statistics vs switching time lead2follow
% and follow2lead
% L.J. Schumacher 22.10.14, 15.08.15

close all
clear

precision = 2; % significant figures for filenames and plot labels etc.

sensingAccuracyValues = [0.1, 0.01];

% load variables from saved collated results
load('manuscripts/VEGF/figures/experiment31conv4collatedResults.mat','xBins','cellDistributions','actualLeaderFraction','lead2followValues','follow2leadValues','numCells','referenceCellDistribution','referenceLeaderFraction','referenceNumCells')

%% check that we run enough repeats of simulations
% plotMuSigmaErr(1:20,squeeze(actualLeaderFraction(2,2,1,:))','<n>',10)

%% plot the results
% for contour plots, take to take the transpose if data is (x,y), to ensure that
% columns are plotted at x-values and rows at y-values (after transpose)
[L2F,F2L] = meshgrid(lead2followValues,follow2leadValues); % for plotting 'patterns' inside contours
nLevels = 10;
plotColors = linspace(0,1,nLevels)'*[251 101 4]/255 + (1 - linspace(0,1,nLevels)')*[113 18 160]/255;
actualLeaderNumber = actualLeaderFraction.*numCells;
axisticks = [1 8 16 24 32 40 48 56];
for defaultFollow = [2]
    for sensAccCtr = 1
        fractionFig = figure;
        leaderFraction = squeeze(mean(actualLeaderFraction(defaultFollow + 1,sensAccCtr,:,:,:),5));
        leaderFractionCOV = squeeze(std(actualLeaderFraction(defaultFollow + 1,sensAccCtr,:,:,:),0,5)...
            ./mean(actualLeaderFraction(defaultFollow + 1,sensAccCtr,:,:,:),5));
        contourf(lead2followValues,follow2leadValues,leaderFraction',linspace(0,0.5,nLevels+1),...
            'EdgeColor','none')

        cb = colorbar; cb.Label.String = 'leader fraction';
        caxis([0 0.5]), colormap(plotColors)
        xlabel('switch time lead -> follow, \tau_{LF} (min)')
        ylabel('switch time follow -> lead, \tau_{FL} (min)')
        set(gca,'xtick',axisticks,'ytick',axisticks)
%         % plot contours showing coefficient of variation
%         hold on
%         [~, contourHandle] = contourf(lead2followValues,follow2leadValues,...
%             leaderFractionCOV',[0.2 0.2],'LineColor',[1 1 1]);
%         hatchHandle = hatchfill2(contourHandle,'single','HatchAngle',45,...
%             'HatchSpacing',10,'Fill','off');
% 
%         [~, contourHandle2] = contourf(lead2followValues,follow2leadValues,...
%             leaderFractionCOV',[0.3 0.3],'LineColor',[1 1 1]/2);
%         hatchHandle2 = hatchfill2(contourHandle2,'single','HatchAngle',-45,...
%             'HatchSpacing',10,'Fill','off');
%         
        numberFig = figure;
        leaderNumber = squeeze(mean(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),5));
        leaderNumberCOV = squeeze(std(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),0,5)...
            ./mean(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),5));
        contourf(lead2followValues,follow2leadValues,leaderNumber',linspace(0,25,nLevels+1),...
            'EdgeColor','none')

        cb = colorbar; cb.Label.String = 'number of leaders';
        caxis([0 25]), colormap(parula(nLevels))
        xlabel('switch time lead -> follow, \tau_{LF} (min)')
        ylabel('switch time follow -> lead, \tau_{FL} (min)')
        set(gca,'xtick',axisticks,'ytick',axisticks)
%         % plot contours showing coefficient of variation
%         hold on
%         [~, contourHandle] = contourf(lead2followValues,follow2leadValues,...
%           leaderNumberCOV',[0.2 0.2],'LineColor',[1 1 1]);
%         hatchHandle = hatchfill2(contourHandle,'single','HatchAngle',45,...
%                         'HatchSpacing',10,'Fill','off');
% 
%         [~, contourHandle2] = contourf(lead2followValues,follow2leadValues,...
%           leaderNumberCOV',[0.3 0.3],'LineColor',[1 1 1]/2);
%         hatchHandle2 = hatchfill2(contourHandle2,'single','HatchAngle',-45,...
%                         'HatchSpacing',10,'Fill','off');
    end
end
%% export figure
exportOptions = struct('Format','eps2',...
    'Width','14.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = '../../Thesis/figures/switchingTimesLeaderFraction';
pos = get(fractionFig,'Position');
set(fractionFig,'PaperUnits','centimeters','Position',pos,'color','none');
exportfig(fractionFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);

filename = '../../Thesis/figures/switchingTimesLeaderNumber';
pos = get(numberFig,'Position');
set(numberFig,'PaperUnits','centimeters','Position',pos,'color','none');
exportfig(numberFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);