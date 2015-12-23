% plot various migration statistics vs switching time lead2follow
% and follow2lead
% L.J. Schumacher 22.10.14

close all
clear

precision = 2; % significant figures for filenames and plot labels etc.
numReps = 20;
plotError = 0;

sensingAccuracyValues = [0.1, 0.01];

% load variables from saved collated results
load('manuscripts/VEGF/figures/experiment31conv4collatedResults.mat','xBins','cellDistributions','actualLeaderFraction','lead2followValues','follow2leadValues','numCells','referenceCellDistribution','referenceLeaderFraction','referenceNumCells')
%% check that we run enough repeats of simulations
% plotMuSigmaErr(1:20,squeeze(numCells(2,2,1,:))','<n>',10)

%% plot the results
% for contour plots, take to take the transpose if data is (x,y), to ensure that
% columns are plotted at x-values and rows at y-values (after transpose)
[L2F,F2L] = meshgrid(lead2followValues,follow2leadValues); % for plotting 'patterns' inside contours
nLevels = 20;
axisticks = [1 8 16 24 32 40 48 56];
defaultFollow = [2];
sensAccCtr = 1;
contourFig = figure;
migrationEfficiency = squeeze(mean(numCells(defaultFollow + 1,sensAccCtr,:,:,:),5))/mean(referenceNumCells(:,sensAccCtr));
migrationCOV = squeeze(std(numCells(defaultFollow + 1,sensAccCtr,:,:,:),0,5)...
    ./mean(numCells(defaultFollow + 1,sensAccCtr,:,:,:),5));
contourf(lead2followValues,follow2leadValues,migrationEfficiency',linspace(0,1,nLevels+1),...
    'EdgeColor','none')

cb = colorbar; cb.Label.String = 'relative migration efficiency';
caxis([0 1]), colormap(parula(nLevels))
xlabel('switch time lead \rightarrow follow, \tau_{LF} (min)')
ylabel('switch time follow \rightarrow lead, \tau_{FL} (min)')
set(gca,'xtick',axisticks,'ytick',axisticks)

% plot contours showing coefficient of variation
hold on
[~, contourHandle] = contourf(lead2followValues,follow2leadValues,migrationCOV',[0.2 0.2],'LineColor',[1 1 1]);
hatchHandle = hatchfill2(contourHandle,'single','HatchAngle',45,...
    'HatchSpacing',10,'Fill','off');

if any(migrationCOV(:)>0.3)
    [~, contourHandle2] = contourf(lead2followValues,follow2leadValues,migrationCOV',[0.3 0.3],'LineColor',[1 1 1]/2);
    hatchHandle2 = hatchfill2(contourHandle2,'single','HatchAngle',-45,...
        'HatchSpacing',10,'Fill','off');
end
%%
ratioFig = figure;
plotHandles = NaN(length(lead2followValues')+1,1);
plotHandles(1) = plot(0,0,'.','Color',[1 1 1]);
hold on
for tauCtr = 1:length(lead2followValues)
    if plotError
        plotHandles(tauCtr+1) = errorbar(lead2followValues./follow2leadValues(tauCtr),migrationEfficiency(tauCtr,:),...
            migrationCOV(tauCtr,:)./sqrt(numReps));
    else
        plotHandles(tauCtr+1) = plot(lead2followValues./follow2leadValues(tauCtr),migrationEfficiency(tauCtr,:));
        plot(lead2followValues./follow2leadValues(tauCtr),migrationCOV(tauCtr,:),...
            '--','Color',get(plotHandles(tauCtr+1),'Color'));
    end
end
plot([1 1],[0 1],'k:')
xlim([0 20])
xlabel('relative switch time \tau_{LF}/\tau_{FL}')
if plotError
    ylabel('migration efficiency, \mu \pm \mu/\sigma')
else
    ylabel('migration efficiency, \mu (-), \mu/\sigma (--)')
end
legendHandle = legend(plotHandles,['\tau_{FL}';mat2cell([num2str(follow2leadValues'),repmat(' min',length(follow2leadValues),1)],ones(length(follow2leadValues),1),6)]);
% text(15,1,'\tau_{FL}')
box on
%% export figure
exportOptions = struct('Format','eps2',...
    'Width','14.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = '../../Thesis/figures/migrationEfficiencySwitchingTimesFull';
set(contourFig,'PaperUnits','centimeters');
exportfig(contourFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);

filename = '../../Thesis/figures/migrationEfficiencySwitchingTimes';
set(ratioFig,'PaperUnits','centimeters');
exportfig(ratioFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);