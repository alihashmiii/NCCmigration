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
% pcolor(lead2followValues,follow2leadValues,migrationEfficiency)
% shading interp
cb = colorbar; cb.Label.String = 'relative migration efficiency';
caxis([0 1]), colormap(parula(nLevels))
xlabel('switch time lead -> follow, \tau_{L->F} (min)')
ylabel('switch time follow -> lead, \tau_{F->L} (min)')
set(gca,'xtick',axisticks,'ytick',axisticks)
% plot contours showing coefficient of variation
hold on
contour(lead2followValues,follow2leadValues,migrationCOV',[0.15 0.15],'Color',[1 1 1]);
contourIdcs = migrationCOV'>=0.15&migrationCOV'<0.2;
scatter(L2F(contourIdcs),F2L(contourIdcs),'.','MarkerEdgeColor',[1 1 1])
%         [~, covBins] = hist(migrationCOV(:),25);
%         contour(lead2followValues,follow2leadValues,migrationCOV',covBins(covBins>=0.15&covBins<0.2),...
%             'Color',[1 1 1], 'LineStyle', ':')
contour(lead2followValues,follow2leadValues,migrationCOV',[0.2 0.2],...
    'Color',[1 1 1]/2, 'LineWidth', 1)
contourIdcs = migrationCOV'>0.2;
scatter(L2F(contourIdcs),F2L(contourIdcs),'.','MarkerEdgeColor',[1 1 1]/2)
%         contour(lead2followValues,follow2leadValues,migrationCOV',covBins(covBins>=0.2),...
%             'Color',[1 1 1]*0.5, 'LineStyle', ':')
%         contour(lead2followValues,follow2leadValues,migrationCOV',[0.25 0.25],...
%             'Color',[1 1 1]*0, 'LineWidth', 1)
%%
ratioFig = figure;
hold on
plotHandles = NaN(size(lead2followValues'));
for tauCtr = 1:length(lead2followValues)
    if plotError
        plotHandles(tauCtr) = errorbar(lead2followValues./follow2leadValues(tauCtr),migrationEfficiency(tauCtr,:),...
            migrationCOV(tauCtr,:)./sqrt(numReps));
    else
        plotHandles(tauCtr) = plot(lead2followValues./follow2leadValues(tauCtr),migrationEfficiency(tauCtr,:));
        plot(lead2followValues./follow2leadValues(tauCtr),migrationCOV(tauCtr,:),...
            '--','Color',get(plotHandles(tauCtr),'Color'));
    end
end
plot([1 1],[0 1],'k:')
xlim([0 20])
xlabel('relative switch time \tau_{L->F}/\tau_{F->L}')
if plotError
    ylabel('migration efficiency, \mu \pm \mu/\sigma')
else
    ylabel('migration efficiency, \mu (-), \mu/\sigma (--)')
end
legendHandle = legend(plotHandles,[num2str(follow2leadValues'), ...
    repmat(' min',length(follow2leadValues),1)]);
text(15,1,'\tau_{F->L}')
%% export figure
exportOptions = struct('Format','eps2',...
    'Width','10.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = 'manuscripts/VEGF/figures/Fig3F';
set(contourFig,'PaperUnits','centimeters');
exportfig(contourFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);

filename = 'manuscripts/VEGF/figures/Fig3G';
set(ratioFig,'PaperUnits','centimeters');
exportfig(ratioFig,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);