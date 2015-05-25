% plot various migration statistics vs switching time lead2follow
% and follow2lead
% L.J. Schumacher 22.10.14

close all
clear

precision = 2; % significant figures for filenames and plot labels etc.

sensingAccuracyValues = [0.1, 0.01];

% load variables from saved collated results
load('manuscripts/VEGF/figures/experiment31conv4collatedResults.mat','xBins','cellDistributions','actualLeaderFraction','lead2followValues','follow2leadValues','numCells','referenceCellDistribution','referenceLeaderFraction','referenceNumCells')
%% check that we run enough repeats of simulations
% plotMuSigmaErr(1:20,squeeze(actualLeaderNumber(2,2,1,:))','<n>',10)

%% plot the results
% for contour plots, take to take the transpose if data is (x,y), to ensure that
[L2F,F2L] = meshgrid(lead2followValues,follow2leadValues); % for plotting 'patterns' inside contours
nLevels = 20;
actualLeaderNumber = actualLeaderFraction.*numCells;
axisticks = [1 8 16 24 32 40 48 56];
for defaultFollow = [0 1 2]
    for sensAccCtr = 1:2
        subplot(3,2,defaultFollow*2 + sensAccCtr)
        leaderFraction = squeeze(mean(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),5));
        leaderFractionCOV = squeeze(std(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),0,5)...
            ./mean(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),5));
        contourf(lead2followValues,follow2leadValues,leaderFraction',linspace(0,40,nLevels+1),...
            'EdgeColor','none')
        % pcolor(lead2followValues,follow2leadValues,leaderFraction)
        % shading interp
        cb = colorbar; cb.Label.String = 'number of leaders';
        caxis([0 40]), colormap(parula(nLevels))
        xlabel('lead -> follow (min)'), ylabel('follow -> lead (min)')
        set(gca,'xtick',axisticks,'ytick',axisticks)
        % plot contours showing coefficient of variation
        hold on
        contour(lead2followValues,follow2leadValues,leaderFractionCOV',[0.15 0.15],...
            'Color',[1 1 1], 'LineWidth', 1);
        contourIdcs = leaderFractionCOV'>=0.15&leaderFractionCOV'<0.2;
        scatter(L2F(contourIdcs),F2L(contourIdcs),'.','MarkerEdgeColor',[1 1 1])
        %         [~, covBins] = hist(leaderFractionCOV(:),20);
        %         contour(lead2followValues,follow2leadValues,leaderFractionCOV',covBins(covBins>=0.15&covBins<0.2),...
        %             'Color',[1 1 1], 'LineStyle', ':')
        contour(lead2followValues,follow2leadValues,leaderFractionCOV',[0.2 0.2],...
            'Color',[1 1 1]/2, 'LineWidth', 1)
        contourIdcs = leaderFractionCOV'>0.2;
        scatter(L2F(contourIdcs),F2L(contourIdcs),'.','MarkerEdgeColor',[1 1 1]/2)
        %         contour(lead2followValues,follow2leadValues,leaderFractionCOV',covBins(covBins>=0.2),...
        %             'Color',[1 1 1]*0.5, 'LineStyle', ':')
        %         contour(lead2followValues,follow2leadValues,leaderFractionCOV',[0.25 0.25],...
        %             'Color',[1 1 1]*0, 'LineWidth', 1)
        if defaultFollow
            if defaultFollow > 1
                title(['default = follow0, acc. = ' num2str(sensingAccuracyValues(sensAccCtr))])
            elseif defaultFollow == 1
                title(['default = follow, acc. = ' num2str(sensingAccuracyValues(sensAccCtr))])
            elseif ~defaultFollow
                title(['default = leader, acc. = ' num2str(sensingAccuracyValues(sensAccCtr))])
            end
        end
    end
end
    %% export figure
    exportOptions = struct('Format','eps2',...
        'Width','18.0',...
        'Color','rgb',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    
    filename = 'manuscripts/VEGF/figures/FigS4N';
    pos = get(gcf,'Position');
    pos(4) = 4/3*pos(3); % adjust height to fraction of width
    set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
    exportfig(gcf,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);