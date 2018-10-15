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
%         figure
        subplot(3,2,defaultFollow*2 + sensAccCtr)
        leaderNumber = squeeze(mean(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),5));
        leaderNumberCOV = squeeze(std(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),0,5)...
            ./mean(actualLeaderNumber(defaultFollow + 1,sensAccCtr,:,:,:),5));
        contourf(lead2followValues,follow2leadValues,leaderNumber',linspace(0,40,nLevels+1),...
            'EdgeColor','none')
        % pcolor(lead2followValues,follow2leadValues,leaderFraction)
        % shading interp
        cb = colorbar; cb.Label.String = 'number of leaders';
        caxis([0 40]), colormap(parula(nLevels))
        xlabel('lead -> follow (min)'), ylabel('follow -> lead (min)')
        set(gca,'xtick',axisticks,'ytick',axisticks)
        % plot contours showing coefficient of variation
        hold on
        [~, contourHandle] = contourf(lead2followValues,follow2leadValues,leaderNumberCOV',[0.2 0.2],'LineColor',[1 1 1]);
        hatchHandle = hatchfill2(contourHandle,'single','HatchAngle',45,...
                        'HatchSpacing',10,'Fill','off');

        [~, contourHandle2] = contourf(lead2followValues,follow2leadValues,leaderNumberCOV',[0.3 0.3],'LineColor',[1 1 1]/2);
        hatchHandle2 = hatchfill2(contourHandle2,'single','HatchAngle',-45,...
                        'HatchSpacing',10,'Fill','off');
        
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