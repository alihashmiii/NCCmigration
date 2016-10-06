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
slowSpeeds = [10 20 40];

fileName = 'exp42_slowEntryDynamicDAN';

plotMarkers = {'-+','-s','-o'};
% to calculate the density profile of cells along the x-direction
dx = 50;
xBins = 0:dx:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(length(slowSpeeds),numReps,length(xBins));

exportOptions = struct('Format','eps2',...
    'Width','13.8',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2,...
    'renderer','painters');

for cntGdn = {'toward', 'parallel'}
    contactGuidance = char(cntGdn);
    for diffus = diffusivities
        sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        figure
        hold on
        for slowSpeedCtr = 1:length(slowSpeeds)
            slowSpeed = slowSpeeds(slowSpeedCtr);
            for repCtr = 1:numReps
                loadInfo = [fileName '/' fileName ...
                    '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
                    '_slowSpeed_' num2str(slowSpeed) '_contactGuidance_' ...
                    contactGuidance '_Run_' num2str(repCtr)];
                load(['results/' loadInfo '.mat'])
                
                % find time index to plot
                timeIdx = find(out.t_save>timeToPlot,1,'first');
                % load cell positions into variables
                cells = out.cells_save{timeIdx}; % all cells

                % calculate migration profile
                cellDistributions(slowSpeedCtr,repCtr,:) = histc(cells(1,:),xBins); % leaders
            end
            
            %% plot migration profile
            % plot lines for migration profile (simplest, if less accurate)
            meanNumCells = squeeze(mean(cellDistributions(slowSpeedCtr,:,:),2));
            stdNumCells = squeeze(std(cellDistributions(slowSpeedCtr,:,:),[],2));
            errorbar(xBins + dx/2,meanNumCells,stdNumCells./sqrt(numReps),...
                plotMarkers{slowSpeedCtr});
        end
        xlabel('x/$\mu$m','interpreter','latex')
        ylabel('number of cells / 50$\,\mu$m','interpreter','latex')
        hl = legend(num2str(slowSpeeds'));
        hl.Title.String = 'v_{DAN} (\mum/hr)';
        hl.Title.FontWeight = 'normal';
        xlim([0 800])
        grid off, set(gca,'Layer','top')
        box on
        %% save figure as .fig file
        savename = ['~/Dropbox/projects/cellMigration/DAN/figures/' fileName ...
            '_D_' num2str(diffus)  '_sensingAcc_' num2str(sensingAccuracy,precision) ...
            '_contactGuidance_' contactGuidance '_migrationProfiles'];
        saveas(gcf,[savename '.fig'])
        
        %% export figure
        pos = get(gcf,'Position');
        pos(4) = 1/2*pos(3); % adjust height to fraction of width
        set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
        exportfig(gcf,[savename '.eps'],exportOptions);
        system(['epstopdf ' savename '.eps']);
    end
end