% load & plot distributions of cells as overlapping bar-charts
% L.J. Schumacher 21.05.14

close all
clear

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;

followerFraction = [1];
sensingAccuracy = 0.1;
needNeighbours = 0;

experiments = {['experiment35/exp35_']; ['experiment31/exp31_']};

% to calculate the density profile of cells along the x-direction
xBins = 0:50:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(length(experiments),numRepeats,3,length(xBins));

% preallocate variables for saving collated results
actualLeaderFraction = NaN(length(experiments),numRepeats);

figure
hold all
for expCtr = length(experiments):-1:1
    for repCtr = 1:numRepeats
        loadInfo = [experiments{expCtr} 'followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
            '_Run_' num2str(repCtr)];
        try % sometimes we get corrupt files, which crashes the script
            load(['results/' loadInfo '.mat'])
        catch
            error(['Could not load results/' loadInfo '.mat'])
        end
        
        % load cell positions into variables
        cells = out.cells_save{end}; % all cells
        numberOfCells = size(cells,2);
        followIdcs = out.cellsFollow{end}(1:numberOfCells);
        attachIdcs = out.attach_save{end}(1:numberOfCells);
        leaders = cells(:,followIdcs==0);
        followers = cells(:,followIdcs==1&attachIdcs~=0);
        losts = cells(:,followIdcs==1&attachIdcs==0);
        
        actualLeaderFraction(expCtr,repCtr) = size(leaders,2)/numberOfCells;

        % calculate migration profile
        cellDistributions(expCtr,repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
        cellDistributions(expCtr,repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
        cellDistributions(expCtr,repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached
    end
    % plot migration profile
    f_L = mean(actualLeaderFraction(expCtr,:));
    plotColor = f_L*[251 101 4]/255 + (1 - f_L)*[113 18 160]/255;
    stairs(xBins,squeeze(mean(sum(cellDistributions(expCtr,:,:,:),3),2)),'Color',plotColor)
    h(expCtr) = bar(xBins + 25,squeeze(mean(sum(cellDistributions(expCtr,:,:,:),3),2)),1,...
        'FaceColor', plotColor, 'EdgeColor', plotColor);
    %     h(followFracCtr) = area(xBins + 25,squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2)),...
    %         'FaceColor', plotColor, 'EdgeColor', plotColor);
%     alpha(get(h(expCtr),'children'),0.5)
end
xlabel('x/\mum')
ylabel('# cells / 50\mum'), 
legend(h, {'lead+'; 'WT'})
ylim([0 16]), xlim([0 800]), set(gca,'YTick',[0 4 8 12 16])
grid on, set(gca,'Layer','top')
box on

%% export figure
exportOptions = struct('Format','eps2',...
    'Width','14.4',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

pos = get(gcf,'Position');
pos(4) = 1/2*pos(3); % adjust height to fraction of width
set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
filename = ['manuscripts/subpopulations/figures/resultsFig2C'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);