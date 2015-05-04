% load & plot distributions of cells as overlapping bar-charts
% L.J. Schumacher 21.05.14

close all
clear
%%
time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

followFracValues = [0, 3/4, 7/8, 15/16, 1];
plotMarkers = {'-x','-.','-s','-+','-o'};
sensingAccuracy = 0.1;
needNeighbours = 0;

% to calculate the density profile of cells along the x-direction
dx = 50;
xBins = 0:dx:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(length(followFracValues),numRepeats,3,length(xBins));

% preallocate variables for saving collated results
actualLeaderFraction = NaN(length(followFracValues),numRepeats);

figure
hold all
for followFracCtr = length(followFracValues):-1:1
    followerFraction = followFracValues(followFracCtr);
    for repCtr = 1:numRepeats
        loadInfo = ['experiment31/exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
            '_Run_' num2str(repCtr)];
        try % sometimes we get corrupt files, which crashes the script
            load(['results/' loadInfo '.mat'])
        catch
            delete(['results/' loadInfo '.mat']) % delete the corrupt file
            experiment31leaderFractionWithoutPlasticity; % recreate the missing results file
            load(['results/' loadInfo '.mat']) % load again
        end
        
        % load cell positions into variables
        cells = out.cells_save{end}; % all cells
        numberOfCells = size(cells,2);
        if isfield(out,'cellsFollow')
            followIdcs = out.cellsFollow{end}(1:numberOfCells);
        else
            followIdcs = out.cellsFollow_save{end}(1:numberOfCells);
        end
        attachIdcs = out.attach_save{end}(1:numberOfCells);
        leaders = cells(:,followIdcs==0);
        followers = cells(:,followIdcs==1&attachIdcs~=0);
        losts = cells(:,followIdcs==1&attachIdcs==0);
        
        actualLeaderFraction(followFracCtr,repCtr) = size(leaders,2)/numberOfCells;

        % calculate migration profile
        numberOfCells = size(out.cells_save{end},2);
        cellDistributions(followFracCtr,repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
        cellDistributions(followFracCtr,repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
        cellDistributions(followFracCtr,repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached
    end
    %% plot migration profile
    f_L = mean(actualLeaderFraction(followFracCtr,:));
    plotColor = f_L*[251 101 4]/255 + (1 - f_L)*[113 18 160]/255;
    % plot lines for migration profile (simplest, if less accurate)
    plot(xBins + dx/2,squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2)),plotMarkers{followFracCtr},'color',plotColor);
%     %   plot as area-graph instead of bar as matlab will sometimes get the alpha wrong for barplots (on mac) 
%     [xs, ys] = stairs(xBins,squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2)));
%     h(followFracCtr) = area(xs,ys,...
%             'FaceColor', plotColor, 'EdgeColor', plotColor);
%     % alternatively plot regions of sem
%     h(followFracCtr) = area([xBins+dx/2, fliplr(xBins+dx/2)],...
%     [squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2))+std(squeeze(sum(cellDistributions(followFracCtr,:,:,:),3)))'/sqrt(numRepeats);...
%             flipud(squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2))-std(squeeze(sum(cellDistributions(followFracCtr,:,:,:),3)))'/sqrt(numRepeats))],...
%     'FaceColor', plotColor, 'EdgeColor', plotColor);
%     alpha(get(h(followFracCtr),'children'),0.5)
end
    
xlabel('x/\mum')
ylabel('# cells / 50\mum'), 
legend(num2str(flipud(mean(actualLeaderFraction,2)),precision),'Location',...
    'NorthEast'); 
text(675,10,'\langle f_L\rangle')
ylim([0 17]), xlim([0 800]), set(gca,'YTick',[0 4 8 12 16])
grid off, set(gca,'Layer','top')
box on
%% save figure as .fig file
filename = 'manuscripts/subpopulations/figures/resultsFig1C';
saveas(gcf,[filename '.fig'])

%% export figure
exportOptions = struct('Format','eps2',...
    'Width','13.8',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

pos = get(gcf,'Position');
pos(4) = 1/2*pos(3); % adjust height to fraction of width
set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
system(['pdfcrop ' filename '.pdf']);