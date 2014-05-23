% load & plot distributions of cells as point-cloud-bar-charts
% L.J. Schumacher 21.05.14

close all
clear

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;

followFracValues = [0, 3/4, 7/8, 15/16, 1];
sensingAccuracy = 0.1;
needNeighbours = 0;

actualLeaderFraction = NaN(length(followFracValues),numRepeats);
distanceMigrated = cell(size(followFracValues));
totalCellNumbers = NaN(size(followFracValues));

scatterFig = figure('Visible','on');
axis equal
hold all
markerSize = 176.7146; % area in points squared
 
%% load and plot data for every run of this parameter combination
for followFracCtr = 1:length(followFracValues)
    followerFraction = followFracValues(followFracCtr);
    yOffset = (followFracCtr - 1)*120;
    for repCtr = 1:numRepeats
                loadInfo = ['experiment31/exp31_followFrac_' num2str(followerFraction,precision) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
                    '_Run_' num2str(repCtr)];
        try % sometimes we get corrupt files, which crashes the script
            load(['results/' loadInfo '.mat'])
        catch
            delete(['results/' loadInfo '.mat']) % delete the corrupt file
            experiment311leaderFractionWithoutPlasticity; % recreate the missing results file
            load(['results/' loadInfo '.mat']) % load again
        end
        % make a scatter plot of all repeats on top of each other
        cells = out.cells_save{end}; % all cells
        numberOfCells = size(cells,2);
        followIdcs = out.cellsFollow{end}(1:numberOfCells);
        attachIdcs = out.attach_save{end}(1:numberOfCells);
        leaders = cells(:,followIdcs==0);
        followers = cells(:,followIdcs==1&attachIdcs~=0);
        losts = cells(:,followIdcs==1&attachIdcs==0);
        scatter(leaders(1,:),leaders(2,:) + yOffset,markerSize,[251 101 4]/255,'fill')
        scatter(followers(1,:),followers(2,:) + yOffset,markerSize,[113 18 160]/255,'fill')
        scatter(losts(1,:),losts(2,:) + yOffset,markerSize,0.5*[1 1 1],'fill')
        
        actualLeaderFraction(followFracCtr,repCtr) = size(leaders,2)/numberOfCells;
        distanceMigrated{followFracCtr} = [distanceMigrated{followFracCtr}, cells(1,:)];
    end
    % record the total number of data points for later use
    totalCellNumbers(followFracCtr) = length(distanceMigrated{followFracCtr});
end
%% boxplot the distribution of distances migrated over all repeats
yTicks = linspace(60,120*length(followFracValues)-60,length(followFracValues));
maxCellNumber = max(totalCellNumbers);
for followFracCtr = 1:length(followFracValues)
    boxplot(distanceMigrated{followFracCtr},'plotstyle','compact','medianstyle','target','orientation','horizontal',...
        'colors',[1 1 1] - totalCellNumbers(followFracCtr)/maxCellNumber,'symbol','x','positions',yTicks(followFracCtr))
end   
xlabel('x/\mum')
ylabel('f_L')
set(gca,'YTick',yTicks,'YTickLabel', num2str(mean(actualLeaderFraction,2),precision))
xlim([0, 800])
ylim([0, max(yTicks) + 60])

% add a colorbar for relative stream density
ch = colorbar;
title(ch,'\rho^{stream}_{relative}')
colormap(linspace(1,0)'*[1 1 1])

%% export figure
exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

pos = get(gcf,'Position');
%     pos(4) = 3/2*pos(3);% adjust height to 3/2 width
set(gcf,'PaperUnits','centimeters','Position',pos);
filename = ['manuscripts/subpopulations/figures/resultsFig1B'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
