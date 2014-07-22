% load & plot distributions of cells as point-cloud-bar-charts
% L.J. Schumacher 21.05.14

close all
clear

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;

sensingAccuracy = 0.1;
needNeighbours = 0;

cellRadius = 7.5;
opacity = 3/numRepeats;

experiments = {'experiment36/exp36'; 'experiment31/exp31_followFrac_1'};

actualLeaderFraction = NaN(length(experiments),numRepeats);
lostFraction = NaN(length(experiments),numRepeats);
distanceMigrated = cell(size(experiments));
totalCellNumbers = NaN(size(experiments));

scatterFig = figure('Visible','on');
axis equal
hold all
markerSize = 4; % area in points squared
 
%% load and plot data for every run of this parameter combination
for expCtr = 1:length(experiments)
    yOffset = (expCtr - 1)*120;
    for repCtr = 1:numRepeats
                loadInfo = [experiments{expCtr} ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours) ...
                    '_Run_' num2str(repCtr)];
        try % sometimes we get corrupt files, which crashes the script
            load(['results/' loadInfo '.mat'])
        catch
            error(['Could not load results/' loadInfo '.mat'])
        end
        % make a scatter plot of all repeats on top of each other
        cells = out.cells_save{end}; % all cells
        numberOfCells = size(cells,2);
        if expCtr==1 % an inconvenient if-statement to deal with a change in naming of a variable - urgh.
            followIdcs = out.cellsFollow_save{end}(1:numberOfCells);
        elseif expCtr==2
            followIdcs = out.cellsFollow{end}(1:numberOfCells);
        end
        attachIdcs = out.attach_save{end}(1:numberOfCells);
        leaders = cells(:,followIdcs==0);
        followers = cells(:,followIdcs==1&attachIdcs~=0);
        losts = cells(:,followIdcs==1&attachIdcs==0);
%         scatter(leaders(1,:),leaders(2,:) + yOffset,markerSize,[251 101 4]/255,'fill')
%         scatter(followers(1,:),followers(2,:) + yOffset,markerSize,[113 18 160]/255,'fill')
%         scatter(losts(1,:),losts(2,:) + yOffset,markerSize,0.5*[1 1 1],'fill')
        plotCircles([leaders(1,:); leaders(2,:) + yOffset],cellRadius,[251 101 4]/255,opacity);
        plotCircles([followers(1,:); followers(2,:) + yOffset],cellRadius,[113 18 160]/255,opacity);
        plotCircles([losts(1,:); losts(2,:) + yOffset],cellRadius,0*[1 1 1],opacity);

        actualLeaderFraction(expCtr,repCtr) = size(leaders,2)/numberOfCells;
        lostFraction(expCtr,repCtr) = size(losts,2)/numberOfCells;
        distanceMigrated{expCtr} = [distanceMigrated{expCtr}, cells(1,:)];
    end
    % record the total number of data points for later use
    totalCellNumbers(expCtr) = length(distanceMigrated{expCtr});
end
%% boxplot the distribution of distances migrated over all repeats
yTicks = linspace(60,120*length(experiments)-60,length(experiments));
maxCellNumber = max(totalCellNumbers);
for expCtr = 1:length(experiments)
    boxplot(distanceMigrated{expCtr},'plotstyle','compact','medianstyle','target','orientation','horizontal',...
        'colors',[1 1 1]*2*mean(lostFraction(expCtr,:)),'symbol','x','positions',yTicks(expCtr))
end   
xlabel('x/\mum')
set(gca,'YTick',yTicks,'YTickLabel', num2str(mean(actualLeaderFraction,2),precision))
xlim([0, 800])
ylim([0, max(yTicks) + 60])

% add a colorbar for relative stream density
ch = colorbar;
title(ch,'stream cohesion')
colormap(linspace(1,0)'*[1 1 1])
set(ch,'YTickLabel',[0.5 0.75 1])

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
set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
filename = ['manuscripts/subpopulations/figures/resultsFigS2A'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
