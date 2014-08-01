% load & plot distributions of cells as point-cloud-bar-charts --
% supplementary figure
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
tstep = 1/60;                   % time step in hours

actualLeaderFraction = NaN(length(followFracValues),numRepeats);
distanceMigrated = cell(size(followFracValues,1),numRepeats);
totalCellNumbers = NaN(size(followFracValues,1),numRepeats);

caCmap = load('cmap_blue2cyan.txt');

scatterFig = figure('Visible','on');
hold all
markerSize = 4; % area in points squared

%% load and plot data for every run of this parameter combination
for followFracCtr = 1:length(followFracValues)
    followerFraction = followFracValues(followFracCtr);
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
         
        actualLeaderFraction(followFracCtr,repCtr) = size(leaders,2)/numberOfCells;
        distanceMigrated{followFracCtr,repCtr} = cells(1,:);
        
        scatter(leaders(1,:),log2(actualLeaderFraction(followFracCtr,repCtr))*ones(size(leaders(2,:))),markerSize,[251 101 4]/255,'fill')
        scatter(followers(1,:),log2(actualLeaderFraction(followFracCtr,repCtr))*ones(size(followers(2,:))),markerSize,[113 18 160]/255,'fill')
        scatter(losts(1,:),log2(actualLeaderFraction(followFracCtr,repCtr))*ones(size(losts(2,:))),markerSize,[0.5 0.5 0.5],'fill')
        
        % record the total number of data points for later use
        totalCellNumbers(followFracCtr,repCtr) = numberOfCells;
    end
end

%% plot the median distances migrated over all repeats
symbols = {'>'; 'v'; '^'; 'd'; 's'};
hsymbols = NaN(size(symbols)); % handle vector for symbols to show in legend
maxCellNumber = max(max(totalCellNumbers));
for followFracCtr = 1:length(followFracValues)
    for repCtr = 1:numRepeats
    plotColor = [1 1 1] - totalCellNumbers(followFracCtr,repCtr)/maxCellNumber;
    hsymbols(followFracCtr) = plot(median(distanceMigrated{followFracCtr,repCtr}),log2(actualLeaderFraction(followFracCtr,repCtr)),...
        symbols{followFracCtr},'MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor,'MarkerSize',5);
    end
end
xlabel('x/\mum')
ylabel('leader fraction f_L')
xlim([0, 950])
reverseOrderfL = 2.^-(5:-1:0)';
ylim(log2([min(reverseOrderfL), 1]))
set(gca,'YDir','reverse','YTick',log2(reverseOrderfL),'YTickLabel', num2str(reverseOrderfL,'%0.2f'))

hlegend = legend(flipud(hsymbols),num2str(18*(1 - fliplr(followFracValues)'),precision),'Location','SouthEast');
set(get(hlegend,'title'),'string','t_{LF}/h')

% add a colorbar for relative stream density
ch = colorbar;
title(ch,'\rho^{stream}_{relative}')
colormap(linspace(1,0)'*[1 1 1])
box on

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
filename = ['manuscripts/subpopulations/figures/resultsFigS1'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);
