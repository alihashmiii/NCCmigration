% load & plot distributions of cells as point-cloud-bar-charts
% L.J. Schumacher 21.05.14

close all
clear

time = 18;
numRepeats = 20;

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;

followFracValues = [0, 3/4, 7/8, 15/16, 1];
sensingAccuracy = 0.1;
needNeighbours = 0;

actualLeaderFraction = NaN(length(followFracValues),numRepeats);

caCmap = load('cmap_blue2cyan.txt');

scatterFig = figure('Visible','on');
axis equal
hold all
markerSize = 1;

for followFracCtr = 1:length(followFracValues)
    followerFraction = followFracValues(followFracCtr);
    yOffset = (followFracCtr - 1)*120;
    %% load and plot data for every run of this parameter combination
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
        scatter(leaders(1,:),leaders(2,:) + yOffset,markerSize,[1 0.5 0],'fill')
        scatter(followers(1,:),followers(2,:) + yOffset,markerSize,[1 0 1],'fill')
        scatter(losts(1,:),losts(2,:) + yOffset,markerSize,[0.5 0.5 0.5],'fill')
        
        actualLeaderFraction(followFracCtr,repCtr) = size(leaders,2)/numberOfCells;
%         
%         scatter(leaders(1,:),actualLeaderFraction(followFracCtr,repCtr)*ones(size(leaders(2,:))),markerSize,[1 0.5 0],'fill')
%         scatter(followers(1,:),actualLeaderFraction(followFracCtr,repCtr)*ones(size(followers(2,:))),markerSize,[1 0 1],'fill')
%         scatter(losts(1,:),actualLeaderFraction(followFracCtr,repCtr)*ones(size(losts(2,:))),markerSize,[0.5 0.5 0.5],'fill')
    end
    
end
xlabel('x/\mum')
ylabel('f_L')
set(gca,'YTick',linspace(60,120*length(followFracValues)-60,length(followFracValues)),'YTickLabel', ...
    num2str(mean(actualLeaderFraction,2),precision))

