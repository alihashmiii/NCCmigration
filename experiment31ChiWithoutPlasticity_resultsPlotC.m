% load & plot distributions of cells as overlapping bar-charts
% L.J. Schumacher 21.05.14

close all
clear
%%
time = 18;
numRepeats = 10;

precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;

followFracValues = [0, 1];
sensingAccuracy = 0.1;

% to calculate the density profile of cells along the x-direction
xBins = 0:50:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(length(followFracValues),numRepeats,3,length(xBins));

% preallocate variables for saving collated results
actualLeaderFraction = NaN(length(followFracValues),numRepeats);
alternativeColors = lines(length(followFracValues));

figure
hold all
for followFracCtr = length(followFracValues):-1:1
    followerFraction = followFracValues(followFracCtr);
    for repCtr = 1:numRepeats
        loadInfo = ['experiment31Chi/exp31_followFrac_' num2str(followerFraction,precision) ...
            '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
        try % sometimes we get corrupt files, which crashes the script
            load(['results/' loadInfo '.mat'])
        catch
            delete(['results/' loadInfo '.mat']) % delete the corrupt file
            experiment31ChiWithoutPlasticity; % recreate the missing results file
            load(['results/' loadInfo '.mat']) % load again
        end
        
        % load cell positions into variables
        cells = out.cells_save{end}; % all cells
        numberOfCells = size(cells,2);
        followIdcs = out.cellsFollow_save{end}(1:numberOfCells);
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
%     plotColor = f_L*[251 101 4]/255 + (1 - f_L)*[113 18 160]/255;
    plotColor = alternativeColors(followFracCtr,:);
    [xs, ys] = stairs(xBins,squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2)));
    %   plot as area-graph instead of bar as matlab will sometimes get the alpha wrong for barplots (on mac) 
    h(followFracCtr) = area(xs,ys,...
            'FaceColor', plotColor, 'EdgeColor', plotColor);
%     % alternatively plot regions of sem
%     h(followFracCtr) = area([xBins+25, fliplr(xBins+25)],...
%     [squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2))+std(squeeze(sum(cellDistributions(followFracCtr,:,:,:),3)))'/sqrt(numRepeats);...
%             flipud(squeeze(mean(sum(cellDistributions(followFracCtr,:,:,:),3),2))-std(squeeze(sum(cellDistributions(followFracCtr,:,:,:),3)))'/sqrt(numRepeats))],...
%     'FaceColor', plotColor, 'EdgeColor', plotColor);
    alpha(get(h(followFracCtr),'children'),0.5)
end
    
xlabel('x/\mum')
ylabel('# cells / 50\mum'), 
hlegend = legend(h,num2str(mean(actualLeaderFraction,2),precision),'Location',...
    'SouthWest');
set(get(hlegend,'title'),'string','<f_L>')
ylim([0 16]), xlim([0 800]), set(gca,'YTick',[0 4 8 12 16])
grid on, set(gca,'Layer','top')
box on
%% save figure as .fig file
filename = 'manuscripts/subpopulations/figures/resultsFigSChiC';
saveas(gcf,[filename '.fig'])

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
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);