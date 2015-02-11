% plot small multiples of migration profiles for switching time lead2follow
% vs follow2lead
% L.J. Schumacher 05.09.14

close all
clear all

time = 30;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
defaultFollowValues = [0 1 2];
lead2followValues = [1 2 4 8 12 16 24 32 40 48 56];
follow2leadValues = [1 2 4 8 12 16 24 32 40 48 56];
sensingAccuracyValues = [0.1, 0.01];
numParamCombinations = length(defaultFollowValues)*length(sensingAccuracyValues)*length(lead2followValues)*length(follow2leadValues);

xBins = 0:50:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(numParamCombinations,numRepeats,3,length(xBins));
referenceCellDistribution = NaN(numRepeats,3,length(xBins));

% preallocate variables for saving collated results
actualLeaderFraction = NaN(length(defaultFollowValues),length(sensingAccuracyValues),length(lead2followValues),length(follow2leadValues),numRepeats);
numCells = NaN(length(defaultFollowValues),length(sensingAccuracyValues),length(lead2followValues),length(follow2leadValues),numRepeats);
referenceLeaderFraction = NaN(numRepeats,1);
referenceNumCells = NaN(numRepeats,1);

paramCtr = 1;

for sensAccCtr = 1:length(sensingAccuracyValues)
    sensingAccuracy = sensingAccuracyValues(sensAccCtr);
    %% load non-switching reference experiments
    for repCtr = 1:numRepeats
        loadInfo = ['experiment31/exp31_followFrac_1_sensingAcc_' num2str(sensingAccuracy)...
            '_needNeighbours_0_Run_' num2str(repCtr)];
        
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
        
        referenceLeaderFraction(repCtr) = size(leaders,2)/numberOfCells;
        referenceNumCells(repCtr) = numberOfCells;
        
        % calculate migration profile
        referenceCellDistribution(repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
        referenceCellDistribution(repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
        referenceCellDistribution(repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached
    end
    for defaultFollow = defaultFollowValues
        figure
        hold on
        %% load data
        for lead2followCtr = 1:length(lead2followValues)
            lead2follow = lead2followValues(lead2followCtr);
            for follow2leadCtr = 1:length(follow2leadValues)
                follow2lead = follow2leadValues(follow2leadCtr);
                numSteps = [lead2follow, follow2lead];
                for repCtr = 1:numRepeats
                    loadInfo = ['experiment31conversion4/exp31'...
                        '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                        '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                    try % sometime we get corrupt files, which crashes the script
                        load(['results/' loadInfo '.mat'])
                    catch
                        delete(['results/' loadInfo '.mat']) % delete the corrupt file
                        experiment31leaderFractionWithConversion4; % recreate the missing results file
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
                    
                    actualLeaderFraction(defaultFollow + 1,sensAccCtr,lead2followCtr,follow2leadCtr,repCtr) = size(leaders,2)/numberOfCells;
                    numCells(defaultFollow + 1,sensAccCtr,lead2followCtr,follow2leadCtr,repCtr) = numberOfCells;
                    
                    % calculate migration profile
                    cellDistributions(paramCtr,repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
                    cellDistributions(paramCtr,repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
                    cellDistributions(paramCtr,repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached
                end
                %% plot migration profile
                f_L = mean(actualLeaderFraction(defaultFollow + 1,sensAccCtr,lead2followCtr,follow2leadCtr,:));
                n_C = mean(numCells(defaultFollow + 1,sensAccCtr,lead2followCtr,follow2leadCtr,:));
                plotColor = f_L*[251 101 4]/255 + (1 - f_L)*[113 18 160]/255;
%                 % plot diagonal 'gridline' to aid the eye
%                 plot(max(xBins)*(lead2followCtr - [1 0]), 16*(follow2leadCtr - [0 1]),'--','Color',[0.5 0.5 0.5])
                % plot reference migration profile
                stairs(xBins + max(xBins)*(lead2followCtr - 1) ... % add x-offset
                    ,squeeze(mean(sum(referenceCellDistribution,2),1)) + 16*(follow2leadCtr - 1),... % add y-offset
                    'color',[0.5 0.5 0.5],'LineWidth',1);
                % plot offset migration profile
                stairs(xBins + max(xBins)*(lead2followCtr - 1) ... % add x-offset
                    ,squeeze(mean(sum(cellDistributions(paramCtr,:,:,:),3),2)) + 16*(follow2leadCtr - 1),... % add y-offset
                    'color',plotColor,'LineWidth',2);
                % add label with mean leader fraction and number of cells
                text(max(xBins)*(lead2followCtr - 9.5/12), 16*(follow2leadCtr - 2/12),...
                    ['$\bar{f_L}$ = ' num2str(f_L,precision)],'FontSize',4,'Interpreter','Latex')
                text(max(xBins)*(lead2followCtr - 7/12), 16*(follow2leadCtr - 4/12),...
                    ['$\bar{n}$ = ' num2str(n_C,precision)],'FontSize',4,'Interpreter','Latex')
                paramCtr = paramCtr + 1;
            end
        end
        grid on
        set(gca,'xtick',max(xBins)*(1:length(lead2followValues)),'xticklabel',num2str(lead2followValues'))
        set(gca,'ytick',16*(1:length(follow2leadValues)),'yticklabel',num2str(follow2leadValues'))
        set(gca,'GridLineStyle','-')
        xlabel('lead to follow (min)')
        ylabel('follow to lead (min)')
        %% export figure
        exportOptions = struct('Format','eps2',...
            'Width','26.0',...
            'Color','rgb',...
            'Resolution',300,...
            'LineWidth',2);
        %             'FontMode','fixed',...
%             'FontSize',8,...
        
        filename = ['manuscripts/VEGF/figures/FigS3_defaultFollow_' num2str(defaultFollow) '_sensAcc_' num2str(sensingAccuracy)];
        pos = get(gcf,'Position');
        % pos(4) = 1/2*pos(3); % adjust height to fraction of width
        set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
        exportfig(gcf,[filename '.eps'],exportOptions);
        system(['epstopdf ' filename '.eps']);
    end
end

%% save data from workspace
save('manuscripts/VEGF/figures/experiment31conv4collatedResults','xBins','cellDistributions','actualLeaderFraction','lead2followValues','follow2leadValues','sensingAccuracyValues','numCells','referenceCellDistribution','referenceLeaderFraction','referenceNumCells')
