% plot small multiples of migration profiles for switching time lead2follow
% vs follow2lead
% L.J. Schumacher 05.09.14

close all
clear all

time = 18;
numRepeats = 20;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
lead2followValues = [2 4 8 16];
follow2leadValues = [2 4 8 16];
numParamCombinations = length(lead2followValues)*length(follow2leadValues); 

xBins = 0:50:800; % bins for counting cell num vs. x profiles
cellDistributions = NaN(numParamCombinations,numRepeats,3,length(xBins));
referenceCellDistribution = NaN(numRepeats,3,length(xBins));

% preallocate variables for saving collated results
actualLeaderFraction = NaN(numParamCombinations,numRepeats);
referenceLeaderFraction = NaN(numRepeats);

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

paramCtr = 1;

for sensingAccuracy = [0.1, 0.01]
    figure
    hold on
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

        % calculate migration profile
        referenceCellDistribution(repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
        referenceCellDistribution(repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
        referenceCellDistribution(repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached
    end
    %% load data
    for lead2followCtr = 1:length(lead2followValues)
        lead2follow = lead2followValues(lead2followCtr);
        for follow2leadCtr = 1:length(follow2leadValues)
            follow2lead = follow2leadValues(follow2leadCtr);
            numSteps = [lead2follow, follow2lead];
            for repCtr = 1:numRepeats
                loadInfo = ['experiment31conversion4/exp31'...
                    '_conversion_' num2str(conversionType) '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_0'...
                    '_Run_' num2str(repCtr)];
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
                
                actualLeaderFraction(paramCtr,repCtr) = size(leaders,2)/numberOfCells;
                
                % calculate migration profile
                numberOfCells = size(out.cells_save{end},2);
                cellDistributions(paramCtr,repCtr,1,:) = histc(out.cells_save{end}(1,out.cellsFollow_save{end}(1:numberOfCells)==0),xBins); % leaders
                cellDistributions(paramCtr,repCtr,2,:) = histc(out.cells_save{end}(1,(out.cellsFollow_save{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)~=0)),xBins); % followers, attached
                cellDistributions(paramCtr,repCtr,3,:) = histc(out.cells_save{end}(1,(out.cellsFollow_save{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)==0)),xBins); % followers, dettached
                
                % calculate migration profile
                numberOfCells = size(out.cells_save{end},2);
                cellDistributions(paramCtr,repCtr,1,:) = histc(leaders(1,:),xBins); % leaders
                cellDistributions(paramCtr,repCtr,2,:) = histc(followers(1,:),xBins); % followers, attached
                cellDistributions(paramCtr,repCtr,3,:) = histc(losts(1,:),xBins); % followers, attached
            end
            %% plot migration profile
            f_L = mean(actualLeaderFraction(paramCtr,:));
            plotColor = f_L*[251 101 4]/255 + (1 - f_L)*[113 18 160]/255;
            % plot diagonal 'gridline' to aid the eye
            plot(max(xBins)*(lead2followCtr - [1 0]), 16*(follow2leadCtr - [0 1]),'--','Color',[0.5 0.5 0.5])
            % plot reference migration profile
            stairs(xBins + max(xBins)*(lead2followCtr - 1) ... % add x-offset
                ,squeeze(mean(sum(referenceCellDistribution,2),1)) + 16*(follow2leadCtr - 1),... % add y-offset
                'color',[0.5 0.5 0.5],'LineWidth',1);
            % plot offset migration profile
            stairs(xBins + max(xBins)*(lead2followCtr - 1) ... % add x-offset
                ,squeeze(mean(sum(cellDistributions(paramCtr,:,:,:),3),2)) + 16*(follow2leadCtr - 1),... % add y-offset
                'color',plotColor,'LineWidth',2);
            % add label with mean leader fraction
            text(max(xBins)*(lead2followCtr - 2/3), 16*(follow2leadCtr - 1/4),['<f_L> = ' num2str(f_L,precision)])
            paramCtr = paramCtr + 1;
        end
    end
    grid on
    set(gca,'xtick',max(xBins)*(1:length(lead2followValues)),'xticklabel',num2str(lead2followValues'))
    set(gca,'ytick',16*(1:length(follow2leadValues)),'yticklabel',num2str(follow2leadValues'))
    set(gca,'GridLineStyle','-')
    xlabel('lead to follow (min)')
    ylabel('follow to lead (min)')
end
