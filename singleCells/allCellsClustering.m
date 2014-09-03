% plot histograms of single cell gene expression data for each cell group

nBins = 16;

bins = linspace(0,15,nBins);
distribution = NaN(nBins,96,6);

loadExpressionData;

% pool data from all cells in each group into a histogram to create
% distribution of expression (for each gene and group)
distribution(:,:,1) = hist(trailblazers16h',bins)/size(trailblazers16h,2);
distribution(:,:,2) = hist(trailblazers24h',bins)/size(trailblazers24h,2);
distribution(:,:,3) = hist(quartile1',bins)/size(quartile1,2);
distribution(:,:,4) = hist(quartile2',bins)/size(quartile2,2);
distribution(:,:,5) = hist(quartile3',bins)/size(quartile3,2);
distribution(:,:,6) = hist(quartile4',bins)/size(quartile4,2);

% concatenate data to cluster all cells together
allCells = [trailblazers16h, trailblazers24h, quartile1, quartile2, ...
    quartile3, quartile4];
% keep track of original cell groupings for later comparison
cellIDs = [ones(1,size(trailblazers16h,2)), 2*ones(1,size(trailblazers24h,2)),...
    3*ones(1,size(quartile1,2)),4*ones(1,size(quartile2,2)),...
    5*ones(1,size(quartile3,2)),6*ones(1,size(quartile4,2))];

%% calculate MI matrix between each of the cells
% - doesn't make so much sense to pool expression from all genes into a
% distribution of (any) expression for each cell, so might not be useful
% nCells = size(allCells,2);
% MI = NaN(nCells);
% for ii = 1:nCells
%     for jj = ii:nCells
%         MI(ii,jj) = mutualInformation(round(allCells(:,ii)),round(allCells(:,jj)));
%         % force symmetry (otherwise accuracy to 1e-15 might not be good
%         % enough for Iclust
%         MI(jj,ii) = MI(ii,jj); 
%     end
% end
% save('MI','MI')
%% calculate pearson correlation matrix btw each of the cells
% i.e. correlation of each cell's profile over all genes
R = corr(allCells);
% force symmetry (otherwise accuracy to 1e-15 might not be good
% enough for Iclust
R = (R + R')/2;        
%% make precision trade-off curve - set number of clusters and inverse temperature for which to do clusterings
numClusts = 2:6;
inverseTs = [10, 15, 20, 25, 30, 35, 40];
%% run Iclust for various numbers of clusters and temperatures
% for numClustCtr = length(numClusts):-1:1
%     for invTempCtr = length(inverseTs):-1:1
%         Cmi(numClustCtr,invTempCtr) = Iclust(MI,1/inverseTs(invTempCtr),numClusts(numClustCtr),[]);
%     end
% end
% save('Cmi','Cmi')
%%
for numClustCtr = length(numClusts):-1:1
    for invTempCtr = length(inverseTs):-1:1
        Cr(numClustCtr,invTempCtr) = Iclust(R,1/inverseTs(invTempCtr),numClusts(numClustCtr),[]);
    end
end
save('Cr','Cr')
%% plot precision trade-off curve to show that information is saturating
figure
hold all
symbols = {'+:';'o:';'^:';'s:';'x:'};
for numClustCtr = 1:length(numClusts)
    plot([Cr(numClustCtr,:).Ici_end], [Cr(numClustCtr,:).avgS_end],symbols{numClustCtr})
end
xlabel('I(C;i) (bits)')
ylabel('<s> (bits)')
legend(num2str(numClusts'),'Location','SouthEast')
%% explore hierarchy
bestParents = NaN(min(numClusts)+1,max(numClusts));
inclusion = NaN(min(numClusts)+1,max(numClusts));
for numClustCtr = 2:length(numClusts)
   numClusters = numClusts(numClustCtr);
   [~, children] = max(Cr(numClustCtr,end).Qc_i,[],2);
   [~, parents] = max(Cr(numClustCtr-1,end).Qc_i,[],2);
   for clustCtr = 1:numClusters
       % find the cluster in current 'generation' that includes most cells of
       % the cluster in next 'generation'
       bestParents(numClustCtr-1,clustCtr) = mode(parents(children==clustCtr));
       % find how many of the cells in the child cluster are also in the best
       % parent cluster, and divide by total number in child cluster
       inclusion(numClustCtr-1,clustCtr) = sum(children==clustCtr&...
           parents==bestParents(numClustCtr-1,clustCtr))...
           ./sum(children==clustCtr);
   end
end

%% plot clusters and their proportional content of original groupings
pieLabels = {'t13','t15','Q1','Q2','Q3','Q4'};
pieColors = [[1 0 0]; [0 1 0]; [0 0 1]; [0.7 0.7 0]; [0.5 0 0.5]; [1 0 1]];

numGroups = size(distribution,3);

for numClustCtr = 1:length(numClusts)
    numClusters = numClusts(numClustCtr);
    [~, clustIDs] = max(Cr(numClustCtr,end).Qc_i,[],2);
    for clustCtr = 1:numClusters
        subplot(numClusts(end)-1,numClusts(end), numClusts(end)*(numClustCtr -1 )...
            + clustCtr)
        h = pie(hist(cellIDs(clustIDs==clustCtr),1:numGroups),pieLabels);
        % change color of slices to reflect original groups
        h = findobj(h,'Type','patch');
        groups = unique(cellIDs(clustIDs==clustCtr));
        for groupCtr = 1:length(groups)
            set(h(groupCtr),'FaceColor',pieColors(groups(groupCtr),:),...
                'EdgeColor','none')
        end
        title(['n = ' num2str(sum(clustIDs==clustCtr))])
    end
end