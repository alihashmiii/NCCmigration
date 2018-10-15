% load data for gene expression, to analyse time scales of change
clear
load integrateAndSwitchGeneExpression.mat

%% normalise data relative to baseline (t -120)
normExpression = meanExpression;...
%     ./repmat(meanExpression(1,:),size(meanExpression,1),1);

%% ignore genes with missing data
% meanExpression is in times by genes, each values average of 2-3
% replicates
% remove all columns with NaN values
nanIdcs = ~any(isnan(normExpression));
cleanExpression = normExpression(:,nanIdcs);
cleanGenes = genes(nanIdcs);

% % remove only columns with all NaN values
% cleanExpression = meanExpression(:,~all(isnan(meanExpression)));

%% split data into T and V branches
offExpression = cleanExpression(2:10,:);
onExpression = cleanExpression(10:end,:);

%% interweave gradients for better clustering?
% offExpression(1:2:15,:) = cleanExpression(2:10,:);
% onExpression(1:2:15,:) = cleanExpression(10:end,:);
% offExpression(2:2:14,:) = diff(cleanExpression(2:10,:))./repmat(diff(timeData(2:10)),1,66);
% onExpression(2:2:14,:) = diff(cleanExpression(10:end,:))./repmat(diff(timeData(2:10)),1,66);

%% calculate pearson correlation matrix btw each of the genes
% i.e. correlation of each genes's profile over all timepoints
Roff = corr(offExpression);
Ron = corr(onExpression);
% force symmetry (otherwise accuracy to 1e-15 might not be good
% enough for Iclust, add offset to make all values positive (doesn't seem
% to affect results, but surpresses warning)
Roff = (Roff + Roff')/2 + 1;   
Ron = (Ron + Ron')/2 + 1;        
%% make precision trade-off curve - set number of clusters and inverse temperature for which to do clusterings
numClusts = 2:4;
inverseTs = [3, 5, 10, 15, 20, 25, 30];
%% run Iclust for various numbers of clusters and temperatures
for numClustCtr = length(numClusts):-1:1
    for invTempCtr = length(inverseTs):-1:1
        Coff(numClustCtr,invTempCtr) = Iclust(Roff,1/inverseTs(invTempCtr),numClusts(numClustCtr),[]);
    end
end
for numClustCtr = length(numClusts):-1:1
    for invTempCtr = length(inverseTs):-1:1
        Con(numClustCtr,invTempCtr) = Iclust(Ron,1/inverseTs(invTempCtr),numClusts(numClustCtr),[]);
    end
end
save('Coff','Coff')
save('Con','Con')
%% plot precision trade-off curve to show that information is saturating
figure
symbols = {'+:';'o:';'^:';'s:';'x:'};
subplot(1,2,1), hold all
for numClustCtr = 1:length(numClusts)
    plot([Coff(numClustCtr,:).Ici_end], [Coff(numClustCtr,:).avgS_end],symbols{numClustCtr})
end
title('noVEGF')
xlabel('I(C;i) (bits)')
ylabel('<s> (bits)')
legend(num2str(numClusts'),'Location','SouthEast')

subplot(1,2,2), hold all
for numClustCtr = 1:length(numClusts)
    plot([Con(numClustCtr,:).Ici_end], [Con(numClustCtr,:).avgS_end],symbols{numClustCtr})
end
title('VEGF')
xlabel('I(C;i) (bits)')
ylabel('<s> (bits)')
legend(num2str(numClusts'),'Location','SouthEast')

%% plot curvers in each cluster
clusterFig = figure;
numClusters = 3;
colors = lines(2*numClusters);
plotOffset = 3;
subplot(1,2,1), hold all
[~, clustIDs] = max(Coff(numClusts==numClusters,end).Qc_i,[],2);
handles = cell(numClusters,1);
for clustCtr = 1:numClusters
    handles{clustCtr} = plot(timeData(2:10),offExpression(1:9,clustIDs==clustCtr) + plotOffset*...
        (clustCtr - 1),'Color',colors(clustCtr,:));
end
CIcounts = hist(clustIDs,1:numClusters);
% lh = legend(flipud(cellfun(@min,handles)),flipud(num2str(CIcounts')));
% set(get(lh,'title'),'string','# genes')
xlabel('t (min)')
xlim(timeData([2, 10])')
ylim([0 11])
title('-VEGF')

subplot(1,2,2), hold all
[~, clustIDs] = max(Con(numClusts==numClusters,end).Qc_i,[],2);
for clustCtr = 1:numClusters
    handles{clustCtr} = plot(timeData(10:end),onExpression(1:9,clustIDs==clustCtr) + plotOffset*...
        (clustCtr - 1),'Color',colors(numClusters + clustCtr,:));
end
CIcounts = hist(clustIDs,1:numClusters);
% lh = legend(flipud(cellfun(@min,handles)),flipud(num2str(CIcounts')));
% set(get(lh,'title'),'string','# genes')
xlabel('t (min)')
xlim(timeData([10, end])')
ylim([0 11])
title('+VEGF')

%% export figures as eps & convert to PDF

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

fileName = 'clusteredTrajectories';

pos = get(clusterFig,'Position');
pos(4) = 2/3*pos(3);% adjust height to 2/3 width
set(clusterFig,'PaperUnits','centimeters','Position',pos);
exportfig(clusterFig,[fileName '.eps'],exportOptions);
system(['epstopdf ' fileName '.eps']);
