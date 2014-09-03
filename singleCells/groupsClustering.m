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

%% make precision trade-off curve - set number of clusters and inverse temperature for which to do clusterings
numClusts = 2:5;
inverseTs = [0.5, 1, 30, 100];
numGroups = size(distribution,3);

%% calculate symmetric KL between each of the groups' distributions
sKL = NaN(numGroups);

pseudocount =  eps; % for regularisation of 0 probabilities
for ii = 1:numGroups
    for jj = ii:numGroups
        sKLtemp = NaN(size(distribution,2),1);
        for gg = 1:size(distribution,2)
        sKLtemp(gg) = ... 
            distribution(:,gg,ii)'*log((distribution(:,gg,ii) + pseudocount)...
            ./(distribution(:,gg,jj) + pseudocount)) + ...
            distribution(:,gg,jj)'*log((distribution(:,gg,jj) + pseudocount)...
            ./(distribution(:,gg,ii) + pseudocount));
        end
        sKL(ii,jj) = mean(sKLtemp);
        sKL(jj,ii) = sKL(ii,jj);
    end
end

clear sKLtemp

%% calculate Jensen-Shannon divergence between each of the groups' distributions
JS = NaN(numGroups);

pseudocount =  eps; % for regularisation of 0 probabilities
for ii = 1:numGroups
    for jj = ii:numGroups
        JStemp = NaN(size(distribution,2),1);
        for gg = 1:size(distribution,2)
        JStemp(gg) = ... 
            distribution(:,gg,ii)'*log((distribution(:,gg,ii) + pseudocount)...
            ./(mean(distribution(:,gg,[ii jj]),3) + pseudocount)) + ...
            distribution(:,gg,jj)'*log((distribution(:,gg,jj) + pseudocount)...
            ./(mean(distribution(:,gg,[ii jj]),3) + pseudocount));
        end
        JS(ii,jj) = mean(JStemp);
        JS(jj,ii) = JS(ii,jj);
    end
end

clear JStemp

%% compare sKL and JS
figure
subplot(1,2,1)
imagesc(sKL), colorbar, title('symmetrised Kulback-Leibler divergence')

subplot(1,2,2)
imagesc(JS), colorbar, title('Jensen-Shannon divergence')

%% I guess we're now minimizing s while maximizing I, since we're using a dissimilarity score (?)
for numClustCtr = length(numClusts):-1:1
    for invTempCtr = length(inverseTs):-1:1
        Ckl(numClustCtr,invTempCtr) = IclustDist(sKL,1/inverseTs(invTempCtr),numClusts(numClustCtr),[]);
    end
end
save('Ckl','Ckl')
%%
for numClustCtr = length(numClusts):-1:1
    for invTempCtr = length(inverseTs):-1:1
        Cjs(numClustCtr,invTempCtr) = IclustDist(JS,1/inverseTs(invTempCtr),numClusts(numClustCtr),[]);
    end
end
save('Cjs','Cjs')
%% plot precision trade-off curve (notice it approaches a lower bound / upper bound for -s now)
figure
hold all
symbols = {'+:';'o:';'^:';'s:';'x:'};
for numClustCtr = 1:length(numClusts)
    plot([Cjs(numClustCtr,:).Ici_end], -[Cjs(numClustCtr,:).avgS_end],symbols{numClustCtr})
end
xlabel('I(C;i) (bits)')
ylabel('<s> (bits)')
legend(num2str(numClusts'),'Location','SouthEast')
%% explore hierarchy
bestParents = NaN(min(numClusts)+1,max(numClusts));
inclusion = NaN(min(numClusts)+1,max(numClusts));
for numClustCtr = 2:length(numClusts)
   numClusters = numClusts(numClustCtr);
   [~, children] = max(Cjs(numClustCtr,end).Qc_i,[],2);
   [~, parents] = max(Cjs(numClustCtr-1,end).Qc_i,[],2);
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

