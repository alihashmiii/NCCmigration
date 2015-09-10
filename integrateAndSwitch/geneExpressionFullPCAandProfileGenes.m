% load data for gene expression, to analyse time scales of change
clear
load integrateAndSwitchGeneExpression.mat

%% build the times by genesXreplicates matrix
sampleMatrices = {B1, B2, P20, P10, T0, T02, T04, T08, T16, T30, T45, T60,...
    T90, V02, V04, V08, V16, V30, V45, V60, V90};
numSamples = length(sampleMatrices);
numGenes = max(cellfun('size',sampleMatrices,2));
numReplicates = cellfun('size',sampleMatrices,1);
fullExpression = NaN(numSamples, max(numReplicates)*numGenes);

% - we should repeat this multiple times after shuffling the replicate
% identity at every time point (but keep it the same for all genes in a sample),
% and build up a distribution of results (as the data are not truly longitudinal
% and the replicate assignment was arbitrary
% - could also try setting all NaN to zero, see what difference it makes?
for sampleCtr = 1:numSamples
    fullExpression(sampleCtr,1:(numReplicates(sampleCtr)*numGenes)) = ...
        sampleMatrices{sampleCtr}(:);
end

%% cut out all-NaN rows and columns

fullExpression = fullExpression(:,~all(isnan(fullExpression)));
fullExpression = fullExpression(~all(isnan(fullExpression')),:);

%% normalise data relative to overall baseline (t -120)
normExpression = fullExpression...
    ./repmat(nanmean(nanmean(fullExpression(1:2,:))),...
    size(fullExpression,1),size(fullExpression,2));

%% split data into T and V branches
offExpression = normExpression(3:13,:);
onExpression = normExpression(13:end,:);
% cut out all-NaN rows and columns
offExpression = offExpression(:,~all(isnan(offExpression)));
onExpression = onExpression(:,~all(isnan(onExpression)));
offExpression = offExpression(~all(isnan(offExpression')),:);
onExpression = onExpression(~all(isnan(onExpression')),:);

%% method 1: probabilistic principle component analysis
% use ppca or pca with algorithm 'als'
% do we need to center data and normalise the stdev?
opt = statset('ppca');
opt.MaxIter = 1e4;
opt.TolFun = 0.01;
[offCoeffs, offPCs, offEigC, ~, offResVar] = ppca(offExpression,3,'Options',opt);
[onCoeffs, onPCs, onEigC, ~, onResVar] = ppca(onExpression,3,'Options',opt);

% plot bar chart of coeffs with gene names?

%% boot-strap - build a distribution of eigenvalues from shuffling the matrix
% numIter = 1e4;
% offEigBS = NaN(numIter,length(offEigC));
% shuffledOffExpression = NaN(size(offExpression));
% onEigBS = NaN(numIter,length(onEigC));
% shuffledOnExpression = NaN(size(onExpression));
% for ii = 1:numIter
    %     for jj = 1:size(offExpression,2) % shuffle matrix entries
    %        shuffledOffExpression(:,jj) = offExpression(randperm(size(offExpression,1)),jj);
    %        shuffledOnExpression(:,jj) = onExpression(randperm(size(onExpression,1)),jj);
    %     end
%     
%     [~, ~, offEigBS(ii,:), ~, ~] = pca(shuffledOffExpression,...
%         'Algorithm','eig','Centered',true);
%     [~, ~, onEigBS(ii,:), ~, ~] = pca(shuffledOnExpression,...
%         'Algorithm','eig','Centered',true);
% end
% %%
% nBins = 200;
% [offEigDist, bins] = hist(offEigBS(:),linspace(0,max(offEigBS(:)),nBins));
% % % correct for bin spacing affecting distribution height...
% offEigDist = offEigDist./numIter./size(offExpression,2)./mean(diff(bins));
% % eigenFigure = figure;
% subplot(2,1,1)
% plot(bins,offEigDist)
% hold on
% % we don't have enough eigenvalues to compare distributions, instead just
% % indicate where on the bootstrapped distribution they lie
% for lambda = offEigC'
%     stem(lambda,offEigDist(find(bins>=lambda,1,'first')),'r')
% end
% xlabel('\lambda'), ylabel('\rho(\lambda)')
% legend('bootstrapped distribution','PC eigenvalues')
% % repeat for second experimental regime
% [onEigDist, bins] = hist(onEigBS(:),linspace(0,max(onEigBS(:)),nBins));
% onEigDist = onEigDist./numIter./size(onExpression,2)./mean(diff(bins));
% subplot(2,1,2)
% plot(bins,onEigDist)
% hold on
% for lambda = onEigC'
%     stem(lambda,onEigDist(find(bins>=lambda,1,'first')),'r')
% end
% xlabel('\lambda'), ylabel('\rho(\lambda)')
% legend('bootstrapped distribution','PC eigenvalues')

%% check distribution of PC components against unit normal distribution
% scale variance to unity
scaledOffCoeffs = offCoeffs./sqrt(repmat(sum(offCoeffs.^2),size(offCoeffs,1),1))*sqrt(size(offCoeffs,1));
scaledOnCoeffs = onCoeffs./sqrt(repmat(sum(onCoeffs.^2),size(onCoeffs,1),1))*sqrt(size(onCoeffs,1));
pOffCoeffs = NaN(size(offPCs,2),1); pOnCoeffs = pOffCoeffs;
% figure, hold on -- too few samples to plot the distribution of components?
for pcCtr = 1:size(offPCs,2)
[~, pOffCoeffs(pcCtr)] = kstest(scaledOffCoeffs(:,pcCtr));
[~, pOnCoeffs(pcCtr)] = kstest(scaledOnCoeffs(:,pcCtr));
end
