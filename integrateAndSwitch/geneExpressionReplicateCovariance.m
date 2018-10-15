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

for sampleCtr = 1:numSamples
    fullExpression(sampleCtr,1:(numReplicates(sampleCtr)*numGenes)) = ...
        sampleMatrices{sampleCtr}(:);
end

%% cut out all-NaN rows and columns

fullExpression = fullExpression(:,~all(isnan(fullExpression)));
fullExpression = fullExpression(~all(isnan(fullExpression')),:);

%% calculate and display correlation (variance-normalised covariance)

corrIS = corrcoef(fullExpression,'rows','pairwise');

imagesc(corrIS), colorbar
colormap(autumn(4))
hold on
spy(isnan(corrIS),'k')