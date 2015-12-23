% normality test for the multiple replicates of gene expressiond ata
% LJS 16.12.2015
close all
clear
load integrateAndSwitchGeneExpression.mat

%% 
sampleMatrices = {B1, B2, P20, P10, T0, T02, T04, T08, T16, T30, T45, T60,...
    T90, V02, V04, V08, V16, V30, V45, V60, V90};
numSamples = length(sampleMatrices);
numGenes = max(cellfun('size',sampleMatrices,2));
numReplicates = cellfun('size',sampleMatrices,1);
pValues = NaN(numSamples,numGenes);

for sampleCtr = 1:numSamples
    for geneCtr = 1:numGenes
       samples = sampleMatrices{sampleCtr}(:,geneCtr);
       if any(~isnan(samples))&&nanstd(samples)~=0
           samples = (samples - nanmean(samples))/nanstd(samples);
       [~, pValues(sampleCtr,geneCtr)] = kstest(samples);
       end
    end
end
%% 
histogram(pValues)
normals = nnz(pValues>0.9);
all = nnz(pValues<=0.9) + normals;
normals/all