% load data for gene expression, to analyse time scales of change
load integrateAndSwitchGeneExpression.mat

% meanExpression is in times by genes, each values average of 2-3
% replicates
% remove all columns with NaN values
cleanedExpression = meanExpression(:,~any(isnan(meanExpression)));

% % remove only columns with all NaN values
% cleanedExpression = meanExpression(:,~all(isnan(meanExpression)));

% principle component analysis
% algorithms eig or svd (or als), Centered false or true
[coeff, score, ~, ~, explained] = pca(cleanedExpression,...
    'Algorithm','eig','Centered',true');

figure
subplot(2,1,1)
plot(score(:,1))
hold all
plot(score(:,2))
set(gca,'XTick',1:size(times,1),'XTickLabel',times)
legend(num2str(explained(1:2)))

subplot(2,1,2)
plot(timeData,score(:,1))
hold all
plot(timeData,score(:,2))
xlabel('time (min)')

