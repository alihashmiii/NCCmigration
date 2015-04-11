% load data for gene expression, to analyse time scales of change
clear, close all
load integrateAndSwitchGeneExpression.mat
plotError = 1;
%% normalise data relative to baseline (t -120) for each gene
normExpression = meanExpression;...
%     ./repmat(meanExpression(1,:),size(meanExpression,1),1);
normError = errorInMean;...
%        ./repmat(meanExpression(1,:),size(meanExpression,1),1);

%% impute the error for 1-replicate means
% take the mean error for each gene and time, and average for each gt-combi
meanErrorEstimate = nanmean(normError,2)*nanmean(normError)/2;
nanIndcs = isnan(normError);
imputedError = normError;
imputedError(nanIndcs) = meanErrorEstimate(nanIndcs);

%% select genes that show significant change, consitently
% find the indices of the relevant genes in the list of all genes
onSigIdcs = false(size(genes));
offSigIdcs = false(size(genes));
preSigIdcs = false(size(genes));

for gene = onSignificants'
    onSigIdcs = onSigIdcs|strcmpi(genes,gene);
end
for gene = offSignificants'
    offSigIdcs = offSigIdcs|strcmpi(genes,gene);
end
for gene = preSignificants'
    preSigIdcs = preSigIdcs|strcmpi(genes,gene);
end
allSigIdcs = onSigIdcs&offSigIdcs&preSigIdcs;

%% method 1: principle component analysis
% algorithms eig or svd (or als), Centered false or true
[offCoeffs, offPCs, offEigC, ~, offVarExplained] = pca(normExpression(2:10,allSigIdcs),...
    'Algorithm','eig','Centered',true,'variableWeights',1./mean(imputedError(2:10,allSigIdcs).^2));
%     'VariableWeights',1./nanmean(normError(2:10,allSigIdcs)./normExpression(2:10,allSigIdcs)),...
%     'Weights',1./nanmean(normError(2:10,allSigIdcs)./normExpression(2:10,allSigIdcs),2));
    
[onCoeffs, onPCs, onEigC, ~, onVarExplained] = pca(normExpression(10:end,allSigIdcs),...
    'Algorithm','eig','Centered',true,'variableWeights',1./mean(imputedError(10:end,allSigIdcs).^2));
%     'VariableWeights',1./nanmean(normError(10:end,allSigIdcs)./normExpression(10:end,allSigIdcs)),...
%     'Weights',1./nanmean(normError(10:end,allSigIdcs)./normExpression(10:end,allSigIdcs),2));

% plot bar chart of coeffs with gene names?
%% boot-strap - build a distribution of eigenvalues from shuffling the matrix
% numIter = 1e4;
% offEigBS = NaN(numIter,length(offEigC));
% offExpression = normExpression(2:10,allSigIdcs);
% shuffledOffExpression = NaN(size(offExpression));
% onEigBS = NaN(numIter,length(onEigC));
% onExpression = normExpression(10:end,allSigIdcs);
% shuffledOnExpression = NaN(size(onExpression));
% for ii = 1:numIter
%     %     for jj = 1:size(offExpression,2) % shuffle matrix entries
%     %        shuffledOffExpression(:,jj) = offExpression(randperm(size(offExpression,1)),jj);
%     %        shuffledOnExpression(:,jj) = onExpression(randperm(size(onExpression,1)),jj);
%     %     end
%     shuffledOffExpression = reshape(offExpression(randperm(numel(offExpression))),size(offExpression));
%     shuffledOnExpression = reshape(onExpression(randperm(numel(onExpression))),size(onExpression));
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


%% smooth data with cubic splines or akima interpolation with zero slope at final timepoint
nSmooth = 1000;
mSmooth = 'akima';
[offSmoothTime, offSmoothPCs] = smoothGeneExpression(timeData(2:10),offPCs,nSmooth,mSmooth,0);
[onSmoothTime, onSmoothPCs] = smoothGeneExpression(timeData(10:end),onPCs,nSmooth,mSmooth,0);

[offSmoothTime, offSmooth] = smoothGeneExpression(timeData(2:10),normExpression(2:10,:),nSmooth,mSmooth,0);
[onSmoothTime, onSmooth] = smoothGeneExpression(timeData(10:end),normExpression(10:end,:),nSmooth,mSmooth,0);

%% plot data and smoothed curves
%% significant genes
sigFig = figure;
subplot(2,1,1)
gh = plot(timeData(2:10),normExpression(2:10,allSigIdcs),'Color',[0.5 0.5 0.5]);
hold on
[weightedMean,~,weightedMeanAverageError] = ...
    weightedStats(normExpression(2:10,allSigIdcs)', sqrt(imputedError(2:10,allSigIdcs)'),'s');
% when providing sigma, the function weights by 1/sigma^2, but if we
% want single power weighting, provide srqt(sigma)
if plotError
    mh = errorbar(timeData(2:10),weightedMean,...
        weightedMeanAverageError,'m','LineWidth',3);
else
    mh = plot(timeData(2:10),weightedMean,'m','LineWidth',3);
end

legend([gh(1),mh],'genes','mean')
set(gca,'ytick',0:1:7,'xtick',0:30:90)
ylim([0 3])
xlim([0 90])
xlabel('time (min)')
ylabel('relative expression')
subplot(2,1,2)
gh = plot(timeData(10:end),normExpression(10:end,allSigIdcs),'Color',[0.5 0.5 0.5]);
hold on
[weightedMean,~,weightedMeanAverageError] = ...
    weightedStats(normExpression(10:end,allSigIdcs)', sqrt(imputedError(10:end,allSigIdcs)'),'s');
% when providing sigma, the function weights by 1/sigma^2, but if we
% want single power weighting, provide srqt(sigma)
if plotError
    mh = errorbar(timeData(10:end),weightedMean,...
        weightedMeanAverageError,'g','LineWidth',3);
else
    mh = plot(timeData(10:end),weightedMean,'g','LineWidth',3);
end
legend([gh(1),mh],'genes','mean')
set(gca,'ytick',0:1:7,'xtick',90:30:180)
xlim([90 180])
ylim([0 3])
xlabel('time (min)')
ylabel('relative expression')

%% smoothed genes
smoothFig = figure;
subplot(2,1,1)
plot(offSmoothTime,offSmooth(allSigIdcs,:),'m')
set(gca,'ytick',0:1:7,'xtick',0:30:90)
ylim([0 4])
xlabel('time (min)')
ylabel('relative expression')
subplot(2,1,2)
plot(onSmoothTime,onSmooth(allSigIdcs,:),'g')
set(gca,'ytick',0:1:7,'xtick',90:30:180)
ylim([0 4])
xlabel('time (min)')
ylabel('relative expression')


%% PCs
PCfig = figure;
nPlots = 1;
plotColor = lines(nPlots);
hLines = NaN(nPlots,2);
subplot(2,1,1), hold on
subplot(2,1,2), hold on
for plotCtr = 1:nPlots
    subplot(2,1,1)
    hLines(plotCtr,1) = plot(offSmoothTime,offSmoothPCs(plotCtr,:),'LineWidth',offVarExplained(plotCtr)/10,...
        'Color',plotColor(plotCtr,:));
    subplot(2,1,2)
    hLines(plotCtr,2) = plot(onSmoothTime,onSmoothPCs(plotCtr,:),'LineWidth',onVarExplained(plotCtr)/10,...
        'Color',plotColor(plotCtr,:));
%     plot(timeData(2:end),PCs(:,plotCtr),'+','LineWidth',varExplained(plotCtr)/20,...
%         'MarkerSize',sqrt(varExplained(plotCtr))*2,'Color',plotColor(plotCtr,:))
end
subplot(2,1,1)
xlim([offSmoothTime(1) offSmoothTime(end)])
% ylim([-4 5]), ylim([-2 4]), box on
set(gca,'ytick',-4:1:5,'xtick',0:30:90)
xlabel('time (min)')
ylabel('relative expression')
hLegend = legend(hLines(:,1),num2str(offVarExplained(1:nPlots),2),'Location','NorthEast');
%set(get(hLegend,'title'),'string',' % var explained')

subplot(2,1,2)
xlim([onSmoothTime(1) onSmoothTime(end)])
% ylim([-4 5]),  ylim([-2 4]), box on
set(gca,'ytick',-4:1:5,'xtick',90:30:180)
xlabel('time (min)')
ylabel('relative expression')
hLegend = legend(hLines(:,2),num2str(onVarExplained(1:nPlots),2),'Location','NorthEast');
%set(get(hLegend,'title'),'string',' % var explained')

%% plot histogram of first response times
% the genes in first response times are intersect(offSig,onSig), excluding
% Hprt, which is a housekeeping gene
timeFig = figure;
binCntrs = [2 4 8 16 30 45 90];
M = hist3(firstResponseTimes,{binCntrs,binCntrs});
imagesc(M')
set(gca,'xticklabel',num2str(binCntrs'),'yticklabel',num2str(binCntrs'),'TickLength',[0 0])
set(gca,'ydir','normal')
cb = colorbar;
cm = colormap(jet(max(max(M))+2));
colormap(cm(2:end,:)) %shift colormap for better readibility of text on blue...
caxis([-0.5 5.5])
set(cb,'ytick',0:5,'TickLength',[0 0])
set(get(cb,'xlabel'),'string','number of genes')
xlabel('first response after VEGF removal (min)')
ylabel('first response after VEGF readdition (min)')

% add gene names (taken from spreadsheet 'I&S comparison 090814'
% highlighted ones are allSigIdcs
highlight = struct('fontweight','bold');
text(1.55,7.35,'Adam10',highlight), text(1.55,7.125,'Cdh11'),
text(1.55,6.875,'Itga9',highlight),text(2.1,6.875,',b5'), text(1.55,6.65,'Pax3')
text(1.55,4.3,'Bmpr1a'), text(1.55,4,'Tgfbr1',highlight)
text(6.55,2.3,'Bmpr1b')
text(5.55,2.3,'Cdh7',highlight)
text(1.55,1.35,'Col2a1',highlight), text(1.55,1.125,'Itgb1',highlight),
text(1.55,0.875,'Nrp2',highlight), text(1.55,0.65,'Robo1',highlight)
text(6.55,3.3,'Cxcl12')
text(1.55,2.3,'Erbb2',highlight), text(1.55,2,'Itga4',highlight)
text(3.55,5.3,'Foxd3',highlight)
text(1.55,6.3,'Itgav'), text(1.55,6,'Tfap2a')
text(1.55,5.3,'Notch1',highlight), text(1.55,5,'Slit1')
text(4.55,7.3,'Pcdh1',highlight)
text(6.55,4.3,'Robo2')
text(6.55,1.3,'Sox10')
text(3.55,6.3,'Slit2')

%% export figures as eps & convert to PDF

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
figs = [sigFig, smoothFig, PCfig, timeFig];
filenames = {'significantGenes','significantGenesSmoothed','PCA','responseTimes'};
for figCtr = 1:length(figs)
    pos = get(figs(figCtr),'Position');
    if figCtr ~= 4
    pos(4) = 1/2*pos(3);% adjust height to fraction of width
    else
       exportOptions.Width = 15;
    end
    set(figs(figCtr),'PaperUnits','centimeters','Position',pos);
    exportfig(figs(figCtr),[filenames{figCtr} '.eps'],exportOptions);
    system(['epstopdf ' filenames{figCtr} '.eps']);
end