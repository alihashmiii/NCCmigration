% plot histograms of single cell gene expression data for each cell group

nBins = 16;

bins = linspace(0,15,nBins);
distribution = NaN(nBins,96,6);

loadExpressionData;

distribution(:,:,1) = hist(trailblazers16h',bins)/size(trailblazers16h,2);
distribution(:,:,2) = hist(trailblazers24h',bins)/size(trailblazers24h,2);
distribution(:,:,3) = hist(quartile1',bins)/size(quartile1,2);
distribution(:,:,4) = hist(quartile2',bins)/size(quartile2,2);
distribution(:,:,5) = hist(quartile3',bins)/size(quartile3,2);
distribution(:,:,6) = hist(quartile4',bins)/size(quartile4,2);

%% plot the average expression for each group
plot(bins,squeeze(mean(distribution,2)))
legend('t16','t24','q1','q2','q3','q4')

%% plot all genes for each group

figure
hold on
nGroups = 6;
colors = lines(nGroups);
for gCtr = 1:nGroups
    subplot(2,3,gCtr)
    plot(bins,distribution(:,:,gCtr),'Color',colors(gCtr,:))
    xlim([1 15])
end

%% make violin plots (uses distributionPlot.m)

%% plot distributions of genes in trailblazer profile for trailblazers at 16h & 24h, and the rest of the stream

