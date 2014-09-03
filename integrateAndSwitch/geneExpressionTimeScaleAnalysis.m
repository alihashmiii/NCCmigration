% load data for gene expression, to analyse time scales of change
clear
load integrateAndSwitchGeneExpression.mat

%% normalise data relative to baseline (t -120)
normExpression = meanExpression...
    ./repmat(meanExpression(1,:),size(meanExpression,1),1);

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

%% method 1: principle component analysis
% algorithms eig or svd (or als), Centered false or true
[offCoeffs, offPCs, ~, ~, offVarExplained] = pca(offExpression,...
    'Algorithm','eig','Centered',true);
[onCoeffs, onPCs, ~, ~, onVarExplained] = pca(onExpression,...
    'Algorithm','eig','Centered',true);

% plot bar chart of coeffs with gene names?

%% method 2: select only those genes that are part of the 16 gene profile
geneProfile = {'aqp1 ';'bambi ';'cdh7 ';'cdh11 ';'cxcr4 ';'ephb3 ';'itgb5 ';'Nedd9 ';'Notch1 ';'Pkp2 ';'tfap2a '};
% this ignores ccr9,ctnnb1,cxcr1&7,snai2 which weren't looked at. tfap2a is
% added in, but for hand2 there are too many missing data

% find the indices of the relevant genes in the list of all genes
geneIdcs = NaN(size(geneProfile));
for geneCtr = 1:length(geneProfile)
    geneIdcs(geneCtr) = find(strcmpi(cleanGenes,geneProfile(geneCtr)));
end
offProfileGenes = offExpression(:,geneIdcs);
onProfileGenes = onExpression(:,geneIdcs);

% compare with plot of all genes?

%% smooth data with cubic splines or akima interpolation with zero slope at final timepoint
nSmooth = 1000;
mSmooth = 'akima';
[offSmoothTime, offSmoothPCs] = smoothGeneExpression(timeData(2:10),offPCs,nSmooth,mSmooth,0);
[onSmoothTime, onSmoothPCs] = smoothGeneExpression(timeData(10:end),onPCs,nSmooth,mSmooth,0);

[offSmoothTime, offSmoothPGs] = smoothGeneExpression(timeData(2:10),offProfileGenes,nSmooth,mSmooth,0);
[onSmoothTime, onSmoothPGs] = smoothGeneExpression(timeData(10:end),onProfileGenes,nSmooth,mSmooth,0);

[offSmoothTime, offSmooth] = smoothGeneExpression(timeData(2:10),offExpression,nSmooth,mSmooth,0);
[onSmoothTime, onSmooth] = smoothGeneExpression(timeData(10:end),onExpression,nSmooth,mSmooth,0);

%% plot data and smoothed curves
%% raw genes
rawFig = figure;
subplot(2,1,1)
plot(timeData(2:10),offExpression)
xlabel('time (min)')
ylabel('relative expression')
title('all genes -VEGF condition')

subplot(2,1,2)
plot(timeData(10:end),onExpression)
xlabel('time (min)')
ylabel('relative expression')
title('all genes +VEGF condition')

%% smoothed genes
smoothFig = figure;
subplot(2,1,1)
plot(offSmoothTime,offSmooth)
hold on
plot(timeData(2:10),offExpression,'+')
xlabel('time (min)')
ylabel('relative expression')
title('all genes -VEGF condition')

subplot(2,1,2)
plot(onSmoothTime,onSmooth)
hold on
plot(timeData(10:end),onExpression,'+')
xlabel('time (min)')
ylabel('relative expression')
title('smoothed genes +VEGF condition')

%% Profile genes
PGfig = figure;
symbols = {'+','o','*','.','x','s','d','^','v','<','>','p','h'};
subplot(2,1,1)
nPlots = size(offProfileGenes,2);
plotColor = lines(nPlots);
hLines = NaN(nPlots,1);
hold on
for plotCtr = nPlots:-1:1
    plot(offSmoothTime,offSmoothPGs(plotCtr,:),'Color',plotColor(plotCtr,:));
    hLines(plotCtr) = plot(timeData(2:10),offProfileGenes(:,plotCtr),...
        symbols{plotCtr},'Color',plotColor(plotCtr,:));
end
xlim([offSmoothTime(1) offSmoothTime(end)])
xlabel('time (min)')
ylabel('relative expression')
title('Trailblazer genes -VEGF condition')
hLegend = legend(hLines,cleanGenes(geneIdcs),'Location','EastOutside');
set(get(hLegend,'title'),'string','gene IDs')

subplot(2,1,2)
nPlots = size(onProfileGenes,2);
plotColor = lines(nPlots);
hLines = NaN(nPlots,1);
hold on
for plotCtr = nPlots:-1:1
    plot(onSmoothTime,onSmoothPGs(plotCtr,:),'Color',plotColor(plotCtr,:));
    hLines(plotCtr) = plot(timeData(10:end),onProfileGenes(:,plotCtr),...
        symbols{plotCtr},'Color',plotColor(plotCtr,:));
end
xlim([onSmoothTime(1) onSmoothTime(end)])
xlabel('time (min)')
ylabel('relative expression')
title('Trailblazer genes +VEGF condition')
hLegend = legend(hLines,cleanGenes(geneIdcs),'Location','EastOutside');
set(get(hLegend,'title'),'string','gene IDs')

%% PCs
PCfig = figure;
subplot(2,1,1)
nPlots = size(offPCs,2);
plotColor = lines(nPlots);
hLines = NaN(nPlots,1);
hold on
for plotCtr = nPlots:-1:1
    hLines(plotCtr) = plot(offSmoothTime,offSmoothPCs(plotCtr,:),'LineWidth',offVarExplained(plotCtr)/10,...
        'Color',plotColor(plotCtr,:));
    plot(timeData(2:10),offPCs(:,plotCtr),'+','LineWidth',offVarExplained(plotCtr)/20,...
        'MarkerSize',sqrt(offVarExplained(plotCtr))*2,'Color',plotColor(plotCtr,:))
end
xlim([offSmoothTime(1) offSmoothTime(end)])
xlabel('time (min)')
ylabel('relative expression')
title('PCA -VEGF condition')
hLegend = legend(hLines,num2str(offVarExplained,2),'Location','EastOutside');
set(get(hLegend,'title'),'string',' % \sigma explained')

subplot(2,1,2)
nPlots = size(onPCs,2);
plotColor = lines(nPlots);
hLines = NaN(nPlots,1);
hold on
for plotCtr = nPlots:-1:1
    hLines(plotCtr) = plot(onSmoothTime,onSmoothPCs(plotCtr,:),'LineWidth',onVarExplained(plotCtr)/10,...
        'Color',plotColor(plotCtr,:));
    plot(timeData(10:end),onPCs(:,plotCtr),'+','LineWidth',onVarExplained(plotCtr)/20,...
        'MarkerSize',sqrt(onVarExplained(plotCtr))*2,'Color',plotColor(plotCtr,:))
end
xlim([onSmoothTime(1) onSmoothTime(end)])
xlabel('time (min)')
ylabel('relative expression')
title('PCA +VEGF condition')
hLegend = legend(hLines,num2str(onVarExplained,2),'Location','EastOutside');
set(get(hLegend,'title'),'string',' % \sigma explained')

%% export figures as eps & convert to PDF

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
figs = [rawFig, smoothFig, PGfig, PCfig];
filenames = {'allGenes','smoothedGenes','profileGenes','PCA'};
for figCtr = 1:4
    pos = get(figs(figCtr),'Position');
    % pos(4) = 3/2*pos(3);% adjust height to 3/2 width
    set(figs(figCtr),'PaperUnits','centimeters','Position',pos);
    exportfig(figs(figCtr),[filenames{figCtr} '.eps'],exportOptions);
    system(['epstopdf ' filenames{figCtr} '.eps']);
end