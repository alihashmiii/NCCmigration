% load data for gene expression, to analyse time scales of change
clear
load integrateAndSwitchGeneExpression.mat

%% normalise data relative to baseline (t -120)
normExpression = meanExpression...
    ./repmat(meanExpression(1,:),size(meanExpression,1),1);
normError = errorInMean...
       ./repmat(meanExpression(1,:),size(meanExpression,1),1);

%% select genes that show significant change from t0/90
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
[offCoeffs, offPCs, ~, ~, offVarExplained] = pca(normExpression(2:10,allSigIdcs),...
    'Algorithm','eig','Centered',true);
[onCoeffs, onPCs, ~, ~, onVarExplained] = pca(normExpression(10:end,allSigIdcs),...
    'Algorithm','eig','Centered',true);

% plot bar chart of coeffs with gene names?
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
plot(timeData(2:10),normExpression(2:10,allSigIdcs),'m')
set(gca,'ytick',0:1:7,'xtick',0:30:90)
ylim([0 4])
xlabel('time (min)')
ylabel('relative expression')
subplot(2,1,2)
plot(timeData(10:end),normExpression(10:end,allSigIdcs),'g')
set(gca,'ytick',0:1:7,'xtick',90:30:180)
ylim([0 4])
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
nPlots = 3;
plotColor = lines(nPlots);
hLines = NaN(nPlots,2);
subplot(2,1,1), hold on
subplot(2,1,2), hold on
for plotCtr = nPlots:-1:1
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
ylim([-4 5]), ylim([-2 4]), box on
set(gca,'ytick',-4:1:5,'xtick',0:30:90)
xlabel('time (min)')
ylabel('relative expression')
hLegend = legend(hLines(:,1),num2str(offVarExplained,2),'Location','EastOutside');
set(get(hLegend,'title'),'string',' % var explained')

subplot(2,1,2)
xlim([onSmoothTime(1) onSmoothTime(end)])
ylim([-4 5]),  ylim([-2 4]), box on
set(gca,'ytick',-4:1:5,'xtick',90:30:180)
xlabel('time (min)')
ylabel('relative expression')
hLegend = legend(hLines(:,2),num2str(onVarExplained,2),'Location','EastOutside');
set(get(hLegend,'title'),'string',' % var explained')

%% plot histogram of first response times
% the genes in first response times are intersect(offSig,onSig), excluding
% Hprt, which is a housekeeping gene
timeFig = figure;
binCntrs = [2 4 8 16 30 45 90];
M = hist3(firstResponseTimes,{binCntrs,binCntrs});
imagesc(M')
set(gca,'xticklabel',num2str(binCntrs'),'yticklabel',num2str(binCntrs'))
set(gca,'ydir','normal')
cb = colorbar;
colormap(jet(max(max(M))+1))
caxis([-0.5 5.5])
set(cb,'ytick',0:5)
set(get(cb,'xlabel'),'string','n')
xlabel('first response after VEGF removal (min)')
ylabel('first response after VEGF readdition (min)')

% add gene names (taken from spreadsheet 'I&S comparison 090814'
text(1.6,7.35,'Adam10'), text(1.6,7.125,'Cdh11'), text(1.6,6.875,'Itga9,-b5'), text(1.6,6.65,'Pax3')
text(1.6,4.3,'Bmpr1a'), text(1.6,4,'Tgfbr1')
text(6.6,2.3,'Bmpr1b')
text(5.6,2.3,'Cdh7')
text(1.6,1.35,'Col2a1'), text(1.6,1.125,'Itgb1'), text(1.6,0.875,'Nrp2'), text(1.6,0.65,'Robo1')
text(6.6,3.3,'Cxcl12')
text(1.6,2.3,'Erbb2'), text(1.6,2,'Itga4')
text(3.6,5.3,'Foxd3')
text(1.6,6.3,'Itgav'), text(1.6,6,'Tfap2a')
text(1.6,5.3,'Notch1'), text(1.6,5,'Slit1')
text(4.6,7.3,'Pcdh1')
text(6.6,4.3,'Robo2')
text(6.6,1.3,'Sox10')
text(3.6,6.3,'Slit2')

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
       exportOptions.Width = 9;
    end
    set(figs(figCtr),'PaperUnits','centimeters','Position',pos);
    exportfig(figs(figCtr),[filenames{figCtr} '.eps'],exportOptions);
    system(['epstopdf ' filenames{figCtr} '.eps']);
end