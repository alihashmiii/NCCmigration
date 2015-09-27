% perform anova analysis for each gene
% issues to-do
% - manova1 & manovacluster might be of interest (?), but currently gives 
% the error that The within-group sum of squares and cross products matrix 
% is singular (even when excluding the reference genes). I think this might
% be due to the fact that we only have three replicates for each group, so
% that within-group sum-of-squares becomes singular for more than three
% genes?

clear, close all
load('integrateAndSwitchGeneExpression.mat')

pThreshold = 0.1;
precision = 2;
%% VEGF removal
TX = [T0; T02; T04; T08; T16; T30; T45; T60; T90];
timesRep = [0 0 0 2 2 2 4 4 4 8 8 16 16 16 30 30 30 45 45 45 60 60 60 90 90 90];
% remove genes with data missing from all timepoints
missingGenes = all(isnan(TX));
TX = TX(:,~missingGenes);
Tgenes = genes(~missingGenes);
numTgenes = length(Tgenes);
TpVals = NaN(numTgenes,1);
offFigure = figure;
for geneCtr = 1:numTgenes
    missingTimes = isnan(TX(:,geneCtr));
    [TpVals(geneCtr),~,Tstats(geneCtr)] = anova1(TX(~missingTimes,geneCtr),timesRep(~missingTimes),'off');
    subplot(11,8,geneCtr)
    boxplot(TX(~missingTimes,geneCtr),timesRep(~missingTimes),'color','r');
    ylim([0 ceil(max(TX(~missingTimes,geneCtr)))]);
    htitle = title([Tgenes{geneCtr} ' p=' num2str(TpVals(geneCtr),precision)]);
    if TpVals(geneCtr)<=pThreshold&&Tstats(geneCtr).df>0 % if the means are different, check which groups are different from the first
        multipleComparison = multcompare(Tstats(geneCtr),'Alpha',pThreshold,'Display','off',...
            'CType','lsd'); %Tukey's least significant difference procedure. This procedure is a simple t-test. It is reasonable if the preliminary test (say, the one-way ANOVA F statistic) shows a significant difference. If it is used unconditionally, it provides no protection against multiple comparisons.
        differentMeanGroups = multipleComparison(multipleComparison(:,1)==1 ...% select comparison with first group
            &multipleComparison(:,6)<=pThreshold,:);% select "significant" differences
        if ~isempty(differentMeanGroups)
            % highlight values on plot
            yrange = ylim;
            hold on
            plot(repmat(differentMeanGroups(:,2),1,2)',yrange,'k:')
        end
    else
        htitle.FontWeight = 'normal';
    end
end

%% VEGF re-addition
VX = [T90; V02; V04; V08; V16; V30; V45; V60; V90];
timesRep = [0 0 0 2 2 2 4 4 8 8 8 16 16 30 30 45 45 45 60 60 90 90 90];
% remove genes with data missing from all timepoints
missingGenes = all(isnan(VX));
VX = VX(:,~missingGenes);
Vgenes = genes(~missingGenes);
numVgenes = length(Vgenes);
VpVals = NaN(numVgenes,1);
onFigure = figure;
for geneCtr = 1:numVgenes
    missingTimes = isnan(VX(:,geneCtr));
    [VpVals(geneCtr),~,Vstats(geneCtr)] = anova1(VX(~missingTimes,geneCtr),timesRep(~missingTimes),'off');
    subplot(11,8,geneCtr)
    boxplot(VX(~missingTimes,geneCtr),timesRep(~missingTimes),'color','b');
    ylim([0 ceil(max(VX(~missingTimes,geneCtr)))]);
    htitle = title([Vgenes{geneCtr} ' p=' num2str(VpVals(geneCtr),precision)]);
    if VpVals(geneCtr)<=pThreshold&&Vstats(geneCtr).df>0 % if the means are different, check which groups are different from the first
        multipleComparison = multcompare(Vstats(geneCtr),'Alpha',pThreshold,'Display','off',...
            'CType','lsd'); %Tukey's least significant difference procedure. This procedure is a simple t-test. It is reasonable if the preliminary test (say, the one-way ANOVA F statistic) shows a significant difference. If it is used unconditionally, it provides no protection against multiple comparisons.
        differentMeanGroups = multipleComparison(multipleComparison(:,1)==1 ...% select comparison with first group
            &multipleComparison(:,6)<=pThreshold,:);% select "significant" differences
        if ~isempty(differentMeanGroups)
            % highlight values on plot
            yrange = ylim;
            hold on
            plot(repmat(differentMeanGroups(:,2),1,2)',yrange,'k:')
        end
    else
        htitle.FontWeight = 'normal';
    end
end

%% export figures
exportOptions = struct('Format','eps2',...
    'Width','19.0',...
    'Height','22',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',4,...
    'LineWidth',1);

filename = ['../../../Thesis/figures/anovaForEachGeneOffCondition'];
set(offFigure,'PaperUnits','centimeters');
exportfig(offFigure,[filename '.eps'],exportOptions);
%system(['epstopdf ' filename '.eps']);
%system(['pdfcrop ' filename '.pdf'])

filename = ['../../../Thesis/figures/anovaForEachGeneOnCondition'];
set(onFigure,'PaperUnits','centimeters');
exportfig(onFigure,[filename '.eps'],exportOptions);
%system(['epstopdf ' filename '.eps']);
%system(['pdfcrop ' filename '.pdf'])