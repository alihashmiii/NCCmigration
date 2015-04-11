% perform anova analysis for each gene
% issues to-do
% - manova1 & manovacluster might be of interest (?), but currently gives the error that The within-group sum of squares and cross products matrix is singular.
clear, close all
load('integrateAndSwitchGeneExpression.mat')

pThreshold = 0.1;
%% VEGF removal
TX = [T0; T02; T04; T08; T16; T30; T45; T60; T90];
timesRep = [0 0 0 2 2 2 4 4 4 8 8 16 16 16 30 30 30 45 45 45 60 60 60 90 90 90];
% remove genes with data missing from all timepoints
missingData = all(isnan(TX));
TX = TX(:,~missingData);
Tgenes = genes(~missingData);
numTgenes = length(Tgenes);
TpVals = NaN(numTgenes,1);
figure
for geneCtr = 1:numTgenes
    missingData = isnan(TX(:,geneCtr));
    [TpVals(geneCtr),~,Tstats(geneCtr)] = anova1(TX(~missingData,geneCtr),timesRep(~missingData),'off');
    subplot(8,11,geneCtr)
    boxplot(TX(~missingData,geneCtr),timesRep(~missingData),'color','r');
    htitle = title([Tgenes{geneCtr} ' p=' num2str(TpVals(geneCtr))]);
    if TpVals(geneCtr)<=pThreshold&&Tstats(geneCtr).df>0 % if the means are different, check which groups are different from the first
        multipleComparison = multcompare(Tstats(geneCtr),'Alpha',pThreshold,'Display','off',...
            'CType','lsd'); %Tukey's least significant difference procedure. This procedure is a simple t-test. It is reasonable if the preliminary test (say, the one-way ANOVA F statistic) shows a significant difference. If it is used unconditionally, it provides no protection against multiple comparisons.
        differentMeanGroups = multipleComparison(multipleComparison(:,1)==1 ...% select comparison with first group
            &multipleComparison(:,6)<=pThreshold,:);% select "significant" differences
        if ~isempty(differentMeanGroups)
            % highlight values on plot
            yrange = ylim;
            hold on
            plot(repmat(differentMeanGroups(:,2),1,2)',yrange,'k--')
        end
    else
        htitle.FontWeight = 'normal';
    end
end

%% VEGF re-addition
VX = [T90; V02; V04; V08; V16; V30; V45; V60; V90];
timesRep = [0 0 0 2 2 2 4 4 8 8 8 16 16 30 30 45 45 45 60 60 90 90 90];
% remove genes with data missing from all timepoints
missingData = all(isnan(VX));
VX = VX(:,~missingData);
Vgenes = genes(~missingData);
numVgenes = length(Vgenes);
VpVals = NaN(numVgenes,1);
figure
for geneCtr = 1:numVgenes
    missingData = isnan(VX(:,geneCtr));
    [VpVals(geneCtr),~,Vstats(geneCtr)] = anova1(VX(~missingData,geneCtr),timesRep(~missingData),'off');
    subplot(8,11,geneCtr)
    boxplot(VX(~missingData,geneCtr),timesRep(~missingData),'color','b');
    htitle = title([Vgenes{geneCtr} ' p=' num2str(VpVals(geneCtr))]);
    if VpVals(geneCtr)<=pThreshold&&Vstats(geneCtr).df>0 % if the means are different, check which groups are different from the first
        multipleComparison = multcompare(Vstats(geneCtr),'Alpha',pThreshold,'Display','off',...
            'CType','lsd'); %Tukey's least significant difference procedure. This procedure is a simple t-test. It is reasonable if the preliminary test (say, the one-way ANOVA F statistic) shows a significant difference. If it is used unconditionally, it provides no protection against multiple comparisons.
        differentMeanGroups = multipleComparison(multipleComparison(:,1)==1 ...% select comparison with first group
            &multipleComparison(:,6)<=pThreshold,:);% select "significant" differences
        if ~isempty(differentMeanGroups)
            % highlight values on plot
            yrange = ylim;
            hold on
            plot(repmat(differentMeanGroups(:,2),1,2)',yrange,'k--')
        end
    else
        htitle.FontWeight = 'normal';
    end
end