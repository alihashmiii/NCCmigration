% load
% L.J. Schumacher 28.10.13

close all
clear

time = 18;
numRepeats = 100;
maxRuns2plot = 20;
numParamCombinations = 8;
% to calculate the density profile of cells and chemoattractant along the x-direction
cellRadius = 7.5;
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
xBins = 0:(cellRadius + filolength):24*(cellRadius + filolength); % bins for counting cell num vs. x profiles
cellDistributions = NaN(numParamCombinations,numRepeats,3,length(xBins));
caDistribution = NaN(numParamCombinations,numRepeats,50);
xlat_save = NaN(50,1);
% preallocate variables for saving collated results
actualLeaderFraction = NaN(numParamCombinations,1);
eatRates = NaN(numParamCombinations,1);
volumeExclusions = NaN(numParamCombinations,1);
standStills = NaN(numParamCombinations,1);
tsteps = NaN(numParamCombinations,1);

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);
precision = 2; % significant figures for filenames and plot labels etc.
paramCtr = 1;
eatRate = 100;
conversionType = 4;
numSteps = 4*5;
standStill = 0;
volumeExclusion = 1;
tstep = 1/4*5/60;
insertStepsValues = [8 16];
followFracValues = [0, 1];
diffValues = [0.1 100];
profilesFig = figure('Visible','off');
for insertStepsCtr = 1:length(insertStepsValues)
    insertEverySteps = insertStepsValues(insertStepsCtr);
    for followFracCtr = 1:length(followFracValues)
        followerFraction = followFracValues(followFracCtr);
        for diffCtr = 1:length(diffValues)
            diffus = diffValues(diffCtr);
            runsFig = figure('Visible','off');
            for repCtr = 1:numRepeats
                loadInfo = ['experiment34/exp34_followFrac_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                    '_diff_' num2str(diffus) '_conversion_' num2str(conversionType) '_numSteps_' num2str(numSteps) ...
                    '_insertSteps_' num2str(insertEverySteps) '_tstep_' num2str(tstep,precision) '_Run_' num2str(repCtr)];
                try % sometimes we get corrupt files, which crashes the script
                    load(['results/' loadInfo '.mat'])
                catch
                    delete(['results/' loadInfo '.mat']) % delete the corrupt file
                    experiment34decreasedSizeWithConversion4; % recreate the missing results file
                    load(['results/' loadInfo '.mat']) % load again
                end
                % make a plot of all repeats
                if repCtr <= maxRuns2plot
                    subplot(min(numRepeats,maxRuns2plot)/2 + 2,2,repCtr+2)
                    make_plot(out.cells_save{end},out.cellsFollow{end},out.xlat_save{end},out.ylat_save{end}, ...
                        out.ca_save{end},out.filopodia_save{end},out.numFilopodia,out.attach_save{end},out.cellRadius,filolength,sensingAccuracy,0,1)
                    title([num2str(size(out.cells_save{end},2)) ' cells, ' num2str(min([size(out.cells_save{end},2) nnz(out.cellsFollow{end}==0)])) ' leaders.'])
                end
                % calculate migration profile
                numberOfCells = size(out.cells_save{end},2);
                cellDistributions(paramCtr,repCtr,1,:) = histc(out.cells_save{end}(1,out.cellsFollow{end}(1:numberOfCells)==0),xBins); % leaders
                cellDistributions(paramCtr,repCtr,2,:) = histc(out.cells_save{end}(1,(out.cellsFollow{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)~=0)),xBins); % followers, attached
                cellDistributions(paramCtr,repCtr,3,:) = histc(out.cells_save{end}(1,(out.cellsFollow{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)==0)),xBins); % followers, attached
                caDistribution(paramCtr,repCtr,:) = mean(out.ca_save{end},2);
                if paramCtr==1, xlat_save = out.xlat_save{end}; end % load the x-coordinated of the CA profile, only once as they're always the same
            end
            % plot migration profile
            subplot(min(numRepeats,maxRuns2plot)/2 + 2,2,[1 2])
            plot_migration_profile
            xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
            legend([H3(3);H3(2);H3(1)],'lead','follow','lost');
            
            % title has parameter values and actual leader fraction
            actualLeaderFraction(paramCtr) = sum(mean(squeeze(cellDistributions(paramCtr,:,1,:)))); % mean number of leader cells
            actualLeaderFraction(paramCtr) = actualLeaderFraction(paramCtr)/(actualLeaderFraction(paramCtr) + sum(sum(mean(squeeze(cellDistributions(paramCtr,:,2:3,:)))))); % divide by mean total number of cells
            title(['Exp3.4: leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', eatRate=' num2str(eatRate) '_D=' num2str(diffus) ', followDefault=' num2str(followerFraction) ', numSteps=' num2str(numSteps) ', insertSteps=' num2str(insertEverySteps)])
            % save plot
            filename = ['results/experiment34/figures/exp34_defaultFollow_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                '_diff_' num2str(diffus) '_numSteps_' num2str(numSteps) ...
                '_insertStep_' num2str(insertEverySteps) '_tstep_' num2str(tstep,precision) '_allRuns.eps'];
            pos = get(runsFig,'Position');
            pos(4) = 3/2*pos(3);% adjust height to 3/2 width
            set(runsFig,'PaperUnits','centimeters','Position',pos);
            exportfig(runsFig,filename,exportOptions);
            system(['epstopdf ' filename]);
            close(runsFig);
            
            set(0,'CurrentFigure',profilesFig);
            subplot(length(followFracValues)+length(insertStepsValues),length(diffValues),length(diffValues)*(followFracCtr - 1) + (length(diffValues) + length(followFracValues))*(insertStepsCtr - 1) + diffCtr)
            plot_migration_profile
            % xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
            % %                     legend([H3;H1;H2],'leaders','followers','chemoattractant');
            
             % title has parameter values and actual leader fraction
            if diffCtr==1,
                title(['Exp3.4: followDefault=' num2str(followerFraction) ', leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', D=' num2str(diffus) ', insertSteps=' num2str(insertEverySteps)])
            else
                title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', D=' num2str(diffus)])
            end
            
            eatRates(paramCtr) = eatRate;
            volumeExclusions(paramCtr) = volumeExclusion;
            standStills(paramCtr) = standStill;
            tsteps(paramCtr) = tstep;
            
            paramCtr = paramCtr + 1;
        end
    end
end
% make one summary plot with a migration profile for each of the 8
% parameter combinations
pos = get(profilesFig,'Position');
pos(4) = 3/2*pos(3);% adjust height to 3/2 width
set(profilesFig,'PaperUnits','centimeters','Position',pos);
filename = ['results/experiment34/figures/exp34_eatRate_' num2str(eatRate) ...
                '_numSteps_' num2str(numSteps) ...
                '_tstep_' num2str(tstep,precision) '_migrationProfiles.eps'];
exportfig(profilesFig,filename,exportOptions);
system(['epstopdf ' filename]);
close(profilesFig);

save('results/experiment34/figures/experiment34collatedResults','xBins','cellDistributions','xlat_save','caDistribution','actualLeaderFraction','eatRates','volumeExclusions','standStills','tsteps')