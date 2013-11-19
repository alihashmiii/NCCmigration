% load
% L.J. Schumacher 28.10.13

close all
clear

time = 18;
numRepeats = 100;
maxRuns2plot = 20;
numParamCombinations = 2;
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
eatRateValues = 100;
conversionType = 4;
num_stepsValues = 4*4;
num_stepsColors = jet(length(num_stepsValues));
standStill = 0;
volumeExclusion = 1;
tstep = 1/4*5/60;

for insert_step = [8 16]
    for followerFraction = [0, 1]
        profilesFig = figure('Visible','off');
        profiles2getherFig = figure('Visible','off');
        for num_stepsCtr = 1:length(num_stepsValues)
            num_steps = num_stepsValues(num_stepsCtr);
            for eatRateCtr = 1:length(eatRateValues)
                eatRate = eatRateValues(eatRateCtr);
                runsFig = figure('Visible','off');
                for repCtr = 1:numRepeats
                    loadInfo = ['experiment34/exp34_followFrac_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                        '_conversion_' num2str(conversionType) '_numSteps_' num2str(num_steps) ...
                        '_insertSteps_' num2str(insert_step) '_tstep_' num2str(tstep,precision) '_Run_' num2str(repCtr)];
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
                            out.ca_save{end},out.filopodia_save{end},out.numFilopodia,out.attach_save{end},out.cellRadius,0,1)
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
                title(['Exp3.4: leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', eatRate=' num2str(eatRate) ', followDefault=' num2str(followerFraction) ', num_steps=' num2str(num_steps) ', insert_step=' num2str(insert_step)])
                % save plot
                filename = ['results/experiment34/figures/exp34_defaultFollow_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                    '_numSteps_' num2str(num_steps) ...
                    '_tstep_' num2str(tstep,precision) '_allRuns.eps'];
                pos = get(runsFig,'Position');
                pos(4) = 3/2*pos(3);% adjust height to 3/2 width
                set(runsFig,'PaperUnits','centimeters','Position',pos);
                exportfig(runsFig,filename,exportOptions);
                system(['epstopdf ' filename]);
                close(runsFig);
                
                set(0,'CurrentFigure',profilesFig);
                subplot(length(num_stepsValues),length(eatRateValues),length(eatRateValues)*(num_stepsCtr - 1) + eatRateCtr)
                plot_migration_profile
                % xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
                % %                     legend([H3;H1;H2],'leaders','followers','chemoattractant');
                
                % title has parameter values and actual leader fraction
                if eatRateCtr==1&&num_stepsCtr==1
                    title(['Exp3.4: followDefault=' num2str(followerFraction) ', leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', insert_step=' num2str(insert_step)])
                elseif num_stepsCtr==1
                    title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', eatRate=' num2str(eatRate)])
                elseif num_stepsCtr==2&&eatRateCtr==1
                    title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', eatRate=' num2str(eatRate)])
                else
                    title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', steps=' num2str(num_steps)])
                end
                
                set(0,'CurrentFigure',profiles2getherFig);
                subplot(length(eatRateValues),1,eatRateCtr)
                if num_stepsCtr==1, hold on, end
                plot(xBins,squeeze(mean(sum(cellDistributions(paramCtr,:,:,:),3),2)),'Color',num_stepsColors(num_stepsCtr,:));
                if num_stepsCtr==length(num_stepsValues)
                    title(['Exp3.4: eatRate=' num2str(eatRate) ', followDefault=' num2str(followerFraction) ', tstep=' num2str(tstep,precision) ', insert_step=' num2str(insert_step)])
                    xlabel('x/\mum'), ylabel('N(cells)'), legend(num2str(num_stepsValues'))
                    ylim([0 10]), xlim([0 735]), set(gca,'YTick',[0 2 4 6 8 10]), grid on
                end
                
                eatRates(paramCtr) = eatRate;
                volumeExclusions(paramCtr) = volumeExclusion;
                standStills(paramCtr) = standStill;
                tsteps(paramCtr) = tstep;
                
                paramCtr = paramCtr + 1;
            end
        end
        % for each default lead/follow (followFrac = 0,1) make a plot with migration profiles for all set leader fractions and
        % eatRates - in total 12 plots with 21 subplots each
        pos = get(profilesFig,'Position');
        pos(4) = 3/2*pos(3);% adjust height to 3/2 width
        set(profilesFig,'PaperUnits','centimeters','Position',pos);
        filename = ['results/experiment34/figures/exp34_defaultFollow_' num2str(followerFraction) ...
            '_tstep_' num2str(tstep,precision) '_migrationProfiles.eps'];
        exportfig(profilesFig,filename,exportOptions);
        system(['epstopdf ' filename]);
        close(profilesFig);
        
        pos = get(profiles2getherFig,'Position');
        pos(4) = 3/2*pos(3);% adjust height to 3/2 width
        set(profiles2getherFig,'PaperUnits','centimeters','Position',pos);
        filename = ['results/experiment34/figures/exp34_defaultFollow_' num2str(followerFraction) ...
            '_tstep_' num2str(tstep,precision) '_migrationProfiles2gether.eps'];
        exportfig(profiles2getherFig,filename,exportOptions);
        system(['epstopdf ' filename]);
        close(profiles2getherFig);
        
        paramCtr = paramCtr + 1;
    end
end

save('results/experiment34/figures/experiment34collatedResults','xBins','cellDistributions','xlat_save','caDistribution','actualLeaderFraction','eatRates','volumeExclusions','standStills','tsteps')