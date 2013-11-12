% load
% L.J. Schumacher 28.10.13

close all
clear

time = 18;
numRepeats = 20;
numParamCombinations = 42;
% to calculate the density profile of cells and chemoattractant along the x-direction
cellRadius = 7.5;
xBins = 0:2*cellRadius:49*2*cellRadius; % check the end length for time step run, or try loading this from the out-file somehow
cellDistributions = NaN(numParamCombinations,numRepeats,3,length(xBins));
caDistribution = NaN(numParamCombinations,numRepeats,length(xBins));
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
eatRateValues = [100, 150, 300];
followFracValues = [0,(10:15)/16];
followFracColors = jet(length(followFracValues));
for volumeExclusion = 1
    for standStill = [1, 0]
        for tstep = 1/4*5/60
            profilesFig = figure('Visible','off');
            profiles2getherFig = figure('Visible','off');
            for followFracCtr = 1:length(followFracValues)
                followerFraction = followFracValues(followFracCtr);
                for eatRateCtr = 1:length(eatRateValues)
                    eatRate = eatRateValues(eatRateCtr);
                    runsFig = figure('Visible','off');
                    for repCtr = 1:numRepeats
                        loadInfo = ['experiment31/exp31_followFrac_' num2str(followerFraction, precision) '_eatRate_' num2str(eatRate) ...
                            '_volumeExclusion_' num2str(volumeExclusion) '_standStill_' num2str(standStill) ...
                            '_tstep_' num2str(tstep, precision) '_Run_' num2str(repCtr)];
                        load(['results/' loadInfo '.mat'])
                        % make a plot of all repeats
                        subplot(numRepeats/2 + 2,2,repCtr+2)
                        make_plot(out.cells_save{end},out.cellsFollow{end},out.xlat_save{end},out.ylat_save{end}, ...
                            out.ca_save{end},out.filopodia_save{end},out.numFilopodia,out.attach_save{end},out.cellRadius,0,1)
                        title([num2str(size(out.cells_save{end},2)) ' cells, ' num2str(min([size(out.cells_save{end},2) nnz(out.cellsFollow{end}==0)])) ' leaders.'])
                        % calculate migration profile
                        numberOfCells = size(out.cells_save{end},2);
                        cellDistributions(paramCtr,repCtr,1,:) = histc(out.cells_save{end}(1,out.cellsFollow{end}(1:numberOfCells)==0),xBins); % leaders
                        cellDistributions(paramCtr,repCtr,2,:) = histc(out.cells_save{end}(1,(out.cellsFollow{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)~=0)),xBins); % followers, attached
                        cellDistributions(paramCtr,repCtr,3,:) = histc(out.cells_save{end}(1,(out.cellsFollow{end}(1:numberOfCells)==1)&(out.attach_save{end}(1:numberOfCells)==0)),xBins); % followers, attached
                        caDistribution(paramCtr,repCtr,:) = sum(out.ca_save{end},2);
                    end
                    % plot migration profile
                    subplot(numRepeats/2 + 2,2,[1 2])
                    plot_migration_profile
                    xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
                    legend([H3(3);H3(2);H3(1)],'lead','follow','lost');
                    
                    % title has parameter values and actual leader fraction
                    actualLeaderFraction(paramCtr) = sum(mean(squeeze(cellDistributions(paramCtr,:,1,:)))); % mean number of leader cells
                    actualLeaderFraction(paramCtr) = actualLeaderFraction(paramCtr)/(actualLeaderFraction(paramCtr) + sum(sum(mean(squeeze(cellDistributions(paramCtr,:,2:3,:)))))); % divide by mean total number of cells
                    title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', eatRate=' num2str(eatRate) ', volExcl=' num2str(volumeExclusion) ', standStill=' num2str(standStill) ', tstep=' num2str(tstep,precision) ])
                    % save plot
                    filename = ['results/experiment31/figures/exp31_followFrac_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                        '_volumeExclusion_' num2str(volumeExclusion) '_standStill_' num2str(standStill) ...
                        '_tstep_' num2str(tstep,precision) '_allRuns.eps'];
                    pos = get(runsFig,'Position');
                    pos(4) = 3/2*pos(3);% adjust height to 3/2 width
                    set(runsFig,'PaperUnits','centimeters','Position',pos);
                    exportfig(runsFig,filename,exportOptions);
                    system(['epstopdf ' filename]);
                    close(runsFig);
                    
                    set(0,'CurrentFigure',profilesFig);
                    subplot(length(followFracValues),length(eatRateValues),length(eatRateValues)*(followFracCtr - 1) + eatRateCtr)
                    plot_migration_profile
                   % xlabel('x/\mum'), ylabel(AX(1),'N(cells)'), ylabel(AX(2),'C(chemoattractant)')
% %                     legend([H3;H1;H2],'leaders','followers','chemoattractant');

                    % title has parameter values and actual leader fraction
                    if eatRateCtr==1&&followFracCtr==1
                        title(['volExcl=' num2str(volumeExclusion) ', standStill=' num2str(standStill) ', tstep=' num2str(tstep,precision) ])
                    else
                        title(['leadFrac=' num2str(actualLeaderFraction(paramCtr),precision) ', eatRate=' num2str(eatRate)])
                    end
                    
                    set(0,'CurrentFigure',profiles2getherFig);
                    subplot(length(eatRateValues),1,eatRateCtr)
                    if followFracCtr==1, hold on, end
                    plot(xBins,squeeze(mean(sum(cellDistributions(paramCtr,:,:,:),3),2)),'Color',followFracColors(followFracCtr,:));
                    if followFracCtr==length(followFracValues)
                        title(['eatRate=' num2str(eatRate) ', volExcl=' num2str(volumeExclusion) ', standStill=' num2str(standStill) ', tstep=' num2str(tstep,precision) ])
                        xlabel('x/\mum'), ylabel('N(cells)'), legend(num2str(actualLeaderFraction((length(eatRateValues):length(eatRateValues):length(eatRateValues)*length(followFracValues)) - length(eatRateValues)*length(followFracValues) + paramCtr),precision))
                        ylim([0 7]), xlim([0 735]), set(gca,'YTick',[0 2 4 6])
                    end
                    
                    eatRates(paramCtr) = eatRate;
                    volumeExclusions(paramCtr) = volumeExclusion;
                    standStills(paramCtr) = standStill;
                    tsteps(paramCtr) = tstep;
                    
                    paramCtr = paramCtr + 1;
                end
            end
            % for each volumeExclusion, standStill and timeStep, make a plot with migration profiles for all set leader fractions and
            % eatRates - in total 12 plots with 21 subplots each
            pos = get(profilesFig,'Position');
            pos(4) = 3/2*pos(3);% adjust height to 3/2 width
            set(profilesFig,'PaperUnits','centimeters','Position',pos);
            filename = ['results/experiment31/figures/exp31_volExcl_' num2str(volumeExclusion) ...
                '_standStill_' num2str(standStill) '_tstep_' num2str(tstep,precision) '_migrationProfiles.eps'];
            exportfig(profilesFig,filename,exportOptions);
            system(['epstopdf ' filename]);
            close(profilesFig);
            
            pos = get(profiles2getherFig,'Position');
            pos(4) = 3/2*pos(3);% adjust height to 3/2 width
            set(profiles2getherFig,'PaperUnits','centimeters','Position',pos);
            filename = ['results/experiment31/figures/exp31_volExcl_' num2str(volumeExclusion) ...
                '_standStill_' num2str(standStill) '_tstep_' num2str(tstep,precision) '_migrationProfiles2gether.eps'];
            exportfig(profiles2getherFig,filename,exportOptions);
            system(['epstopdf ' filename]);
            close(profiles2getherFig);
            
            paramCtr = paramCtr + 1;
        end
        paramCtr = paramCtr + 1;
    end
    paramCtr = paramCtr + 1;
end

save('results/experiment31/figures/experiment31collatedResults','xBins','cellDistributions','caDistribution','actualLeaderFraction','eatRates','volumeExclusions','standStills','tsteps')