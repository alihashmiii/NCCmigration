% plot small multiples of migration profiles for switching time lead2follow
% vs follow2lead
% L.J. Schumacher 05.09.14

close all
clear all

makeFrames = 1;
makeAllMovie = 1;
caCmap = load('cmap_blue2cyan.txt');

time = 18;

precision = 2; % significant figures for filenames and plot labels etc.

conversionType = 4;
defaultFollowValues = [2];
switchingTimes = [4 8];
sensingAccuracyValues = [0.1 0.01];

for sensAccCtr = 1:length(sensingAccuracyValues)
    sensingAccuracy = sensingAccuracyValues(sensAccCtr);
    for defaultFollow = defaultFollowValues
        %% load data
        for switchTimeCtr = 1:length(switchingTimes)
            lead2follow = switchingTimes(switchTimeCtr);
            follow2lead = switchingTimes(switchTimeCtr);
            numSteps = [lead2follow, follow2lead];
            for repCtr = 1:2
                loadInfo = ['experiment31conversion4/exp31'...
                    '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                    '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                try % sometime we get corrupt files, which crashes the script
                    load(['results/' loadInfo '.mat'])
                catch
                    delete(['results/' loadInfo '.mat']) % delete the corrupt file
                    experiment31leaderFractionWithConversion4; % recreate the missing results file
                    load(['results/' loadInfo '.mat']) % load again
                end
                % load results from the output structure into
                % variables
                load_results
                saveInfo = ['exp31'...
                    '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                    '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                %% make frames %%
                if makeFrames==1
                    make_frames
                    disp('made frames')
                    framesFig = open(['avi_mat/frames/',saveInfo,'.fig']);
                    %% export figure
                    exportOptions = struct('Format','eps2',...
                        'Width','18.0',...
                        'Color','rgb',...
                        'Resolution',300,...
                        'FontMode','fixed',...
                        'FontSize',10,...
                        'LineWidth',2);
                    
                    filename = ['results/experiment31conversion4/figures/Frames_' saveInfo];
                    set(framesFig,'PaperUnits','centimeters');
                    exportfig(framesFig,[filename '.eps'],exportOptions);
                    system(['epstopdf ' filename '.eps']);
                    system(['rm -f ' filename '.eps']);
                    close(framesFig)
                end
                %% make movies %%
                %%% make cells+ca movie (allmovie.avi)%%%
                if makeAllMovie==1
                    make_all_movie_hidden(saveInfo)
                end
            end
        end
        
        
    end
end