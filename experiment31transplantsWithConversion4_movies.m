% plot migration profiles for simulated VEGF transplants.
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
lead2follow = [8];
follow2lead = [8];
sensingAccuracyValues = [0.1];
experiments = [11 12];
numParamCombinations = length(defaultFollowValues)*length(sensingAccuracyValues)...
    *length(experiments);


xBins = 0:50:800; % bins for counting cell num vs. x profiles

for defaultFollow = defaultFollowValues
    for sensAccCtr = 1:length(sensingAccuracyValues)
        sensingAccuracy = sensingAccuracyValues(sensAccCtr);
        for expCtr = 1:length(experiments)
            experiment = experiments(expCtr);
            
            %% load data
            numSteps = [lead2follow, follow2lead];
            for repCtr = 2
                if makeAllMovie==1
                % load the relevant reference simulation
                loadInfo = ['experiment31conversion4/exp31'...
                    '_conversion_4_defaultFollow_' num2str(defaultFollow) ...
                    '_numSteps_' num2str(lead2follow) '_' num2str(follow2lead) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                load(['results/' loadInfo '.mat'])
                % load results from the output structure into
                % variables
                load_results
                saveInfo = ['exp' num2str(experiment) ...
                    '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                    '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                keepFrames=1;
                stopTime = 12;
                make_all_movie_hidden
                startAtFrame = frameCtr;
                end
                % load the transplant simulation
                loadInfo = ['experiment31transplants/exp' num2str(experiment) ...
                    '_conversion_' num2str(conversionType) '_defaultFollow_' num2str(defaultFollow) ...
                    '_numSteps_' num2str(numSteps(1)) '_' num2str(numSteps(2)) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                load(['results/' loadInfo '.mat'])
               
                % load results from the output structure into
                % variables
                load_results
                saveInfo = ['exp' num2str(experiment) ...
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
                    
                    filename = ['results/experiment31transplants/figures/Frames_' saveInfo];
                    set(framesFig,'PaperUnits','centimeters');
                    exportfig(framesFig,[filename '.eps'],exportOptions);
                    system(['epstopdf ' filename '.eps']);
                    system(['rm -f ' filename '.eps']);
                    close(framesFig)
                end
                %% make movies %%
                %%% make cells+ca movie (allmovie.avi)%%%
                if makeAllMovie==1
                    keepFrames = 1;
                    stopTime = 25;
                    make_all_movie_hidden
                    startAtFrame = 1;
                end
            end
        end
    end
end