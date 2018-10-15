% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width

clear
close all

makeMoves = 1;
makeFrames = 1;
makeAllMovie = 1;

experiment = 39;
precision = 2;

diffusivities = 1;%[1 10 100];
speeds = [40];
insertEveryStepsValues = [10]; % corresponding to 4 and 2 (attempted) cell insertions per hour at tstep = 1min

for cntGdn = {'parallel'}
    result.contactGuidance = char(cntGdn);
    for diffus = diffusivities
        result.diffus = diffus;
        result.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        for insertEverySteps = insertEveryStepsValues
            result.insertEverySteps = insertEverySteps;
            for speed = speeds
                result.leadSpeed = speed;
                result.followSpeed = speed;
                result.loadInfo = ['exp39_widerDomainStripe' ...
                    '_D_' num2str(result.diffus) '_sensingAcc_' num2str(result.sensingAccuracy,precision) ...
                        '_speed_' num2str(result.leadSpeed) '_insertEvry_' num2str(insertEverySteps) ...
                    '_contactGuidance_' result.contactGuidance '_Run_1'];
                load(['results/exp39_widerDomainStripe/' result.loadInfo '.mat'])
                saveInfo = out.saveInfo(25:end); % had the folder name repeated in the saveInfo
                make_frames(saveInfo)
                keepFrames =1;
                make_all_movie_hidden(saveInfo,[],keepFrames)
            end
        end
    end
end