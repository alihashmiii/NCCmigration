% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width
% has a zone of reduced cell speeds at the start (grows with tissue)
% DAN zone near domain entrance slows down cells, increases then decreases
% with time linearly

clear
close all

makeMoves = 1;
makeFrames = 1;
makeAllMovie = 1;
keepFrames = 0;

experiment = 42;
precision = 2;

diffusivities = [1];
slowSpeeds = [30 10];
insertEveryStepsValues = [10]; % corresponding to 4 and 2 (attempted) cell insertions per hour at tstep = 1min

fileName = 'exp42_slowEntryDynamicDAN';
numReps = 10;
for repCtr = 1:2
    for cntGdn = {'parallel'}
        result.contactGuidance = char(cntGdn);
        for diffus = diffusivities
            result.diffus = diffus;
            result.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
            for insertEverySteps = insertEveryStepsValues
                result.insertEverySteps = insertEverySteps;
                for slowSpeed = slowSpeeds
                    result.slowSpeed = slowSpeed;
                    result.loadInfo = [fileName '_D_' num2str(result.diffus) ...
                        '_sensingAcc_' num2str(result.sensingAccuracy,precision) ...
                        '_slowSpeed_' num2str(result.slowSpeed) ...
                        '_insertEvry_' num2str(result.insertEverySteps)...
                        '_contactGuidance_' result.contactGuidance '_Run_' num2str(repCtr)];
                    load(['results/' fileName '/' result.loadInfo '.mat'])
                    saveInfo = out.saveInfo(27:end); % had the folder name repeated in the saveInfo
                    make_frames(saveInfo)
                    keepFrames=1;
                    make_all_movie_hidden(saveInfo,[],keepFrames)
                end
            end
        end
    end
end