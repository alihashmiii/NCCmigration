function experiment33increasedSizeWithConversion4
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 18;
input.conversionType = 4;
numReps = 20;

precision = 2; % significant figures for filenames and plot labels etc.
input.volumeExclusion = 1;
input.standStill = 0;
input.tstep = 1/4*5/60;
eatRate = 1000;

for insertEverySteps = [1 2]
    input.insertEverySteps = insertEverySteps;
    for followerFraction = [0, 1] % determines which is the default behaviour of cells, before switching
        input.followerFraction = followerFraction; 
            for diffus = [0.1 100]
                input.diffus = diffus;
                for numSteps = 4*4
                    input.numSteps = numSteps;
                    for repCtr = 1:numReps
                        input.saveInfo = ['experiment33/exp33_followFrac_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                            '_diff_' num2str(diffus) '_conversion_' num2str(input.conversionType) '_numSteps_' num2str(numSteps) ...
                            '_insertSteps_' num2str(input.insertEverySteps) '_tstep_' num2str(input.tstep,precision) '_Run_' num2str(repCtr)];
                        if isempty(dir(['results/' input.saveInfo '_running.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
                            rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                            CA6(input,0);
                            %                         quit
                        end
                    end
                end
            end
    end
end