% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 18;
input.conversionType = 4;
numReps = 100;

precision = 2; % significant figures for filenames and plot labels etc.
input.volumeExclusion = 1;
input.standStill = 0;
input.tstep = 1/4*5/60;

for followerFraction = [0, 1] % determines which is the default behaviour of cells, before switching
    input.followerFraction = followerFraction;
    for eatRate = [100, 300, 900]
        input.eatRate = eatRate;
        for numSteps = 4*[1 2 3 4 5 6 7]
            input.numSteps = numSteps;
            for repCtr = 1:numReps
                input.saveInfo = ['experiment31conversion4/exp31_followFrac_' num2str(followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                    '_conversion_' num2str(input.conversionType) '_numSteps_' num2str(numSteps) ...
                    '_tstep_' num2str(input.tstep,precision) '_Run_' num2str(repCtr)];
                if isempty(dir(['results/' input.saveInfo '_running.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
                    rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                    CA6(input,0);
                    %quit
                end
            end
        end
    end
end