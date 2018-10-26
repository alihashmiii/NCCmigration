function experiment31contStatesDiffus

input.time = 18;
numReps = 40;
input.saveEvery = 10;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')
sensingAccuracyUnscaled = 0.1;
guidanceModes = {'choice','combination'};

for guidanceMode = guidanceModes
    input.guidanceMode = guidanceMode{1};
    for diffus = [0.1, 1, 10, 100, 1e3, 1e4, 10^(4.5), 1e5]
        input.diffus = diffus;
        sensingAccuracy = sensingAccuracyUnscaled*sqrt(0.1/diffus); % sensing accuracy scales with diffusivity
        input.sensingAccuracy = sensingAccuracy;
        for repCtr = 1:numReps
            input.saveInfo = ['experiment31contStates_diffus/exp31' ...
                '_contStates_' input.guidanceMode '_D_' num2str(diffus,precision) ...
                '_sensingAcc_' num2str(sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            if isempty(dir(['../results/' input.saveInfo '_running*.mat']))...
                    &&isempty(dir(['../results/' input.saveInfo '.mat']))
                rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                CA8(input,0);
            end
        end
    end
end