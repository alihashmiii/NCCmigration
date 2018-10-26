function experiment31contStatesEat

input.time = 18;
numReps = 40;
input.saveEvery = 10;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')
input.sensingAccuracy = 0.1;
guidanceModes = {'choice','combination'};

for guidanceMode = guidanceModes
    input.guidanceMode = guidanceMode{1};
    for eatRate = [10, 25, 50, 75, 100, 125, 250, 500, 1000, 1500, 2000]
        input.eatRate = eatRate;
        for repCtr = 1:numReps
            input.saveInfo = ['experiment31contStates_eat/exp31' ...
                '_contStates_' input.guidanceMode '_eat_' num2str(eatRate,precision) ...
                '_sensingAcc_' num2str(input.sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            if isempty(dir(['../results/' input.saveInfo '_running*.mat']))...
                    &&isempty(dir(['../results/' input.saveInfo '.mat']))
                rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                CA8(input,0);
            end
        end
    end
end