function experiment31contStatesCombination

input.time = 18;
numReps = 40;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')
input.guidanceMode = 'combination';
input.saveEvery = 10;

for sensingAccuracy = [0.1, 0.01]
    input.sensingAccuracy = sensingAccuracy;
    for repCtr = 1:numReps
        input.saveInfo = ['experiment31contStates_combination/exp31' ...
            '_contStates_' input.guidanceMode  ...
            '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
        if isempty(dir(['../results/' input.saveInfo '_running*.mat']))&&...
                isempty(dir(['../results/' input.saveInfo '.mat']))
            rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
            CA8(input,0);
        end
    end
end