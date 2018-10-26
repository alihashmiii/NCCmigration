function experiment31contStatesNeedNbrs

input.time = 18;
numReps = 40;
input.saveEvery = 10;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')
input.sensingAccuracy = 0.1;
guidanceModes = {'choice','combination'};

for guidanceMode = guidanceModes
    input.guidanceMode = guidanceMode{1};
    for needNeighbours = [0:4]
        input.needNeighbours = needNeighbours;
        for repCtr = 1:numReps
            input.saveInfo = ['experiment31contStates_needNbrs/exp31' ...
                '_contStates_' input.guidanceMode '_needNbrs_' num2str(needNeighbours,precision) ...
                '_sensingAcc_' num2str(input.sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            if isempty(dir(['../results/' input.saveInfo '_running*.mat']))...
                    &&isempty(dir(['../results/' input.saveInfo '.mat']))
                rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                CA8(input,0);
            end
        end
    end
end