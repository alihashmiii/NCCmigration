function experiment31contStatesPstayattached

input.time = 18;
numReps = 40;
input.saveEvery = 10;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')
input.sensingAccuracy = 0.1;
guidanceModes = {'choice','combination'};

for guidanceMode = guidanceModes
    input.guidanceMode = guidanceMode{1};
    for p_stayattached = [0.8, 0.85, 0.9, 0.95, 0.98, 0.99, 1]
        input.p_stayattached = p_stayattached;
        for repCtr = 1:numReps
            input.saveInfo = ['experiment31contStates_Psa/exp31' ...
                '_contStates_' input.guidanceMode '_Psa_' num2str(p_stayattached,precision) ...
                '_sensingAcc_' num2str(input.sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            if isempty(dir(['../results/' input.saveInfo '_running*.mat']))...
                    &&isempty(dir(['../results/' input.saveInfo '.mat']))
                rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                CA8(input,0);
            end
        end
    end
end