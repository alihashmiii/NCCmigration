function experiment31contStatesChi

input.time = 18;
numReps = 40;
input.saveEvery = 10;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')
input.sensingAccuracy = 0.1;
guidanceModes = {'choice','combination'};

for guidanceMode = guidanceModes
    input.guidanceMode = guidanceMode{1};
    for chi = [1e-2, 1e-3, 1e-4, 1e-5, 10^(-5.5), 1e-6]
        input.chi = chi;
        for repCtr = 1:numReps
            input.saveInfo = ['experiment31contStates_chi/exp31' ...
                '_contStates_' input.guidanceMode '_chi_' num2str(chi,precision) ...
                '_sensingAcc_' num2str(input.sensingAccuracy,precision) '_Run_' num2str(repCtr)];
            if isempty(dir(['../results/' input.saveInfo '_running*.mat']))...
                    &&isempty(dir(['../results/' input.saveInfo '.mat']))
                rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                CA8(input,0);
            end
        end
    end
end