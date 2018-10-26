function experiment31leaderFractionWithConversion4repeat
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 30;
input.conversionType = 4;
numReps = 20;

precision = 2; % significant figures for filenames and plot labels etc.
addpath('../')

for defaultFollow = [0 1 2]
    input.followerFraction = defaultFollow;
    for sensingAccuracy = [0.1, 0.01]
        input.sensingAccuracy = sensingAccuracy;
        for lead2follow = [8]
            for follow2lead = [8]
                input.numSteps = [lead2follow, follow2lead];
                for repCtr = 1:numReps
                    input.saveInfo = ['experiment31conversion4repeat/exp31' ...
                        '_conversion_' num2str(input.conversionType) '_defaultFollow_' num2str(input.followerFraction) ...
                        '_numSteps_' num2str(input.numSteps(1)) '_' num2str(input.numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                    if isempty(dir(['../results/' input.saveInfo '_running*.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
                        rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                        CA6(input,0);
                    end
                end
            end
        end
    end
end