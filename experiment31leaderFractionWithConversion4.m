function experiment31leaderFractionWithConversion4
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 18;
input.conversionType = 4;
numReps = 20;

precision = 2; % significant figures for filenames and plot labels etc.
input.standStill = 0;
input.tstep = 5/4/60;
eatRate = 1000;
defaultFollow = 1;
input.followerFraction = defaultFollow;

for sensingAccuracy = [0.1, 0.01]
    input.sensingAccuracy = sensingAccuracy;
    for needNeighbours = [0, 1, 2]
        input.needNeighbours = needNeighbours;
        for lead2follow = [4 8]
            for follow2lead = [8 16 30]
                input.numSteps = [lead2follow, follow2lead];
                for repCtr = 1:numReps
                    input.saveInfo = ['experiment31conversion4/exp31_followFrac_' num2str(input.followerFraction,precision) '_eatRate_' num2str(eatRate) ...
                        '_conversion_' num2str(input.conversionType) '_numSteps_' num2str(input.numSteps(1)) '_' num2str(input.numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_needNeighbours_' num2str(needNeighbours)...
                        '_tstep_' num2str(input.tstep,precision) '_Run_' num2str(repCtr)];
                    if isempty(dir(['results/' input.saveInfo '_running.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
                        rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                        CA6(input,0);
                        % check if anyone else is logged into the
                        % current machine. If so, quit. If not,
                        % keep looping and run next job.
                        [~, logins_check] = system('/mi/libexec/check-logins | grep -v schumacher');
                        if size(logins_check,1)>0
                            quit
                        end
                    end
                end
            end
        end
    end
end