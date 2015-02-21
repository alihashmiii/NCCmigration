function experiment31halfFiloWithConversion4
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 18;
numReps = 20;
input.conversionType = 4;
cellRadius = 7.5;
input.filolength = cellRadius + 10;
precision = 2; % significant figures for filenames and plot labels etc.

for defaultFollow = [2]
    input.followerFraction = defaultFollow;
    for sensingAccuracy = 0.1
        input.sensingAccuracy = sensingAccuracy;
        for lead2follow = [4]
            for follow2lead = [4]
                input.numSteps = [lead2follow, follow2lead];
                for repCtr = 1:numReps
                    input.saveInfo = ['experiment31halfFilo/exp31'...
                        '_conversion_' num2str(input.conversionType) '_defaultFollow_' num2str(input.followerFraction) ...
                        '_numSteps_' num2str(input.numSteps(1)) '_' num2str(input.numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
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