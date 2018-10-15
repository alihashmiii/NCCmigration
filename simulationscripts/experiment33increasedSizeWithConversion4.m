function experiment33increasedSizeWithConversion4
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13, 14.08.2015

input.time = 18;
input.conversionType = 4;
numReps = 20;

precision = 2; % significant figures for filenames and plot labels etc.

for insertEverySteps = [1 3]
    input.insertEverySteps = insertEverySteps;
    for defaultFollow = [1 2] % determines which is the default behaviour of cells, before switching
        input.followerFraction = defaultFollow;
        for sensingAccuracy = [0.1, 0.01]
            input.sensingAccuracy = sensingAccuracy;
            for switchingTime = [4 8]
                input.numSteps = [switchingTime, switchingTime];
                for repCtr = 1:numReps
                    input.saveInfo = ['experiment33/exp33'...
                        '_conversion_' num2str(input.conversionType) '_defaultFollow_' num2str(input.followerFraction) ...
                        '_numSteps_' num2str(input.numSteps(1)) '_' num2str(input.numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_insertSteps_' num2str(input.insertEverySteps)...
                        '_Run_' num2str(repCtr)];
                    if isempty(dir(['results/' input.saveInfo '_running*.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
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

