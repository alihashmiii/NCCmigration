function experiment31ChiWithoutPlasticity
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 18;
numReps = 20;
input.chi = 1;
precision = 2; % significant figures for filenames and plot labels etc.

for followerFraction = [0, 1]
    input.followerFraction = followerFraction;
        for sensingAccuracy = 0.1
            input.sensingAccuracy = sensingAccuracy;
            for repCtr = 1:numReps
                input.saveInfo = ['experiment31Chi/exp31_followFrac_' num2str(followerFraction,precision) ...
                    '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                if isempty(dir(['results/' input.saveInfo '_running.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
                    rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                    CA6(input,0);
                    % check if anyone else is logged into the
                    % current machine. If so, quit. If not,
                    % keep looping and run next job.
                    [~, logins_check] = system('/mi/libexec/check-logins | grep -v schumacher');
                    if size(logins_check,1)>0
                        %                         quit
                    end
                end
            end
        end
end