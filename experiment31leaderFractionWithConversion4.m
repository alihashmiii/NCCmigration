b   b   function experiment31leaderFractionWithConversion4
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 30;
input.conversionType = 4;
numReps = 20;

precision = 2; % significant figures for filenames and plot labels etc.

for defaultFollow = [0 1 2]
    input.followerFraction = defaultFollow;
    for sensingAccuracy = [0.1, 0.01]
        input.sensingAccuracy = sensingAccuracy;
        for lead2follow = [1 2 4 8 12 16 24 32 40 48 56]
            for follow2lead = [1 2 4 8 12 16 24 32 40 48 56]
                input.numSteps = [lead2follow, follow2lead];
                for repCtr = 1:numReps
                    input.saveInfo = ['experiment31conversion4/exp31' ...
                        '_conversion_' num2str(input.conversionType) '_defaultFollow_' num2str(input.followerFraction) ...
                        '_numSteps_' num2str(input.numSteps(1)) '_' num2str(input.numSteps(2)) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
% %                     if isempty(dir(['results/' input.saveInfo '_running.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
% %                         rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
% %                         CA6(input,0);
% %                         % check if anyone else is logged into the
% %                         % current machine. If so, quit. If not,
% %                         % keep looping and run next job.
% %                         [~, logins_check] = system('/mi/libexec/check-logins | grep -v schumacher');
% %                         if size(logins_check,1)>0
% %                             quit
% %                         end
% %                     end
                                % delete some of the saved variables that aren't needed
            % here to save disk space
            if ~isempty(dir(['results/' input.saveInfo '.mat']))
                load(['results/',input.saveInfo,'.mat'])
                out.moved = out.moved(:,1:length(out.cells_save{end}));
                out.happiness = out.happiness(:,1:length(out.cells_save{end}));
                for timeCtr = 1:length(cellsFollow_save)
                    out.cellsFollow_save{timeCtr} = ...
                        out.cellsFollow_save{timeCtr}(1:size(out.cells_save{timeCtr},2));
                    out.attach_save{timeCtr} = ...
                        out.attach_save{timeCtr}(1:size(out.cells_save{timeCtr},2));
                end
                save(['results/',input.saveInfo,'.mat'],'out')
                clear out
            end
                end
            end
        end
    end
end