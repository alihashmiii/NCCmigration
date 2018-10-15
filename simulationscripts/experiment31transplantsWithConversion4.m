function experiment31transplantsWithConversion4
% run the next job in line of a bunch
% called by scan-submit shell script
% looks which job is the next one to do,runs this one job, then quits
% L.J. Schumacher 28.10.13

input.time = 18 - 6;
input.conversionType = 4;
numReps = 20;
input.transplantTime = 12; % time at which CA-production will be locally increased

precision = 2; % significant figures for filenames and plot labels etc.

for experiment = [12 11]
    for defaultFollow = [0 1 2]
        input.followerFraction = defaultFollow;
        for sensingAccuracy = [0.1, 0.01]
            input.sensingAccuracy = sensingAccuracy;
            for switchingTime = [4 8]
                input.numSteps = [switchingTime, switchingTime];
                for repCtr = 1:numReps
                    input.saveInfo = ['experiment31transplants/exp' num2str(experiment) ...
                        '_conversion_' num2str(input.conversionType) '_defaultFollow_' num2str(input.followerFraction) ...
                        '_numSteps_' num2str(switchingTime) '_' num2str(switchingTime) ...
                        '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                    if isempty(dir(['results/' input.saveInfo '_running*.mat']))&&isempty(dir(['results/' input.saveInfo '.mat']))
                        % load the relevant reference simulation
                        loadInfo = ['experiment31conversion4/exp31'...
                            '_conversion_4_defaultFollow_' num2str(defaultFollow) ...
                            '_numSteps_' num2str(switchingTime) '_' num2str(switchingTime) ...
                            '_sensingAcc_' num2str(sensingAccuracy) '_Run_' num2str(repCtr)];
                        load(['results/' loadInfo '.mat'])
                        % transfer relevant outputs to inputs
                        timeIdx = find(out.t_save>=input.transplantTime,1,'first');
                        input.t_0 = double(out.t_save(timeIdx));
                        input.cells = double(out.cells_save{timeIdx});
                        input.xlat_new = double(out.xlat_save{timeIdx});
                        input.cellsFollow = out.cellsFollow_save{timeIdx};
                        input.attach = out.attach_save{timeIdx};
                        input.happiness = double(out.happiness(timeIdx,1:size(input.cells,2)));
                        % adapt the chemoattractant distribution to transplant
                        switch experiment
                            case 12 %VEGF transplant back half
                                transplantXLocation = 0; % left edge of square region in which CA-production will be increased
                                xindcs =  (input.xlat_new>=transplantXLocation)&(input.xlat_new<=(transplantXLocation + 1/8*max(input.xlat_new)));
                                yindcs = out.ylat_save{timeIdx}<=out.param.domainHeight/2;
                            case 13 %VEGF transplant middle half
                                transplantXLocation = 70; % left edge of square region in which CA-production will be increased
                                xindcs =  (input.xlat_new>=transplantXLocation)&(input.xlat_new<=(transplantXLocation + 1/8*max(input.xlat_new)));
                                yindcs = out.ylat_save{timeIdx}<=out.param.domainHeight/2;
                            case 11 %VEGF transplant back edge
                                transplantXLocation = 0; % left edge of square region in which CA-production will be increased
                                xindcs =  (input.xlat_new>=transplantXLocation)&(input.xlat_new<=(transplantXLocation + 1/8*max(input.xlat_new)));
                                yindcs = out.ylat_save{timeIdx}<=out.param.domainHeight/20;
                            case 14 %increased VEGF at far (right-most) edge
                                transplantXLocation = 425; % left edge of square region in which CA-production will be increased
                                xindcs =  (input.xlat_new>=transplantXLocation)&(input.xlat_new<=max(input.xlat_new));
                                yindcs = true(size(out.ylate_save{timeIdx}));
                            otherwise
                                input.transplantXLocation = NaN; xindcs = []; yindcs = [];
                        end
                        % set region of transplant to high CA
                        transplant = implicit_heat2D(double(xindcs)*double(yindcs)',...
                            100,mean(diff(out.xlat_save{timeIdx})),mean(diff(out.ylat_save{timeIdx})),0.1,1); % we need some smoothing to reduce the sharp boundaries
                        input.ca_new = double(out.ca_save{timeIdx}) + transplant;
                        input.ca_new = min(input.ca_new,1);
                        
                        rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                        CA6(input,experiment);
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