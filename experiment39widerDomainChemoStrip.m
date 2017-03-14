% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width

clear
close all

in.time = 24;
in.conversionType = 4;
in.domainHeight = 360;

experiment = 39;
precision = 2;

% setting the initial CA outside the migratory path to zero in case of a wider domain
x = linspace(0,300,64);
y = linspace(0,in.domainHeight,32);
ca = zeros(length(x),length(y));
middleStripe = y<=(in.domainHeight/2 + 60)&y>=(in.domainHeight/2 - 60);
ca(:,middleStripe) = 1;
in.ca_new = implicit_heat2D(ca,100,mean(diff(x)),mean(diff(y)),0.1,1); % we need some smoothing to reduce the sharp boundaries

diffusivities = [1];
speeds = [40];
insertEveryStepsValues = [10]; % corresponding to 4 and 2 (attempted) cell insertions per hour at tstep = 1min
numReps = 10;

for repCtr = 1:numReps
    for cntGdn = {'parallel'}
        in.contactGuidance = char(cntGdn);
        for diffus = diffusivities
            in.diffus = diffus;
            in.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
            for insertEverySteps = insertEveryStepsValues
                in.insertEverySteps = insertEverySteps;
                for speed = speeds
                    in.leadSpeed = speed;
                    in.followSpeed = speed;
                    in.saveInfo = ['exp39_widerDomainStripe/exp39_widerDomainStripe' ...
                        '_D_' num2str(in.diffus)  '_sensingAcc_' num2str(in.sensingAccuracy,precision) ...
                        '_speed_' num2str(in.leadSpeed) '_insertEvry_' num2str(insertEverySteps) ...
                        '_contactGuidance_' in.contactGuidance '_Run_' num2str(repCtr)];
                    CA6(in,experiment)
                end
            end
        end
    end
end