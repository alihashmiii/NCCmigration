% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width
% has a zone of reduced cell speeds at the start (grows with tissue)
% DAN zone near domain entrance slows down cells, increases then decreases
% with time linearly, but is broken down as cells move through it

clear
close all

in.time = 24;
in.conversionType = 4;
in.domainHeight = 360;

experiment = 44;
precision = 2;

% setting the initial CA outside the migratory path to zero in case of a wider domain
x = linspace(0,300,64);
y = linspace(0,in.domainHeight,32);
ca = zeros(length(x),length(y));
middleStripe = y<=(in.domainHeight/2 + 60)&y>=(in.domainHeight/2 - 60);
ca(:,middleStripe) = 1;
in.ca_new = implicit_heat2D(ca,100,mean(diff(x)),mean(diff(y)),0.1,1); % we need some smoothing to reduce the sharp boundaries

diffusivities = [1];
slowSpeeds = [10, 30, 40];
in.leadSpeed = 40;
in.followSpeed = 40;
fileName = 'exp44_slowEntryDynamicTunnelingDAN';
numReps = 10;
for repCtr = 1:numReps
    for cntGdn = {'parallel','toward'}
        in.contactGuidance = char(cntGdn);
        for diffus = diffusivities
            in.diffus = diffus;
            in.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
            for slowSpeed = slowSpeeds
                in.slowSpeed = slowSpeed;
                in.saveInfo = [fileName '/' fileName ...
                    '_D_' num2str(in.diffus)  '_sensingAcc_' num2str(in.sensingAccuracy,precision) ...
                    '_slowSpeed_' num2str(in.slowSpeed) '_contactGuidance_' in.contactGuidance  '_Run_' num2str(repCtr)];
                if isempty(dir(['results/' in.saveInfo '_running*.mat']))&&isempty(dir(['results/' in.saveInfo '.mat']))
                    rng('shuffle'); % shuffle random number sequences to not repeat result from previous matlab sessions
                    CA6(in,experiment)
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