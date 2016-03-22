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

diffusivities = [1, 10, 100];
speeds = [10, 20, 40];

for cntGdn = ['toward', 'parallel']
    in.contactGuidance = cntGdn;
    for diffus = diffusivities
        in.diffus = diffus;
        in.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        for speed = speeds
            in.leadSpeed = speed;
            in.followSpeed = speed;
            in.saveInfo = ['exp39_widerDomainStripe/exp39_widerDomainStripe_' ...
                '_D_' num2str(in.diffus)  '_sensingAcc_' num2str(in.sensingAccuracy,precision) ...
                '_speed_' num2str(in.leadSpeed) '_contactGuidance_' cntGdn];
            CA6(in,experiment)
        end
    end
end