% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width
% has a zone of reduced cell speeds at the start (grows with tissue & dilutes)

clear
close all

in.time = 24;
in.conversionType = 4;
in.domainHeight = 360;

experiment = 41;
precision = 2;

% setting the initial CA outside the migratory path to zero in case of a wider domain
x = linspace(0,300,64);
y = linspace(0,in.domainHeight,32);
ca = zeros(length(x),length(y));
middleStripe = y<=(in.domainHeight/2 + 60)&y>=(in.domainHeight/2 - 60);
ca(:,middleStripe) = 1;
in.ca_new = implicit_heat2D(ca,100,mean(diff(x)),mean(diff(y)),0.1,1); % we need some smoothing to reduce the sharp boundaries

diffusivities = [1, 10, 100];
slowSpeeds = [5, 10, 20];
in.leadSpeed = 40;
in.followSpeed = 40;
fileName = 'exp41_slowEntryDilute';
for cntGdn = {'toward', 'parallel'}
    in.contactGuidance = char(cntGdn);
    for diffus = diffusivities
        in.diffus = diffus;
        in.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        for slowSpeed = slowSpeeds
            in.slowSpeed = slowSpeed;
            in.saveInfo = [fileName '/' fileName ...
                '_D_' num2str(in.diffus)  '_sensingAcc_' num2str(in.sensingAccuracy,precision) ...
                '_slowSpeed_' num2str(in.slowSpeed) '_contactGuidance_' in.contactGuidance];
            CA6(in,experiment)
        end
    end
end