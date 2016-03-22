% simulate neural crest cell migration on wider comain with vegf production
% only in the middle

clear
close all

% issues/to-dos:
% - change file name automatically
% - vary cell speed in simulations

input.time = 24;
input.conversionType = 4;
input.domainHeight = 360;
input.initYFrac = 1/3;

experiment = 38;

% setting the initial CA outside the migratory path to zero in case of a wider domain
x = linspace(0,300,64);
y = linspace(0,input.domainHeight,32);
ca = zeros(length(x),length(y));
middleStripe = y<=(input.domainHeight/2 + 60)&y>=(input.domainHeight/2 - 60);
ca(:,middleStripe) = 1;
input.ca_new = implicit_heat2D(ca,100,mean(diff(x)),mean(diff(y)),0.1,1); % we need some smoothing to reduce the sharp boundaries

input.eatRate = 1000;
input.diffus = 0.1;
input.sensingAccuracy = 0.01; % sens acc scales with 1/sqrt(diffus)
input.leadSpeed = 40;
input.followSpeed = 50;

input.saveInfo = ['exp39_widerDomainStripe_eatRate' num2str(input.eatRate

CA6(input,experiment)