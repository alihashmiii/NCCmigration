% parameter sweep for Dyson model CA6
% perturbed parameters around reference parameter set
% plot results and calculate elasticities
% LJSchumacher 29.07.2013

clear all
close all

% these parameters are not sweeped, but needed to calculate some reference values
cellRadius = 7.5;              % radius in um (= 7.5um)

% reference parameter set -- these should be the same as defautls in CA6.m
leadSpeed = 41.6;                     % speed of the leader cells in mu/h
followSpeed = 49.9;                 % speed of the follower cells in mu/h
numFilopodia = [3,1];  % the number of filopodia for lead cells and follower cells
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
diffus = 0.1;%252e3;    % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
chi = 0.0001;                  % chemoattractant production term (usually 0.0001)
eatRate = 100;                      % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
eatWidth = cellRadius;         % width of eating chemoattractant, equivalent to gaussian sigma
followerFraction = 6/8;        % proportion of cells that are followers (0<=follow_per<=1)
tstep = 5/2/60;                   % time step in hours

dateString = '2013_10_10'; % date of simulation that is to be loaded - needs to be adapted if ran on multiple days. Use wildcard (*) for time.
    
collatedResults = load(['results/parameterSweeps/' dateString '-collatedResults']);
saveInfo = collatedResults.saveInfo;
numCells = collatedResults.numCells;
numRepeats = size(numCells,2); % number of runs per parameter combination, to gather stats
numPerturbations = size(numCells,1); % how many parameter combinations have been run altogether, not incl. repeats
meanDirectionality = collatedResults.meanDirectionality; % mean directionality for leaders and followers
meanSpeed = collatedResults.meanSpeed; % mean speed for leaders and followers
xBins = 0:2*cellRadius:49*2*cellRadius; % check the end length for time step run, or try loading this from the out-file somehow
cellDistributions = collatedResults.cellDistributions;
caDistribution = collatedResults.caDistribution;

speedFactor = 1.25;
filoLengthFactor = 1.5;
followerFractions = [5/8, 7/8];
diffusFactor = 252e4;
chiFactor = 10;
eatRateFactor = 10;
eatWidthFactor = 1.6;
tStepFactor = 2;

deltaSpeed = speedFactor - 1/speedFactor;
deltaFiloLength = filoLengthFactor - 1/filoLengthFactor;
deltaFollowerFraction = 7/6 - 5/6;
deltaDiffus = diffusFactor - 1/diffusFactor;
deltaChi = chiFactor - 1/chiFactor;
deltaEatRate = eatRateFactor - 1/eatRateFactor;
deltaEatWidth = eatWidthFactor - 1/eatWidthFactor;
deltaTstep = tStepFactor - 1/tStepFactor;

%% analyse effects of parameter perturbation on cell number
meanNumCells = mean(numCells,2); % average over multiple runs (of same perturbation)
outputQuantity = meanNumCells;

% parameters with experimental reference values
elasticityLeadSpeed = (outputQuantity(3) - outputQuantity(2))/outputQuantity(1)/deltaSpeed;
elasticityFollowSpeed = (outputQuantity(5) - outputQuantity(4))/outputQuantity(1)/deltaSpeed;
elasticityFiloLength = (outputQuantity(9) - outputQuantity(8))/outputQuantity(1)/deltaFiloLength;
elasticityFollowerFraction = (outputQuantity(11) - outputQuantity(10))/outputQuantity(1)/deltaFollowerFraction;

% parameters without experimental reference values
elasticityDiffus = (outputQuantity(13) - outputQuantity(12))/outputQuantity(1)/deltaDiffus;
elasticityChi = (outputQuantity(15) - outputQuantity(14))/outputQuantity(1)/deltaChi;
elasticityEatRate = (outputQuantity(17) - outputQuantity(16))/outputQuantity(1)/deltaEatRate;
elasticityEatWidth = (outputQuantity(19) - outputQuantity(18))/outputQuantity(1)/deltaEatWidth;
elasticityTstep = (outputQuantity(21) - outputQuantity(20))/outputQuantity(1)/deltaTstep;

labels = {'reference', ['leadSpeed /' num2str(speedFactor)], ['leadSpeed x' num2str(speedFactor)],...
    ['followSpeed /' num2str(speedFactor)],['followSpeed x' num2str(speedFactor)], ...%'numFilopodia [2;1]', 'numFilopodia [6;2]', ...
    ['filoLength /' num2str(filoLengthFactor)], ['filoLength x' num2str(filoLengthFactor)], ...
    ['followerFraction ' num2str(followerFractions(1))], ['followerFraction ' num2str(followerFractions(2))], ...
    ['diffusion /' num2str(diffusFactor)], ['diffusion x' num2str(diffusFactor)],...
    ['CA production /' num2str(chiFactor)], ['CA production x' num2str(chiFactor)], ...
    ['eatRate /' num2str(eatRateFactor)], ['eatRate x' num2str(eatRateFactor)],...
    ['eatWidth /' num2str(eatWidthFactor)], ['eatWidth x' num2str(eatWidthFactor)] ...%,'timeStep', ['e = ' num2str(elasticityTstep)]...
    };

figCellNumber = figure;
boxplot(fliplr(numCells([1:5, 8:19],:)'),'labels',fliplr(labels),'orientation','horizontal')
text(50,15.5,['e = ' num2str(elasticityLeadSpeed)]);
text(50,13.5,['e = ' num2str(elasticityFollowSpeed)]);
text(50,11,['e = ' num2str(elasticityFiloLength)]);
text(50,10,['e = ' num2str(elasticityFollowerFraction)]);
text(50,8,['e = ' num2str(elasticityDiffus)]);
text(50,5.5,['e = ' num2str(elasticityChi)]);
text(50,3.5,['e = ' num2str(elasticityEatRate)]);
text(50,1.5,['e = ' num2str(elasticityEatWidth)]);
title('effects of parameter perturbation on cell number')
xlabel('# of cells at end of simulation')

%% analyse effects of parameter perturbation on leader directionality
meanMeanDirectionality = squeeze(mean(meanDirectionality,2)); % average over multiple runs (of same perturbation)
leaderMMDirectionality = meanMeanDirectionality(:,1);

outputQuantity = leaderMMDirectionality;

% parameters with experimental reference values
elasticityLeadSpeed = (outputQuantity(3) - outputQuantity(2))/outputQuantity(1)/deltaSpeed;
elasticityFollowSpeed = (outputQuantity(5) - outputQuantity(4))/outputQuantity(1)/deltaSpeed;
elasticityFiloLength = (outputQuantity(9) - outputQuantity(8))/outputQuantity(1)/deltaFiloLength;
elasticityFollowerFraction = (outputQuantity(11) - outputQuantity(10))/outputQuantity(1)/deltaFollowerFraction;

% parameters without experimental reference values
elasticityDiffus = (outputQuantity(13) - outputQuantity(12))/outputQuantity(1)/deltaDiffus;
elasticityChi = (outputQuantity(15) - outputQuantity(14))/outputQuantity(1)/deltaChi;
elasticityEatRate = (outputQuantity(17) - outputQuantity(16))/outputQuantity(1)/deltaEatRate;
elasticityEatWidth = (outputQuantity(19) - outputQuantity(18))/outputQuantity(1)/deltaEatWidth;
elasticityTstep = (outputQuantity(21) - outputQuantity(20))/outputQuantity(1)/deltaTstep;

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [2;1]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)], ...
    'timeStep', ['e = ' num2str(elasticityTstep)]};

figLeaderDirectionality = figure;
boxplot(fliplr(squeeze(meanDirectionality(:,:,1))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on leader directionality')
xlabel('directionality of leader cells (mean per simulation)')

%% analyse effects of parameter perturbation on follower directionality
followerMMDirectionality = meanMeanDirectionality(:,2);

outputQuantity = followerMMDirectionality;

% parameters with experimental reference values
elasticityLeadSpeed = (outputQuantity(3) - outputQuantity(2))/outputQuantity(1)/deltaSpeed;
elasticityFollowSpeed = (outputQuantity(5) - outputQuantity(4))/outputQuantity(1)/deltaSpeed;
elasticityFiloLength = (outputQuantity(9) - outputQuantity(8))/outputQuantity(1)/deltaFiloLength;
elasticityFollowerFraction = (outputQuantity(11) - outputQuantity(10))/outputQuantity(1)/deltaFollowerFraction;

% parameters without experimental reference values
elasticityDiffus = (outputQuantity(13) - outputQuantity(12))/outputQuantity(1)/deltaDiffus;
elasticityChi = (outputQuantity(15) - outputQuantity(14))/outputQuantity(1)/deltaChi;
elasticityEatRate = (outputQuantity(17) - outputQuantity(16))/outputQuantity(1)/deltaEatRate;
elasticityEatWidth = (outputQuantity(19) - outputQuantity(18))/outputQuantity(1)/deltaEatWidth;
elasticityTstep = (outputQuantity(21) - outputQuantity(20))/outputQuantity(1)/deltaTstep;

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [2;1]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)], ...
    'timeStep', ['e = ' num2str(elasticityTstep)]};

figFollowerDirectionality = figure;
boxplot(fliplr(squeeze(meanDirectionality(:,:,2))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on follower directionality')
xlabel('directionality of follower cells (mean per simulation)')

%% analyse effects of parameter perturbation on leader speed
meanMeanSpeed = squeeze(mean(meanSpeed,2)); % average over multiple runs (of same perturbation)
leaderMMSpeed = meanMeanSpeed(:,1);

outputQuantity = leaderMMSpeed;

% parameters with experimental reference values
elasticityLeadSpeed = (outputQuantity(3) - outputQuantity(2))/outputQuantity(1)/deltaSpeed;
elasticityFollowSpeed = (outputQuantity(5) - outputQuantity(4))/outputQuantity(1)/deltaSpeed;
elasticityFiloLength = (outputQuantity(9) - outputQuantity(8))/outputQuantity(1)/deltaFiloLength;
elasticityFollowerFraction = (outputQuantity(11) - outputQuantity(10))/outputQuantity(1)/deltaFollowerFraction;

% parameters without experimental reference values
elasticityDiffus = (outputQuantity(13) - outputQuantity(12))/outputQuantity(1)/deltaDiffus;
elasticityChi = (outputQuantity(15) - outputQuantity(14))/outputQuantity(1)/deltaChi;
elasticityEatRate = (outputQuantity(17) - outputQuantity(16))/outputQuantity(1)/deltaEatRate;
elasticityEatWidth = (outputQuantity(19) - outputQuantity(18))/outputQuantity(1)/deltaEatWidth;
elasticityTstep = (outputQuantity(21) - outputQuantity(20))/outputQuantity(1)/deltaTstep;

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [2;1]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)], ...
    'timeStep', ['e = ' num2str(elasticityTstep)]};

figLeaderSpeed = figure;
boxplot(fliplr(squeeze(meanSpeed(:,:,1))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on effective leader speed')
xlabel('speed of leader cells (\mum/hr, mean per simulation)')

%% analyse effects of parameter perturbation on follower speed
followerMMSpeed = meanMeanSpeed(:,2);

outputQuantity = followerMMSpeed;

% parameters with experimental reference values
elasticityLeadSpeed = (outputQuantity(3) - outputQuantity(2))/outputQuantity(1)/deltaSpeed;
elasticityFollowSpeed = (outputQuantity(5) - outputQuantity(4))/outputQuantity(1)/deltaSpeed;
elasticityFiloLength = (outputQuantity(9) - outputQuantity(8))/outputQuantity(1)/deltaFiloLength;
elasticityFollowerFraction = (outputQuantity(11) - outputQuantity(10))/outputQuantity(1)/deltaFollowerFraction;

% parameters without experimental reference values
elasticityDiffus = (outputQuantity(13) - outputQuantity(12))/outputQuantity(1)/deltaDiffus;
elasticityChi = (outputQuantity(15) - outputQuantity(14))/outputQuantity(1)/deltaChi;
elasticityEatRate = (outputQuantity(17) - outputQuantity(16))/outputQuantity(1)/deltaEatRate;
elasticityEatWidth = (outputQuantity(19) - outputQuantity(18))/outputQuantity(1)/deltaEatWidth;
elasticityTstep = (outputQuantity(21) - outputQuantity(20))/outputQuantity(1)/deltaTstep;

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [2;1]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)], ...
    'timeStep', ['e = ' num2str(elasticityTstep)]};

figFollowerSpeed = figure;
boxplot(fliplr(squeeze(meanSpeed(:,:,2))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on effective follower speed')
xlabel('speed of follower cells (\mum/hr, mean per simulation)')


%% migration profile along x

figxProfile = figure;
plotOrder = [10, 11, 11];
followerFractions = [5/8, 6/8, 7/8];
for plotCtr = 1:3
subplot(3,1,plotCtr)
[AX,H1,H2] = plotyy(xBins,mean(squeeze(cellDistributions(plotOrder(plotCtr),:,2,:))),...
    xBins,mean(squeeze(caDistribution(plotOrder(plotCtr),:,:)))/numRepeats);
hold(AX(1))
H3 = plot(AX(1),xBins,mean(squeeze(cellDistributions(plotOrder(plotCtr),:,1,:))),'b-');
set(H1,'LineStyle','--','Color','b')
set(H2,'LineStyle',':','Color',[0 0.5 0])
xlabel('x/\mum')
ylabel(AX(1),'P(cells)')
ylabel(AX(2),'C(chemoattractant)')
ylim(AX(1),[0, 0.075])
xlim(AX(1),[0 735]), xlim(AX(2),[0 735])
title(['follower fraction ' num2str(followerFractions(plotCtr))])
end
legend('leaders','followers','chemoattractant')

%% for exporting as pdf
% % set(gcf, 'PaperPositionMode', 'auto'); % set size to as seen on screen
% % print(gcf, '-r0','../../TransferOfStatus/ThesisProposal/figures/velocityCorrelation.eps', '-depsc2');
% % or ...
% set(figCellNumber,'PaperUnits','centimeters');
% set(figLeaderDirectionality,'PaperUnits','centimeters');
% set(figFollowerDirectionality,'PaperUnits','centimeters');
% set(figLeaderSpeed,'PaperUnits','centimeters');
% set(figFollowerSpeed,'PaperUnits','centimeters');
set(figxProfile,'PaperUnits','centimeters');
% 
exportOptions = struct('Format','eps2',...
    'Width','15.0',...
    'Color','cmyk',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',8,...
    'LineWidth',2);

% filename = 'results/tmp/statsCellNumber5.eps';
% exportfig(figCellNumber,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsLeaderDirectionality5.eps';
% exportfig(figLeaderDirectionality,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsFollowerDirectionality5.eps';
% exportfig(figFollowerDirectionality,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsLeaderSpeed5.eps';
% exportfig(figLeaderSpeed,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsFollowerSpeed5.eps';
% exportfig(figFollowerSpeed,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
filename = 'results/tmp/migrationProfile5.eps';
exportfig(figxProfile,filename,exportOptions);
system(['epstopdf ' filename]);