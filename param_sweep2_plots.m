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
numFilopodia = [2,2];  % the number of filopodia for lead cells and follower cells
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
diffus = 1;%252e3;    % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
chi = 0.0001;                  % chemoattractant production term (usually 0.0001)
eatRate = 1;                      % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
eatWidth = cellRadius;         % width of eating chemoattractant, equivalent to gaussian sigma
followerFraction = 0.7;        % proportion of cells that are followers (0<=follow_per<=1)

dateString = '2013_07_3*'; % date of simulation that is to be loaded - needs to be adapted if ran on multiple days. Use wildcard (*) for time.
    
collatedResults = load(['results/parameterSweeps/' dateString '-collatedResults']);
saveInfo = collatedResults.saveInfo;
numCells = collatedResults.numCells;
numRepeats = size(numCells,2); % number of runs per parameter combination, to gather stats
numPerturbations = size(numCells,1); % how many parameter combinations have been run altogether, not incl. repeats
meanDirectionality = collatedResults.meanDirectionality; % mean directionality for leaders and followers
meanSpeed = collatedResults.meanSpeed; % mean speed for leaders and followers

%% analyse effects of parameter perturbation on cell number
meanNumCells = mean(numCells,2); % average over multiple runs (of same perturbation)

% parameters with experimental reference values
elasticityLeadSpeed = (meanNumCells(3) - meanNumCells(2))/meanNumCells(1)/(1.1 - 1/1.1);
elasticityFollowSpeed = (meanNumCells(5) - meanNumCells(4))/meanNumCells(1)/(1.1 - 1/1.1);
elasticityFiloLength = (meanNumCells(11) - meanNumCells(10))/meanNumCells(1)/(1.1 - 1/1.1);
elasticityFollowerFraction = (meanNumCells(13) - meanNumCells(12))/meanNumCells(1)/(1.1 - 1/1.1);

% parameters without experimental reference values
elasticityDiffus = (meanNumCells(15) - meanNumCells(14))/meanNumCells(1)/(10 - 1/10);
elasticityChi = (meanNumCells(17) - meanNumCells(16))/meanNumCells(1)/(10 - 1/10);
elasticityEatRate = (meanNumCells(19) - meanNumCells(18))/meanNumCells(1)/(10 - 1/10);
elasticityEatWidth = (meanNumCells(21) - meanNumCells(20))/meanNumCells(1)/(1.1 - 1/1.1);

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [3;2]', 'numFilopodia [4;2]', 'numFilopodia [5;2]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)]};
figure;
boxplot(fliplr(numCells'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on cell number')
xlabel('# of cells at end of simulation')

%% analyse effects of parameter perturbation on leader directionality
meanMeanDirectionality = squeeze(mean(meanDirectionality,2)); % average over multiple runs (of same perturbation)
leaderMMDirectionality = meanMeanDirectionality(:,1);

% parameters with experimental reference values
elasticityLeadSpeed = (leaderMMDirectionality(3) - leaderMMDirectionality(2))/leaderMMDirectionality(1)/(1.1 - 1/1.1);
elasticityFollowSpeed = (leaderMMDirectionality(5) - leaderMMDirectionality(4))/leaderMMDirectionality(1)/(1.1 - 1/1.1);
elasticityFiloLength = (leaderMMDirectionality(11) - leaderMMDirectionality(10))/leaderMMDirectionality(1)/(1.1 - 1/1.1);
elasticityFollowerFraction = (leaderMMDirectionality(13) - leaderMMDirectionality(12))/leaderMMDirectionality(1)/(1.1 - 1/1.1);

% parameters without experimental reference values
elasticityDiffus = (leaderMMDirectionality(15) - leaderMMDirectionality(14))/leaderMMDirectionality(1)/(10 - 1/10);
elasticityChi = (leaderMMDirectionality(17) - leaderMMDirectionality(16))/leaderMMDirectionality(1)/(10 - 1/10);
elasticityEatRate = (leaderMMDirectionality(19) - leaderMMDirectionality(18))/leaderMMDirectionality(1)/(10 - 1/10);
elasticityEatWidth = (leaderMMDirectionality(21) - leaderMMDirectionality(20))/leaderMMDirectionality(1)/(1.1 - 1/1.1);

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [3;2]', 'numFilopodia [4;2]', 'numFilopodia [5;2]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)]};
figure;
boxplot(fliplr(squeeze(meanDirectionality(:,:,1))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on leader directionality')
xlabel('directionality of leader cells (mean per simulation)')

%% analyse effects of parameter perturbation on follower directionality
followerMMDirectionality = meanMeanDirectionality(:,2);

% parameters with experimental reference values
elasticityLeadSpeed = (followerMMDirectionality(3) - followerMMDirectionality(2))/followerMMDirectionality(1)/(1.1 - 1/1.1);
elasticityFollowSpeed = (followerMMDirectionality(5) - followerMMDirectionality(4))/followerMMDirectionality(1)/(1.1 - 1/1.1);
elasticityFiloLength = (followerMMDirectionality(11) - followerMMDirectionality(10))/followerMMDirectionality(1)/(1.1 - 1/1.1);
elasticityFollowerFraction = (followerMMDirectionality(13) - followerMMDirectionality(12))/followerMMDirectionality(1)/(1.1 - 1/1.1);

% parameters without experimental reference values
elasticityDiffus = (followerMMDirectionality(15) - followerMMDirectionality(14))/followerMMDirectionality(1)/(10 - 1/10);
elasticityChi = (followerMMDirectionality(17) - followerMMDirectionality(16))/followerMMDirectionality(1)/(10 - 1/10);
elasticityEatRate = (followerMMDirectionality(19) - followerMMDirectionality(18))/followerMMDirectionality(1)/(10 - 1/10);
elasticityEatWidth = (followerMMDirectionality(21) - followerMMDirectionality(20))/followerMMDirectionality(1)/(1.1 - 1/1.1);

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [3;2]', 'numFilopodia [4;2]', 'numFilopodia [5;2]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)]};
figure;
boxplot(fliplr(squeeze(meanDirectionality(:,:,2))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on follower directionality')
xlabel('directionality of follower cells (mean per simulation)')

%% analyse effects of parameter perturbation on leader speed
meanMeanSpeed = squeeze(mean(meanSpeed,2)); % average over multiple runs (of same perturbation)
leaderMMSpeed = meanMeanSpeed(:,1);

% parameters with experimental reference values
elasticityLeadSpeed = (leaderMMSpeed(3) - leaderMMSpeed(2))/leaderMMSpeed(1)/(1.1 - 1/1.1);
elasticityFollowSpeed = (leaderMMSpeed(5) - leaderMMSpeed(4))/leaderMMSpeed(1)/(1.1 - 1/1.1);
elasticityFiloLength = (leaderMMSpeed(11) - leaderMMSpeed(10))/leaderMMSpeed(1)/(1.1 - 1/1.1);
elasticityFollowerFraction = (leaderMMSpeed(13) - leaderMMSpeed(12))/leaderMMSpeed(1)/(1.1 - 1/1.1);

% parameters without experimental reference values
elasticityDiffus = (leaderMMSpeed(15) - leaderMMSpeed(14))/leaderMMSpeed(1)/(10 - 1/10);
elasticityChi = (leaderMMSpeed(17) - leaderMMSpeed(16))/leaderMMSpeed(1)/(10 - 1/10);
elasticityEatRate = (leaderMMSpeed(19) - leaderMMSpeed(18))/leaderMMSpeed(1)/(10 - 1/10);
elasticityEatWidth = (leaderMMSpeed(21) - leaderMMSpeed(20))/leaderMMSpeed(1)/(1.1 - 1/1.1);

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [3;2]', 'numFilopodia [4;2]', 'numFilopodia [5;2]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)]};
figure;
boxplot(fliplr(squeeze(meanSpeed(:,:,1))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on effective leader speed')
xlabel('speed of leader cells (\mum/hr, mean per simulation)')

%% analyse effects of parameter perturbation on follower speed
followerMMSpeed = meanMeanSpeed(:,2);

% parameters with experimental reference values
elasticityLeadSpeed = (followerMMSpeed(3) - followerMMSpeed(2))/followerMMSpeed(1)/(1.1 - 1/1.1);
elasticityFollowSpeed = (followerMMSpeed(5) - followerMMSpeed(4))/followerMMSpeed(1)/(1.1 - 1/1.1);
elasticityFiloLength = (followerMMSpeed(11) - followerMMSpeed(10))/followerMMSpeed(1)/(1.1 - 1/1.1);
elasticityFollowerFraction = (followerMMSpeed(13) - followerMMSpeed(12))/followerMMSpeed(1)/(1.1 - 1/1.1);

% parameters without experimental reference values
elasticityDiffus = (followerMMSpeed(15) - followerMMSpeed(14))/followerMMSpeed(1)/(10 - 1/10);
elasticityChi = (followerMMSpeed(17) - followerMMSpeed(16))/followerMMSpeed(1)/(10 - 1/10);
elasticityEatRate = (followerMMSpeed(19) - followerMMSpeed(18))/followerMMSpeed(1)/(10 - 1/10);
elasticityEatWidth = (followerMMSpeed(21) - followerMMSpeed(20))/followerMMSpeed(1)/(1.1 - 1/1.1);

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [3;2]', 'numFilopodia [4;2]', 'numFilopodia [5;2]', 'numFilopodia [6;2]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)]};
figure;
boxplot(fliplr(squeeze(meanSpeed(:,:,2))'),'labels',fliplr(labels),'orientation','horizontal')
title('effects of parameter perturbation on effective follower speed')
xlabel('speed of follower cells (\mum/hr, mean per simulation)')

