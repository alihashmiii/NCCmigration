% parameter sweep for Dyson model CA6
% perturbed parameters around reference parameter set
% plot results and calculate elasticities
% LJSchumacher 29.07.2013

clear all

% these parameters are not sweeped, but needed to calculate some reference
% values
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

dateString = '2013_06_29'; % date of simulation that is to be loaded - needs to be adapted if ran on multiple days. Use wildcard (*) for time.
    
collatedResults = load(['results/parameterSweeps/' dateString '-collatedResults']);
saveInfo = collatedResults.saveInfo;
numCells = collatedResults.numCells;
numRepeats = size(numCells,2); % number of runs per parameter combination, to gather stats
numPerturbations = size(numCells,1); % how many parameter combinations have been run altogether, not incl. repeats

meanNumCells = mean(numCells,2);

% parameters with experimental reference values
elasticityLeadSpeed = (meanNumCells(3) - meanNumCells(2))/meanNumCells(1)/(1.1 - 1/1.1);
elasticityFollowSpeed = (meanNumCells(5) - meanNumCells(4))/meanNumCells(1)/(1.1 - 1/1.1);
elasticityFiloLength = (meanNumCells(11) - meanNumCells(10))/meanNumCells(1)/(1.1 - 1/1.1);
elasticityFollowerFraction = (meanNumCells(13) - meanNumCells(12))/meanNumCells(1)/(1.1 - 1/1.1);

% parameters without experimental reference values
elasticityDiffus = (meanNumCells(15) - meanNumCells(14))/meanNumCells(1)/(10 - 1/10);
elasticityChi = (meanNumCells(17) - meanNumCells(16))/meanNumCells(1)/(10 - 1/10);
elasticityEatRate = (meanNumCells(19) - meanNumCells(18))/meanNumCells(1)/(10 - 1/10);
elasticityEatWidth = (meanNumCells(21) - meanNumCells(20))/meanNumCells(1)/(2 - 1/2);

labels = {'reference', 'leadSpeed', ['e = ' num2str(elasticityLeadSpeed)], 'followSpeed',['e = '  num2str(elasticityFollowSpeed)], ...
    'numFilopodia [3;2]', 'numFilopodia [2;1]', 'numFilopodia [1;1]', 'numFilopodia [3;3]', ...
    'filoLength', ['e = ' num2str(elasticityFiloLength)], 'followerFraction', ['e = ' num2str(elasticityFollowerFraction)], ...
    'diffusion', ['e = ' num2str(elasticityDiffus)], 'CA production', ['e = ' num2str(elasticityChi)], ...
    'eatRate', ['e = ' num2str(elasticityEatRate)], 'eatWidth', ['e = ' num2str(elasticityEatWidth)]};
figure;
boxplot(fliplr(numCells'),'labels',fliplr(labels),'orientation','horizontal')
xlabel('# of cells at end of simulation')
