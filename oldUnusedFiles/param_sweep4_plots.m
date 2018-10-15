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
numFilopodia = [6,2];  % the number of filopodia for lead cells and follower cells
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
diffus = 1;%252e3;    % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
chi = 0.0001;                  % chemoattractant production term (usually 0.0001)
eatRate = 1;                      % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
eatWidth = cellRadius;         % width of eating chemoattractant, equivalent to gaussian sigma
followerFraction = 0.7;        % proportion of cells that are followers (0<=follow_per<=1)

dateString = '2013_08_02'; % date of simulation that is to be loaded - needs to be adapted if ran on multiple days. Use wildcard (*) for time.
    
collatedResults = load(['results/parameterSweeps/' dateString '-reference100repeats']);
saveInfo = collatedResults.saveInfo;
numCells = collatedResults.numCells;
numRepeats = size(numCells,2); % number of runs per parameter combination, to gather stats
numPerturbations = size(numCells,1); % how many parameter combinations have been run altogether, not incl. repeats
meanDirectionality = collatedResults.meanDirectionality; % mean directionality for leaders and followers
meanSpeed = collatedResults.meanSpeed; % mean speed for leaders and followers
xBins = 0:2*cellRadius:685.4; % check the end length for time step run, or try loading this from the out-file somehow
cellDistributions = squeeze(collatedResults.cellDistributions);

%% analyse effects of parameter perturbation on cell number
repeats = 10:numRepeats;
nResample = 1;

figCellNumber = figure;

plotMuSigmaErr(repeats,numCells,'number of cells at end of simulation',nResample)

%% analyse effects of parameter perturbation on leader directionality
meanDirectionality = squeeze(meanDirectionality);

figLeaderDirectionality = figure;

plotMuSigmaErr(repeats,meanDirectionality(:,1),'leader directionality',nResample)

%% analyse effects of parameter perturbation on follower directionality

figFollowerDirectionality = figure;

plotMuSigmaErr(repeats,meanDirectionality(:,2),'follower directionality',nResample)

%% analyse effects of parameter perturbation on leader speed
meanSpeed = squeeze(meanSpeed);

figLeaderSpeed = figure;

plotMuSigmaErr(repeats,meanSpeed(:,1),'effective leader speed',nResample)

%% analyse effects of parameter perturbation on follower speed

figFollowerSpeed = figure;

plotMuSigmaErr(repeats,meanSpeed(:,2),'effective follower speed',nResample)

%% migration profile along x

figxProfile = figure;
stairs(xBins,mean(squeeze(cellDistributions(:,1,:))),'k-')
hold on
stairs(xBins,mean(squeeze(cellDistributions(:,2,:))),'k--')
legend('leaders','followers')
xlabel('x/\mum')
ylabel('P')
%% for exporting as pdf
% % set(gcf, 'PaperPositionMode', 'auto'); % set size to as seen on screen
% % print(gcf, '-r0','../../TransferOfStatus/ThesisProposal/figures/velocityCorrelation.eps', '-depsc2');
% % or ...
% set(figCellNumber,'PaperUnits','centimeters');
% set(figLeaderDirectionality,'PaperUnits','centimeters');
% set(figFollowerDirectionality,'PaperUnits','centimeters');
% set(figLeaderSpeed,'PaperUnits','centimeters');
% set(figFollowerSpeed,'PaperUnits','centimeters');
% set(figxProfile,'PaperUnits','centimeters');
% 
% exportOptions = struct('Format','eps2',...
%     'Width','10.0',...
%     'Color','cmyk',...
%     'Resolution',300,...
%     'FontMode','fixed',...
%     'FontSize',8,...
%     'LineWidth',1);
% 
% filename = 'results/tmp/statsCellNumber4.eps';
% exportfig(figCellNumber,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsLeaderDirectionality4.eps';
% exportfig(figLeaderDirectionality,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsFollowerDirectionality4.eps';
% exportfig(figFollowerDirectionality,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsLeaderSpeed4.eps';
% exportfig(figLeaderSpeed,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/statsFollowerSpeed4.eps';
% exportfig(figFollowerSpeed,filename,exportOptions);
% system(['epstopdf ' filename]);
% 
% filename = 'results/tmp/migrationProfile.eps';
% exportfig(figxProfile,filename,exportOptions);
% system(['epstopdf ' filename]);