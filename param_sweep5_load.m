% parameter sweep for Dyson model CA6
% perturbed parameters around reference parameter set
% load results for plotting
% LJSchumacher 29.07.2013

clear all

% these parameters are not sweeped, but needed to calculate some reference
% values
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

numRepeats = 20; % number of runs per parameter combination, to gather stats
numPerturbations = 21; % how many parameter combinations have been run altogether, not incl. repeats
numCells = NaN(numPerturbations,numRepeats);
saveInfo = cell(numPerturbations,1);
meanDirectionality = NaN(numPerturbations,numRepeats,2);
meanSpeed = NaN(numPerturbations,numRepeats,2);
% to calculate the density profile of cells and chemoattractant along the x-direction
xBins = 0:2*cellRadius:49*2*cellRadius; % check the end length for time step run, or try loading this from the out-file somehow
cellDistributions = NaN(numPerturbations,numRepeats,2,length(xBins));
caDistribution = NaN(numPerturbations,numRepeats,length(xBins));

dateString = '2013_10_10'; % date of simulation that is to be loaded - needs to be adapted if ran on multiple days. Use wildcard (*) for time.
dateString2 = '2013_10_1*'; % if the simulations were repeated on a different date

for repCtr = 1:numRepeats
    for ptbCtr = 1:numPerturbations
        switch ptbCtr
            case 1
                % load results for reference set
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {2, 3}
                % perturb parameters with experimental reference values
                perturbFactor = 1.25;
                
                if ptbCtr == 2, newLeadSpeed = leadSpeed/perturbFactor; else newLeadSpeed = leadSpeed*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(newLeadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {4, 5}
                perturbFactor = 1.25;
                
                if ptbCtr ==4, newFollowSpeed = followSpeed/perturbFactor; else newFollowSpeed = followSpeed*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(newFollowSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {6, 7}
                if ptbCtr == 6, newnumFilopodia = [2; 1]; else newnumFilopodia = [6; 2]; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(newnumFilopodia(1)) '_nFiloFollow_' num2str(newnumFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {8, 9}
                perturbFactor = 1.5;
                
                if ptbCtr == 8, newFilolength = filolength/perturbFactor; else newFilolength = filolength*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(newFilolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {10, 11}
                perturbFactor = 1.5;
                
                if ptbCtr == 10, newFollowerFraction = 5/8; else newFollowerFraction = 7/8; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(newFollowerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {12, 13}
                % perturb parameters without experimental reference values
                perturbFactor = 252e4;
                
                if ptbCtr == 12, newDiffus = diffus/perturbFactor; else newDiffus =  diffus*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(newDiffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {14, 15}
                perturbFactor = 10;
                
                if ptbCtr == 14, newChi = chi/perturbFactor; else newChi = chi*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(newChi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {16, 17}
                perturbFactor = 10;
                
                if ptbCtr == 16, newEatRate = eatRate/perturbFactor; else newEatRate = eatRate*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(newEatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {18, 19}
                perturbFactor = 1.6;
                
                if ptbCtr == 18, newEatWidth = eatWidth/perturbFactor; else newEatWidth = eatWidth*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(newEatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(tstep) '_Run_' num2str(repCtr) '.mat']);
                
            case {20, 21}
                perturbFactor = 2;
                if ptbCtr == 20, newTstep = tstep/perturbFactor; else newTstep =  tstep*perturbFactor; end
                loadInfo = dir(['results/parameterSweeps/' dateString '-allResults/' dateString2 '*' ...
                    '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                    '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                    '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                    '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                    '_tstep_' num2str(newTstep) '_Run_' num2str(repCtr) '.mat']);
        end
        
        load(['results/parameterSweeps/' dateString '-allResults/' loadInfo.name], 'out');
        
        % extract cell numbers
        numberOfCells = size(out.cells_save{end},2);
        numCells(ptbCtr,repCtr) =  numberOfCells;
        
        % to calculate directionality and speed, first restructure cell
        % position data on per cell basis (t;x;y;follower?)
        cellPositions = cell(numberOfCells,1);
        for timeCtr = 1:2:size(out.cells_save,1) % only take every other time step to match experimental time resolution
            for cellCtr = 1:size(out.cells_save{timeCtr},2)
                cellPositions{cellCtr} = [cellPositions{cellCtr} [out.t_save(timeCtr + 1); out.cells_save{timeCtr}(:,cellCtr); out.cellsFollow{timeCtr}(cellCtr)]];
            end
        end
        
        % construct distribution of cells along x, for leaders and followers
        cellDistributions(ptbCtr,repCtr,1,:) = histc(out.cells_save{end}(1,out.cellsFollow{end}(1:numberOfCells)==0),xBins)./numberOfCells; % leaders
        cellDistributions(ptbCtr,repCtr,2,:) = histc(out.cells_save{end}(1,out.cellsFollow{end}(1:numberOfCells)==1),xBins)./numberOfCells; % leaders
        
        caDistribution(ptbCtr,repCtr,:) = sum(out.ca_save{end},2);
        
        % calculate directionality and speed, needs to be adapted if phenotype switching is enabled
        directionality = NaN(numberOfCells,2);
        effectiveSpeed = NaN(numberOfCells,2);
        for cellCtr = 1:numberOfCells
            totalPath = sum(sqrt((cellPositions{cellCtr}(2,2:end) - cellPositions{cellCtr}(2,1:end-1)).^2 ... % x^2
                + (cellPositions{cellCtr}(3,2:end) - cellPositions{cellCtr}(3,1:end-1)).^2)); % y^2
            straightPath = sqrt((cellPositions{cellCtr}(2,end) - cellPositions{cellCtr}(2,1)).^2 ... % x^2
                + (cellPositions{cellCtr}(3,end) - cellPositions{cellCtr}(3,1)).^2); % y^2
            directionality(cellCtr,:) = [cellPositions{cellCtr}(4,1) ...% follower?
                straightPath/totalPath];
            effectiveSpeed(cellCtr,:) = [cellPositions{cellCtr}(4,1) ...% follower?
                totalPath/(cellPositions{cellCtr}(1,end) -  cellPositions{cellCtr}(1,1))];
        end
        % when taking averages disregard cells with too few data whos
        % directionality of speed may be NaN
        meanDirectionality(ptbCtr,repCtr,1) = mean(directionality(directionality(:,1)==0 & ~isnan(directionality(:,2)),2)); % mean leader directionality
        meanDirectionality(ptbCtr,repCtr,2) = mean(directionality(directionality(:,1)==1 & ~isnan(directionality(:,2)),2)); % mean follower directionality
        meanSpeed(ptbCtr,repCtr,1) = mean(effectiveSpeed(effectiveSpeed(:,1)==0 & ~isnan(effectiveSpeed(:,2)),2)); % mean leader effectiveSpeed
        meanSpeed(ptbCtr,repCtr,2) = mean(effectiveSpeed(effectiveSpeed(:,1)==1 & ~isnan(effectiveSpeed(:,2)),2)); % mean follower effectiveSpeed
        
    end
    
    % only once, save parameter details
    if repCtr == 1, saveInfo{1} = loadInfo.name; end
    
end
save(['results/parameterSweeps/' dateString '-collatedResults'], 'numCells', 'saveInfo', 'meanDirectionality', 'meanSpeed','cellDistributions','caDistribution')
