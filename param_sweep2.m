% parameter sweep for Dyson model CA6
% perturb parameters around reference parameter set
% save results for later plotting
% LJSchumacher 29.06.2013

clear all
close all

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

input.time = 18; % duration of each run, simulation time in hours
numRepeats = 20; % number of runs per parameter combination, to gather stats

dateString = '2013_07_30'; % for the folder name to save all results in, even if simulations actually run over multiple days or are repeated on another day

totalTime = tic;
for repCtr = 1:numRepeats
    
    % run model for reference set
    input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
        '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
        '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
        '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
        '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
        '_Run_' num2str(repCtr)];
    CA6(input,0);

    % perturb parameters with experimental reference values
    perturbFactor = 1.1;
    
    for newLeadSpeed = [leadSpeed/perturbFactor, leadSpeed*perturbFactor]
        input.leadSpeed = newLeadSpeed;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(newLeadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'leadSpeed');
    
    for newFollowSpeed = [followSpeed/perturbFactor, followSpeed*perturbFactor]
        input.followSpeed = newFollowSpeed;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(newFollowSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'followSpeed');
    
    for newnumFilopodia = [[3; 2], [4; 2], [5; 2], [6; 2]]
        input.numFilopodia = newnumFilopodia';
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(newnumFilopodia(1)) '_nFiloFollow_' num2str(newnumFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'numFilopodia');
    
    for newFilolength = [filolength/perturbFactor, filolength*perturbFactor]
        input.filolength = newFilolength;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(newFilolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'filolength');
    
    for newFollowerFraction = [followerFraction/perturbFactor, followerFraction*perturbFactor]
        input.followerFraction = newFollowerFraction;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(newFollowerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'followerFraction');
    
    % perturb parameters without experimental reference values
    perturbFactor = 10;
    
    for newDiffus = [diffus/perturbFactor, diffus*perturbFactor]
        input.diffus = newDiffus;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(newDiffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'diffus');
    
    for newChi = [chi/perturbFactor, chi*perturbFactor]
        input.chi = newChi;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(newChi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'chi');
    
    for newEatRate = [eatRate/perturbFactor, eatRate*perturbFactor]
        input.eatRate = newEatRate;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(newEatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'eatRate');
    
    perturbFactor = 1.1;
    for newEatWidth = [eatWidth/perturbFactor, eatWidth*perturbFactor]
        input.eatWidth = newEatWidth;
        input.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(newEatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_Run_' num2str(repCtr)];
        CA6(input,0);
    end
    input = rmfield(input,'eatWidth');
    
end
toc(totalTime);