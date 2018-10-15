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
numFilopodia = [3,1];  % the number of filopodia for lead cells and follower cells
filolength = cellRadius + 9*2;   % filopodial length (um) (measured from cell centre -- LJS). The average filopodial length found in experiment was 9mu, here I may be choosing a higher effective value to account for interfilopodial contact -- LJS
diffus = 0.1;%252e3;    % chemoattractant diffusivity (in (mu)^2/h?), for VEGF diffusing in the matrix this should probably be around 7e-11m^2/s = 252e3(mu)^2/h, for membrane bound VEGF unknown/near zero -- LJS
chi = 0.0001;                  % chemoattractant production term (usually 0.0001)
eatRate = 100;                      % chemoattractant consumption rate, usually 0.045 -- need to check this -- LJS
eatWidth = cellRadius;         % width of eating chemoattractant, equivalent to gaussian sigma
followerFraction = 6/8;        % proportion of cells that are followers (0<=follow_per<=1)

tstep = 5/2/60;                   % time step in hours

numRepeats = 20; % number of runs per parameter combination, to gather stats
input = cell(numRepeats,1);

dateString = '2013_10_10'; % for the folder name to save all results in, even if simulations actually run over multiple days or are repeated on another day

totalTime = tic;
for repCtr = 1:(numRepeats)
    input{repCtr}.time = 18; % duration of each run, simulation time in hours
    
    % run model for reference set
    if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done, if previous sweeps have been aborted
        input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
            '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
            '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
            '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
            '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
            '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
        CA6(input{repCtr},0);
    end
    
    % perturb parameters with experimental reference values
    perturbFactor = 1.25;
    
    for newLeadSpeed = [leadSpeed/perturbFactor, leadSpeed*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*leadSpeed_' num2str(newLeadSpeed) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.leadSpeed = newLeadSpeed;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(newLeadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'leadSpeed');
        end
    end
    
    for newFollowSpeed = [followSpeed/perturbFactor, followSpeed*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*followSpeed_' num2str(newFollowSpeed) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.followSpeed = newFollowSpeed;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(newFollowSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'followSpeed');
        end
    end
    
    for newnumFilopodia = [[2; 1], [6; 2]]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*nFiloLead_' num2str(newnumFilopodia(1)) '_nFiloFollow_' num2str(newnumFilopodia(2)) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.numFilopodia = newnumFilopodia';
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(newnumFilopodia(1)) '_nFiloFollow_' num2str(newnumFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'numFilopodia');
        end
    end
    
    perturbFactor = 1.5;
    
    for newFilolength = [filolength/perturbFactor, filolength*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*filolength_' num2str(newFilolength) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.filolength = newFilolength;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(newFilolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'filolength');
        end
    end
    
    for newFollowerFraction = [5/8, 7/8]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*followerFraction_' num2str(newFollowerFraction) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.followerFraction = newFollowerFraction;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(newFollowerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'followerFraction');
        end
    end
    
    % perturb parameters without experimental reference values
    perturbFactor = 252e4;
    
    for newDiffus = [diffus/perturbFactor, diffus*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*diffus_' num2str(newDiffus) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.diffus = newDiffus;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(newDiffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'diffus');
        end
    end
    
    perturbFactor = 10;
    
    for newChi = [chi/perturbFactor, chi*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*chi_' num2str(newChi) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.chi = newChi;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(newChi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'chi');
        end
    end
    
    for newEatRate = [eatRate/perturbFactor, eatRate*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*eatRate_' num2str(newEatRate) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.eatRate = newEatRate;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(newEatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'eatRate');
        end
    end
    
    perturbFactor = 1.6;
    
    for newEatWidth = [eatWidth/perturbFactor,  eatWidth*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*eatWidth_' num2str(newEatWidth) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.eatWidth = newEatWidth;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(newEatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(tstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'eatWidth');
        end
    end
    
    perturbFactor = 2;

    for newTstep = [tstep/perturbFactor, tstep*perturbFactor]
        if isempty(dir(['results/parameterSweeps/' dateString '-allResults/*tstep_' num2str(newTstep) '_*Run_' num2str(repCtr) '.mat'])) % check if this run hasn't been done
            input{repCtr}.tstep = newTstep;
            input{repCtr}.saveInfo = ['parameterSweeps/' dateString '-allResults/' datestr(now,'yyyy_mm_dd-HH_MM') ...
                '-leadSpeed_' num2str(leadSpeed) '_followSpeed_' num2str(followSpeed) ...
                '_nFiloLead_' num2str(numFilopodia(1)) '_nFiloFollow_' num2str(numFilopodia(2)) ...
                '_filolength_' num2str(filolength) '_diffus_' num2str(diffus) '_chi_' num2str(chi) ...
                '_eatRate_' num2str(eatRate) '_eatWidth_' num2str(eatWidth) '_followerFraction_' num2str(followerFraction) ...
                '_tstep_' num2str(newTstep) '_Run_' num2str(repCtr)];
            CA6(input{repCtr},0);
            input{repCtr} = rmfield(input{repCtr},'tstep');
        end
    end
end
toc(totalTime);