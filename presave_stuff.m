%% save parameters and stuff %%
if isstruct(in)&&ismember('saveInfo',fields(in))
    saveInfo = in.saveInfo;
else
    if conversionType~=0
        if conversionType==4
            saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),'_steps_',num2str(numSteps(1)),'_',num2str(numSteps(2)),...
                '_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
        else
            saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),'_steps_',num2str(numSteps),...
                '_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
        end
    else
        saveInfo = [datestr(now,'yyyy_mm_dd-HH_MM'),'_foll_',num2str(followerFraction,2),'_convert_',mat2str(conversionType),...
            '_eatRate_',num2str(eatRate),'_diff_',num2str(diffus)];
    end
end
% don't overwrite existing file
if isempty(dir(['/mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'.mat'])) % check if this run hasn't been done, if previous sweeps have been aborted
    
    out.domainLengths = domainLengths;
    out.saveInfo = saveInfo;
    out.numTsteps = numTsteps;
    out.growingDomain = growingDomain;
    out.followerFraction = followerFraction;
    out.tstep = tstep;
    out.cellRadius = cellRadius;
    out.domainHeight = domainHeight;
    out.filolength = filolength;
    out.param = param; %param = [Linf, a, diffus, eatWidth, growingDomain, initialDomainLength, makeChemoattractant, chi, domainHeight, zeroBC, insert, tstep, t_start, eatRate, numSteps, numDirections];
    
    % these are parameters we might sweep
    out.leadSpeed = leadSpeed;
    out.followSpeed = followSpeed;
    out.numFilopodia = numFilopodia;
    
    out.growingDomain = growingDomain;
    out.followerFraction = followerFraction;
    out.divide_cells = divide_cells;
    out.experiment = experiment;
    
    save(['/mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'_running.mat'],'out')
    fprintf(['created results file at /mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'_running.mat \n'])
else
    fprintf(['error in creating results file: /mi/share/scratch/schumacher/Dropbox/DPhil/DysonModel/all_vers2/results/',saveInfo,'.mat already exists \n'])
end