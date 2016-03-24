% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width

clear
close all

makeMoves = 1;
makeFrames = 1;
makeAllMovie = 1;

experiment = 41;
precision = 2;

diffusivities = [1 10 100];
slowSpeeds = [5 10 20];

fileName = 'exp41_slowEntryDilute';

for cntGdn = {'toward', 'parallel'}
    result.contactGuidance = char(cntGdn);
    for diffus = diffusivities
        result.diffus = diffus;
        result.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        for slowSpeed = slowSpeeds
            result.slowSpeed = slowSpeed;
            result.loadInfo = [fileName '_D_' num2str(result.diffus) ...
                '_sensingAcc_' num2str(result.sensingAccuracy,precision) ...
                '_slowSpeed_' num2str(result.slowSpeed) '_contactGuidance_' result.contactGuidance];
            load(['results/' fileName '/' result.loadInfo '.mat'])
            load_results
            make_frames
            saveInfo = out.saveInfo(23:end); % had the folder name repeated in the saveInfo
            make_all_movie_hidden
        end
    end
end