% simulate neural crest cell migration on wider comain with vegf production
% only in the middle, but cells get inserted at full width

clear
close all

makeMoves = 1;
makeFrames = 1;
makeAllMovie = 1;

experiment = 39;
precision = 2;

% setting the initial CA outside the migratory path to zero in case of a wider domain


eatRates = [100 1000];
diffusivities = [1 10 100];
speeds = [20 40];

for eatRate = eatRates
    result.eatRate = eatRate;
    for diffus = diffusivities
        result.diffus = diffus;
        result.sensingAccuracy = 0.1/sqrt(diffus/0.1); % sens acc scales with 1/sqrt(diffus)
        for speed = speeds
            result.leadSpeed = speed;
            result.followSpeed = speed;
            result.loadInfo = ['exp39_widerDomainStripe_eat_' num2str(result.eatRate) ...
                '_D_' num2str(result.diffus) '_sensingAcc_' num2str(result.sensingAccuracy,precision) ...
                '_speed_' num2str(result.leadSpeed)];
            load(['results/exp39_widerDomainStripe/' result.loadInfo '.mat'])
            load_results
            saveInfo = out.saveInfo(25:end); % had the folder name repeated in the saveInfo
            make_frames
            make_all_movie_hidden
        end
    end
end