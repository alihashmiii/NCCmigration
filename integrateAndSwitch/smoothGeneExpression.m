function [smoothTime, smooth] = smoothGeneExpression(time,samples,N,method,slope)
% smooth data with cubic splines or akima interpolation, with optional zero
% slope at last timepoint

smoothTime = linspace(time(1),time(end),N);
smooth = NaN(size(samples,2),N);

if slope==0
    % specify zero slope at right end by duplicating last value at machine
    % epsilon before last time point (no zero slope at left end as conditions
    % discontinuous, i.e., an impulse)
    time = [time(1:end-1); time(end) - eps(time(end)); time(end)];
    samples = samples([1:end end],:);
end

if strcmp(method,'cubic')
    smooth = spline(time,samples',smoothTime);
elseif strcmp(method,'akima')
    for ctr = 1:size(samples,2)
        smooth(ctr,:) = akimai(time,samples(:,ctr),smoothTime);
    end
else
    error(['method ' method ' not recognised']);
end

end

