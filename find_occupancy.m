function occupancy = find_occupancy(lat_x, lat_y, query_x, query_y, tol)
% find occupancy of lattice lat_x, lat_y by off-lattice query points within
% a tolerance (or size) tol
nx = length(lat_x);
ny = length(lat_y);
[X,Y] = meshgrid(lat_x, lat_y);
nX = numel(X);
nY = numel(Y);
nQ = length(query_x);
X = X(:);
Y = Y(:);
% broadcasting using tony's trick
occupancy = (X(:,ones(nQ,1)) - query_x(ones(nX,1),:)).^2 + ...
    (Y(:,ones(nQ,1)) - query_y(ones(nY,1),:)).^2 <= tol^2;
occupancy = any(occupancy,2); % don't mind who is occupying where
occupancy = reshape(occupancy,ny,nx);
end

