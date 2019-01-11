function [ dan ] = tunneling_solve(cells, xsave, ysave, param, dan_prev,cellRadius)
%sets lattice of DAN concentration to zero where cells are
%   L.J. Schumacher

% find points at which cells occupancy with the dan lattice (within a cell
% radius)
xcell = cells(1,:);
ycell = cells(2,:);

% xsave and ysave are the coordinates of the chemoattractant lattice
% points. dan has the same y-coordinates, but only makes up a 1/3 of the x
% range
xrange = xsave(xsave>=0&xsave<=max(xsave)/3)*param.initialDomainLength; % rescale the x-coordinates from stationary domain of unit+ length to units of microns-- LJS
occupancy_indcs = find_occupancy(xrange,ysave,xcell,ycell,cellRadius);

dan = dan_prev;
dan(occupancy_indcs) = 0;

end