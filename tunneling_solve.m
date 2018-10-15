function [ dan ] = tunneling_solve(cells, xsave, ysave, param, dan_prev,t,cellRadius)
%sets lattice of DAN concentration to zero where cells are
%   quick and dirty for now

if (param.growingDomain==1)
    [~, L, ~] = domain_growth([],t-param.tstep,param.tstep,param.Linf,param.a,param.initialDomainLength,param.t_s);
else
    L = param.initialDomainLength;
end

% find points at which cells occupancy with the dan lattice (within a cell
% radius)
xcell = cells(1,:);
ycell = cells(2,:);

% xsave and ysave are the coordinates of the chemoattractant lattice
% points. dan has the same y-coordinates, but only makes up a 1/3 of the x
% range
xrange = xsave(xsave>=0&xsave<=max(xsave)/3)*param.initialDomainLength; % rescale the x-coordinates from stationary domain of unit+ length to units of mu-- LJS
occupancy_indcs = find_occupancy(xrange,ysave,xcell,ycell,cellRadius);

dan = dan_prev;
dan(occupancy_indcs) = 0;

end