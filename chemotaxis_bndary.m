% using d03ra to solve a pde in a rectangular domain
% defines boundary condition 
% for use with chemotaxis_solve.m
% Based on Louise Dyson D.Phil project
% modified by L.J. Schumacher

function [res] = chemotaxis_bndary(npts, npde, t, x, y, u, ut, ux, uy, nbpts, lbnd, res)
global param % we need to use more parameters that used by the d03ra syntax
% and using global variables is much faster than saving & loading from disk -- LJS

domainHeight = param.domainHeight;
zeroBC = param.zeroBC;

if isunix==1
    nbpts = int32(nbpts); % the documentation calls for this parameter to be int32, but commenting this out doesn't seem to make a difference -- LJS
end

if zeroBC==1
    res(lbnd,1) = u(lbnd,1);
else
    % No Flux Boundary conditions
    tolx = 1/sqrt(double(npts))/2 - x02aj(); % we need to account not only for machine precision, but also lattice discretisation -- LJS
    % at x = 0 or x = 1
    xBndIdcs = lbnd(abs(x(lbnd)) <= tolx | abs(x(lbnd) - 1) <= tolx);
    res(xBndIdcs,1) = ux(xBndIdcs,1);
    % at y = 0 or y = domainHeight
    toly = domainHeight/sqrt(double(npts))/2 - x02aj(); % we need to account not only for machine precision, but also lattice discretisation -- LJS
    yBndIdcs = lbnd(abs(y(lbnd)) <= toly | abs(y(lbnd) - domainHeight) <= toly);
    res(yBndIdcs,1) = uy(yBndIdcs,1);
end
end