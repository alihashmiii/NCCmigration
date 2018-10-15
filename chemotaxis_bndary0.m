% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% defines boundary condition G = 0 at x=xmin, xmax, y=ymin, ymax by
% rewriting the solution at those points
% for use with chemotaxis_solve.m

function [res] = chemotaxis_bndary0(npts, npde, t, x, y, u, ut, ux, uy, nbpts, lbnd, res)
global param % using global variables is much faster than saving & loading from disk -- LJS

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
    % at x = 1
    xLBndIdcs = lbnd(abs(x(lbnd) - 1) <= tolx);
    res(xLBndIdcs,1) = ux(xLBndIdcs,1);
    % at x = 0
    x0BndIdcs = lbnd(abs(x(lbnd)) <= tolx);
    res(x0BndIdcs,1) = ux(x0BndIdcs,1);
    % at y = 0 or y = domainHeight
    toly = domainHeight/sqrt(double(npts))/2 - x02aj(); % we need to account not only for machine precision, but also lattice discretisation -- LJS
    yBndIdcs = lbnd(abs(y(lbnd)) <= toly | abs(y(lbnd) - domainHeight) <= toly);
    res(yBndIdcs,1) = uy(yBndIdcs,1);
end
end