% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% defines boundary condition G = 0 at x=xmin, xmax, y=ymin, ymax by
% rewriting the solution at those points
% for use with chemotaxis_solve.m

function [res] = chemotaxis_bndary(npts, npde, t, x, y, u, ut, ux, uy, nbpts, lbnd, res)
global param % using global variables is much faster than saving & loading from disk -- LJS
% load avi_mat/param
% load avi_mat/current_domainLength
domainHeight = param(9);
zero_bc = param(10);

if isunix==1
    nbpts = int32(nbpts);
end

if zero_bc==1
    res(lbnd,1) = u(lbnd,1);
else
    % No Flux Boundary conditions
    tol = x02aj();
    % at x = 0 or x = 1
    xBndIdcs = lbnd(abs(x(lbnd)) <= tol | abs(x(lbnd) - 1) <= tol);
    res(xBndIdcs,1) = ux(xBndIdcs,1);
    % at y = 0 or y = domainHeight
    yBndIdcs = lbnd(abs(y(lbnd)) <= tol | abs(y(lbnd) - domainHeight) <= tol);
    res(yBndIdcs,1) = uy(yBndIdcs,1);
end
end