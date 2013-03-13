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
domainWidth = param(6);

if isunix==1
    nbpts = int32(nbpts);
end

if zero_bc==1
    % zero boundary conditions
    for idks = 1:nbpts
        j = lbnd(idks);
        res(j,1) = u(j,1);
    end
else
    % No Flux Boundary conditions
    tol = x02aj();
%     tol = 20*tol;
    for idks = 1:nbpts
        j = lbnd(idks);
        % at x = 0
        if (abs(x(j)) <= tol)
            res(j,1) = ux(j,1);
            % at y = 0
        elseif (abs(y(j)) <= tol)
            res(j,1) = uy(j,1);
            % at y = domainHeight
        elseif(abs(y(j)-domainHeight) <= tol)
            res(j,1) = uy(j,1);
%             at x = domainLength
        else
            res(j,1) = ux(j,1);
        end
    end
end