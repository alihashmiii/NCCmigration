% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% defines initial conditions
% for use with chemotaxis_solve.m

function [u] = chemotaxis_pdeiv(npts, npde, t, x, y)

global param % using global variables is much faster than saving & loading from disk -- LJS
% load avi_mat/param
initialDomainLength = param.initialDomainLength;
domainHeight = param.domainHeight;
zeroBC = param.zeroBC;
insert = param.insert;
growingDomain= param.growingDomain;

if insert==1&&~isempty(param.ca_new) %% will have to check these when doing tissue transplantations -- LJS
    if length(x)~=length(param.ca_new(:))
        disp('interpolating initial conditions for refined grid')
        ylat = linspace(0,domainHeight,size(param.ca_new,2));
        xlat = linspace(0,1,size(param.ca_new,1));
        [Xlat,Ylat] = meshgrid(xlat,ylat);
        u = interp2(Xlat,Ylat,param.ca_new',x,y)';
    else
        disp('inserting initial conditions for chemotaxis solver')
        u = reshape(param.ca_new,1,length(x));
    end
else
    if zeroBC == 1
        if growingDomain==1
            steepy = 0.1; % the higher steep is the steeper the gradient at the edge
            steepx = 0.15;
        else
            steepy = 0.03;
            steepx = 0.2;
        end
        z1 = 1./(1+exp(-steepy*x*initialDomainLength))-1/2;
        z2 = 1./(1+exp(steepy*(-initialDomainLength + x*initialDomainLength)))-1/2;
        z3 = 1./(1+exp(steepx*(-y)))-1/2;
        z4 = 1./(1+exp(steepx*(-domainHeight+y)))-1/2;
    else % no flux boundary condition -- LJS
        z1 = ones(size(x))/2;
        z2 = ones(size(x))/2;
        z3 = ones(size(y))/2;
        z4 = ones(size(y))/2;
    end
    u = min(min(min(z1,z2),z3),z4)/max(min(min(min(z1,z2),z3),z4));
end