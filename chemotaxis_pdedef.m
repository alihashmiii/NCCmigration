% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% defines pde "res=0"
% for use with chemotaxis_solve.m

function [res] = chemotaxis_pdedef(npts, npde, t, x, y, u, ut, ux, uy, uxx, uxy, uyy)
global param cells % using global variables is much faster than saving & loading from disk -- LJS
% load avi_mat/cells
% load avi_mat/param
xcell = cells(1,:);
ycell = cells(2,:);

Linf = param(1);
a = param(2);
diffus = param(3);
eatWidth = param(4);
growingDomain = param(5);
initialDomainLength = param(6);
make_chemoattractant = param(7);
chi = param(8);
tstep = param(12);
t_start = param(13);
eatRate = param(14);

if (growingDomain==1)
    out = domain_growth([],t-tstep,tstep,Linf,a,initialDomainLength,[],t_start);
    L = out.domainLength;
    Ldiff = out.Ldiff;
else
    L = initialDomainLength;
    Ldiff =0;
end
xcell(xcell>0) = xcell(xcell>0)/L; % rescale the cells x-coordinates to stationary domain of unit length -- LJS
if make_chemoattractant~=1
    chi=0;
end
% the following equation should describe the chemoattractant evolution on
% the stationary domain [0,1]. Therefore, introduce factors of L-- LJS
eatTerm = zeros(size(x));
for ctr = 1:length(xcell)
    eatTerm = eatTerm + eatRate*u/(eatWidth^2*2*pi).*exp(-1/2/eatWidth^2.*((x - xcell(ctr)).^2*L^2 +(y - ycell(ctr)).^2));      
end
res = ut -(diffus.*(1./L^2.*uxx + uyy)...
          - eatTerm...
          + chi.*u.*(1-u));   % res = ut-f(u) means ut=f(u)
res(x>0) = res(x>0) + Ldiff/L*u(x>0);
% max(abs(res-res_check))
% x = linspace(0,1100,50);
% y = linspace(0,250,50);
% xcell = 500;
% ycell = 125;
% eatRate = 0.1;
% e = 100;
% z = eatRate.*exp(-e.*((x./max(x)-xcell./max(x)).^2+(y./max(y)-ycell/max(y)).^2));
