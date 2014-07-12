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

Linf = param.Linf;
a = param.a;
diffus = param.diffus;
eatWidth = param.eatWidth;
growingDomain = param.growingDomain;
initialDomainLength = param.initialDomainLength;
makeChemoattractant = param.makeChemoattractant;
chi = param.chi;
domainHeight = param.domainHeight;
tstep = param.tstep;
t_start = param.t_start;
eatRate = param.eatRate;

experiment = param.experiment;
transplantTime = param.transplantTime;
transplantXLocation = param.transplantXLocation;
secondaryChi = param.secondaryChi;

if (growingDomain==1)
    [~, L, Ldiff] = domain_growth([],t-tstep,tstep,Linf,a,initialDomainLength,t_start);
else
    L = initialDomainLength;
    Ldiff =0;
end
xcell(xcell>0) = xcell(xcell>0)/L; % rescale the cells x-coordinates to stationary domain of unit length -- LJS
if makeChemoattractant~=1
    chi=0;
end
% % the following equation should describe the chemoattractant evolution on
% % the stationary domain [0,1]. Therefore, introduce factors of L -- LJS
% eatTerm = zeros(size(x));
% for ctr = 1:length(xcell)
%     eatTerm = eatTerm + eatRate*u/(eatWidth^2*2*pi).*exp(-1/2/eatWidth^2.*((x - xcell(ctr)).^2*L^2 +(y - ycell(ctr)).^2));      
% end
% % this following vectorisation, using "Tony's trick" for indexing, is
% faster -- LJS
bigEat = eatRate*u(:,ones(length(xcell),1))/(eatWidth^2*2*pi).*exp(-1/2/eatWidth^2.*(...
    (x(:,ones(length(xcell),1)) - xcell(ones(npts,1),:)).^2*L^2 + (y(:,ones(length(xcell),1)) - ycell(ones(npts,1),:)).^2)); %tony's trick
eatTerm = sum(bigEat,2);
res = ut -(diffus.*(1./L^2.*uxx + uyy)...
    - eatTerm...
    + chi.*u.*(1-u));   % res = ut-f(u) means ut=f(u)
if (experiment==11||experiment==12||experiment==13)&&t>=transplantTime
    if (experiment==12||experiment==13)
        indcs =  (x>=transplantXLocation/L)&(x<=(transplantXLocation/L + 1/8))&(y<=domainHeight/2);
    elseif experiment==11
        indcs =  (x>=transplantXLocation/L)&(x<=(transplantXLocation/L + 1/8))&(y<=domainHeight/20);
    end
    res(indcs) = res(indcs) - secondaryChi.*u(indcs).*(1 - u(indcs)); % increased CA production
end
res(x>0) = res(x>0) + Ldiff/L*u(x>0); % dilution through tissue growth -- LJS