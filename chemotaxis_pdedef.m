% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% defines pde "res=0"
% for use with chemotaxis_solve.m

function [res] = chemotaxis_pdedef(npts, npde, t, x, y, u, ut, ux, uy, uxx, uxy, uyy)
load avi_mat/cells
load avi_mat/param
xcell = cells(1,:);
ycell = cells(2,:);

Linf = param(1);
a = param(2);
diffus = param(3);
e = param(4);
growingDomain = param(5);
domainWidth = param(6);
make_chemoattractant = param(7);
chi = param(8);
tstep = param(12);
t_start = param(13);
eatRate = param(14);

if (growingDomain==1)
    out = domain_growth([],t-tstep,tstep,Linf,a,domainWidth,[],t_start);
    L = out.domainLength;
    Ldiff = out.Ldiff;
else
    L = domainWidth;
    Ldiff =0;
end
xcell(xcell>0) = xcell(xcell>0)/L; % rescale the cells x-coordinates to stationary domain of unit length -- LJS
if make_chemoattractant~=1
    chi=0;
end
% the following equation should describe the chemoattractant evolution on
% the stationary domain [0,1]. Therefore, introduce factors of L-- LJS
eatTerm = eatRate*u*ones(size(xcell))/(e*sqrt(pi)).*exp(-1/e^2.*((x*ones(size(xcell))-ones(size(x))*xcell).^2*L^2 +(y*ones(size(ycell))-ones(size(y))*ycell).^2));      
res = ut -(diffus.*(1./L^2.*uxx + uyy)...
          -sum(eatTerm,2)...
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

%% to make CA_eat.fig - USE make_initial_conditions.m
% load avi_mat/multiplier
% if multiplier==0
%     figure
%     x = x(1:50);
%     y = y(1:50:end);
%     multiplier = zeros(length(y),length(x));
%     for i=1:length(x)
%         for j=1:length(y)
%             multiplier(j,i) = sum(eatRate*exp(-e.*((x(i)-xcell).^2+(y(j)-ycell).^2)));
%         end
%     end
%     x = x*L/domainWidth;
%     saved.x = x;
%     saved.y =y;
%     saved.eatRate = multiplier;
%     save avi_mat/saved saved
%     % figure('visible','on');
%     surf(x,y,multiplier)
%     xlim([0,1100])
%     ylim([0,120])
% %    set(gca,'ZTick',0:0.02:0.1,'ZTickLabel',{'0.00'; '0.02';'0.04';'0.06';'0.08';'0.10';}); 
% 
%     saveas(gcf,'CA_eat.fig')
%     close all
%     multiplier = 1;
%     save avi_mat/multiplier multiplier
%     disp('made CA_eat')
%     % pause
%     % return
% end

% %%% old formulation %%%
% for i = 1:npts
%     res_old(i,1) = ut(i,1) -( diff.*(1/L^2*uxx(i,1)+uyy(i,1))...
%         -sum(u(i,1)*exp(-d.*((x(i)/L*ones(size(xcell))-xcell/L).^2+(y(i)/L*ones(size(ycell))-ycell/L).^2)))...
%         -L_diff/L*u(i,1)); % <-- non-apical growth  +L_diff/L*u(i,1); <-- apical growth
% %     res_one(i,1) = sum(u(i,1)*exp(-d.*((x(i)/L*ones(size(xcell))-xcell/L).^2+(y(i)/L*ones(size(ycell))-ycell/L).^2)));
% end
