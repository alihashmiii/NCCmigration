% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% defines initial conditions
% for use with chemotaxis_solve.m

function [u] = chemotaxis_pdeiv(npts, npde, t, x, y)
% u = ones(npts,npde);
% u = sin(pi.*x./1100);

load avi_mat/param
domainWidth = param(6);
height = param(9);
zero_bc = param(10);
insert = param(11);
growingDomain= param(5);
if insert==1
    disp('inserting')
    load avi_mat/ca_new
    load avi_mat/xlat_new
    load avi_mat/ylat_new
    %
    %     load avi_mat/xlat_save
    %     load avi_mat/ylat_save
    %     load avi_mat/ca_save
    %     load avi_mat/after_diffusion
    %     xpt = find(x==x(1),2,'first');
    %     new_x = x(1:xpt(2)-1)
    %     new_y = y(1:length(new_x):end)
    %     figure
    %     surf(xlat_new,ylat_new,ca_new')
    %     ca_new = interp2(xlat_new, ylat_new,ca_new',x,y);
    %     ca_new = ca_new';
    %     size(ca_new)
    %     x
    %     for i=1:length(x)
    %         u(i) = interp2(after_diffusion.x,after_diffusion.y,after_diffusion.ca',x(i),y(i));
    %         u(i) = interp2(xlat_new,ylat_new,ca_new',x(i),y(i));
    % %         u(i) = interp2(xlat_save{120},ylat_save{120},ca_save{120}',x(i),y(i));
    %     end
    % size(ca_new)
    % size(xlat_new)
    % size(ylat_new)
    % size(x)
    % size(y)
%     L_inf = param(1);
%     a = param(2);
%     domainWidth = param(6);
%     t_start = param(13);
%     tstep = param(12);
%     
%     out = domain_growth([],t-tstep,tstep,L_inf,a,domainWidth,[],t_start);
%     L = out.domainLength;
%     for i=1:length(x)
%         u(i) = interp2(xlat_new*domainWidth/L,ylat_new,ca_new',x(i),y(i));
%     end
        u = reshape(ca_new,[1,length(x)]);
%         x(1:length(xlat_new))
%         xlat_new*domainWidth/L
    %     y
    %     u = ca_new;
else
    if zero_bc == 1
        if growingDomain==1
            steepy = 0.1; % the higher steep is the steeper the gradient at the edge
            steepx = 0.15;
        else
            steepy = 0.03;
            steepx = 0.2;
        end
        z1 = 1./(1+exp(-steepy*x))-1/2;
        z2 = 1./(1+exp(steepy*(-domainWidth+x)))-1/2;
        z3 = 1./(1+exp(steepx*(-y)))-1/2;
        z4 = 1./(1+exp(steepx*(-height+y)))-1/2;
    else % no flux boundary condition -- LJS
        z1 = ones(size(x))/2;
        z2 = ones(size(x))/2;
        z3 = ones(size(y))/2;
        z4 = ones(size(y))/2;
    end
    u = min(min(min(z1,z2),z3),z4)/max(min(min(min(z1,z2),z3),z4));
end

% save_initial = input('save initial?:');
% if save_initial==1
%     xpts = x;
%     ypts = y;
%     plotxpts = xpts(1:find(xpts(2:end)<xpts(1:end-1),1,'first'));
%     plotypts = ypts((ypts(2:end)-ypts(1:end-1)>0));
%     plotypts(end+1) = ypts(end);
%
%     plotsol = reshape(u,[length(plotxpts),length(plotypts)]);
%     surf(plotxpts,plotypts,plotsol)
%     xlabel('x')
%     ylabel('y')
%     zlabel('Chemoattractant concentration')
%     saveas(gcf,['avi_mat/initial_conditions.fig'])
% end