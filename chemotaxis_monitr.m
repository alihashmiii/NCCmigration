% Louise Dyson D.Phil project chemotaxis_solve, 14/10/09
% using d03ra to solve a pde in a rectangular domain
% checking functions etc.
% for use with chemotaxis_solve.m

function [ierr] = chemotaxis_monitr(a, t, b, c, tlast, d, ngpts, xpts, ypts, lsol, sol, ierr)
global param plotsol xsave ysave % using global variables is much faster than saving & loading from disk -- LJS

if isunix==1
    ngpts = int32(ngpts);
end

if tlast
    level = 1;
    npts=ngpts(level);
    if isunix==1
        ipsol=int32(lsol(level));
    else
        ipsol=lsol(level);
    end
    %load avi_mat/param
    growingDomain = param(5);
    tstep = param(12);
    
    if growingDomain==1
        L_inf = param(1);
        a = param(2);
        initialDomainLength = param(6);
        t_start = param(13);
        
        out = domain_growth([],t-tstep,tstep,L_inf,a,initialDomainLength,[],t_start);
        L = out.domainLength;

        k = sum(ngpts(1:level-1));
        xpts = xpts(k+1:k+npts);
        ypts = ypts(k+1:k+npts);
%         xsave = xpts(1:find(xpts(2:end)-xpts(1:end-1)<0,1,'first')).*L/initialDomainLength;
        xsave = xpts(1:find(xpts(2:end)-xpts(1:end-1)<0,1,'first'));
        xsave(xsave>=0) = xsave(xsave>=0).*L/initialDomainLength;
        ysave = ypts(1:length(xsave):end);
    else
        k = sum(ngpts(1:level-1));
        xpts = xpts(k+1:k+npts);
        ypts = ypts(k+1:k+npts);
        xsave = xpts(1:find(xpts(2:end)-xpts(1:end-1)<0,1,'first'));
        ysave = ypts(1:length(xsave):end);
    end
    
    plotsol = zeros(1,npts);
    for idks=1:npts
        plotsol(idks) = sol(ipsol+idks);
    end;
    plotsol = reshape(plotsol,[length(xsave),length(ysave)]);
%     surf(xsave,ysave,plotsol)
%     pause
%     save avi_mat/plotsol plotsol
%     save avi_mat/xsave xsave
%     save avi_mat/ysave ysave
end
end