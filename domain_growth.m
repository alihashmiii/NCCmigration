%% This function takes in the current time, position of the cells and
%% parameters relating to the domain growth
%% speed (including 'width', the initial length of the domain) etc.
%% The output is the next position of the cells, end of domain and Ldiff

function [cellsX, r_next, Ldiff] = domain_growth(cellsX,t,tstep,Linf,a,initialDomainLength,t_s)

Linf = Linf/initialDomainLength;

r = (Linf*exp(Linf*a*(t-t_s))./(Linf-1+exp(Linf*a*(t-t_s))) ...
    + 1-Linf*exp(Linf*a*(-t_s))/(Linf-1+exp(Linf*a*(-t_s))))*initialDomainLength;
t = t + tstep;
r_next = (Linf*exp(Linf*a*(t-t_s))./(Linf-1+exp(Linf*a*(t-t_s))) ...
    + 1-Linf*exp(Linf*a*(-t_s))/(Linf-1+exp(Linf*a*(-t_s))))*initialDomainLength;
Ldiff = Linf^2*a*exp(Linf*a*(t-t_s))*initialDomainLength*(Linf-1)/(Linf-1+exp(Linf*a*(t-t_s))).^2;

if ~isempty(cellsX)
    cellsX(cellsX>=0) = cellsX(cellsX>=0)*r_next/r;
end