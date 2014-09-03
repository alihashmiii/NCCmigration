%% A function to find the least squares regression values for Linf and a
%% and then plot the resulting domain growth
% function [] = make_domain_plot()

% end_time = input('Input length of simulation (usually 24hrs): ');
% tstep = input('Input step size (usually 0.05hrs): ') ;
close all
clear all

end_time = 30;
param.tstep = 1/60;                   % time step in hours

numTsteps = floor((end_time)/param.tstep);     % number of time steps
% Linf = 1100;                            % end domain length
% a = 0.05;                                % how fast the domain grows (?)
% initialDomainLength = 300;                % initial width of the domain with growth (um)
close all
t = 0:param.tstep:end_time;

%% Least Squares Regression to find Linf and a
lengths_data    % gives lengths and time
time = time;
a_values=0.055:0.005:0.12; %
Linf_values = 750:10:1000;
L0_values = 250:5:350;
t_start_values = -30:20;
sum_error_squared = NaN(length(a_values),length(t_start_values),length(L0_values),length(Linf_values));
for l = 1:length(L0_values)
    L0 = L0_values(l);
    for m = 1:length(Linf_values)
        Linf = Linf_values(m);
        for i = 1:length(a_values)
            a = a_values(i);
            for j = 1:length(t_start_values)
                t_start = t_start_values(j);
                [~, estimatedLengths, ~] = domain_growth([],time,0,Linf,a,L0,t_start);
                sum_error_squared(i,j,l,m) = sum((estimatedLengths' - lengths).^2);
            end
        end
    end
end
tmp = find(sum_error_squared==min(min(min(min(sum_error_squared)))));
[ind(1),ind(2),ind(3),ind(4)] = ind2sub(size(sum_error_squared),tmp);

figure
subplot(2,2,1)
plot(a_values,sum_error_squared(:,ind(2),ind(3),ind(4)))
title('\alpha')
subplot(2,2,2)
plot(t_start_values,sum_error_squared(ind(1),:,ind(3),ind(4)))
title('t_s')
subplot(2,2,3)
plot(Linf_values,reshape(sum_error_squared(ind(1),ind(2),ind(3),:),length(Linf_values),1))
title('L_{\infty}')
subplot(2,2,4)
plot(L0_values,reshape(sum_error_squared(ind(1),ind(2),:,ind(4)),length(L0_values),1))
title('L_0')

figure
a = a_values(ind(1))
t_start = t_start_values(ind(2))
Linf = Linf_values(ind(4))
L0 = L0_values(ind(3))
 %% domain growth %%
[~, domainLengths, ~] = domain_growth([],t,0,Linf,a,L0,t_start);
hold on
plot(t,domainLengths,'MarkerSize',500)
% lengths_data
plot(time,lengths,'r.')
hold off

xlabel('time (hours)')
ylabel('domain length ({\mu}m)')
legend('fitted curve','experimental data')