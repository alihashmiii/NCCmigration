function plotMuSigmaErr(repeats,quantityOI,titleString)
% plots the mean of a quantity of interest on one y-scale, standard
% deviation and error in mean on another y-scale

% initialise empty variables
meanQOI = [];
stdQOI = [];

% calculate mean, std and error in mean
for repCtr = 1:length(repeats)
    meanQOI = [meanQOI, mean(quantityOI(1:repeats(repCtr)))];
    stdQOI = [stdQOI, std(quantityOI(1:repeats(repCtr)))];
end
errMeanQOI = stdQOI./sqrt(repeats);

% plot the results
[AX,H1,H2] = plotyy(repeats,meanQOI,repeats,stdQOI);
hold(AX(2))
plot(AX(2),repeats,errMeanQOI,'r-')
set(H1,'LineStyle','--')
set(H2,'LineStyle','-.')

% annotate the plot
set(get(AX(1),'Ylabel'),'String','mean') 
set(get(AX(2),'Ylabel'),'String','deviation and error')
xlabel('# of repeats')
legend('\mu','\sigma','error in mean','Location','NorthWest')
title(titleString)

end