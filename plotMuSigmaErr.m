function plotMuSigmaErr(repeats,quantityOI,titleString,nResample)
% plots the mean of a quantity of interest on one y-scale, standard
% deviation and error in mean on another y-scale
% L.J. Schumacher 19.08.2013

% initialise empty variables
meanQOI = NaN(nResample,length(repeats));
stdQOI = NaN(nResample,length(repeats));

% ensure quantityOI is a row vector
quantityOI = quantityOI(:)';

% randomly resample order of quantityOI, if wanted
if nResample >1
    for reSampleCtr = 2:nResample
        quantityOI = [quantityOI; randsample(quantityOI(1,:),length(quantityOI(1,:)))];
    end
end

% calculate mean, std and error in mean
for repCtr = 1:length(repeats)
    meanQOI(:,repCtr) = mean(quantityOI(:,1:repeats(repCtr)),2);
    stdQOI(:,repCtr) = std(quantityOI(:,1:repeats(repCtr)),0,2);
end
errMeanQOI = stdQOI./sqrt(repmat(repeats,nResample,1));

% plot the results
[AX,H1,H2] = plotyy(repeats,meanQOI,repeats,stdQOI);
hold(AX(2)) % to plot more graphs on the same axis
H3 = plot(AX(2),repeats,errMeanQOI,'r-');
set(H1,'LineStyle','--','Color','b')
set(H2,'LineStyle','-.','Color',[0 0.5 0])

% annotate the plot
set(get(AX(1),'Ylabel'),'String','mean') 
set(get(AX(2),'Ylabel'),'String','deviation and error')
xlabel('# of repeats')
legend([H1(1), H2(1), H3(1)], {'\mu','\sigma','error in mean'},'Location','NorthWest')
title(titleString)

% make sure the second axes has range including zero
yLim2 = ylim(AX(2));
ylim(AX(2),[0 yLim2(2)])
set(AX(2),'YTick',linspace(0, yLim2(2), 3))
end