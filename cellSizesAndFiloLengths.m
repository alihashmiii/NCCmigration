% draw a figure to plot cell size and filopodia length distributions
load('cellSizesInclFilopodia.mat')
load('filopodiaInclBranches.mat')

filopodiaAll = [filopodia; filopodia1; filopodia2];

%% plot results
subplot(1,2,1)
h1 = histogram(cellSizes);
h1.BinWidth = 5;
set(gca,'Ytick',0:5)
hold on
plot(mean(cellSizes)*[1 1], [0 6], 'r-')
plot(median(cellSizes)*[1 1], [0 6], 'r--')
xlabel('s/\mum')
ylabel('n')
ylim([0 5])
title({'A: cell size (including filopodia)';''})

subplot(1,2,2)
h2 = histogram(filopodiaAll);
h2.BinWidth = 2.5;
hold on
plot(mean(filopodiaAll)*[1 1], [0 70], 'r-')
plot(median(filopodiaAll)*[1 1], [0 70], 'r--')
xlabel('l/\mum')
ylabel('n')
title({'B: filopodial length'; '(measured from cell body/branching point)';''})

%% export figure

exportOptions = struct('Format','eps2',...
    'Width','18.0',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

filename = ['../../Thesis/figures/cellSizeAndFiloMeasurements'];
pos = get(gcf,'Position');
pos(4) = 1/2*pos(3); % adjust height to fraction of width
set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);