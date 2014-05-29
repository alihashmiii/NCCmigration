% load & plot example simulation runs
% L.J. Schumacher 28.05.14

load results/2014_05_09-17_04_foll_1_convert_0_eatRate_1000_diff_0.1.mat
load_results

caCmap = load('cmap_blue2cyan.txt');

make_frames

open(['avi_mat/frames/frames3',saveInfo,'.fig'])

%% export figure
exportOptions = struct('Format','eps2',...
    'Width','14.4',...
    'Color','rgb',...
    'Resolution',300,...
    'FontMode','fixed',...
    'FontSize',10,...
    'LineWidth',2);

pos = get(gcf,'Position');
% pos(4) = 1/2*pos(3); % adjust height to fraction of width
set(gcf,'PaperUnits','centimeters','Position',pos,'color','none');
filename = ['manuscripts/subpopulations/figures/schematicsFigE'];
exportfig(gcf,[filename '.eps'],exportOptions);
system(['epstopdf ' filename '.eps']);