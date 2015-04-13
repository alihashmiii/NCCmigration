
% load & plot example simulation runs
% L.J. Schumacher 28.05.14

clear all
close all
load results/2014_05_09-17_04_foll_1_convert_0_eatRate_1000_diff_0.1.mat
load_results

caCmap = load('cmap_blue2cyan.txt');

% modified from make_frames
showColorbar = 0;
frames2show = [1, find(t_save<=12,1,'last'), find(t_save<=24,1,'last')];
heights = [3, 3.7, 7.2];
for frameCtr=1:length(frames2show)
    timeCtr = frames2show(frameCtr);
    figure
%     subplot(length(frames2show),1,frameCtr)    
    make_plot(cells_save{timeCtr},cellsFollow_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},t_save(timeCtr-1),ca_save{timeCtr},filopodia_save{timeCtr},...
        numFilopodia,attach_save{timeCtr},cellRadius,filolength,sensingAccuracy,showColorbar,caCmap,0)
    
    title(['t = ',mat2str(t_save(timeCtr)),' hours'])
    xlabel('x/\mum'), ylabel('y/\mum') % -- LJS
    xLim = xlim;
    set(gca,'XTick',[0:200:xLim(2)])
%     xlim([0 1000])
    hcb = colorbar;
    set(hcb,'ytick',[0 1])
    set(get(hcb,'title'),'string','CA')
    % export figure
    exportOptions = struct('Format','eps2',...
        'Height',heights(frameCtr),...
        'Color','rgb',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    
filename = ['manuscripts/subpopulations/figures/schematicsFigE' num2str(frameCtr)];
    exportfig(gcf,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);
end