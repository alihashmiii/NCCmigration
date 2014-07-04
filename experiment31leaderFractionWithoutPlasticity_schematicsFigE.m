
% load & plot example simulation runs
% L.J. Schumacher 28.05.14

clear all
close all
load results/2014_05_09-17_04_foll_1_convert_0_eatRate_1000_diff_0.1.mat
load_results

caCmap = load('cmap_blue2cyan.txt');

% modified from make_frames
showColorbar = 1;
frames2show = [1, find(t_save<=12,1,'last'), find(t_save<=24,1,'last')];
for frameCtr=1:length(frames2show)
    timeCtr = frames2show(frameCtr);
    subplot(length(frames2show),1,frameCtr)    
    make_plot(cells_save{timeCtr},cellsFollow_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},ca_save{timeCtr},filopodia_save{timeCtr},...
        numFilopodia,attach_save{timeCtr},cellRadius,filolength,sensingAccuracy,showColorbar,caCmap,0)
    
    title(['time t = ',mat2str(t_save(timeCtr)),' hours with ' num2str(size(cells_save{timeCtr},2))...
        ' cells and ' num2str(min([size(cells_save{timeCtr},2) nnz(cellsFollow_save{timeCtr}==0)])) ' leaders.'])
    xlabel('x/\mum'), ylabel('y/\mum') % -- LJS
    xlim([0 1000])
end
    %% export figure
    exportOptions = struct('Format','eps2',...
        'Width',num2str(14.4),...
        'Color','rgb',...
        'Resolution',300,...
        'FontMode','fixed',...
        'FontSize',10,...
        'LineWidth',2);
    
filename = ['manuscripts/subpopulations/figures/schematicsFigE'];
    exportfig(gcf,[filename '.eps'],exportOptions);
    system(['epstopdf ' filename '.eps']);