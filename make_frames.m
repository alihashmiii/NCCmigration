% close all
framesFig = figure;

showColorbar = 1;
if experiment==0
    frames2show = [1, find(t_save<=12,1,'last'), find(t_save<=24,1,'last')] 
else
    frames2show = [find(t_save<=12,1,'last'), find(t_save<=24,1,'last')]
end

for frameCtr=1:length(frames2show)
    timeCtr = frames2show(frameCtr);
    subplot(length(frames2show),1,frameCtr)    
    make_plot(cells_save{timeCtr},cellsFollow_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},ca_save{timeCtr},filopodia_save{timeCtr},...
        numFilopodia,attach_save{timeCtr},cellRadius,filolength,sensingAccuracy,showColorbar,caCmap,0)
    if growingDomain==0
        xlim([0,300])
    end
    title(['time t = ',mat2str(t_save(timeCtr)),' hours with ' num2str(size(cells_save{timeCtr},2))...
        ' cells and ' num2str(min([size(cells_save{timeCtr},2) nnz(cellsFollow_save{timeCtr}==0)])) ' leaders.'])
    xlabel('x/\mum'), ylabel('y/\mum') % -- LJS
end

saveas(framesFig,['avi_mat/frames/frames3',saveInfo,'.fig'])
close(framesFig)