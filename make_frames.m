% close all
figure

edge = 1;
if experiment==0
    frames = [1, find(t_save<=12,1,'last'), find(t_save<=24,1,'last')] 
else
    frames = [find(t_save<=12,1,'last'), find(t_save<=24,1,'last')]
end
lengths = cellfun(@length,cells_save);
lengths_change = find(lengths(2:end)-lengths(1:end-1)>2);
lengths(lengths_change);
lengths(lengths_change+1);
for j=1:length(frames)
    k = frames(j);
    subplot(length(frames),1,j)    
    make_plot(cells_save{k},cellsFollow_save{k},xlat_save{k},ylat_save{k},ca_save{k},filopodia_save{k},numFilopodia,attach_save{k},cellRadius,edge,experiment)
    if growingDomain==0
        xlim([0,300])
    end
    title(['Cell migration at time t = ',mat2str(t_save(k)),' hours with ' num2str(size(cells_save{k},2)) ' cells and ' num2str(min([size(cells_save{k},2) nnz(cellsFollow_save{k}==0)])) ' leaders.'])
    xlabel('x (\mum)'), ylabel('y (\mum)') % -- LJS
end

saveas(gcf,['avi_mat/frames/frames3',saveInfo,'.fig'])