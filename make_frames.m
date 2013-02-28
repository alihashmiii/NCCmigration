% close all
figure

edge = 1;
if experiment==0
%     frames = [1, find(t_save==12), find(t_save==24)]
    frames = [1, find(t_save==12), find(t_save==24)] %, find(t_save==36), length(t_save)-1]
else
    frames = [find(t_save==12), find(t_save==24)]
end
lengths = cellfun(@length,cells_save);
lengths_change = find(lengths(2:end)-lengths(1:end-1)>2);
lengths(lengths_change)
lengths(lengths_change+1)
for j=1:length(frames)
    k = frames(j);
    subplot(length(frames),1,j)    
    if (experiment==4)||(experiment==5)
        make_plot(cells_save{k},cellsFollow_save{k},xlat_save{k},ylat_save{k},ca_save{k},filopodia_save{k},num_filopodia,attach_save{k},cellRadius,edge,barrier(k),experiment)
    else
        make_plot(cells_save{k},cellsFollow_save{k},xlat_save{k},ylat_save{k},ca_save{k},filopodia_save{k},num_filopodia,attach_save{k},cellRadius,edge,[],experiment)
    end
    if growingDomain==0
        xlim([0,300])
    end
    title(['Cell invasion at time t = ',mat2str(t_save(k)),' hours'])
    
end

% one colorbar for everything
% if experiment==6
%     L = get(colorbar,'YLim');
%     M = colorbar;
%     set(M,'YTick',L,'YTickLabel',[{0},1]);
%     set(M, 'Position', [.8314 .164 .0381 .71]);
%     set(M,'YTick',L,'YTickLabel',[{0},{1}]);
%     
%     h = findobj(gcf,'Type','axes','Tag','');
%     for i=1:length(h)
%         pos=get(h(i), 'Position');
%         set(h(i), 'Position', [pos(1) pos(2) .6626 pos(4)])
%     end
% end

saveas(gcf,['avi_mat/frames/frames3',save_info,'.fig'])