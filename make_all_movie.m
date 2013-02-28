% Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/10/09
% Edited on 06/01/10
%
% creates a movie of the cells moving with yellow followers and blue leaders
% on a background of chemoattractant
% saves to /scratch/dysonl/avi_mat

close all
clear mex
figure('units','normalized','outerposition',[0 0 1 0.5],'visible','off'); % fullscreen and not visible
comp = input('all_movie with compression? (1 for yes): ');
if comp==1&&isunix==1
    aviobj = VideoWriter(['avi_mat/allmovie/allmovie_',save_info,'compressed.avi'],'Motion JPEG AVI');
elseif comp==1
    aviobj = VideoWriter(['avi_mat/allmovie/allmovie_',save_info,'compressed.avi'],'MPEG-4'); % "Compressed MPEG-4 file with H.264 encoding (Windows 7 systems only)" -- LJS
else
    aviobj = VideoWriter(['avi_mat/allmovie/allmovie_',save_info,'uncompressed.avi'],'Uncompressed AVI');
end
aviobj.FrameRate = 15;    % frames per second - fewer frames will make the movie slower
set(gcf,'Renderer','zbuffer')
open(aviobj);
for k=1:tsteps
    disp(['step ',mat2str(k),' of ',mat2str(tsteps)])
    make_plot(cells_save{k},cellsFollow_save{k},xlat_save{k},ylat_save{k},ca_save{k},filopodia_save{k},num_filopodia,attach_save{k},cellRadius,0,barrier(k),experiment)
    
    axis image
    xlim([min(xlat_save{end}),1100])
    ylim([-50,120+50])
    T_t = get(gca,'TightInset');
    T = get(gca,'Position');
    set(gca,'position',[T_t(1)+0.05 T_t(2) 1-T_t(1)-T_t(3)-0.17 1-T_t(2)-T_t(4)]);
    
    %% Insert title here
    %     title([mat2str(follow_perc),' followers without growing domain and ',mat2str(tsteps),' time steps of ',mat2str(tstep),' hours at time = ', mat2str(t_save(k))])
    if experiment>0
        title(['Cell invasion simulating experiment ',mat2str(experiment)])
    elseif follow_perc==0
        title('Cell invasion without followers')
    else
        title('Cell invasion with followers')
    end
    writeVideo(aviobj,getframe(gca));
end
close(aviobj);

disp(['saved to avi_mat/allmovie/allmovie_',save_info,'.avi'])