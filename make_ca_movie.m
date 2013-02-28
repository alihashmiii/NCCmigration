close all
clear mex
fig = figure('units','normalized','outerposition',[0 0 1 1],'visible','off');
comp = input('ca_movie with compression? (1 for yes): ');
if comp==1&&isunix==1
    aviobj = VideoWriter(['avi_mat/camovie/camovie_',save_info,'compressed.avi'],'Motion JPEG AVI');
elseif comp==1
    aviobj = VideoWriter(['avi_mat/camovie/camovie_',save_info,'compressed.avi'],'MPEG-4'); % "Compressed MPEG-4 file with H.264 encoding (Windows 7 systems only)" -- LJS
else
    aviobj = VideoWriter(['avi_mat/camovie/camovie_',save_info,'uncompressed.avi'],'Uncompressed AVI');
end
aviobj.FrameRate = 20;    % frames per second - fewer frames will make the movie slower
set(gcf,'Renderer','zbuffer')
open(aviobj);
for k=1:tsteps
    surf(xlat_save{k},ylat_save{k},ca_save{k}')
    zlim([0,1])
    xlim([0,max(domainLength)])
    ylim([0,height])
    
    %% insert title here
    if experiment>0
        title(['Cell invasion simulating experiment ',mat2str(experiment)])
    elseif follow_perc==0
        title('Cell invasion without followers')
    else
        title('Cell invasion with followers')
    end
    writeVideo(aviobj,getframe(fig));
end
close(fig)
close(aviobj);
