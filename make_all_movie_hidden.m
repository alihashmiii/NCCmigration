% saves images and compiles movie externally - LJS
% saves to /avi_mat

% close all
  
movieFig = figure('units','points','outerposition',[0 0 1600 430],'position',[1 1 1599 363],...
    'PaperPositionMode','manual','PaperUnits','points','PaperPosition',[0 0 1600 430],'visible','off');
filePath = ['avi_mat/allmovie/',saveInfo,'/'];
mkdir(filePath)
frameRate = 30;    % frames per second - fewer frames will make the movie slower
skip = round(5/60/param.tstep); % skip more time steps for faster saving
frameCtr = 1;
% % set(gcf,'Renderer','zbuffer')

caCmap = load('cmap_blue2cyan.txt');

minx = min([min(xlat_save{end}) 0]);
for timeCtr=1:skip:numTsteps
    disp(['step ',mat2str(timeCtr),' of ',mat2str(numTsteps)])
    make_plot(cells_save{timeCtr},cellsFollow_save{timeCtr},xlat_save{timeCtr},...
        ylat_save{timeCtr},ca_save{timeCtr},filopodia_save{timeCtr},numFilopodia,attach_save{timeCtr},cellRadius,filolength,sensingAccuracy,1,caCmap,0);
    
    xlim([min(xlat_save{timeCtr}), 1000])
    
    %% Insert title here
    title(['Cell migration at time = ' num2str(t_save(timeCtr),'%2.2f') ' hours with '  ...
        num2str(size(cells_save{timeCtr},2)) ' cells and ' num2str(min([size(cells_save{timeCtr},2) ...
        nnz(cellsFollow_save{timeCtr}==0)])) ' leaders.'],'FontSize',23)
    % resize figure to tighter margins
    T_t = get(gca,'TightInset');
    T = get(gca,'Position');
    set(gca,'position',[T_t(1)+0.04 T_t(2) 1-T_t(1)-T_t(3)-0.1 1-T_t(2)-T_t(4)]);
    
%     axis image
    print(movieFig,'-r0','-dtiff',[filePath,int2str(frameCtr),'.tiff'])
    frameCtr = frameCtr + 1;
end
close(movieFig)
%% using external video compiling
disp('compiling the movie now');

%help on this to be found here:
%https://libav.org/avconv.html

system(['avconv -i ', filePath, '%d.tiff -r ', int2str(frameRate), ' -qscale 10 -y ',filePath,saveInfo,'.mp4']);
% To have a constant quality (but a variable bitrate), use the option 
% ?-qscale n? when ?n? is between 1 (excellent quality) and 31 (worst quality). 
% delete the individual frames
system(['rm -f ' filePath '*.tiff']);
disp(['saved to ',filePath])
