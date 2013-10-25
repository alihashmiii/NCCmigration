% Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/10/09
% Edited on 06/01/10
%
% creates a movie of the cells moving with yellow followers and blue leaders
% on a background of chemoattractant
% saves to /avi_mat

% close all

plotHandle = figure('units','points','outerposition',[0 0 1600 430],'position',[1 1 1599 363],...
    'PaperPositionMode','manual','PaperUnits','points','PaperPosition',[0 0 1600 430],'visible','off');
filePath = ['avi_mat/allmovie/',saveInfo,'/'];
mkdir(filePath)
frameRate = 30;    % frames per second - fewer frames will make the movie slower
skip = 2; % skip more time steps for faster saving
frameCtr = 1;
% % set(gcf,'Renderer','zbuffer')
minx = min([min(xlat_save{end}) 0]);
for timeCtr=1:skip:numTsteps
    disp(['step ',mat2str(timeCtr),' of ',mat2str(numTsteps)])
    make_plot(cells_save{timeCtr},cellsFollow_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},ca_save{timeCtr},filopodia_save{timeCtr},numFilopodia,attach_save{timeCtr},cellRadius,0,experiment);
    
%     xlim([minx,1100]) %this can cause unpleasant jittering in movies
%     ylim([-50,120+50])
    T_t = get(gca,'TightInset');
    T = get(gca,'Position');
    set(gca,'position',[T_t(1)+0.05 T_t(2) 1-T_t(1)-T_t(3)-0.17 1-T_t(2)-T_t(4)]);
    
    %% Insert title here
    title(['Cell migration at time = ' num2str(t_save(timeCtr),'%2.2f') ' hours with '  num2str(size(cells_save{timeCtr},2)) ' cells and ' num2str(min([size(cells_save{k},2) nnz(cellsFollow_save{k}==0)])) ' leaders.'])
    end
%     axis image
    print(plotHandle,'-dpng',[filePath,int2str(frameCtr),'.png'])
    frameCtr = frameCtr + 1;
end
% using external video compiling as recommended by Jochen Kursawe (JK)
% compiling the movie with external shell command... -- JK
disp('compiling the movie now');

%help on this to be found here (ffmpeg and avconv work nearly the same way,
%the latter is the new version of the first):
%http://ffmpeg.org/trac/ffmpeg/wiki/x264EncodingGuide -- JK

system(['avconv -r ', int2str(frameRate), ' -i ', filePath,...
    '%d.png -c:v libx264 -crf 35 -loglevel quiet -y ',filePath,saveInfo,'.mp4']);

% delete the individual frames
system(['rm -f ' filePath '*.png']);
disp(['saved to ',filePath])