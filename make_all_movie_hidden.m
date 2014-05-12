% Louise Dyson D.Phil project CA program describing the migration of cranial neural crest cells, 22/10/09
% Edited on 06/01/10
%
% creates a movie of the cells moving with yellow followers and blue leaders
% on a background of chemoattractant
% saves to /avi_mat

% close all
  
movieFig = figure('units','points','outerposition',[0 0 1600 430],'position',[1 1 1599 363],...
    'PaperPositionMode','manual','PaperUnits','points','PaperPosition',[0 0 1600 430],'visible','off');
filePath = ['avi_mat/allmovie/',saveInfo,'/'];
mkdir(filePath)
frameRate = 30;    % frames per second - fewer frames will make the movie slower
skip = round(5/60/tstep); % skip more time steps for faster saving
frameCtr = 1;
% % set(gcf,'Renderer','zbuffer')

minx = min([min(xlat_save{end}) 0]);
for timeCtr=1:skip:numTsteps
    disp(['step ',mat2str(timeCtr),' of ',mat2str(numTsteps)])
    make_plot(cells_save{timeCtr},cellsFollow_save{timeCtr},xlat_save{timeCtr},ylat_save{timeCtr},ca_save{timeCtr},filopodia_save{timeCtr},numFilopodia,attach_save{timeCtr},cellRadius,filolength,sensingAccuracy,1,caCmap,0);
    
%     xlim([minx,1100]) %this can cause unpleasant jittering in movies
%     ylim([-50,120+50])
    T_t = get(gca,'TightInset');
    T = get(gca,'Position');
    set(gca,'position',[T_t(1)+0.05 T_t(2) 1-T_t(1)-T_t(3)-0.17 1-T_t(2)-T_t(4)]);
    
    %% Insert title here
    title(['Cell migration at time = ' num2str(t_save(timeCtr),'%2.2f') ' hours with '  ...
        num2str(size(cells_save{timeCtr},2)) ' cells and ' num2str(min([size(cells_save{timeCtr},2) ...
        nnz(cellsFollow_save{timeCtr}==0)])) ' leaders.'],'FontSize',23)
    
%     axis image
    print(movieFig,'-dtiff',[filePath,int2str(frameCtr),'.tiff'])
    frameCtr = frameCtr + 1;
end
% using external video compiling as recommended by Jochen Kursawe (JK)
% compiling the movie with external shell command... -- JK
disp('compiling the movie now');

%help on this to be found here (ffmpeg and avconv work nearly the same way,
%the latter is the new version of the first):
%http://ffmpeg.org/trac/ffmpeg/wiki/x264EncodingGuide -- JK

system(['avconv -r ', int2str(frameRate), ' -i ', filePath, '%d.tiff -c:v libx264 -crf 23 -loglevel quiet -y ',filePath,saveInfo,'.mp4']);
%Choose a CRF value
%The range of the quantizer scale is 0-51: where 0 is lossless, 23 is default, and 51 is worst possible. A lower value is a higher quality and a subjectively sane range is 18-28. Consider 18 to be visually lossless or nearly so: it should look the same or nearly the same as the input but it isn't technically lossless.
%The range is exponential, so increasing the CRF value +6 is roughly half the bitrate while -6 is roughly twice the bitrate. General usage is to choose the highest CRF value that still provides an acceptable quality. If the output looks good, then try a higher value and if it looks bad then choose a lower value.
%Note: The CRF quantizer scale mentioned on this page only applies to 8-bit x264 (10-bit x264 quantizer scale is 0-63). You can see what you are using with x264 --help listed under Output bit depth. 8-bit is more common among distributors.

% delete the individual frames
% system(['rm -f ' filePath '*.tiff']);
disp(['saved to ',filePath])
close(movieFig)