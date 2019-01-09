ind_a=32000;
ind_b=40350;
%%------------------
outputVideo = VideoWriter(fullfile('images','decode_sep27_2d.avi'));
outputVideo.FrameRate = 5;
open(outputVideo);
for ii = ind_a:ind_b
   if exist(['images\result_' num2str(ii) '.png'], 'file')
    img = imread(['images\result_' num2str(ii) '.png']);
    writeVideo(outputVideo,img)
   end
end
close(outputVideo);
