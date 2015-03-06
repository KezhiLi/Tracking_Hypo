
workingDir = 'S:\Kezhi\Video_Sample\02_2015';

filename = 'S:\GeckoVideo\Compressed\20150216\CaptureTest_90pc_Ch1_16022015_174636.hdf5';

if ~exist('infoFile', 'var')
    infoFile =  h5info(filename);
end
 
chunkSize = [2048 2048 1];

writerObj = VideoWriter('Sample_Video\try_2702.avi');
open(writerObj);

 
image = h5read(filename, '/mask', [1,1,1], chunkSize);
[~, rect] = imcrop(image);
figure,
for frame_number = 1:2:2000;
    image = h5read(filename, '/mask', [1,1,frame_number], chunkSize);
    img_rect = imcrop(image, rect);
    img_rect(img_rect<1)=195;
    imshow(img_rect);
    
    frame = getframe;
    writeVideo(writerObj,frame);
end
close(writerObj);