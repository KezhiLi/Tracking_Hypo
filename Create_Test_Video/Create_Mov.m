function Create_Mov(masks, fname, fps, workingDir)

% fname = ['test_',date,'.avi' ];

% video rate
% fps = 25;

% workingDir = 'C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_17\Sample_Video\hdf5';

outputVideo = VideoWriter(fullfile(workingDir,fname));
outputVideo.FrameRate = fps;
open(outputVideo)

size_1 = size(masks,1);
size_2 = size(masks,2);
size_3 = size(masks,3);

bg = round(max(max(masks(:,:,1)))*1.1);

for ii = 1:size_3
   img1 = masks(:,:,ii);
   img1(img1<1) = bg; 
   %img = img1'; 
   img = img1; 
   writeVideo(outputVideo,img)
end

close(outputVideo)