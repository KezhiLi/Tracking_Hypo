function save_crt_fra(filename, k)
% function to save current frame to the objective file
% Input: filename: a text word corresponding the objective file name
%        k: current frame index
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 24/02/2015 

 frame = getframe(1);
 im = frame2im(frame);
 [imind,cm] = rgb2ind(im,256);
 if k == 3;
     imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
 else
     imwrite(imind,cm,filename,'gif','WriteMode','append');
 end