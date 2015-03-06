function mov = save_crt_fra(filename, k, fps)
% function to save current frame to the objective file
% Input: filename: a text word corresponding the objective file name
%        k: current frame index
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 24/02/2015 
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

 frame = getframe(1);
 im = frame2im(frame);
 [imind,cm] = rgb2ind(im,256);
 if k == 3;
     imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
 else
     imwrite(imind,cm,filename,'gif','WriteMode','append','delaytime',1/fps);
 end
 
 mov = im2frame(im);