function [worm_body_x, worm_body_y] = shape2lin(worm_body, tt, tts)
% transform points in the body to a linear line
%
%
%
% 29/01/2015 Kezhi Li, MRC, Imperial 
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

worm_body_x =  round(interp1(tt,worm_body(:,1),tts,'linear'));
worm_body_y =  round(interp1(tt,worm_body(:,2),tts,'linear'));
        
worm_body_x(worm_body_x<1)=1;
worm_body_y(worm_body_y<1)=1;

