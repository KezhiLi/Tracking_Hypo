function [ang_hypo, mid_start, omg] = ang2ang_1st_narrow(ang, ang_head_i,ang_tail_i ,omg)
% function used to generate angle hypothesis for specific model
%
% Input:   ang: an M*1 vector, the original angle 
%          ang_head_i: a scalar, a parameter regarding head angle variance   
%          ang_tail_i: a scalar, a parameter regarding tail angle variance
%          omg: a scalar, the first derivative of angle, obtained from last iteration 
% Output:  ang_hypo: an M*1 vector, the hypothesis angle
%          mid_start: a scalar, the index of the 1st point of the middle segment
%          omg: a scalar, the current omega, which is the the first derivative of angle
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 10/04/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.


% calculate the length of angle vector
len_ang = length(ang);

if len_ang<1
    len_ang
end

% Find the middle two points, set them as mid_start and mid_end
mid_start = round(len_ang/2);
mid_end = mid_start+1;



% These two are designed functions for variances of tha angle change
distr =  (exp(-0.5*(len_ang-(1:len_ang))))';   % -0.9
distr(end) = distr(end)*0.6;

ang_chan = omg * distr;
ang_hypo_temp = ang + ang_chan;

if rand(1)>0.2  % 80% probability to change angles very near head/tail with large magnitues
    ang_head_chg = ang_head_i * [0.8,1.6]' * rand(1);
    ang_tail_chg = ang_tail_i * [1.6, 0.8]' * rand(1);
    ang_hypo_temp(1:2) = ang_hypo_temp(1:2) + ang_tail_chg;
    ang_hypo_temp(end-1:end) = ang_hypo_temp(end-1:end) + ang_head_chg;
else  % 20% probability to change angles near head/tail with large magnitues
    ang_head_chg = ang_head_i * [0.4,0.9,1.4,1.9]' * rand(1);
    ang_tail_chg = ang_tail_i * [1.9,1.4, 0.9,0.4]' * rand(1);    
    ang_hypo_temp(1:4) = ang_hypo_temp(1:4) + ang_tail_chg;
    ang_hypo_temp(end-3:end) = ang_hypo_temp(end-3:end) + ang_head_chg;    
end


omg = omg + ang_head_i * 0.2;


ang_hypo = zeros(size(ang_hypo_temp));
ang_hypo(mid_start) = ang_hypo_temp(mid_start);

    

% set the max magnitude of angle variance bound, be careful about the case when angles
% larger tha pi or smaller than -pi
thre = 1.4;
for ii = mid_start:-1:2;
    if abs(ang_hypo_temp(ii-1)-ang_hypo(ii))<thre||abs(abs(ang_hypo_temp(ii-1)-ang_hypo(ii))-2*pi)<thre;
        ang_hypo(ii-1) = ang_hypo_temp(ii-1);
    elseif  abs(abs(ang_hypo_temp(ii-1)-ang_hypo(ii))-2*pi)<=pi;
        angg1 = ang_hypo(ii)-sign(ang_hypo_temp(ii-1)-ang_hypo(ii))*thre;
        flag1 = abs(angg1)-pi >0;
        ang_hypo(ii-1) = (-1)^flag1*(2*pi*flag1+(-1)^flag1*abs(angg1))*sign(angg1);
    else
        ang_hypo(ii-1) = ang_hypo(ii)+sign(ang_hypo_temp(ii-1)-ang_hypo(ii))*thre;
    end
end
for ii = mid_start:len_ang-1;
    if abs(ang_hypo_temp(ii+1)-ang_hypo(ii))<thre||abs(abs(ang_hypo_temp(ii+1)-ang_hypo(ii))-2*pi)<thre;
        ang_hypo(ii+1) = ang_hypo_temp(ii+1);
    elseif  abs(abs(ang_hypo_temp(ii+1)-ang_hypo(ii))-2*pi)<=pi;
        angg2 = ang_hypo(ii)-sign(ang_hypo_temp(ii+1)-ang_hypo(ii))*thre;
        flag2 = abs(angg2)-pi >0;
        ang_hypo(ii+1) = (-1)^flag2*(2*pi*flag2+(-1)^flag2*abs(angg2))*sign(angg2);
    else
        ang_hypo(ii+1) = ang_hypo(ii)+sign(ang_hypo_temp(ii+1)-ang_hypo(ii))*thre;
    end
end 












