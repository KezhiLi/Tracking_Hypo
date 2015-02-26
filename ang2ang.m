function [ang_hypo,ang_chan, mid_start, speed, omg,len_short] = ang2ang(ang, omg,jj, speed,len_short)
% calculate the hypothesis angles based on previous angles
%
% Input: angle: an M *1 vector
%        var_ang: a scalar to control the variance of the angle change
% Output: ang_hypo: an M * 1 vector, the hypothesis angle created
%         ang_chan: an M * 1 vector, the angle change, comparing to the input change
%         mid_start: the index of middle starting point in pt_ske
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 14/01/2014

% calculate the length of angle vector
len_ang = length(ang);

if len_ang<1
    len_ang
end

% Find the middle two points, set them as mid_start and mid_end
mid_start = round(len_ang/2);
mid_end = mid_start+1;

% % These two are designed functions for variances of tha angle change
% % parameter: 0.8, 7, etc are set accroding the len_ang =17 or 18, and can be changed
% distr_ind_1half = 1./(mid_end-0.1-[1:mid_start]);
% distr_ind_2half = 1./(len_ang-mid_end+3-[1:(len_ang-mid_end+1)]);
% 
% % some adjustment about the distribution, to make head and tail not so
% % vibrant
% %distr = exp([fliplr(distr_ind_1half),distr_ind_2half]')-1;
% distr = [fliplr(distr_ind_1half),distr_ind_2half]'.*[fliplr(distr_ind_1half),distr_ind_2half]';
% distr(1) = distr(1)/2;
% distr = flipud(distr);

distr =  (exp(-0.9*(len_ang-(1:len_ang))))';
distr(end) = distr(end)*0.6;



% calculate the change of angle and output angle vector
% % 
%     if jj<24
%         ang_hypo_temp =ang;
%         if jj<5
%             ang_hypo_temp(1:jj) = ang_hypo_temp(1:jj)+0.4; %0.8
%         elseif jj<9&&jj>4
%             ang_hypo_temp(1:(jj-4)) = ang_hypo_temp(1:(jj-4))-0.4;
%         elseif jj<13&&jj>8
%             ang_hypo_temp((end-(jj-9)):end) = ang_hypo_temp((end-(jj-9)):end)+0.4;
%         elseif jj<17&&jj>12
%             ang_hypo_temp((end-(jj-13)):end) = ang_hypo_temp((end-(jj-13)):end)-0.4;
%         elseif jj == 17;
%             ang_hypo_temp(1:4) = ang_hypo_temp(1:4)+[0.8,0.6,0.4,0.2]';
%         elseif jj==18;
%             ang_hypo_temp(1:4) = ang_hypo_temp(1:4)-[0.8,0.6,0.4,0.2]';
%         elseif jj==19;
%             ang_hypo_temp(end-3:end) = ang_hypo_temp(end-3:end)+[0.2,0.4,0.6,0.8]'; %0.3,0.6,0.9,1.2
%         elseif jj==20;       
%             ang_hypo_temp(end-3:end) = ang_hypo_temp(end-3:end)-[0.2,0.4,0.6,0.8]';
%         elseif jj==21
%             ang_hypo_temp(1) = ang_hypo_temp(1)-1.5*sign(rand(1));
%         else
%             ang_hypo_temp(end) = ang_hypo_temp(end)-1.5*sign(rand(1));
%         end
%     else
%         ang_chan = var_ang*distr.*randn(len_ang,1);
%         ang_hypo_temp = ang + ang_chan;
%     end


    if jj<23
        ang_hypo_temp =ang;
        if jj<4
            ang_hypo_temp(1:jj) = ang_hypo_temp(1:jj)+0.4; %0.8
        elseif jj<7&&jj>3
            ang_hypo_temp(1:(jj-3)) = ang_hypo_temp(1:(jj-3))-0.4;
        elseif jj<10&&jj>6
            ang_hypo_temp((end-(jj-7)):end) = ang_hypo_temp((end-(jj-7)):end)+0.4;
            omg = omg + 0.2 + (jj-7)*0.07;
        elseif jj<13&&jj>9
            ang_hypo_temp((end-(jj-10)):end) = ang_hypo_temp((end-(jj-10)):end)-0.4;
            omg = omg - 0.2 - (jj-10)*0.07;
        elseif jj == 13;
            ang_hypo_temp(1:4) = ang_hypo_temp(1:4)+[0.6,0.45,0.3,0.15]';
        elseif jj==14;
            ang_hypo_temp(1:4) = ang_hypo_temp(1:4)-[0.6,0.45,0.3,0.15]';
        elseif jj==15;
            ang_hypo_temp(end-3:end) = ang_hypo_temp(end-3:end)+[0.15,0.3,0.45,0.6]'; %0.3,0.6,0.9,1.2
            omg = omg + 0.4;
        elseif jj==16;       
            ang_hypo_temp(end-3:end) = ang_hypo_temp(end-3:end)-[0.15,0.3,0.45,0.6]';
            omg = omg - 0.4;
        elseif jj==17
            ang_hypo_temp(1:2) = ang_hypo_temp(1:2)+[1.0,0.6]';
        elseif jj==18
            ang_hypo_temp(1:2) = ang_hypo_temp(1:2)-[1.0,0.6]';
        elseif jj==19
            ang_hypo_temp(end-1:end) = ang_hypo_temp(end-1:end)+[0.6,1.0]';
            omg = omg + 0.7;
        elseif jj==20
            ang_hypo_temp(end-1:end) = ang_hypo_temp(end-1:end)-[0.6,1.0]';
            omg = omg - 0.7;
        elseif jj==21
            ang_hypo_temp(1) = ang_hypo_temp(1)+sign(randn(1))*1.4;
        else
            ang_hypo_temp(end) = ang_hypo_temp(end)+sign(randn(1))*1.4;
        end
    elseif jj<41;
        omg = omg+randn(1)*0.2;
        ang_chan = omg * distr;
        ang_hypo_temp = ang + ang_chan;
    elseif jj<47;
        ang_hypo_temp = ang;
    elseif jj == 47;
        ang_hypo_temp = mid_pt_chg(ang,mid_start,1);
    elseif jj == 48;
        ang_hypo_temp = mid_pt_chg(ang,mid_start,-1);
    else
        ang_hypo_temp = ang;
    end

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


ang_chan = ang_hypo-ang;


% ii=1:16;xx=1./(17-ii);yy=1./(20-ii); figure, plot(1:32,[(flipud(yy'))',xx],'ro')