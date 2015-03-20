function [ang_hypo,mid_start,jj_4] = ang2ang_2nd(ang, jj, scenario)
% 
% 
% 
% 
% 
% 
% default scenario = [4,5,7,2,12];


% calculate the length of angle vector
len_ang = length(ang);

if len_ang<1
    len_ang
end

num_pt = 6;

% Find the middle two points, set them as mid_start and mid_end
mid_start = round(len_ang/2);
mid_end = mid_start+1;

jj_1 = scenario(1)*num_pt+1;
jj_2 = jj_1 + scenario(2);
jj_3 = jj_2 + scenario(3);
jj_4 = jj_3 + scenario(4);
jj_5 = jj_4 + scenario(5);

ang_hypo_temp = ang;
    
var_case1 = 0.2; 
var_case2 = 0.2;
var_case3 = 0.3;
var_case4 = 0.2;
var_case5 = 0.2;

if jj < jj_1
    % case 1
    jj_case1_pt = floor((jj-1)/scenario(1)); %  1 ~ num_pt
    jj_case1_idx =  mod(jj,scenario(1));       %  1 ~ scenario(1) 
    if jj_case1_pt == 1 || jj_case1_pt == 2
        ang_hypo_temp(jj_case1_pt) = ang_hypo_temp(jj_case1_pt) + randn(1) * var_case1;
    else 
        idx_case1 = len_ang-(num_pt-jj_case1_pt);
        ang_hypo_temp(idx_case1) = ang_hypo_temp(idx_case1) + randn(1) * var_case1;
    end
elseif jj < jj_2
    % case 2
    % generate case2 index between 2 ~ len_ang-2
    case2_idx = floor(rand(1)*(len_ang-3)+2); 
    ang_hypo_temp(case2_idx: case2_idx+1) = ang_hypo_temp(case2_idx: case2_idx+1) + sign(randn(1))*[var_case2; -var_case2];
elseif jj < jj_3
    % case 3
    % generate case3 index between 3 ~ len_ang-5
    case3_idx = floor(rand(1)*(len_ang-7)+3);
    % generate case3 index between 2 ~ 4
    case3_interval = floor(rand(1)*3)+2;
    ang_hypo_temp(case3_idx: case3_idx+case3_interval) = ang_hypo_temp(case3_idx: case3_idx+case3_interval) ...
        + sign(randn(1))*[var_case3; zeros(case3_interval-1,1); -var_case3];
elseif jj < jj_4
    % case 4
    if jj == jj_4-2
        ang_hypo_temp = mid_pt_chg(ang_hypo_temp,mid_start,1,var_case4);
    else % jj == jj_4-1
        ang_hypo_temp = mid_pt_chg(ang_hypo_temp,mid_start,-1,var_case4);
    end
else
    % case 5
    % case5_pt = 0 or 1
    case5_pt = floor((jj - jj_4)/(scenario(5)/2));
    if case5_pt == 0 % adjust tail
         ang_hypo_temp(1:2) = ang_hypo_temp(1:2) + randn(1) * var_case5 * [1.5 ; 0.8];
    else  % adjust head
         ang_hypo_temp(end-1:end) = ang_hypo_temp(end-1:end) + randn(1) * var_case5 * [1.5 ; 0.8];
    end
end

ang_hypo = zeros(size(ang_hypo_temp));
ang_hypo(mid_start) = ang_hypo_temp(mid_start);

% set the max magnitude of angle variance bound, be careful about the case when angles
% larger tha pi or smaller than -pi
thre = 1.3;
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













