function XX = Hypo_1st(Npop_particles, sub_num, var_speed, var_len, X, seg_len, len_max, len_min, estimate_shift)
% function to calculate the rough hypothese as the 1st layer of the 2-lay
% hypothesis tracking algoirthm.
%
% Input:
%       Npop_particles: a scalar, the number of particles(hypotheses) are saved after each iteration
%              sub_num: a scalar, the number of sub-particles generated for each particle saved
%            var_speed: a scalar, the first derivative of the step size of worm velocity
%                        (pixels/second), pointing from tail to head
%              var_len: a scalar, the variance of worm length during length-changing
%                       hypothesis
%                    X: 1 * Npop_particle cell, the saved Npop_particle hypotheses
%           Frenet_k_1: only use its format: .xy, .T, .N
%                       .xy: points on skeleton from tail to head
%                       .T : frenet tangent direction
%                       .N : frenet perpendicular direction
%              seg_len: a scalar, the number of points in each segment of the skeleton
%              len_max: the upper bound of the length of the worm
%              len_min: the lower bound of the length of the worm
% Output:
%                   XX: Npop_particle * sub_num cell, all hypotheses
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 30/03/2015
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.


% calculate hypothesis middle point
% variance along the tangent direction
var_ske_direc = 0; %3;
% random noise
var_direc = 0; %3;

% the location variance of the middle point of the worm spline
var_mid = 0.25;

% loop for every particle, generate sub_num sub-particles, in which
for ii=1:Npop_particles;
    %% calcaulte extension curves to head and tail
    
    % The length of points on skeleton
    m_fre_pt = size(X{1,1}.xy,1);
    
    % It is the state of k-1 time
    Frenet_k_1 = X{ii};
    
    % 't' and 'ts' are used for spline/interp1 to extend current skeleton to
    % head and tail direction, respectively
    t = 1: m_fre_pt;
    ts = -2:1/seg_len:m_fre_pt+3;
    
    %% the extension curve, which is a combination of spline and linear curves.
    %(spline seems weird sometimes, particularly the extension part)
    
    % linear interpolation
    linear_xy = interp1(t,Frenet_k_1.xy,ts,'linear','extrap');
    % spline interpolation
    spline_xy = interp1(t,Frenet_k_1.xy,ts,'spline','extrap');
    
    % extension part
    ext_xy = 0.8*linear_xy + 0.2*spline_xy;
    % inner part
    inn_xy = 0.2*linear_xy + 0.8*spline_xy;
    
    % the interpolation obtained
    pt_sp_xy = [ext_xy(1:3*seg_len,:);inn_xy(3*seg_len+1:end-3*seg_len,:);ext_xy(end-3*seg_len+1:end,:)];
    
    % the number '3' comes from 'ts'. The curve extends 3 segments to both head and
    % tail directions

    head_loc = pt_sp_xy(seg_len*3+1,:);
    tail_loc = pt_sp_xy(end-seg_len*3,:);
    
    % skeleton to angle, only use the length of each segments/vector
    [pt_ske, vec_len, angle_ske] = ske2ang(Frenet_k_1.xy, 1);
    % record the sum of segments lengths
    sum_len(ii) = sum(vec_len);
    
    %% parameters
    
    % the moving velocity
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);
    
    %% keep the same speed
    speed_ori = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
    
    % modify the speed_ori accordingly based on estimate_shift
    shift = -estimate_shift;
    if shift>0.5 && speed_ori<shift*1
        speed_ori = shift*1;
    elseif shift<-0.5 && speed_ori>shift*1
        speed_ori = shift*1;
    end
    
    % parameters prepared for 'sub_num' sub-particles
    speed = ones(sub_num,1)*speed_ori;
    omg = ones(sub_num,1) * X{ii}.omg;
    
    % New skeleton prediction
    ske_pred_speed_ori = ske_prediction(speed_ori, pt_sp_xy, flipud(pt_sp_xy),head_loc, tail_loc, seg_len);
    
    ske_pred_1 = ske_pred_speed_ori;
    % interpolation to keep the same number of samples
    len_ske_pred = size(ske_pred_1,1);
    len_ske_old = size(inn_xy(3*seg_len+1:end-3*seg_len,:),1);
    
    [length_ske_pred_1, dis_ske_1] = pt_len(ske_pred_1);
    if len_ske_pred == len_ske_old && (max(dis_ske_1)/min(dis_ske_1)<2)
        ske_pred =  ske_pred_1;
    else
        ske_pred = re_interp(len_ske_pred, len_ske_old, ske_pred_1);
        % this is slower
        %ske_pred_tmp1 = curvspace(ske_pred_1,len_ske_old);
    end 

    % calculate the angle of each point
    [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred, seg_len);
    
    
    
    %% generate sub-hypotheses
    for jj = 1:sub_num;
        if jj<76 || (jj==sub_num-1) || (jj==sub_num)
            if jj<76
                vel_i = floor((jj-1)/15)-2;
                speed(jj) = speed_ori + vel_i * var_speed;
                % New skeleton prediction
                ske_pred_1 = ske_prediction(speed(jj), pt_sp_xy, flipud(pt_sp_xy),head_loc, tail_loc, seg_len);
            elseif jj==sub_num-1
                ske_pred_1 = pt_sp_xy(1*seg_len+1:end-5*seg_len,:);
            elseif jj==sub_num
                ske_pred_1 = pt_sp_xy(5*seg_len+1:end-1*seg_len,:);
            end
            % interpolation to keep the same number of samples
            len_ske_pred = size(ske_pred_1,1);
            len_ske_old = size(inn_xy(3*seg_len+1:end-3*seg_len,:),1);
            if len_ske_pred == len_ske_old && (max(dis_ske_1)/min(dis_ske_1)<2)
                ske_pred_mv =  ske_pred_1;
            else
                ske_pred_mv = re_interp(len_ske_pred, len_ske_old, ske_pred_1);
                % this is slower
                % ske_pred_tmp1 = curvspace(ske_pred_1,len_ske_old);
            end

            % calculate the angle of each point
            [pt_ske_pred_mv, vec_len_pred_mv, angle_ske_pred_1] = ske2ang(ske_pred_mv, seg_len);
            
            if jj<76
                ang_head_i = floor((mod(jj-1,15))/3)-2;  % ang_head_i = -2, -1, 0, 1, 2
                ang_tail_i = floor(mod(jj-1,3))-1;       % ang_tail_i = -1, 0, 1
            else % (jj==sub_num-1) || (jj==sub_num)
                ang_head_i = 0;
                ang_tail_i = 0;
            end
            
            % calculate the angle hypothesis. Two different functions for first
            % half and second half of the worm. The mid segment is mid_start to
            % mid_end.
            [ang_hypo,mid_start, omg(jj)] = ang2ang_1st(angle_ske_pred_1, ang_head_i,ang_tail_i ,omg(jj),jj);
            
        elseif jj>75 && jj<94
            
            vel_i = 0;
            ang_tail_i = 0;
            
            speed(jj) = speed_ori + vel_i * var_speed;
            
            if  jj<93
                ang_head_i = floor((jj-76)/6)-1;   % ang_head_i = -1, 0, 1
            else
                ang_head_i = rand(1) * 2 - 1;       % ang_head_i = [-1, 1]
            end
            
            % calculate the angle hypothesis. Two different functions for first
            % half and second half of the worm. The mid segment is mid_start to
            % mid_end.
            [ang_hypo,mid_start, omg(jj)] = ang2ang_1st(angle_ske_pred, ang_head_i,ang_tail_i ,omg(jj),jj);
            
            jjj = mod((jj-75),6);
            vec_len_pred_2 = vec_len_pred;
            para_1_2 = 0.9;  % 90% probability excute the first command, 10% probability else 
            
            len_stp_siz = rand(1);
            if jjj == 1  % head segment length plus
                if rand(1)<para_1_2
                    vec_len_pred_2(end) = vec_len_pred_2(end)+var_len*len_stp_siz;
                    speed(jj) = speed(jj) + len_stp_siz;
                else
                    vec_len_pred_2(end-1:end) = vec_len_pred_2(end-1:end)+[var_len*rand(1);var_len*rand(1)];
                end
            elseif jjj == 2  % head segment length minus
                if rand(1)<para_1_2
                    vec_len_pred_2(end) = max(vec_len_pred_2(end)-var_len*len_stp_siz,1);
                    speed(jj) = speed(jj) - len_stp_siz;
                else
                    vec_len_pred_2(end-1:end) = max(vec_len_pred_2(end-1:end)-[var_len*rand(1);var_len*rand(1)],[1;1]);
                end
            elseif jjj == 3   % tail segment length plus
                if rand(1)<para_1_2
                    vec_len_pred_2(1) = vec_len_pred_2(1)+var_len*rand(1);
                    speed(jj) = speed(jj) + len_stp_siz;
                else
                    vec_len_pred_2(1:2) = vec_len_pred_2(1:2)+[var_len*rand(1);var_len*rand(1)];
                end
            elseif jjj == 4   % tail segment length minus
                if rand(1)<para_1_2
                    vec_len_pred_2(1) = max(vec_len_pred_2(1)-var_len*rand(1),1);
                    speed(jj) = speed(jj) - len_stp_siz;
                else
                    vec_len_pred_2(1:2) = max(vec_len_pred_2(1:2)-[var_len*rand(1);var_len*rand(1)],[1;1]);
                end
            elseif jjj == 5   % all segment length plus
                if rand(1)<para_1_2
                    vec_len_pred_2 = vec_len_pred_2 + (var_len*rand(1))/length(vec_len_pred_2);
                else
                    vec_len_pred_2 = vec_len_pred_2 + (2*var_len*rand(1))/length(vec_len_pred_2);
                end
            else             % all segment length minus
                if rand(1)<para_1_2
                    vec_len_pred_2 = max(vec_len_pred_2 - (2*var_len*rand(1))/length(vec_len_pred_2),1);
                else
                    vec_len_pred_2 = max(vec_len_pred_2 - var_len*rand(1)/length(vec_len_pred_2),1);
                end
            end
        elseif jj>93 && jj<96
            % jj=93,96 are for middle point shift, pi/2 means swing 90
            % degree to the perpendicular direction
            [X_mid, Y_mid] =  pol2cart(ang_hypo(mid_start)+pi/2,0.2*ave_vec_len_pred);
            if jj==94
                % new coordinates of the middle point
                mid_hypo = pt_ske_pred(mid_start,:)+[X_mid, Y_mid];
                % the angle change due to the middle point change
                ang_hypo = mid_pt_chg(angle_ske_pred,mid_start,1, var_mid);
            else
                % new coordinates of the middle point, twist to the other side
                mid_hypo = pt_ske_pred(mid_start,:)-[X_mid, Y_mid];
                % the angle change due to the middle point change
                ang_hypo = mid_pt_chg(angle_ske_pred,mid_start,-1, var_mid);
            end
        elseif jj>95 && jj<max(99,sub_num-5)
            vel_i = 0;
            speed(jj) = speed_ori + vel_i * 4;
            
            head_or_tail = randn(1)>0;
            if head_or_tail
                ang_head_i = sign(randn(1))*2;
                ang_tail_i = 0;
            else
                ang_head_i = 0;
                ang_tail_i = sign(randn(1))*2;
            end
            % omg(jj) = omg(jj) + rand(1) * 0.2;
            
            % calculate the angle hypothesis. The angles can be changed
            % here are only first/last two points near tail/head.
            [ang_hypo,mid_start, omg(jj)] = ang2ang_1st_narrow(angle_ske_pred, ang_head_i,ang_tail_i ,omg(jj));
        else
            omg(jj) = omg(jj) + rand(1) * 0.2;
            [ang_hypo,mid_start, omg(jj)] = ang2ang_1st(angle_ske_pred, ang_head_i,ang_tail_i ,omg(jj),jj);
        end
        
        if jj~=94&& jj ~=95
            if jj < 76 || (jj==sub_num-1) || (jj==sub_num)
                % function to calculate new middle point
                mid_hypo = mid2mid(ske_pred_mv, pt_ske_pred_mv, mid_start, var_ske_direc, var_direc);
            else
                % function to calculate new middle point
                mid_hypo = mid2mid(ske_pred, pt_ske_pred, mid_start, var_ske_direc, var_direc);
            end
        end
        

        if jj<76
            vec_len_pred_now = vec_len_pred_mv;
        else
            vec_len_pred_now = vec_len_pred_2;
        end
        
        [vec_len_pred_adjust, ave_vec_len_pred] = bound_vec(vec_len_pred_now, len_min, len_max);
        
        % calculate hypothesis skeleton according to angles
        ske_hypo = ang2ske(ang_hypo, vec_len_pred_adjust, mid_hypo, mid_start);
        
        
        % calculate frenet vector 'T' and 'N'
        [XX{ii,jj}.T,XX{ii,jj}.N] = frenet_TN(ske_hypo(:,1),ske_hypo(:,2));
        
        XX{ii,jj}.xy = ske_hypo;
        if norm(XX{ii,jj}.T(1:2))~= 0
            XX{ii,jj}.vel = -speed(jj)*XX{ii,jj}.T(round(end/2),1:2)/norm(XX{ii,jj}.T(round(end/2),1:2));
        else
            XX{ii,jj}.vel = X{ii}.vel;
        end
        
        % omega is the radial swing velocity
        XX{ii,jj}.omg = omg(jj);
        
    end % end of jj
    
end   % end of ii
