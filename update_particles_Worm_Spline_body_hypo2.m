function [XX, ske_pred, ske_hypo] = update_particles_Worm_Spline_body_hypo2(Npop_particles, var_speed, Xstd_vec, X, sub_num, Frenet_k_1,seg_len)

N = Npop_particles*sub_num;
%seg_len = 8;

for ii=1:Npop_particles;
    
    m_fre_pt = size(X{1,1}.xy,1);
    
    Frenet_k_1.xy = X{ii}.xy;
    
    t = 1: m_fre_pt;
    ts = -2:1/seg_len:m_fre_pt+3;
    
    pt_sp_x = 0.5*spline(t,Frenet_k_1.xy(:,1),ts) + 0.5*interp1(t,Frenet_k_1.xy(:,1),ts,'linear','extrap');
    pt_sp_y = 0.5*spline(t,Frenet_k_1.xy(:,2),ts) + 0.5*interp1(t,Frenet_k_1.xy(:,2),ts,'linear','extrap');
    
    head_loc = [pt_sp_x(25),pt_sp_y(25)];
    tail_loc = [pt_sp_x(end-24),pt_sp_y(end-24)];
    
    [pt_ske, vec_len, angle_ske] = ske2ang(Frenet_k_1.xy, 1);
    sum_len(ii) = sum(vec_len);

    %% eigenworm part
    centr_curv = [pt_sp_x(25:end-24)',pt_sp_y(25:end-24)'];
    
    % calculate the angle of each point
    [pt_centr, centr_vec_len, centr_angle] = ske2ang(centr_curv, 1); 

        load eigenWorms.mat
        eigen_num = 5;
    [projectedAmps, fit_cuv_NonSaLeng] = eigenWormProject_NonSaLeng(eigenWorms, centr_angle,eigen_num);

    %%


    
    % Frenet Transform
    %[TT_k,NN_k,B_k,k_fre,t_fre,Frenet_k_1.T,Frenet_k_1.N] = frenet(Frenet_k_1.xy(:,1),Frenet_k_1.xy(:,2));
    [Frenet_k_1.T,Frenet_k_1.N] = frenet_TN(Frenet_k_1.xy(:,1),Frenet_k_1.xy(:,2));

% 
%     t = 1: m_fre_pt;
%     ts = -2:1/seg_len:m_fre_pt;
%     ts_flp = -2:1/seg_len:m_fre_pt;
%     
%     % A combination prediction of spline and linear on both sides 
%     % 0.3 and 0.7 are weightes that can be changed adaptively
%     pt_sp_x = 0.3*spline(t,Frenet_k_1.xy(:,1),ts) + 0.7*interp1(t,Frenet_k_1.xy(:,1),ts,'linear','extrap');
%     pt_sp_y = 0.3*spline(t,Frenet_k_1.xy(:,2),ts) + 0.7*interp1(t,Frenet_k_1.xy(:,2),ts,'linear','extrap');
%    
%     pt_sp_x_flp = 0.3*spline(t,Frenet_k_1.xy_flp(:,1),ts_flp) + 0.7*interp1(t,Frenet_k_1.xy_flp(:,1),ts_flp,'linear','extrap');
%     pt_sp_y_flp = 0.3*spline(t,Frenet_k_1.xy_flp(:,2),ts_flp) + 0.7*interp1(t,Frenet_k_1.xy_flp(:,2),ts_flp,'linear','extrap');
%     
    % draw two splines
%     hold on
%     plot(pt_sp_x,pt_sp_y,'m'),
%     plot(pt_sp_x_flp,pt_sp_y_flp,'g'),
%     hold off

%     pt_sp_x_mid = pt_sp_x(round(length(pt_sp_x)/2));
%     pt_sp_y_mid = pt_sp_y(round(length(pt_sp_x)/2));
%     % the distance between the mid point and other points on the skeleton
%     dis = [pt_sp_x;pt_sp_y] - kron(ones(1,length(pt_sp_x)),[pt_sp_x_mid;pt_sp_y_mid]);
%     % absolute value of the distance
%     dis_norm = sqrt(dis(1,:).^2+dis(2,:).^2);
    
    % the moving velocity 
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);

    speed = zeros(sub_num,1);
    
%     Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))')
%     
%     Vel_norm
%     
        var_ang = 0.01; % 1.5;
        
    for jj = 1:sub_num;
    
        if jj>20; %21
    speed(jj) = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))')+ var_speed* randn(1);
        else 
            speed(jj) = 0.5*Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
        end
        
        if jj==sub_num-1
            speed(jj) = 8;
        elseif jj ==sub_num
            speed(jj) = -8;
        end
        
     if jj<21
        if jj<4
            speed(jj) = speed(jj)+1;
        elseif jj<7&&jj>3
            speed(jj) = speed(jj)+1;
        elseif jj<10&&jj>6
            speed(jj) = speed(jj)-1;
        elseif jj<13&&jj>9
            speed(jj) = speed(jj)-1;
        elseif jj == 13;
            speed(jj) = speed(jj)+1;
        elseif jj==14;
            speed(jj) = speed(jj)+1;
        elseif jj==15;
            speed(jj) = speed(jj)-1;
        elseif jj==16;       
            speed(jj) = speed(jj)-1;
        elseif jj==17
            speed(jj) = speed(jj)+3;
        elseif jj==18
            speed(jj) = speed(jj)+3;
        elseif jj==19
            speed(jj) = speed(jj)-3;
        else
            speed(jj) = speed(jj)-3;
        end
    end
    
    % predict the next time skeleton based on previous skeleton
    %ske_pred = ske_prediction(speed(jj), [pt_sp_x',pt_sp_y'], [pt_sp_x_flp',pt_sp_y_flp']);
    ske_pred = ske_prediction_1(speed(jj), [pt_sp_x',pt_sp_y'], [flipud(pt_sp_x'),flipud(pt_sp_y')],head_loc, tail_loc);
    
    % calculate the angle of each point
    [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred, seg_len);

    sum_len_pred = sum(vec_len_pred)+randn(1);
    len_short = sum_len(ii) - sum_len_pred;
    
    if sum_len_pred <92
        len_short = len_short + (92 - sum_len_pred);
    elseif sum_len_pred >105
        len_short = len_short + (105 - sum_len_pred);
    elseif ii>19&&mod(ii,5)==0
        len_short = len_short +2*sign(rand(1));
    end
    
    vec_len_pred = vec_len_pred + len_short/(length(vec_len_pred));
    ave_vec_len_pred = sum_len_pred/(length(vec_len_pred));
    vec_len_pred_adjust = 0.8*vec_len_pred +0.2*ave_vec_len_pred;
    sum(vec_len_pred_adjust)
    
    % calculate the angle hypothesis. Two different functions for first
    % half and second half of the worm. The mid segment is mid_start to
    % mid_end. 

    [ang_hypo,ang_chan, mid_start, speed(jj)] = ang2ang(angle_ske_pred, var_ang,jj,speed(jj));


    % calculate hypothesis middle point
    var_ske_direc = 1; %3; % acturally duplicate to speed error
    var_direc = 1; %3;
    mid_hypo = mid2mid(ske_pred, pt_ske_pred, mid_start, var_ske_direc, var_direc);
    
    % calculate hypothesis skeleton according to angles
    ske_hypo = ang2ske(ang_hypo, vec_len_pred_adjust, mid_hypo, mid_start);

    
    [XX{ii,jj}.T,XX{ii,jj}.N] = frenet_TN(ske_hypo(:,1),ske_hypo(:,2));
    
    XX{ii,jj}.xy = ske_hypo;
    if norm(XX{ii,jj}.T(1:2))~= 0
        %XX{ii,jj}.vel = speed(jj)*XX{ii,jj}.T(1:2)/norm(XX{ii,jj}.T(1:2)); 
        XX{ii,jj}.vel = -speed(jj)*XX{ii,jj}.T(round(end/2),1:2)/norm(XX{ii,jj}.T(round(end/2),1:2));
    else 
        XX{ii,jj}.vel = X{ii}.vel;
    end      
    
    end
    
    
%    br_mtx = kron(eye(size(X,2)),ones(1,sub_num));
%    XX(1:2,:) = X(1:2,:)*br_mtx + Xstd_pos * randn(2, N);


% hold on
% plot(mean(XX(1,:)),mean(XX(2,:)),'ro')
% hold off


end


