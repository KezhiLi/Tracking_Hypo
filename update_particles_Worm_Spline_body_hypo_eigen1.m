function [XX, ske_pred, ske_hypo] = update_particles_Worm_Spline_body_hypo_eigen1(Npop_particles, var_speed, X, sub_num, Frenet_k_1,seg_len)

N = Npop_particles*sub_num;

eigen_num = 4;
        
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

    %    load eigenWorms.mat

    %[projectedAmps, fit_cuv_NonSaLeng] = eigenWormProject_NonSaLeng(eigenWorms, centr_angle,eigen_num, centr_vec_len);
        %[projectedAmps, fit_cuv_NonSaLeng] = eigenWormProject_NonSaLeng(eigenWorms, centr_angle,eigen_num);

 %   [Frenet_k_1.T,Frenet_k_1.N] = frenet_TN(Frenet_k_1.xy(:,1),Frenet_k_1.xy(:,2));

    
    % the moving velocity 
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);

    speed = zeros(sub_num,1);
    omg = ones(sub_num,1) * X{ii}.omg;

    
        
    for jj = 1:sub_num;
        
        % radial magnitude 

    
        if jj>22&&jj<41; %21
            speed(jj) = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))')+ var_speed* randn(1);
        elseif jj>40&&jj<47
            speed(jj) = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
        else 
            speed(jj) = 0.5*Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
        end
        
        if jj==sub_num-1;
            speed(jj) = 8;
        elseif jj ==sub_num
            speed(jj) = -8;
        end

%     if jj > 40&& jj < 40+3^eigen_num+1;
%                 jjj = jj - 40;  % from 41
%         mid_start = round(eigen_num/2);
%         mid_hypo = pt_centr(mid_start,:) + Vel;
%         alpha_new = alpha2alpha(projectedAmps, jjj);
%         ang_alp = alpha2angle(alpha_new, eigenWorms,eigen_num, centr_vec_len);
%         centr_curv_hypo = ang2ske(ang_alp', centr_vec_len, mid_hypo, mid_start);
%         ske_hypo = centr_curv_hypo(round(1:(end/(m_fre_pt)):end),:);
%     else
    % predict the next time skeleton based on previous skeleton
    %ske_pred = ske_prediction(speed(jj), [pt_sp_x',pt_sp_y'], [pt_sp_x_flp',pt_sp_y_flp']);
    ske_pred = ske_prediction_1(speed(jj), [pt_sp_x',pt_sp_y'], [flipud(pt_sp_x'),flipud(pt_sp_y')],head_loc, tail_loc);
    
    % calculate the angle of each point
    [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred, seg_len);

    sum_len_pred = sum(vec_len_pred)+randn(1);
    len_short = sum_len(ii) - sum_len_pred;

    
    % calculate the angle hypothesis. Two different functions for first
    % half and second half of the worm. The mid segment is mid_start to
    % mid_end. 

    [ang_hypo,ang_chan, mid_start, speed(jj), omg(jj),len_short] = ang2ang_acc(angle_ske_pred, omg(jj),jj,speed(jj),len_short);

    if sum_len_pred <91
        len_short = len_short + (91 - sum_len_pred);
    elseif sum_len_pred >104
        len_short = len_short + (104 - sum_len_pred);
%     elseif ii>19&&mod(ii,5)==0
%         len_short = len_short +2*sign(rand(1));
    end
    
    vec_len_pred = vec_len_pred + len_short/(length(vec_len_pred));
    ave_vec_len_pred = sum_len_pred/(length(vec_len_pred));
    vec_len_pred_adjust = 0.8*vec_len_pred +0.2*ave_vec_len_pred;
    if jj>40 && jj<47
        if jj ==41
            vec_len_pred_adjust(end) = vec_len_pred_adjust(end)+4;
        elseif jj ==42
            vec_len_pred_adjust(end) = max(vec_len_pred_adjust(end)-4,1);
        elseif jj ==43
            vec_len_pred_adjust(1) = vec_len_pred_adjust(1)+4;
        elseif jj == 44
            vec_len_pred_adjust(1) = max(vec_len_pred_adjust(1)-4,1);
        elseif jj == 45
            vec_len_pred_adjust = vec_len_pred_adjust + 4/length(vec_len_pred_adjust);
        else
            vec_len_pred_adjust = max(vec_len_pred_adjust - 4/length(vec_len_pred_adjust),1);  
        end
    end
    
    sum(vec_len_pred_adjust)
    
    
    % calculate hypothesis middle point
    var_ske_direc = 1; %3; % acturally duplicate to speed error
    var_direc = 0; %3;
    
    if jj==47 || jj ==48;
        [X_mid, Y_mid] =  pol2cart(ang_hypo(mid_start)+pi/2,0.2*ave_vec_len_pred);
        if jj==47
            mid_hypo = pt_ske_pred(mid_start,:)+[X_mid, Y_mid];
        else
            mid_hypo = pt_ske_pred(mid_start,:)-[X_mid, Y_mid];
        end
    else
        mid_hypo = mid2mid(ske_pred, pt_ske_pred, mid_start, var_ske_direc, var_direc);
    end
    
    
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
    
        XX{ii,jj}.omg = omg(jj); 

    end
    
     end   
    
%    br_mtx = kron(eye(size(X,2)),ones(1,sub_num));
%    XX(1:2,:) = X(1:2,:)*br_mtx + Xstd_pos * randn(2, N);


% hold on
% plot(mean(XX(1,:)),mean(XX(2,:)),'ro')
% hold off


end


