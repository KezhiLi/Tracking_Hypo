function [XX, ske_pred, ske_hypo] = update_particles_Worm_Spline_body_hypo(Npop_particles, var_speed, Xstd_vec, X, sub_num, Frenet_k_1)

N = Npop_particles*sub_num;
seg_len = 10;

for ii=1:Npop_particles;
    
    m_fre_pt = size(X{1,1}.xy,1);
    
    Frenet_k_1.xy = X{ii}.xy;
    
    t = 1: m_fre_pt;
    ts = -2:1/seg_len:m_fre_pt+3;
    
    pt_sp_x = 0.3*spline(t,Frenet_k_1.xy(:,1),ts) + 0.7*interp1(t,Frenet_k_1.xy(:,1),ts,'linear','extrap');
    pt_sp_y = 0.3*spline(t,Frenet_k_1.xy(:,2),ts) + 0.7*interp1(t,Frenet_k_1.xy(:,2),ts,'linear','extrap');
    
    head_loc = [pt_sp_x(31),pt_sp_y(31)];
    tail_loc = [pt_sp_x(end-30),pt_sp_y(end-30)];
    
    
%     m_fre_pt = size(X{1,1}.xy,1);
%     t = 1:m_fre_pt ;
%     ts = 1:1/seg_len:m_fre_pt;
%     
%     pt_sp_x = spline(t,X{ii}.xy(:,1),ts);
%     pt_sp_y = spline(t,X{ii}.xy(:,2),ts);
%     
%     ske_spline = [pt_sp_x',pt_sp_y'];
% 
% %     %  Seg skeleton from tail to head
% %     Frenet_k_1.xy = ske_spline([end:-seg_len:2,1],:);
% %     % Seg skeleton from head to tail
% %     Frenet_k_1.xy_flp = ske_spline([1:seg_len:end-1,end],:);  
%     %  Seg skeleton from tail to head
%     Frenet_k_1.xy = ske_spline([end:-seg_len:2,1],:);
%     % Seg skeleton from head to tail
%     Frenet_k_1.xy_flp = ske_spline([1:seg_len:end-1,end],:);  

    
    % Frenet Transform
    [TT_k,NN_k,B_k,k_fre,t_fre,Frenet_k_1.T,Frenet_k_1.N] = frenet(Frenet_k_1.xy(:,1),Frenet_k_1.xy(:,2));

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

    pt_sp_x_mid = pt_sp_x(round(length(pt_sp_x)/2));
    pt_sp_y_mid = pt_sp_y(round(length(pt_sp_x)/2));
    
    % the distance between the mid point and other points on the skeleton
    dis = [pt_sp_x;pt_sp_y] - kron(ones(1,length(pt_sp_x)),[pt_sp_x_mid;pt_sp_y_mid]);
    % absolute value of the distance
    dis_norm = sqrt(dis(1,:).^2+dis(2,:).^2);
    % the moving velocity 
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);

    speed = zeros(sub_num,1);
    
    Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))')
    
    Vel_norm
    
    
    for jj = 1:sub_num;
        
    speed(jj) = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))')+ var_speed* randn(1);
    
    % predict the next time skeleton based on previous skeleton
    %ske_pred = ske_prediction(speed(jj), [pt_sp_x',pt_sp_y'], [pt_sp_x_flp',pt_sp_y_flp']);
    ske_pred = ske_prediction_1(speed(jj), [pt_sp_x',pt_sp_y'], [flipud(pt_sp_x'),flipud(pt_sp_y')],head_loc, tail_loc);
    
    % calculate the angle of each point
    [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred, seg_len);

    % calculate the angle hypothesis. Two different functions for first
    % half and second half of the worm. The mid segment is mid_start to
    % mid_end. 
    var_ang = 1; % 1.5;
    [ang_hypo,ang_chan, mid_start] = ang2ang(angle_ske_pred, var_ang);
    
    % calculate hypothesis middle point
    var_ske_direc = 5; %3; % acturally duplicate to speed error
    var_direc = 3; %3;
    mid_hypo = mid2mid(ske_pred, pt_ske_pred, mid_start, var_ske_direc, var_direc);
    
    % calculate hypothesis skeleton according to angles
    % still wondering why 0.1 needs to be added here
    ske_hypo = ang2ske(ang_hypo, vec_len_pred+0.1, mid_hypo, mid_start);
    
    
    [TT_k_jj,NN_k_jj,B_k_jj,k_fre_jj,t_fre_jj,XX{ii,jj}.T,XX{ii,jj}.N] = frenet(ske_hypo(:,1),ske_hypo(:,2));
    
    XX{ii,jj}.xy = ske_hypo;
    if norm(XX{ii,jj}.T(1:2))~= 0
        XX{ii,jj}.vel = speed(jj)*XX{ii,jj}.T(1:2)/norm(XX{ii,jj}.T(1:2));
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


