function XX = Hypo_1st(Npop_particles, sub_num, var_speed, X, seg_len, len_max, len_min)
% 
% 
% 
% 
% 
% 

        % calculate hypothesis middle point
        % variance along the tangent direction
        var_ske_direc = 0; %3; 
        % random noise 
        var_direc = 0; %3;
        

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
    
    % the extension curve, which is a combination of spline and linear
    % curves. (spline only seems weird sometimes)
    linear_x = interp1(t,Frenet_k_1.xy(:,1),ts,'linear','extrap');
    linear_y = interp1(t,Frenet_k_1.xy(:,2),ts,'linear','extrap');
    spline_x = spline(t,Frenet_k_1.xy(:,1),ts);
    spline_y = spline(t,Frenet_k_1.xy(:,2),ts);
    ext_x = 0.8*linear_x + 0.2*spline_x;
    ext_y = 0.8*linear_y + 0.2*spline_y;
    inn_x = 0.2*linear_x + 0.8*spline_x;
    inn_y = 0.2*linear_y + 0.8*spline_y;    
    
    pt_sp_x = [ext_x(1:3*seg_len),inn_x(3*seg_len+1:end-3*seg_len),ext_x(end-3*seg_len+1:end)];
    pt_sp_y = [ext_y(1:3*seg_len),inn_y(3*seg_len+1:end-3*seg_len),ext_y(end-3*seg_len+1:end)];
    %pt_sp_x = 0.5*spline(t,Frenet_k_1.xy(:,1),ts) + 0.5*interp1(t,Frenet_k_1.xy(:,1),ts,'linear','extrap');
    %pt_sp_y = 0.5*spline(t,Frenet_k_1.xy(:,2),ts) + 0.5*interp1(t,Frenet_k_1.xy(:,2),ts,'linear','extrap');
    
    % the number '3' comes from 'ts'. The curve extends 3 segments to both head and
    % tail directions
    head_loc = [pt_sp_x(seg_len*3+1),pt_sp_y(seg_len*3+1)];
    tail_loc = [pt_sp_x(end-seg_len*3),pt_sp_y(end-seg_len*3)];
    
    % skeleton to angle, only use the length of each segments/vector
    [pt_ske, vec_len, angle_ske] = ske2ang(Frenet_k_1.xy, 1);
    % record the sum of segments lengths
    sum_len(ii) = sum(vec_len);

    %% parameters
    
    % the moving velocity 
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);

    % parameters prepared for 'sub_num' sub-particles
    speed = zeros(sub_num,1);
    omg = ones(sub_num,1) * X{ii}.omg;

    %% keep the same speed  
    speed_ori = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
    
    % New skeleton prediction
    ske_pred_speed_ori = ske_prediction(speed_ori, [pt_sp_x',pt_sp_y'], [flipud(pt_sp_x'),flipud(pt_sp_y')],head_loc, tail_loc, seg_len);

    ske_pred_1 = ske_pred_speed_ori;
    % interpolation to keep the same number of samples
            len_ske_pred = size(ske_pred_1,1);
            len_ske_old = length(inn_x(3*seg_len+1:end-3*seg_len));
            if len_ske_pred == len_ske_old
                 ske_pred =  ske_pred_1;
            else
                 pp = 1:len_ske_pred;
                 qq = 1:((len_ske_pred-1)/(len_ske_old-1)):len_ske_pred;
                 ske_pred_x0 = interp1(pp,ske_pred_1(:,1),qq,'spline');
                 ske_pred_y0 = interp1(pp,ske_pred_1(:,2),qq,'spline');
                 ske_pred_2 = [ske_pred_x0',ske_pred_y0'];

                 qq2 = zeros(len_ske_old,1);
                 [length_ske_pred, dis_ske] = pt_len(ske_pred_2);
                 reci_dis_ske = [0;1./dis_ske];
                 sum_reci_dis_ske = sum(reci_dis_ske);
                 for q = 1:len_ske_old;
                      qq2(q) = (len_ske_old-1)*(sum(reci_dis_ske(1:q)))/sum_reci_dis_ske+1;
                 end
                 ske_pred_x1 = interp1(1:len_ske_old,ske_pred_x0,qq2,'spline');
                 ske_pred_y1 = interp1(1:len_ske_old,ske_pred_y0,qq2,'spline');
                 ske_pred = [ske_pred_x1,ske_pred_y1];
            end
                   
            size(ske_pred,1)
                   % calculate the angle of each point
            [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred, seg_len);
            

            
    %% generate sub-hypotheses    
    for jj = 1:sub_num;
        
        
    if jj<76;
        vel_i = floor((jj-1)/15)-2;
                   speed(jj) = speed_ori + vel_i * 4;
          
                   % New skeleton prediction
                   ske_pred_1 = ske_prediction(speed(jj), [pt_sp_x',pt_sp_y'], [flipud(pt_sp_x'),flipud(pt_sp_y')],head_loc, tail_loc, seg_len);

                   % interpolation to keep the same number of samples
                   len_ske_pred = size(ske_pred_1,1);
                   len_ske_old = length(inn_x(3*seg_len+1:end-3*seg_len));
                   if len_ske_pred == len_ske_old
                       ske_pred =  ske_pred_1;
                   else
                       pp = 1:len_ske_pred;
                       qq = 1:((len_ske_pred-1)/(len_ske_old-1)):len_ske_pred;
                       ske_pred_x0 = interp1(pp,ske_pred_1(:,1),qq,'spline');
                       ske_pred_y0 = interp1(pp,ske_pred_1(:,2),qq,'spline');
                       ske_pred_2 = [ske_pred_x0',ske_pred_y0'];

                       qq2 = zeros(len_ske_old,1);
                       [length_ske_pred, dis_ske] = pt_len(ske_pred_2);
                       reci_dis_ske = [0;1./dis_ske];
                       sum_reci_dis_ske = sum(reci_dis_ske);
                       for q = 1:len_ske_old;
                           qq2(q) = (len_ske_old-1)*(sum(reci_dis_ske(1:q)))/sum_reci_dis_ske+1;
                       end
                       ske_pred_x1 = interp1(1:len_ske_old,ske_pred_x0,qq2,'spline');
                       ske_pred_y1 = interp1(1:len_ske_old,ske_pred_y0,qq2,'spline');
                       ske_pred = [ske_pred_x1,ske_pred_y1];
                   end
                   
                   size(ske_pred,1)
                   % calculate the angle of each point
                   [pt_ske_pred, vec_len_pred, angle_ske_pred_1] = ske2ang(ske_pred, seg_len);
                   
                    ang_head_i = floor((mod(jj-1,15))/3)-2;
                    ang_tail_i = floor(mod(jj-1,3))-1;


                   % calculate the angle hypothesis. Two different functions for first
                   % half and second half of the worm. The mid segment is mid_start to
                   % mid_end. 
                   [ang_hypo,mid_start, omg(jj)] = ang2ang_1st(angle_ske_pred_1, ang_head_i,ang_tail_i ,omg(jj),jj);
               
                               

    elseif jj>75 && jj<94
            
            vel_i = 0;
            ang_tail_i = 0; 
            
            speed(jj) = speed_ori + vel_i * 4;
                               
            if  jj<93
                 ang_head_i = floor((jj-76)/6)-1;
            else
                 ang_head_i = rand(1) * 2 - 1;
            end
                     
            % calculate the angle hypothesis. Two different functions for first
            % half and second half of the worm. The mid segment is mid_start to
            % mid_end. 
            [ang_hypo,mid_start, omg(jj)] = ang2ang_1st(angle_ske_pred, ang_head_i,ang_tail_i ,omg(jj),jj);
                                                   
            jjj = mod((jj-75),6); 
            if jjj == 1  % head segment length plus
                ang_hypo(end) = ang_hypo(end)+4;
            elseif jjj == 2  % head segment length minus
                ang_hypo(end) = max(ang_hypo(end)-4,1);
            elseif jjj == 3   % tail segment length plus
                ang_hypo(1) = ang_hypo(1)+4;
            elseif jjj == 4   % tail segment length minus
                ang_hypo(1) = max(ang_hypo(1)-4,1);
            elseif jjj == 5   % all segment length plus
                ang_hypo = ang_hypo + 4/length(ang_hypo);
            else             % all segment length minus
                ang_hypo = max(ang_hypo - 4/length(ang_hypo),1);  
            end      
    elseif jj>93 && jj<96
        % jj=93,96 are for middle point shift
            [X_mid, Y_mid] =  pol2cart(ang_hypo(mid_start)+pi/2,0.2*ave_vec_len_pred);
            if jj==94
                mid_hypo = pt_ske_pred(mid_start,:)+[X_mid, Y_mid];
                ang_hypo = mid_pt_chg(angle_ske_pred,mid_start,1);
            else
                mid_hypo = pt_ske_pred(mid_start,:)-[X_mid, Y_mid];
                ang_hypo = mid_pt_chg(angle_ske_pred,mid_start,-1);
            end
    else
            vel_i = 0;
            ang_tail_i = 0; 
            
            speed(jj) = speed_ori + vel_i * 4;
                               
            omg(jj) = omg(jj) + rand(1) * 0.2;         
            % calculate the angle hypothesis. Two different functions for first
            % half and second half of the worm. The mid segment is mid_start to
            % mid_end. 
            [ang_hypo,mid_start, omg(jj)] = ang2ang_1st(angle_ske_pred, ang_head_i,ang_tail_i ,omg(jj),jj);
            
    end
        
        if jj~=94&& jj ~=95
            % function to calculate new middle point
            mid_hypo = mid2mid(ske_pred, pt_ske_pred, mid_start, var_ske_direc, var_direc);
        end
        
        
           % sum lengths prediction
           sum_len_pred = sum(vec_len_pred)+randn(1);
           % The total length change
           %         len_short = sum_len(ii) - sum_len_pred;
           len_short = 0;
                            
                            



        % set a threshold of the skeleton lengths
           if sum_len_pred <len_min
                 len_short = len_short + (len_min - sum_len_pred);
           elseif sum_len_pred >len_max
                 len_short = len_short + (len_max - sum_len_pred);
           end
               
     
        % Compensate the length change to the last 3 segments near tail
        vec_len_pred = vec_len_pred + len_short/length(vec_len_pred);
%         % Adjust the length of each segement, make them with more uniform
%         % lengths
          ave_vec_len_pred = sum_len_pred/(length(vec_len_pred));
%         vec_len_pred_adjust = 0.8*vec_len_pred +0.2*ave_vec_len_pred;
          vec_len_pred_adjust = vec_len_pred;
          
         % indication use, show the total length after adjustment 
        sum(vec_len_pred_adjust)
          
          

        



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
