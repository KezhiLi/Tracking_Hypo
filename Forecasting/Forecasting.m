function [XX] = Forecasting(Npop_particles, sub_num, var_speed, X, seg_len, len_max, len_min)
%% The function generates hypotheses to forecast the next worm position
% Input: 
%       Npop_particles: a scalar, the number of particles(hypotheses) are saved after each iteration
%              sub_num: a scalar, the number of sub-particles generated for each particle
%       saved
%            var_speed: a scalar, the first derivative of the worm velocity
%       (pixels/second), pointing from tail to head
%                    X: 1 * Npop_particle cell, the saved Npop_particle hypotheses  
%           Frenet_k_1: only use its format: .xy, .T, .N
%                       .xy: points on skeleton from tail to head
%                       .T : frenet tangent direction
%                       .N : frenet perpendicular direction
%              seg_len: a scalar, the length (pixels) of each segment of
%       the skeleton
% Output: 
%                   XX: Npop_particle * sub_num cell, all hypotheses
%       
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 23/02/2015     
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%% Code

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
    pt_sp_x = 0.5*spline(t,Frenet_k_1.xy(:,1),ts) + 0.5*interp1(t,Frenet_k_1.xy(:,1),ts,'linear','extrap');
    pt_sp_y = 0.5*spline(t,Frenet_k_1.xy(:,2),ts) + 0.5*interp1(t,Frenet_k_1.xy(:,2),ts,'linear','extrap');
    
    % the number '3' comes from 'ts'. The curve extends 3 segments to both head and
    % tail directions
    head_loc = [pt_sp_x(seg_len*3+1),pt_sp_y(seg_len*3+1)];
    tail_loc = [pt_sp_x(end-seg_len*3),pt_sp_y(end-seg_len*3)];
    
    % skeleton to angle, only use the length of each segments/vector
    [pt_ske, vec_len, angle_ske] = ske2ang(Frenet_k_1.xy, 1);
    % record the sum of segments lengths
    sum_len(ii) = sum(vec_len);

    %% eigenworm part
%     centr_curv = [pt_sp_x(25:end-24)',pt_sp_y(25:end-24)'];
%     
%     % calculate the angle of each point
%     [pt_centr, centr_vec_len, centr_angle] = ske2ang(centr_curv, 1); 
% 
%         load eigenWorms.mat
%     [projectedAmps, fit_cuv_NonSaLeng] = eigenWormProject_NonSaLeng(eigenWorms, centr_angle,eigen_num, centr_vec_len);
%        [Frenet_k_1.T,Frenet_k_1.N] = frenet_TN(Frenet_k_1.xy(:,1),Frenet_k_1.xy(:,2));

    %% parameters
    
    % the moving velocity 
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);

    % parameters prepared for 'sub_num' sub-particles
    speed = zeros(sub_num,1);
    omg = ones(sub_num,1) * X{ii}.omg;

    
    %% generate sub-hypotheses    
    for jj = 1:sub_num;
        
        % radial magnitude 

        % jj = 23~40 are randomly generated radial ones
        if jj>22&&jj<41; 
            % new speed is equal to original speed plus a random change 
            speed(jj) = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))')+ var_speed* randn(1);
        % jj = 41~46 are length modification    
        elseif jj>40&&jj<47
            % new speed is equal to original speed
            speed(jj) = Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
        % jj = 1~22 are artificially direction change
        else 
            % new speed is equal to 0.8 of the original speed
            speed(jj) = 0.8*Vel_norm*sign(-Vel*(Frenet_k_1.T(round(end/2),1:2))');
        end
        
        % The last two sub-hypothese are dramatic speed increase or
        % decrease
        if jj==sub_num-1;
            speed(jj) = 8;
        elseif jj ==sub_num
            speed(jj) = -8;
        end
        
        % New skeleton prediction
        ske_pred = ske_prediction(speed(jj), [pt_sp_x',pt_sp_y'], [flipud(pt_sp_x'),flipud(pt_sp_y')],head_loc, tail_loc);

        % calculate the angle of each point
        [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred, seg_len);
        
        % sum lengths prediction
        sum_len_pred = sum(vec_len_pred)+randn(1);
        % The total length change
        len_short = sum_len(ii) - sum_len_pred;


        % calculate the angle hypothesis. Two different functions for first
        % half and second half of the worm. The mid segment is mid_start to
        % mid_end. 
        [ang_hypo,ang_chan, mid_start, speed(jj), omg(jj),len_short] = ang2ang_acc(angle_ske_pred, omg(jj),jj,speed(jj),len_short);

        % set a threshold of the skeleton lengths
        if sum_len_pred <len_min
            len_short = len_short + (len_min - sum_len_pred);
        elseif sum_len_pred >len_max
            len_short = len_short + (len_max - sum_len_pred);
        end

        % Compensate the length change to each segment
        vec_len_pred = vec_len_pred + len_short/(length(vec_len_pred));
        % Adjust the length of each segement, make them with more uniform
        % lengths
        ave_vec_len_pred = sum_len_pred/(length(vec_len_pred));
        vec_len_pred_adjust = 0.8*vec_len_pred +0.2*ave_vec_len_pred;
        
        %  jj=41~46: artificially change the total length
        if jj>40 && jj<47
            if jj ==41 % head segment length plus
                vec_len_pred_adjust(end) = vec_len_pred_adjust(end)+4;
            elseif jj ==42 % head segment length minus
                vec_len_pred_adjust(end) = max(vec_len_pred_adjust(end)-4,1);
            elseif jj ==43  % tail segment length plus
                vec_len_pred_adjust(1) = vec_len_pred_adjust(1)+4;
            elseif jj == 44  % tail segment length minus
                vec_len_pred_adjust(1) = max(vec_len_pred_adjust(1)-4,1);
            elseif jj == 45  % all segment length plus
                vec_len_pred_adjust = vec_len_pred_adjust + 4/length(vec_len_pred_adjust);
            else             % all segment length minus
                vec_len_pred_adjust = max(vec_len_pred_adjust - 4/length(vec_len_pred_adjust),1);  
            end
        end
        
        % indication use, show the total length after adjustment 
        sum(vec_len_pred_adjust)


        % calculate hypothesis middle point
        % variance along the tangent direction
        var_ske_direc = 1; %3; 
        % random noise 
        var_direc = 0; %3;
        
        % jj=47,48 are for middle point shift
        if jj==47 || jj ==48;
            [X_mid, Y_mid] =  pol2cart(ang_hypo(mid_start)+pi/2,0.2*ave_vec_len_pred);
            if jj==47
                mid_hypo = pt_ske_pred(mid_start,:)+[X_mid, Y_mid];
            else
                mid_hypo = pt_ske_pred(mid_start,:)-[X_mid, Y_mid];
            end
        else
            % function to calculate new middle point
            mid_hypo = mid2mid(ske_pred, pt_ske_pred, mid_start, var_ske_direc, var_direc);
        end


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



