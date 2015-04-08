function XX = Hypo_2nd(Npop_particles, sub_num, var_speed, X, seg_len, len_max, len_min)
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
        
        
var_len_tail = 1;
var_len_head = 1;

scenario = [4,5,7,2,12];

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
    ts = 1:1/seg_len:m_fre_pt;
    
    % the extension curve, which is a combination of spline and linear
    % curves. (spline only seems weird sometimes)

    pt_sp_xy = interp1(t,Frenet_k_1.xy,ts,'spline');
   
    
    % the number '3' comes from 'ts'. The curve extends 3 segments to both head and
    % tail directions
    head_loc = pt_sp_xy(1,:);
    tail_loc = pt_sp_xy(end,:);

    %% parameters
    
    % the moving velocity 
    Vel = X{ii}.vel;
    Vel_norm = norm(Vel);

    % parameters prepared for 'sub_num' sub-particles
    speed = zeros(sub_num,1);
    omg = ones(sub_num,1) * X{ii}.omg;

    % do not consider velocity 
    speed_2nd = 0;
    
    % New skeleton prediction
    ske_pred_speed_ori = ske_prediction(speed_2nd, pt_sp_xy, flipud(pt_sp_xy),head_loc, tail_loc, seg_len);

    ske_pred_1 = ske_pred_speed_ori;
   
    % calculate the angle of each point
    [pt_ske_pred, vec_len_pred, angle_ske_pred] = ske2ang(ske_pred_1, seg_len);
            
          
    %% generate sub-hypotheses    
    for jj = 1:sub_num;
        
    
            % calculate the angle hypothesis. Two different functions for first
            % half and second half of the worm. The mid segment is mid_start to
            % mid_end. 
            [ang_hypo,mid_start, jj_4] = ang2ang_2nd(angle_ske_pred, jj, scenario);
          
            mid_hypo = mid2mid(ske_pred_1, pt_ske_pred, mid_start, var_ske_direc, var_direc);
            
            if jj > jj_4 - 3
                if jj == jj-2
                    [X_mid, Y_mid] =  pol2cart(ang_hypo(mid_start)+pi/2,0.2*ave_vec_len_pred);
                    mid_hypo = pt_ske_pred(mid_start,:)+[X_mid, Y_mid];
                elseif jj == jj-1
                    [X_mid, Y_mid] =  pol2cart(ang_hypo(mid_start)+pi/2,0.2*ave_vec_len_pred);
                    ang_hypo = mid_pt_chg(angle_ske_pred,mid_start,-1);
                elseif jj < jj_4 + scenario(5)/2
                    vec_len_pred(1:2) = max( vec_len_pred(1:2) + rand(1) * var_len_tail * [1;1],[1;1]);
                else 
                    vec_len_pred(end-1:end) = max( vec_len_pred(end-1:end) + rand(1) * var_len_head * [1;1],[1;1]);
                end
            end
        
           [vec_len_pred_adjust, ave_vec_len_pred] = bound_vec(vec_len_pred, len_min, len_max);
          
        % calculate hypothesis skeleton according to angles
        ske_hypo = ang2ske(ang_hypo, vec_len_pred_adjust, mid_hypo, mid_start);
        
        %% compute results
        
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
        XX{ii,jj}.D = 0;

    end % end of jj

end   % end of ii
