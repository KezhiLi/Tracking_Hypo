%% Parameters

addpath(genpath('C:\Kezhi\MyCode!!!\.'));
addpath(genpath('C:\Kezhi\Software\SegWorm-master\SegWorm-master\.'));

fid = figure;
filename = 'testworm_17Feb15.gif';

F_update = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];

Npop_particles = 40;
sub_num = 60;

Xstd_rgb = 200;
Xstd_pos = 5;
Xstd_vec = 2;
var_speed = 6;

Xrgb_trgt = [51; 255; 51];

width = 

%% Loading Movie

vr = VideoReader('P:\before_09_02_15\Kezhi\MyCode\Sample\video\sample2\MyVideo\animation_short.avi');
load Frenet_Pt2;

% The resolution of each frame
Npix_resolution = [ vr.Height  vr.Width];
% The total number of frames
Nfrm_movie = floor(vr.Duration * vr.FrameRate);

%% Object Tracking by Particle Filter
seg_len =10;
X = create_particles_hypo(Npix_resolution, Npop_particles, Frenet_Pt{2},seg_len);

hf_ske_index = zeros(1,Nfrm_movie);

ske_pred = [];
hypo_num = Npop_particles*sub_num;  

for k = 3:Nfrm_movie
    
%     if k>4
%     hold on
%     plot(ske_pred(:,1),ske_pred(:,2),'LineWidth',2,'Color',[0.1 0.5 0.1]),
%     plot(ske_hypo(:,1),ske_hypo(:,2),'LineWidth',2,'Color',[0.3 0.3 0.8]),
%     hold off
%     end

    % Getting Image
    Y_k = read(vr, k);
    
    ske_pred_old = ske_pred;
    
    % Forecasting
    [XX, ske_pred, ske_hypo] = update_particles_Worm_Spline_body_hypo_eigen1(Npop_particles, var_speed, X, sub_num, Frenet_Pt{k-1},seg_len);

    
       %pause(3)
     
       
    % Calculating Log Likelihood
%     hf_ske_index(k) = round(size(Frenet_Pt{k}.xy,1)/2);
%     L = calc_log_likelihood_Worm_body(Xstd_rgb, Xrgb_trgt, XX(1:2, :), Y_k, Frenet_Pt{k}.xy(hf_ske_index(k),:));
    hf_ske_index(k) = round(size(Frenet_Pt{k}.xy,1)/2);
    %L = calc_log_likelihood_Worm_body_hypo(Xstd_rgb, Xrgb_trgt, XX, Y_k);
    [L, C_k] = calc_log_likelihood_Worm_body_hypo2(Xstd_rgb, Xrgb_trgt, XX, Y_k, width);
    
    % Resampling
    X  = resample_particles_Worm_body_hypo(XX, L, sub_num);
    
%     % Calculating Log Likelihood
%      L_log = calc_log_likelihood_Worm_body_hypo(Xstd_rgb, Xrgb_trgt, XX, Y_k);
%      L = exp(L_log - max(L_log));
%      Q = L / sum(L, 2);

    % Showing Image
    %show_particles(X, Y_k); 
    %show_state_estimated(X, Y_k);
    
    show_state_estimated_Worm_body_hypo(X, Y_k);

    %worm_shape = ske2shape(Frenet_Pt{k}.xy, Frenet_Pt{k}.N, 10, -0.14);
 
 hold on
     fold = 10;
    %figure(fid);
    line(Frenet_Pt{k}.xy(:,1),Frenet_Pt{k}.xy(:,2)), hold on
   quiver(Frenet_Pt{k}.xy(hf_ske_index(k),1),Frenet_Pt{k}.xy(hf_ske_index(k),2),Frenet_Pt{k}.T(hf_ske_index(k),1),Frenet_Pt{k}.T(hf_ske_index(k),2),0.4*fold,'color','m','linewidth',2)
   quiver(Frenet_Pt{k}.xy(hf_ske_index(k),1),Frenet_Pt{k}.xy(hf_ske_index(k),2),Frenet_Pt{k}.N(hf_ske_index(k),1),Frenet_Pt{k}.N(hf_ske_index(k),2),5*fold,'color','c','linewidth',2)
    
   % draw predicted skeleton based on previous skeleton
   % plot(worm_shape(:,1),worm_shape(:,2),'LineWidth',2,'Color',[0.2 0.2 0.3]);
   
 hold off

 drawnow
 
 frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if k == 3;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
      end


   
end

