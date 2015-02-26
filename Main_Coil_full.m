%% Parameters

%addpath(genpath('P:\Kezhi\MyCode\.'));
addpath(genpath('C:\Kezhi\Software\SegWorm-master\SegWorm-master\.'));
addpath(genpath('C:\Kezhi\MyCode!!!\EigenWorm\.'));


fid = figure;
filename = 'testworm_18Feb15-coil-full-0(3.5-10-50).gif';

F_update = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];

Npop_particles = 10;
sub_num = 50;

Xstd_rgb = 40; %150
Xstd_pos = 4; %5
alpha_list =[];
var_speed = 2; %5

Xrgb_trgt = [51; 255; 51];

width = 3.5;


%% Loading Movie

vr = VideoReader('C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Test_Coil\video\Vedeo_coil_full.avi');
load Frenet_Coil;

% The resolution of each frame
Npix_resolution = [ vr.Height  vr.Width];
% The total number of frames
Nfrm_movie = floor(vr.Duration * vr.FrameRate);

inn_result = zeros(Nfrm_movie,1);

%% Object Tracking by Particle Filter
seg_len =8;
X = create_particles_hypo(Npix_resolution, Npop_particles, Frenet_Pt_full{2},seg_len);

hf_ske_index = zeros(1,Nfrm_movie);

ske_pred = [];
hypo_num = Npop_particles*sub_num;  

worm_show = X{1}.xy;

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
    [XX, ske_pred, ske_hypo] = update_particles_Worm_Spline_body_hypo_eigen1(Npop_particles, var_speed, X, sub_num, Frenet_Pt_full{k-1},seg_len);

    
       %pause(3)
     if mod(k,15)==0&&k>5;
         k
     end
       
     k
     
    % Calculating Log Likelihood
%     hf_ske_index(k) = round(size(Frenet_Pt{k}.xy,1)/2);
%     L = calc_log_likelihood_Worm_body(Xstd_rgb, Xrgb_trgt, XX(1:2, :), Y_k, Frenet_Pt{k}.xy(hf_ske_index(k),:));
    hf_ske_index(k) = round(size(Frenet_Pt_full{k}.xy,1)/2);
    [L, C_k] = calc_log_likelihood_Worm_body_hypo2(Xstd_rgb, Xrgb_trgt, XX, Y_k, width);
    
    
    % Resampling
    X  = resample_particles_Worm_body_hypo(XX, L, sub_num);
    
    %% double layer hypothese 
    
        % Forecasting
    [XX, ske_pred, ske_hypo] = update_particles_Worm_Spline_body_hypo_eigen1(Npop_particles, var_speed, X, sub_num, Frenet_Pt_full{k-1},seg_len);

    
     %pause(3)
     if mod(k,15)==0&&k>5;
         k
     end
       
     k
     
    % Calculating Log Likelihood
    hf_ske_index(k) = round(size(Frenet_Pt_full{k}.xy,1)/2);
    [L, C_k] = calc_log_likelihood_Worm_body_hypo2(Xstd_rgb, Xrgb_trgt, XX, Y_k, width);
    
    
    % Resampling
    X  = resample_particles_Worm_body_hypo(XX, L, sub_num);
    
    %%
    
    tt =3;
    [worm_show, X, inn_result(k)] = calculate_estimated_Worm_body_hypo(X, worm_show,  C_k, width,tt);
    show_better_estimated_Worm_body_hypo(worm_show, Y_k, width);
    
    % show_state_estimated_Worm_body_hypo(X, Y_k, width);

 
 hold on
     fold = 10;
    %figure(fid);
    line(Frenet_Pt_full{k}.xy(:,1),Frenet_Pt_full{k}.xy(:,2)), hold on
   % quiver(Frenet_Pt{k}.xy(hf_ske_index(k),1),Frenet_Pt{k}.xy(hf_ske_index(k),2),Frenet_Pt{k}.T(hf_ske_index(k),1),Frenet_Pt{k}.T(hf_ske_index(k),2),0.4*fold,'color','m','linewidth',2)
   % quiver(Frenet_Pt{k}.xy(hf_ske_index(k),1),Frenet_Pt{k}.xy(hf_ske_index(k),2),Frenet_Pt{k}.N(hf_ske_index(k),1),Frenet_Pt{k}.N(hf_ske_index(k),2),5*fold,'color','c','linewidth',2)
    
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

