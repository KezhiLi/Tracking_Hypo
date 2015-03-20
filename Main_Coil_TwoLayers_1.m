% The main functin to track a single worm
%
% In the current version, the inputs are:
%   vr: the input video file route (please change the location accordingly)
%   Frenet_Coil: a cell format data, has 
%               .xy: an n*2 vector, the skeleton of the hypothese worm
%               .T:  an n*3 vector, the T vector in frenet frame of each point on skeleton
%               .N:  an n*3 vector, the N vector in frenet frame of each point on skeleton
%               This is a file we generated in advance. The useful part is
%               its data format, and the skeleton points initially. 
%
% Output: a figure shown realtime,       
%         and a tracking video saved as 'filename'
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 23/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%% Parameters

% please add the folder name here
addpath(genpath('C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_10\.'));

% the file location to save current tracking video
% filename = 'results\testworm_1(3.5-5-50)19-Mar15.gif';
filename = 'results\testworm_20Mar15-0(3.5-20-100).gif';

fname = ['results\testworm_1(3.5-20-100)',date,'.avi' ];

% the number of particles(hypotheses) are saved after each iteration
N_particles = 10;  % 10
% the number of sub-particles generated in each iteration
%sub_num = 50;  % 50
sub_num_1 = 100;
sub_num_2 = 50

% the length (pixels) of each segment of the skeleton (this value relates to Frenet_Coil)
seg_len = 8;  % 8 

% the proximated various of the image (0~255)
Xstd_rgb = 60; % 40
% the first derivative of the worm velocity (pixels/second)
var_speed = 5; % 2

% the width of the worm (pixels= width *2)
width = 3.5; % 3.5

% video rate
fps = 10; 

%% Loading Movie
% the input video
%vr = VideoReader('\Sample_Video\Video_coil.avi');
vr = VideoReader('\Sample_Video\Video_09-Mar-2015.avi');
% the initial state (skeleton, frenent N,T, etc)
%load Frenet_Coil;
load Frenet_1903.mat;

% The resolution of each frame
Npix_resolution = [ vr.Height  vr.Width];
% The total number of frames
Nfrm_movie = floor(vr.Duration * vr.FrameRate);
% The inner correlation between the real image and the predicted image
inn_result = zeros(Nfrm_movie,1);

Y_1 = read(vr, 1);
Y_2 = read(vr, 2);

%% Object Tracking by Particle Filter

% Initial predicted worm
X = create_particles_hypo(N_particles, Frenet_Pt{2});

% The worm's spine shown in figure
worm_show = [];
%worm_show{1} = X{1}.xy;
worm_show = X{1}.xy;

% length max, min
len_max = 88;   % 100
len_min = 70;    % 80

% k represent the index of frame image
for k = 3:Nfrm_movie
    
    % Getting Image
    Y_k = read(vr, k);
    
    
  
    %% 1st layer
    %XX = Forecasting(N_particles, sub_num, var_speed, X, seg_len, len_max, len_min);
    
    XX = Hypo_1st(N_particles, sub_num_1, var_speed, X, seg_len, len_max, len_min);
    
    % Calculating Log Likelihood
    [L, C_k] = calc_log_likelihood_Worm(Xstd_rgb, XX, Y_k, width, seg_len);
      
    % Resampling
    X  = resample_particles_Worm(XX, L);
    
    %% 2nd layer
    XX = Hypo_2nd(N_particles, sub_num_2, var_speed, X, seg_len, len_max, len_min);
    
        % Indication purpose 
    k
    if mod(k,20)==0
        k
    end
    
    % Calculating Log Likelihood
    [L, C_k] = calc_log_likelihood_Worm(Xstd_rgb, XX, Y_k, width, seg_len);
      
    % Resampling
    X  = resample_particles_Worm(XX, L);
    
    
    
    hold on 

    % Weighted averaging best tt result to obtain the worm_show
    tt =3;
    ind = 1;
    [worm_show, X, inn_result(k)] = calculate_estimated_Worm(X, worm_show,  C_k, width,tt, seg_len, ind);
    
    show_Worm(worm_show, Y_k, width, seg_len, ind);
    
    for ii = 1:5;
        ind = ii
    % Show the estimated worm body (worm_show)
        show_Worm(X{ind}.xy, Y_k, width, seg_len, ind);
    
        drawnow
    end
    
    hold off
    % Save the figure shown as a frame of the output video 
    mov(k-2) = save_crt_fra(filename,k, fps);

end

% movie2avi(mov, fname, 'compression', 'None', 'fps', fps);

