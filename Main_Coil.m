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
addpath(genpath('C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo\.'));

% the file location to save current tracking video
filename = 'results\testworm_02Mar15-0(3.5-10-50).gif';
% filename = 'testworm_26Feb15-0(3.5-10-50).gif';

% the number of particles(hypotheses) are saved after each iteration
N_particles = 10;  % 10
% the number of sub-particles generated in each iteration
sub_num = 50;  % 50

% the length (pixels) of each segment of the skeleton (this value relates to Frenet_Coil)
seg_len = 8;  % 8 

% the proximated various of the image (0~255)
Xstd_rgb = 60; % 40
% the first derivative of the worm velocity (pixels/second)
var_speed = 2; % 2

% the width of the worm (pixels= width *2)
width = 3.5; % 3.5

%% Loading Movie
% the input video
vr = VideoReader('\Sample_Video\Video_2702.avi');
% the initial state (skeleton, frenent N,T, etc)
%load Frenet_Coil;
load Frenet_2702.mat;

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
worm_show = X{1}.xy;

% k represent the index of frame image
for k = 3:Nfrm_movie
    
    % Getting Image
    Y_k = read(vr, k);
    
    % Forecasting
    [XX] = Forecasting(N_particles, sub_num, var_speed, X, seg_len);
    
    % Indication purpose 
    k

    % Calculating Log Likelihood
    [L, C_k] = calc_log_likelihood_Worm(Xstd_rgb, XX, Y_k, width, seg_len);
      
    % Resampling
    X  = resample_particles_Worm(XX, L);
    
    % Weighted averaging best tt result to obtain the worm_show
    tt =3;
    [worm_show, X, inn_result(k)] = calculate_estimated_Worm(X, worm_show,  C_k, width,tt, seg_len);
    
    % Show the estimated worm body (worm_show)
    show_Worm(worm_show, Y_k, width, seg_len);

    drawnow
    % Save the figure shown as a frame of the output video 
    save_crt_fra(filename,k);

end
