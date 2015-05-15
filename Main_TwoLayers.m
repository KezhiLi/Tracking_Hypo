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
% 30/03/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

%% Parameters

% please add the folder name here
addpath(genpath('C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_16\.'));

% the file location to save current tracking video
% filename = 'results\testworm_1(3.5-5-50)19-Mar15.gif';
filename = 'results\testworm5_15May15-0(4-50-120).gif';
fname = ['results\testworm5_0(4-50-120)',date,'.avi' ];

%% Loading Movie
% the input video
vr = VideoReader('\Sample_Video\Video_Test6.avi');
%vr = VideoReader('\Sample_Video\Video_Test2.avi');
%vr = VideoReader('\Sample_Video\Video_coil.avi');
%vr = VideoReader('\Sample_Video\Video_09-Mar-2015.avi');

% the initial state (skeleton, frenent N,T, etc)
load .\Data_source\Frenet_1405.mat
%load .\Data_source\Frenet_2704.mat
%load .\Data_source\Frenet_2304(2).mat
%load Frenet_0904.mat;
%load Frenet_Coil;
%load Frenet_Pt_full;
%load Frenet_1903.mat;

% The resolution of each frame
Npix_resolution = [ vr.Height  vr.Width];
% The total number of frames
Nfrm_movie = floor(vr.Duration * vr.FrameRate);

%% Parameters
% the number of particles(hypotheses) are saved after each iteration
N_particles = 10;  % 10
% the number of sub-particles generated in each iteration
%sub_num = 50;  % 50
sub_num_1 = 120;
sub_num_2 = 50;

% the length (pixels) of each segment of the skeleton (this value relates to Frenet_Coil)
seg_len = 8;  % 8 

% the estimated variance of the image (0~255)
Xstd_rgb =  60; % 40  % 75
% the first derivative of the worm velocity (pixels/second)
var_speed = 5; % 5
var_len   = 10;

% the half width of the worm (pixels= width *2)
width = 3; % 3.5      Frenet_1903.mat: 3;  Frenet_Coil: 3.5;
para_thre = 0.90;   % coil: 0.80  normal: 0.92  % 0.99

% length max, min    
% Frenet_1903.mat: (88,70); Frenet_Coil.mat: (105,85);
len_max = 95;   % 105        
len_min = 70;    % 85
size_blk = round((len_max+len_min)/12); 

% video rate
fps = 10; 

%% Initial Setup
% The inner correlation between the real image and the predicted image
inn_result = zeros(Nfrm_movie,1);

Y_1 = read(vr, 1); % first frame
Y_2 = read(vr, 2); % second frame

h_Y = size(Y_1,1);
w_Y = size(Y_1,2);

%% Object Tracking by Particle Filter

% Initial predicted worm
X = create_particles_hypo(N_particles*2, Frenet_Pt{2});

% The worm's spine shown in figure
worm_show = [];
%worm_show{1} = X{1}.xy;
worm_show = X{1}.xy;

N1 = 1;

% k represent the index of frame image
for k = 3:Nfrm_movie   % 3:Nfrm_movie
    
    % Getting Image
    Y_k = read(vr, k);
    

    %% 1st layer
    %XX = Forecasting(N_particles, sub_num, var_speed, X, seg_len, len_max, len_min);
    
    XX1 = Hypo_1st(N_particles, sub_num_1, var_speed, var_len, X, seg_len, len_max, len_min);
    
    % Indication purpose 
    k
    if mod(k,10)==0
         k
    end
    
    % Calculating Log Likelihood
    [L, C_k, II] = calc_log_likelihood_Worm_1st(Xstd_rgb, XX1, Y_k, width, seg_len, para_thre);
      
    % Resampling
    X1  = resample_particles_Worm(XX1, L, N1);
    
    %% 2nd layer
    XX2 = Hypo_2nd(N1*N_particles, sub_num_2, var_speed, X1, seg_len, len_max, len_min);
    

    % Calculating Log Likelihood
    [L, C_k, XX2] = calc_log_likelihood_Worm_2nd_2(Xstd_rgb, XX2, C_k, II, width, seg_len, size_blk, len_min);
      
    % Resampling
    X  = resample_particles_Worm(XX2, L, 1);
    
    % hold on 

    % Weighted averaging best tt result to obtain the worm_show
    tt =3;
    ind = 0;
    [worm_show, X, inn_result(k)] = calculate_estimated_Worm(X, worm_show,  C_k, width,tt, seg_len, ind+1, size_blk);
    
    show_Worm(worm_show, Y_k, width, seg_len, ind, inn_result(k),1);
    
    for ii = 5:-1:1;
        ind = ii;
    % Show the estimated worm body (worm_show)
        show_Worm(X{ind}.xy, Y_k, width, seg_len, ind, X{ind}.D, 7-ind, inn_result);
    end
    
    C_k_outline = bwperim(C_k);
    [N, Vertices, Lines, Vertices_orignal] = curvature_N_areainput(C_k_outline, seg_len*2);
    hold on, plot(Vertices(:,1),h_Y-Vertices(:,2)+1, 'LineWidth',1.2,'Color',[0 0.5 1]);
    
    drawnow
    
    hold off
    
    pt_len(X{1}.xy)
    
    % Save the figure shown as a frame of the output video 
    mov(k-2) = save_crt_fra(filename,k, fps);

end

% movie2avi(mov, fname, 'compression', 'None', 'fps', fps);

