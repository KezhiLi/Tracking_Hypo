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

clear 
clc

%% File Info
hdf5_file = '3955';  % 6544
path = 'C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_21\';

% please add the folder name here
addpath(genpath([path,'.']));

% video rate
fps = 25; 

para_thre_ini = 1;   % coil: 0.80  normal: 0.92 % 2:0.88  %5: 1 %6: 0.90

% the length (pixels) of each segment of the skeleton (this value relates to Frenet_Coil)
seg_len = 8;  % 8 

% the initial state (skeleton, frenent N,T, etc)
%load .\Data_source\frenet(6544)_1906.mat
frenet_name = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ];  % eg. date = 19-Jun-2015
frenet_name_all = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ]; 


frenet_dir=[path,'.\Data_source\Frenet(',hdf5_file,')*.mat'];
frenet_prev_all=dir(frenet_dir);
frenet_name_prev = frenet_prev_all.name;


if exist(frenet_name)
    disp('hdf5 file exists');
    load(frenet_name);
elseif exist(frenet_name_prev)
    disp('hdf5 similar file exists');
    load(frenet_name_prev);
else
    read_hdf5_func(hdf5_file, path, fps, para_thre_ini, seg_len);
    load(frenet_name);
end


%% Loading Movie
% the input video
input_fname = ['test_', hdf5_file, '.avi'];
vr = VideoReader(['\Sample_Video\hdf5\',input_fname]);

% the file location to save current tracking video
filename = ['results\test',hdf5_file,'-',date,'-0.gif'];
fname = ['results\test',hdf5_file,'-',date,'-0.avi' ];

% The resolution of each frame
Npix_resolution = [ vr.Height  vr.Width];
% The total number of frames
Nfrm_movie = floor(vr.Duration * vr.FrameRate);

%% Parameters
% the number of particles(hypotheses) are saved after each iteration
N_particles = 10;  % 10
% the number of sub-particles generated in each iteration
sub_num_1 = 120; 
sub_num_2 = 50;

% the estimated variance of the image (0~255)
Xstd_rgb =  50; % 60 % 40  % 75
% the first derivative of the worm velocity (pixels/second)
var_speed = 4; % 5
var_len   = 8; % 10

% length max, min    
% Frenet_1903.mat: (88,70); Frenet_Coil.mat: (105,85);
len_max = 120;   % 105  95      
len_min = 80;    % 85
size_blk = round((len_max+len_min)/12); 

% the half width of the worm (pixels= width *2)
width_ini = 5; % 3.5      Frenet_1903.mat: 3;  Frenet_Coil: 3.5;

%% Initial Setup
% The inner correlation between the real image and the predicted image
inn_result = zeros(Nfrm_movie,1);

Y_1 = read(vr, 1); % first frame
Y_2 = read(vr, 2); % second frame
Y_k_gray = 255 - rgb2gray(Y_2);

h_Y = size(Y_1,1);
w_Y = size(Y_1,2);
num_pixel = h_Y * w_Y;

% Initial predicted worm !!!!
%X = create_particles_hypo(N_particles*2, Frenet_Pt{2}.xy_flp);
X = create_particles_hypo(N_particles*2, Frenet_Pt{2}.xy);

% para_thre = para_thre_ini;
% width = width_ini;
% % select a proper para_thre/width adaptively 
[para_thre, width, D_mtx] = sel_prop_width(X, Y_k_gray, width_ini, seg_len,  para_thre_ini);


%% Object Tracking by Particle Filter

% texture initilization
texture{1}=ske2tex(X{1}.xy, width, Y_k_gray);
%texture{2}=ske2tex(X{1}.xy, width, Y_2);
texture_newY{1} = texture{1};
%texture_newY{2} = texture{2};

% initialize texture matrix
texture_mtx(1,1:9)=zeros(1,9);
%texture_mtx(2,1:6)=zeros(1,6);

% texture shift range
num1 = 9;
mid_num = (num1+1)/2;

% The worm's spine shown in figure
worm_show = [];
%worm_show{1} = X{1}.xy;
worm_show = X{1}.xy;

% minifold of number of hypo in 2nd layer
N1 = 1;

gap_frame = 3;
X_old_1_xy = X{1}.xy;

% key parameter to decide if there is missing frames
jump1 = 1;

% k represent the index of frame image
for k = 3:Nfrm_movie   % 3:Nfrm_movie
    
    
    % Getting Image
    Y_k = read(vr, k);
    Y_k_gray_old = Y_k_gray;
    Y_k_gray = 255 - rgb2gray(Y_k);

    % num_diff_points1 = sum(sum((double(abs(double(Y_k_gray-mean(mean(Y_k_gray))) - double(Y_k_gray_old - mean(mean(Y_k_gray_old))))) > 30)))
    
    %% if video fails for several frames
    if jump1 > 0
        if num_diff_points1 >  num_pixel*0.02
            jump1 = 0;
            Y_k_sure = Y_k_gray_old;
            k_sure = k;
            
            % texture matching 
            k_5_1 = mod(k,gap_frame);
            k_5_fold = floor(k/gap_frame)+1; 
            if k_5_1 == 0
                texture_mtx(k_5_fold,:) = texture_mtx(k_5_fold-1,:);
                texture{k_5_fold} = texture{k_5_fold-1};
            end       
        end
    else
        num_diff_points2 = sum(sum((double(abs(double(Y_k_gray-mean(mean(Y_k_gray))) - double(Y_k_sure - mean(mean(Y_k_sure))))) > 30)))
        if num_diff_points2 >  num_pixel*0.02 % && inn_result(k-1) > 400
            jump1 = 0;
            
            % texture matching 
            k_5_1 = mod(k,gap_frame);
            k_5_fold = floor(k/gap_frame)+1; 
            if k_5_1 == 0
                texture_mtx(k_5_fold,:) = texture_mtx(k_5_fold-1,:);
                texture{k_5_fold} = texture{k_5_fold-1};
            end       
        else
            jump1 = 1;
        end
    end
    
    %% if video fails for this frame
    if jump1 > 0 

        % central point shift
        X = pt_shift_comp(X, CMs, k);

        %% texture matching 
        k_5_1 = mod(k,gap_frame);
        k_5_fold = floor(k/gap_frame)+1;        
        if k_5_1 == 0
            texture_output = texture_bef(k_5_fold, num1, mid_num, X_old_1_xy, width, Y_k_gray, texture);
            texture_newY{k_5_fold} = texture_output{1};
            texture_mtx(k_5_fold,1:2) = texture_output{2};
        end 


        %% 1st layer
        XX1 = Hypo_1st(N_particles, sub_num_1, var_speed, var_len, X, seg_len, len_max, len_min, texture_mtx(k_5_fold,3));

        % Indication purpose 
        k
        if mod(k,39)==0
             k
        end

        % Calculating Log Likelihood
        [L, C_k, II] = calc_log_likelihood_Worm_1st(Xstd_rgb, XX1, Y_k_gray, width, seg_len, para_thre, len_min);

        % Resampling
        X1  = resample_particles_Worm(XX1, L, N1);

        %% 2nd layer
        XX2 = Hypo_2nd(N1*N_particles, sub_num_2, var_speed, X1, seg_len, len_max, len_min);


        % Calculating Log Likelihood
        [L, C_k, XX2] = calc_log_likelihood_Worm_2nd_2(Xstd_rgb, XX2, C_k, II, width, seg_len, size_blk);

        % Resampling
        X  = resample_particles_Worm(XX2, L, 1);

    else
        % reset jump1   
        jump1 = -1;
    end
    
    % hold on 

    % Weighted averaging best tt result to obtain the worm_show
    tt =3;
    ind = 0;
    [worm_show, X, inn_result(k)] = calculate_estimated_Worm(X, worm_show,  C_k, width,tt, seg_len, ind+1, size_blk);
    
    %% update texture information
    if k_5_1 == 0 && jump1 > 0    
        texture_output2 = texture_aft(k_5_fold, X{1}.xy, width, Y_k_gray, texture, texture_newY, num1, mid_num);   
        texture = texture_output2{1};
        texture_mtx = texture_output2{2};
        X_old_1_xy = X{1}.xy;
    end
    
    %% show sessions
    show_Worm(worm_show, Y_k, width, seg_len, ind, inn_result(k),1);
    
    for ii = 5:-1:1;
        ind = ii;
    % Show the estimated worm body (worm_show)
        show_Worm(X{ind}.xy, Y_k, width, seg_len, ind, X{ind}.D, 7-ind, inn_result);
    end
    

    % draw blue estimated countour
    C_k_outline2 = bwboundaries(C_k);
    C_k_outline2 = C_k_outline2{1};
    hold on, plot(C_k_outline2(:,2)+1,C_k_outline2(:,1), 'LineWidth',1.2,'Color',[0 0.5 1]);
    
    drawnow
    
    hold off
    
    pt_len(X{1}.xy)
    
    % Save the figure shown as a frame of the output video 
    %mov(k-2) = save_crt_fra(filename,k, fps);
    
end

% movie2avi(mov, fname, 'compression', 'None', 'fps', fps);

