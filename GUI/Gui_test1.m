function varargout = Gui_test1(varargin)
% GUI_TEST1 MATLAB code for Gui_test1.fig
%      GUI_TEST1, by itself, creates a new GUI_TEST1 or raises the existing
%      singleton*.
%
%      H = GUI_TEST1 returns the handle to a new GUI_TEST1 or the handle to
%      the existing singleton*.
%
%      GUI_TEST1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TEST1.M with the given input arguments.
%
%      GUI_TEST1('Property','Value',...) creates a new GUI_TEST1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gui_test1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gui_test1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gui_test1

% Last Modified by GUIDE v2.5 22-Jun-2015 16:20:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gui_test1_OpeningFcn, ...
                   'gui_OutputFcn',  @Gui_test1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end




% End initialization code - DO NOT EDIT


% --- Executes just before Gui_test1 is made visible.
function Gui_test1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gui_test1 (see VARARGIN)

% Choose default command line output for Gui_test1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Gui_test1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Gui_test1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = get(handles.edit1,'String');
hdf5_file = get(handles.edit3,'String');
fps =  str2num(get(handles.edit4,'String'));
para_thre_ini = str2num(get(handles.edit5,'String')); 
seg_len = str2num(get(handles.edit6,'String'));

% please add the folder name here
addpath(genpath([path,'.']));


% the initial state (skeleton, frenent N,T, etc)
%load .\Data_source\frenet(6544)_1906.mat
frenet_name = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ];  % eg. date = 19-Jun-2015
frenet_name_all = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ]; 

frenet_dir=[path,'.\Data_source\Frenet(',hdf5_file,')*.mat'];
frenet_prev_all=dir(frenet_dir);


if exist(frenet_name)
    disp('hdf5 file exists');
    load(frenet_name);
elseif exist(frenet_name_all)
    frenet_name_prev = frenet_prev_all.name;
    disp('hdf5 similar file exists');
    load(frenet_name_prev);
else
    read_hdf5_func_GUI(hdf5_file, path, fps, para_thre_ini, seg_len);
    load(frenet_name);
end







% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

path = get(handles.edit1,'String');
hdf5_file = get(handles.edit3,'String');
fps =  str2num(get(handles.edit4,'String'));
para_thre_ini = str2num(get(handles.edit5,'String')); 
seg_len = str2num(get(handles.edit6,'String'));

% please add the folder name here
addpath(genpath([path,'.']));

% the initial state (skeleton, frenent N,T, etc)
%load .\Data_source\frenet(6544)_1906.mat
frenet_name = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ];  % eg. date = 19-Jun-2015
frenet_name_all = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ]; 

frenet_dir=[path,'.\Data_source\Frenet(',hdf5_file,')*.mat'];
frenet_prev_all=dir(frenet_dir);


if exist(frenet_name)
    disp('hdf5 file exists');
    load(frenet_name);
elseif exist(frenet_name_all)
    frenet_name_prev = frenet_prev_all.name;
    disp('hdf5 similar file exists');
    load(frenet_name_prev);
else
    disp('please initialize first');
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

    num_diff_points1 = sum(sum((double(abs(double(Y_k_gray-mean(mean(Y_k_gray))) - double(Y_k_gray_old - mean(mean(Y_k_gray_old))))) > 30)))
    
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
    show_Worm_GUI(worm_show, Y_k, width, seg_len, ind, inn_result(k),1, handles);
    
    for ii = 5:-1:1;
        ind = ii;
    % Show the estimated worm body (worm_show)
        show_Worm_GUI(X{ind}.xy, Y_k, width, seg_len, ind, X{ind}.D, 7-ind, inn_result, handles);
    end
    

    % draw blue estimated countour
    C_k_outline2 = bwboundaries(C_k);
    C_k_outline2 = C_k_outline2{1};
    hold on, plot(C_k_outline2(:,2)+1,C_k_outline2(:,1), 'LineWidth',1.2,'Color',[0 0.5 1]);
    
    drawnow
    
    hold off
    
    pt_len(X{1}.xy)
    
    % Save the figure shown as a frame of the output video 
   % mov(k-2) = save_crt_fra(filename,k, fps);
    
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% clipboard('copy',get(handles.editRel,'string'))


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
