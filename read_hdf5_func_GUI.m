function read_hdf5_func_GUI(hdf5_file, path, fps, thre, samp_step)
% 
% 
% 
% 
% 
% 

%%
% movies_file = '/Volumes/behavgenom$/GeckoVideo/Invidual_videos/20150512/Capture_Ch1_12052015_194303/';

% kezhi_file = 'C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_20\Sample_Video\hdf5\worm_6544.hdf5';

%hdf5_file = '3955';  % 6544
%path = 'C:\Kezhi\MyCode!!!\Tracking\PF_Video_EN_Worm_Kezhi\PF_Video_EN\Tracking_Hypo_20\';
hdf5_path = [path,'Sample_Video\hdf5\worm_',hdf5_file,'.hdf5'];

%select a frame
%fps = 25;
frame = fps;

% threshold
%thre = 0.95;  % 0.96

% The sample step of points on skeleton. Change it accordingly in the Main
% function of video processing.
%samp_step = 8;

% read key data
masks = h5read(hdf5_path, '/masks');
frames = h5read(hdf5_path,'/frames');
CMs = h5read(hdf5_path, '/CMs');

%plot movement of the center of mass
% figure
% plot(CMs(1,:), CMs(2,:))

%%

% please add the folder name here
addpath(genpath([path,'.']));

% Read frames from the video
Y_1 = masks(:,:,1);
Y_2 = masks(:,:,2);
Y_3 = masks(:,:,3);

% Initialize the Frenet data (including points on skeleton, T vector, N vector, etc.)
Frenet_Pt{1}.xy = {};
Frenet_Pt{1}.xy_flp = {};
Frenet_Pt{1}.T = {};
Frenet_Pt{1}.N = {};

I = Y_1;
% I = rgb2gray(Y_1);
%I = rgb;
[mm,nn]= size(I);
%figure, imshow(I), title('Original')
% II = 255 - I;   
% 
% se = strel('disk', 1);
% II = imclose(II, se);
% %figure
% %imshow(II), title('Opening-closing (Ioc II)')
% 
%  III = imclose(II, se);
% % figure
% % imshow(III), title('Opening-closing (Ioc III)')

%figure, hist(double(I),256);
level = graythresh(I)*thre;  % 0.96
BW = im2bw(I,level);
% figure, imshow(BW)

reBW = BW;
% reBW = -im2bw(I,level)+1;
%  figure, imshow(reBW)

log_reBW = logical(reBW);

% fill the holds inside the worm body
BWdfill = imfill(log_reBW, 'holes');
% figure, imshow(BWdfill);
% title('binary image with filled holes');

BW2 = bwareaopen(BWdfill, 50); 
% figure, imshow(BW2);
% title('remove small areas');

% draw the contour of the worm
BWoutline = bwperim(BW2);
% figure, imshow(BWoutline);

%[NN, curv, Verticles, Lines] = curvature(BWoutline);
[NN, curv, Vertices, Lines, Vertices_old] = curvature_areainput_GUI(BW2);


% find the points with maximum curvature
[curv_norm, curv_point1, curv_point2] = Curv_Pt(Vertices, curv);

% draw the magnitudes of curvatures
% figure, plot(curv_norm)

% The matrix will be shown
Segout = gray2rgb(I);
Segout(BWoutline) = 1;

% The indexes of head the tail in the matrix shown
V_x=round([Vertices(curv_point1,1),Vertices(curv_point2,1)]);
V_y=mm-round([Vertices(curv_point1,2),Vertices(curv_point2,2)]);

% Draw the area of 5*5 near points of head and tail with green colour
for ii = V_y(1)-2:V_y(1)+2;
    for jj = V_x(1)-2:V_x(1)+2;
        Segout(ii,jj,:) = [0.2,0.2,1];
        %Segout(ii+2,jj+2,:) = [0.2,0.2,1];
    end
end
for ii = V_y(2)-2:V_y(2)+2;
    for jj = V_x(2)-2:V_x(2)+2;
        Segout(ii,jj,:) = [0.2,0.2,1];
        %Segout(ii+2,jj-5,:) = [0.2,0.2,1];
    end
end

% debug use
%hold on, plot([Verticles(curv_point1,1),Verticles(curv_point2,1)],[mm-Verticles(curv_point1,2),mm-Verticles(curv_point2,2)],'bo','MarkerSize', 10);

% The worm is roughly divided into 24 segments of musculature (i.e., hinges
% that represent degrees of freedom) on each side. Therefore, 48 segments
% around a 2-D contour.
sWormSegs = 24;
cWormSegs = 2 * sWormSegs;

% Clean up the worm's contour.
%contour = cleanWorm(Vertices, size(Vertices, 1) / cWormSegs);
contour = Vertices;

% Compute the contour's local high/low-frequency curvature.
% Note: worm body muscles are arranged and innervated as staggered pairs.
% Therefore, 2 segments have one theoretical degree of freedom (i.e. one
% approximation of a hinge). In the head, muscles are innervated
% individually. Therefore, we sample the worm head's curvature at twice the
% frequency of its body.
% Note 2: we ignore Nyquist sampling theorem (sampling at twice the
% frequency) since the worm's cuticle constrains its mobility and practical
% degrees of freedom.
cCCLengths = circComputeChainCodeLengths(contour);
wormSegLength = (cCCLengths(1) + cCCLengths(end)) / cWormSegs;
hfAngleEdgeLength = wormSegLength;
hfCAngles = circCurvature(contour, hfAngleEdgeLength, cCCLengths);
lfAngleEdgeLength = 2 * hfAngleEdgeLength;
lfCAngles = circCurvature(contour, lfAngleEdgeLength, cCCLengths);

% Compute the contour's local high/low-frequency curvature maxima.
%[mhfCMaxP mhfCMaxI] = maxPeaksCircDist(mhfCAngles, hfAngleEdgeLength, ...
%    cCCLengths);
[lfCMaxP lfCMaxI] = maxPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
    cCCLengths);


 % Compute the contour's local low-frequency curvature minima.
 [lfCMinP lfCMinI] = minPeaksCircDist(lfCAngles, lfAngleEdgeLength, ...
        cCCLengths);
 
 headI = curv_point1;
 tailI = curv_point2;
 % Compute the worm's skeleton.
 [skeleton cWidths] = linearSkeleton(headI, tailI, lfCMinP, lfCMinI, ...
        lfCMaxP, lfCMaxI, contour, wormSegLength, cCCLengths);


 % Measure the skeleton's chain code length.
 sCCLengths = computeChainCodeLengths(skeleton);
 sLength = sCCLengths(end);

 % Measure the skeleton angles (curvature).
 lfAngleEdgeLength = sCCLengths(end) * 2 / sWormSegs;
 sAngles = curvature(skeleton, lfAngleEdgeLength, sCCLengths);
 
 % reverse skeleton to skeleton_sw to make it good for plot
 skeleton_sw = [mm-skeleton(:,2),skeleton(:,1)];
 % find the central point of the skeleton
 hf_ske_sw = skeleton_sw(round(length(skeleton)/2),:);
   
 % To show the skeleton
 Img_ske = zeros(mm,nn);
 skeleton_sw_1d = sub2ind(size(Img_ske),  skeleton_sw(:,1),  skeleton_sw(:,2));
 % as the red colour
 Img_ske(skeleton_sw_1d) = 1;
 Segout(logical(Img_ske)) = 0.8;
 
%  % draw the central point of skeleton as a green point
%  Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,2) = 1;
%  Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,1) = 0.2;
%  Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,3) = 0.2;

k=1;

%  Seg skeleton from tail to head
Frenet_Pt{k}.xy = skeleton([end:-samp_step:1,1],:);
% Seg skeleton from head to tail
Frenet_Pt{k}.xy_flp = skeleton([1:samp_step:end,end],:);  

%% Frenet Transform
% adjust the axis to make it convenient to plot
Frenet_Pt{k}.xy = [Frenet_Pt{k}.xy(:,1), mm-Frenet_Pt{k}.xy(:,2)];
Frenet_Pt{k}.xy_flp = [Frenet_Pt{k}.xy_flp(:,1), mm-Frenet_Pt{k}.xy_flp(:,2)];
    
% if the detected skeleton change oritation, change it back 
if k>=3
    est_tail = 2*Frenet_Pt{k-1}.xy(1,:)-Frenet_Pt{k-2}.xy(1,:);
    if  abs(abs(Frenet_Pt{k-1}.xy(end,:)-est_tail))<abs(abs(Frenet_Pt{k-1}.xy(1,:)-est_tail))
        temp_pt = Frenet_Pt{k}.xy;
        Frenet_Pt{k}.xy = Frenet_Pt{k}.xy_flp;
        Frenet_Pt{k}.xy_flp = temp_pt;
    end
end
    
% Frenet Transform
[TT,NN,B,k_fre,t_fre,Frenet_Pt{k}.T,Frenet_Pt{k}.N] = frenet(Frenet_Pt{k}.xy(:,1),Frenet_Pt{k}.xy(:,2));

% reverse skeleton to skeleton_sw to make it good for plot
skeleton_subsamp_sw = [Frenet_Pt{k}.xy(:,1),Frenet_Pt{k}.xy(:,2)];
% find the central point of the skeleton
hf_ske_index = round(length(skeleton_subsamp_sw)/2);
hf_ske_subsamp_sw = skeleton_subsamp_sw(hf_ske_index,:);

 % draw the central point of skeleton as a green point
 hf_ske_sw = [Frenet_Pt{k}.xy(hf_ske_index,2), Frenet_Pt{k}.xy(hf_ske_index,1)];
 Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,2) = 1;
 Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,1) = 0.2;
 Segout(hf_ske_sw(1)-1:hf_ske_sw(1)+1,hf_ske_sw(2)-1:hf_ske_sw(2)+1,3) = 0.2;
    
%% Show head, tail and the tangent and perpendicular vector at the central worm
fold = 10;

% figure;
imshow(Segout);
line(Frenet_Pt{k}.xy(:,1),Frenet_Pt{k}.xy(:,2)), hold on
quiver(Frenet_Pt{k}.xy(hf_ske_index,1),Frenet_Pt{k}.xy(hf_ske_index,2),Frenet_Pt{k}.T(hf_ske_index,1),Frenet_Pt{k}.T(hf_ske_index,2),0.3*fold,'color','y')
quiver(Frenet_Pt{k}.xy(hf_ske_index,1),Frenet_Pt{k}.xy(hf_ske_index,2),Frenet_Pt{k}.N(hf_ske_index,1),Frenet_Pt{k}.N(hf_ske_index,2),20*fold,'color','g')

%% Save the Frenet for future use
Frenet_Pt{2}=Frenet_Pt{1};

% Save data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
frenet_name = ['.\Data_source\Frenet(',hdf5_file,')-',date,'.mat' ];

%save .\Data_source\frenet(6544)_1906 Frenet_Pt CMs
save(frenet_name, 'Frenet_Pt', 'CMs','Segout');

%fname = ['test_',date,'.avi' ]; %%%% Please change them accordingly!!!!
% Create AVI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% fname = ['test_', hdf5_file, '.avi'];
% workingDir = [path,'Sample_Video\hdf5'];

% Create_Mov(masks, fname, frame, workingDir);



