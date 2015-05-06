function [L,C, II] = calc_log_likelihood_Worm_1st(Xstd_rgb, XX, Y, width, seg_len,  para_thre)

% function to calculate the log likelihood of hypotheses
% Input:    Xstd_rgb: a scalar, the various of the image set in advance
%                 XX: an Npop_particles * sub_num cell, record all hypotheses data
%                  Y: a matrix, current frame captured
%              width: a scalar, the worm's width
% Output:  L:  a vector record the log likelihood
%          C:  a matrix, the image after thresholding
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College London,
% 24/02/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

% size of image
Npix_h = size(Y, 1);
Npix_w = size(Y, 2);

% number of paricles and sub-particles
Npop_particles = size(XX,1);
sub_num = size(XX,2);

% the total number of particles
N = Npop_particles*sub_num;

% likelihood vector
L = zeros(1,N);

% overlap penalty
op = 30;

%% image processing to thresholding
Y = rgb2gray(Y);
II = 255 - Y;

% imclose
se = strel('disk', 1);
II = imclose(II, se);

%figure, hist(double(I),256); 

% we set an artificial parameter 0.0.8 here % coil: 0.9  normal: 0.8
level = graythresh(II)* para_thre;
BW = im2bw(II,level);
%figure, imshow(BW)

reBW = BW;
%reBW = imerode(BW,se);  %reBW = BW;
% reBW = -im2bw(I,level)+1;
% figure, imshow(reBW)

log_reBW = logical(reBW);

BWdfill = log_reBW;
%BWdfill = imfill(log_reBW, 'holes');
%figure, imshow(BWdfill);
%title('binary image with filled holes');

C = bwareaopen(BWdfill, 50); 
%figure, imshow(BW2);
%title('remove small areas');


%% likelihood equation

A = -log(sqrt(2 * pi) * Xstd_rgb);
B = - 0.5 / (Xstd_rgb.^2);

    % 't','ts' for contour
    m_fre_pt = size(XX{1,1}.xy,1);
    t = 1:2*m_fre_pt ;
    ts = 1:1/(seg_len*2):2*m_fre_pt;
    
    % 'tt','tts' for all points inside the contour
    div = floor(width + 1);
    total_seg = (div+1)*2-1;
    tt = 1:total_seg*m_fre_pt ;
    tts = 1:1/(seg_len*2):total_seg*m_fre_pt;
    
    img_ratio = sum(sum(C))/700; 

%% Compare new frame and hypotheses

area_hypo_ori = logical(zeros(Npix_h,0)*zeros(0,Npix_w));

max_diff=1;

for ii = 1:Npop_particles;
    for jj = 1:sub_num;

        area_hypo = area_hypo_ori;
 
        % skeleton to points on contour/points in body
        [worm_contour1,worm_body] = ske2shape(XX{ii,jj}.xy, XX{ii,jj}.N, width, -0.40);
        % points on contour to contour
        [worm_shape_x, worm_shape_y] = shape2curv(worm_contour1, 1, t,ts,size(area_hypo));  % 0.5
        % points in body to area in body     
        [worm_body_x, worm_body_y] = shape2curv(worm_body, 1, tt, tts,size(area_hypo));     % 0.7
        % indexes of all points on contour 
        area_hypo_1d = sub2ind(size(area_hypo),   worm_shape_y , worm_shape_x );
        % indexes of all points inside the contour (body)        
        area_hypo_1d_inside = sub2ind(size(area_hypo),   worm_body_y , worm_body_x );
        
        % number of points on the skeleton
        ske_num = (m_fre_pt-1) * (2*seg_len) + 1;
        % curve of skeleton, it was obtained from the shape_x and shape_y
        ske_curv = [ worm_shape_x(end-ske_num+1:end), worm_shape_y(end-ske_num+1:end)];
        % compare the first 1/3 segment and last 1/3 segment part of the worm, save the
        % number of points on the skeleton that are repeated in both
        % segement. The number is used in calculating the difference. 
        third_ske_num = round(ske_num/3);
        tail_1third_curv = ske_curv(1:third_ske_num,:);
        head_1third_curv = ske_curv(end-third_ske_num+1:end,:);
        % comparing to 'intersect', the following method is faster to count
        % the number of duplicated points
        tail_1third_sorted = sort(tail_1third_curv);
        head_1third_sorted = sort(head_1third_curv);
        dup_points = tail_1third_sorted(ismember(tail_1third_sorted,head_1third_sorted,'rows'),:);
        
        if size(dup_points,1)>1
            dup_points
        end
  
        I = (min(worm_shape_y) >= 1 & max(worm_shape_y) <= Npix_h);
        J = (min(worm_shape_x) >= 1 & max(worm_shape_x) <= Npix_w);

        if I && J

            area_hypo(area_hypo_1d) =  logical(1);

            area_hypo(area_hypo_1d_inside) = logical(1);


            %area_hypo1 = imclose(area_hypo, se);
            log_reBW_hypo = area_hypo;

            % cancel it to save time
            % BWdfill_hypo = imfill(log_reBW_hypo,  fliplr(worm_contour1(round(size(worm_contour1,1)/2),:))  );

            %figure, imshow(BWdfill);
            %title('binary image with filled holes');

            BW2_hypo = log_reBW_hypo;
            %BW2_hypo = bwareaopen(log_reBW_hypo, 50); 
            %figure, imshow(BW2);
            %title('remove small areas');            

            % Calculate the difference, the key step
            
            diff_abs = imabsdiff(C,BW2_hypo);
            %diff_abs_op = imopen(diff_abs,se);
            diff_abs_op = imerode(diff_abs,se);
            D = (sum(sum(diff_abs+diff_abs_op*3))+size(dup_points,1)*op)/img_ratio;  
            %D = (sum(sum(diff_abs+diff_abs_op*3)/2)+size(dup_points,1)*op)/img_ratio;  
            %D = (sum(sum(diff_abs))+size(dup_points,1)*op)/img_ratio;  
            
            D2 = D^2;
            
            XX{ii,jj}.D = D;
            
            % Calculate the likelihood
            L((ii-1)*sub_num+jj) =  A + B * D2;
        else
            %L(k) = 0;
            %L(k) = -Inf;
            L((ii-1)*sub_num+jj) =  -Inf;
        end
    end
end
