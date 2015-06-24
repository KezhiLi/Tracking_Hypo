function [para_prop, width_prop, D_adjust] = sel_prop_width(X, Y_gray, width, seg_len,  para_thre)
% select proper width value based on current frame and skeleton
%
% Input: X: a cell, X{1},X{2}... represent hypotheses 
%        Y_gray: a matrix, current frame in gray, usually in Uint8 format
%        width: a scalar, estimated width, test width will be near this
%        value
%        seg_len: a scalar, length of each segment of skeleton
%        para_thre: a scalar, estimated threshold parameter, near 1   
%
% Output: para_prop: a scalar, proper threshold parameter estimated  
%         width_prop: a scalar, proper width estimated   
% 
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 19/06/2015
% You will not remove any copyright or other notices from the Software; 
% you must reproduce all copyright notices and other proprietary 
% notices on any copies of the Software.

width_default = 4.5; 
para_default = 1;

Large_num = 10000;

% size of image
Npix_h = size(Y_gray, 1);
Npix_w = size(Y_gray, 2);

% test range for threshold parameter
num_para_test = 7;
% test parameter start from 
para_start = para_thre - (num_para_test-1)/2*0.05;
% test range for width
num_width_test = 7;
% test width start from 
width_start = width - (num_width_test-1)/2*0.5;

D = zeros(num_para_test,num_width_test);

%% image processing to thresholding

% imclose
se = strel('disk', 1);
II = imclose(Y_gray, se);

% initialize zero matrix, to save processing time
area_hypo_ori = logical(zeros(Npix_h,0)*zeros(0,Npix_w));

% test all parameter
for ii = 1:num_para_test;
    % convert ii to corresponding parameter value
    para_test = para_start + (ii-1)*0.05;
    % we set an artificial parameter 0.0.8 here % coil: 0.9  normal: 0.8
    level = graythresh(II)* para_test;
    BW = im2bw(II,level);
    %figure, imshow(BW)

    reBW = BW;
    %reBW = imerode(BW,se);  %reBW = BW;
    % reBW = -im2bw(I,level)+1;
    % figure, imshow(reBW)

    BWdfill =  logical(reBW);

    C = bwareaopen(BWdfill, 50); 
    %figure, imshow(BW2);
    %title('remove small areas');

    % 't','ts' for contour
    m_fre_pt = size(X{1}.xy,1);
    t = 1:2*m_fre_pt ;
    ts = 1:1/(seg_len*2):2*m_fre_pt;


    %% Compare new frame and hypotheses

    % test width
    for jj = 1:num_width_test;
            % convert jj to correspoinding test width value
            width_test = width_start+0.5*(jj-1);

            % 'tt','tts' for all points inside the contour
            div = floor(width_test + 1);
            total_seg = (div+1)*2-1;
            tt = 1:total_seg*m_fre_pt ;
            tts = 1:1/(seg_len*2):total_seg*m_fre_pt;

            area_hypo = area_hypo_ori;

            % skeleton to points on contour/points in body
            [worm_contour1,worm_body] = ske2shape(X{1}.xy, X{1}.N, width_test, -0.40);
            % points on contour to contour
            [worm_shape_x, worm_shape_y] = shape2curv(worm_contour1, 1, t,ts,size(area_hypo));  % 0.5
            % points in body to area in body     
            [worm_body_x, worm_body_y] = shape2curv(worm_body, 1, tt, tts,size(area_hypo));     % 0.7
            % indexes of all points on contour 
            area_hypo_1d = sub2ind(size(area_hypo),   worm_shape_y , worm_shape_x );
            % indexes of all points inside the contour (body)        
            area_hypo_1d_inside = sub2ind(size(area_hypo),   worm_body_y , worm_body_x );

            I = (min(worm_shape_y) >= 1 & max(worm_shape_y) <= Npix_h);
            J = (min(worm_shape_x) >= 1 & max(worm_shape_x) <= Npix_w);

            if I && J
                area_hypo(area_hypo_1d) =  logical(1);

                area_hypo(area_hypo_1d_inside) = logical(1);

                BW2_hypo = area_hypo;

                diff_abs = imabsdiff(C,BW2_hypo);

                diff_abs_op = imerode(diff_abs,se);
                D(ii,jj) = sum(sum(diff_abs+diff_abs_op))/sum(area_hypo(:));

            else
                 D(ii,jj) = Large_num;
            end
    end
end

% add more weights according to the distance between estimated para/width with the tested para/width, 
% in which the central point is added by 0, then adjecent points are added
% more weights. In this case we are willing to use the estimated
% para/width, except for the fact that other para/width values are much
% better.
[xx,yy] = ndgrid(1:num_para_test,1:num_width_test);
xx = abs(xx - round((num_para_test+1)/2));
yy = abs(yy - round((num_width_test+1)/2));
D_adjust = D + (xx+yy)*0.03;

% find the minimal D element and find its location
[D_min, D_min_pos] = min(D_adjust(:));
[D_min_row, D_min_col] = ind2sub(size(D_adjust),D_min_pos);

% convert the location of matrix to parameter and widht value
para_prop =  para_start + (D_min_row-1)*0.05;
width_prop = width_start + 0.5*(D_min_col-1);

