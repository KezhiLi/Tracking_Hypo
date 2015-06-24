function texture = ske2tex(xy, width, Y_k_gray)
% calculate texture information based on video frame and skeleton
% The result is an N*2 matrix, in which each row represents a circle whose central point is on the skeleton. 
% First column represents the mean of points in the circle, and second
% column respents the variance. 
% 
% Input: xy: coordinates of points on skeleton
%        width: the width of the worm
%        Y_k_gray: matrix of frame Y_k   
% Output: texture: N*2 matrix, with 1st column means and 2nd column
% variance
%
% Copyrighit: author: Kezhi Li, CSC, MRC, Imperial College, London
% 16/06/2015
% You will not remove any copyright or other notices from the Software;
% you must reproduce all copyright notices and other proprietary
% notices on any copies of the Software.


% number of points on skeleton
num_pt = size(xy,1);

% Y_k_gray = 255 - rgb2gray(Y_k);

% neglect first 2 and last 2 points
neglect_pt = 2; 

num_step = (num_pt-2*neglect_pt)*2-1;

% new points set on skeleton
pt_on_ske = zeros(num_step,2);

% copy points on skeleton to the new points set
for ii = 1:num_pt-neglect_pt*2;
    pt_on_ske(ii*2-1,:) = xy(ii+neglect_pt,:);
end

% insert points between original adjacent points
for ii = 2:2:size(pt_on_ske,1);
    pt_on_ske(ii,:) = (pt_on_ske(ii-1,:)+pt_on_ske(ii+1,:))/2;
end


stat_pt = zeros(num_step,2);


for jj = 1:num_step;
    % create a circle mask
    mask_cir = create_cir(pt_on_ske(jj,2), pt_on_ske(jj,1), width, size(Y_k_gray));
    % circle in a full image
    cir_full = double(Y_k_gray).* mask_cir;
    % select circle with pixels > 0
    cir = cir_full(cir_full>0);
    % calculate mean of points in circle
    stat_pt(jj,1) = mean2(cir);
    % calculate variance of points in circle
    stat_pt(jj,2) = std2(cir);
end

%output statistical property 
texture = stat_pt;



